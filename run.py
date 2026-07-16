import os
import pandas as pd
import shutil
import argparse
import subprocess
import time
from Bio import SeqIO
from pandas.core.frame import DataFrame

parser = argparse.ArgumentParser(description='PairAMG: host-context-aware interpretation on roles of viral auxiliary metabolic genes\n')
parser.add_argument('--input_viral', '-viral', type=str, help='root of viral genomes (.fasta/ .fa) or viral KO lists (.txt)', default='')
parser.add_argument('--input_host', '-host', type=str, help='root of host genomes (.fasta/ .fa) or host KO lists (.txt)', default='')
parser.add_argument('--input_viral_host_link_file', '-link', type=str, help='path of viral-host link file (.csv)', default='')
parser.add_argument('--database', '-d', type=str, help='database directory', default='./database/')
parser.add_argument('--mode', '-m', type=int, help='mode of running PairAMG (0: end-to-end 1: only viral AMG identification 2: only viral AMG function interpretation', default=0)
parser.add_argument('--input_type', '-t', type=int, help='type of input (0: require purification on viral AMG candidates 1: without purification on viral AMG candidates ', default=0)
parser.add_argument('--output', '-o', type=str, help='output directory', default='./result/')

inputs = parser.parse_args()
viral_root = inputs.input_viral
host_root = inputs.input_host
link_path = inputs.input_viral_host_link_file
root_database = inputs.database
mode = inputs.mode
input_type = inputs.input_type
result_root = inputs.output

if not os.path.exists(result_root):
    os.makedirs(result_root)

'''
========================================================
               Required functions
========================================================
'''
def step_split(str_module):
    step = []
    list_module = str_module.split(' ')
    temp = ''
    i = 0
    while i < len(list_module):
        if '(' in list_module[i] or ')' in list_module[i]:
            temp = list_module[i]
            j = 1
            while temp.count('(') != temp.count(')'):
                temp += ' '
                temp += list_module[i + j]
                j += 1
            step.append(temp)
            i += j
        else:
            step.append(list_module[i])
            i += 1
    return step


def substep_split(str_module):
    step = []
    list_module = str_module.split(',')
    temp = ''
    i = 0
    while i < len(list_module):
        if '(' in list_module[i] or ')' in list_module[i]:
            temp = list_module[i]
            j = 1
            while temp.count('(') != temp.count(')'):
                temp += ','
                temp += list_module[i + j]
                j += 1
            step.append(temp)
            i += j
        else:
            step.append(list_module[i])
            i += 1
    return step


def find_match_branch(step, KOs):
    flag = False
    branches = step[1:-1].split(',')
    for branch in branches:
        if ' ' in branch:
            keys = branch.split(' ')
            subflag = True
            for key in keys:
                if key not in KOs:
                    subflag = False
            if subflag:
                flag = True
        else:
            if '-' in branch:
                if len(branch.split('-')[-1])>6 and '+' in branch.split('-')[-1]:
                    branch = branch.split('-')[0] + branch.split('-')[-1][6:]
                else:
                    branch = branch.split('-')[0]
            if '+' in branch:
                if find_match_combination(branch, KOs):
                    flag = True
            else:
                if branch in KOs:
                    flag = True
    return flag


def find_match_combination(step, KOs):
    if '(' in step:
        keys = step.split('+')
        temp = []
        i = 0
        while i<len(keys):
            if '(' in keys[i]:
                if keys[i].count('(') == 1 and keys[i].count(')') == 1:
                    temp.append(keys[i])
                    i += 1
                else:
                    temp.append(keys[i]+'+'+keys[i+1])
                    i += 2
            else:
                temp.append(keys[i])
                i += 1
        keys = temp
        for key in keys:
            if key.startswith('('):
                if not find_match_branch(key, KOs):
                    return False
            else:
                if '-' in key:
                    if len(key.split('-')[-1])>6 and '+' in key.split('-')[-1]:
                        key = key.split('-')[0] + key.split('-')[-1][6:]
                    else:
                        key = key.split('-')[0]
                if key not in KOs:
                    return False
        return True

    else:
        keys = step.split('+')
        for key in keys:
            if '-' in key:
                if len(key.split('-')[-1])>6 and '+' in key.split('-')[-1]:
                    key = key.split('-')[0] + key.split('-')[-1][6:]
                else:
                    key = key.split('-')[0]
            if key not in KOs:
                return False
        return True


def find_match_step(step, KOs):
    if step.startswith('('):
        substeps = substep_split(step[1:-1])
        flag = False
        for substep in substeps:
            if substep.startswith('('):
                if substep.count('(')==1 and substep.count(')')==1:
                    if ' ' in substep:
                        subflag = True
                        keys = step_split(substep)
                        for key in keys:
                            if key.startswith('('):
                                if not find_match_branch(key, KOs):
                                    subflag = False
                            else:
                                if key not in KOs:
                                    subflag = False
                        if subflag:
                            flag = True
                    else:
                        if ',' in substep:
                            if find_match_branch(substep, KOs):
                                flag = True
                        else:
                            if '-' in substep:
                                substep = substep.replace('(', '')
                                substep = substep.replace(')', '')
                                if len(substep.split('-')[-1])>6 and '+' in substep.split('-')[-1]:
                                    substep = substep.split('-')[0] + substep.split('-')[-1][6:]
                                else:
                                    substep = substep.split('-')[0]
                                if substep in KOs:
                                    flag = True
                            else:
                                keys = step_split(substep[1:-1])
                                subflag = True
                                for key in keys:
                                    if key not in KOs:
                                        subflag = False
                                if subflag:
                                    flag = True
                elif substep.count('(') > 1 and substep.count(')') > 1:
                    index_split = substep.find(' ')
                    if index_split>0 and substep[index_split-1]==')' and substep[index_split+1]=='(':
                        if substep[0]=='(' and substep[1]=='(' and substep[-1]==')' and substep[-2]==')':
                            keys = step_split(substep[1:-1])
                        else:
                            keys = step_split(substep)
                        subflag = True
                        for key in keys:
                            if ' ' in key:
                                subsubflag = True
                                subkeys = step_split(key)
                                for subkey in subkeys:
                                    if subkey.startswith('('):
                                        if not find_match_branch(subkey, KOs):
                                            subsubflag = False
                                    else:
                                        if subkey not in KOs:
                                            subsubflag = False
                                if not subsubflag:
                                    subflag = False
                                    
                            if key.startswith('('):
                                if not find_match_branch(key, KOs):
                                    subflag = False
                            else:
                                if key not in KOs:
                                    subflag = False
                        if subflag:
                            flag = True
                    else:
                        if substep.split(')')[0].count('(')>1:
                            keys = substep_split(substep[1:-1])
                        else:
                            keys = substep_split(substep)
                        subflag = False
                        for key in keys:
                            if ' ' in key:
                                subsubflag = True
                                subkeys = step_split(key)
                                for subkey in subkeys:
                                    if subkey.startswith('('):
                                        if not find_match_branch(subkey, KOs):
                                            subsubflag = False
                                    else:
                                        if subkey not in KOs:
                                            subsubflag = False
                                if subsubflag:
                                    subflag = True
                            if '+' in key:
                                subsubflag = True
                                subkeys = key.split('+')
                                temp = []
                                i = 0
                                while i < len(subkeys):
                                    if '(' in subkeys[i]:
                                        if subkeys[i].count('(') == 1 and subkeys[i].count(')') == 1:
                                            temp.append(subkeys[i])
                                            i += 1
                                        else:
                                            temp.append(subkeys[i] + '+' + subkeys[i + 1])
                                            i += 2
                                    else:
                                        temp.append(subkeys[i])
                                        i += 1
                                subkeys = temp
                                for subkey in subkeys:
                                    if subkey.startswith('('):
                                        if not find_match_branch(subkey, KOs):
                                            subsubflag = False
                                    else:
                                        if subkey not in KOs:
                                            subsubflag = False
                                if subsubflag:
                                    subflag = True
                            else:
                                if key.startswith('('):
                                    if find_match_branch(key, KOs):
                                        subflag = True
                                else:
                                    if key in KOs:
                                        subflag = True
                        if subflag:
                            flag = True
            else:
                if ' ' in substep:
                    subflag = True
                    keys = step_split(substep)
                    for key in keys:
                        if key.startswith('('):
                            if not find_match_branch(key, KOs):
                                subflag = False
                        else:
                            if key not in KOs:
                                subflag = False
                    if subflag:
                        flag = True
                else:
                    if '+' in substep:
                        if find_match_combination(substep, KOs):
                            flag = True
                    else:
                        if '-' in substep:
                            substep = substep.replace('(', '')
                            substep = substep.replace(')', '')
                            if len(substep.split('-')[-1])>6 and '+' in substep.split('-')[-1]:
                                substep = substep.split('-')[0] + substep.split('-')[-1][6:]
                            else:
                                substep = substep.split('-')[0]
                        if substep in KOs:
                            flag = True
        return flag
    else:
        if '+' in step:
            return find_match_combination(step, KOs)
        else:
            if step.startswith('-'):
                return True
            if step in KOs:
                return True
            else:
                return False


'''
========================================================
                Checking required files
========================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Initialization of PairAMG') 

if not os.path.exists(root_database + 'KEGG_sequences.faa'):
    print('Error: Missing KEGG gene database')
    exit(1)
if not os.path.exists(root_database + 'cazyme_protease.fasta'):
    print('Error: Missing cazyme_protease database')
    exit(1)
if not os.path.exists(root_database + 'vogdb.faa'):
    print('Error: Missing vogdb database')
    exit(1)
if not os.path.exists(root_database + 'virus.fasta'):
    print('Error: Missing NCBI Refseq Virus protein database')
    exit(1)
if not os.path.exists(root_database + 'metabolism.txt'):
    print('Error: Missing KO function category file')
    exit(1)
if not os.path.exists(root_database + 'modules.txt'):
    print('Error: Missing metabolic pathway annotation file')
    exit(1)
if not os.path.exists(root_database + 'sequences_KO.txt'):
    print('Error: Missing KEGG gene annotation file')
    exit(1)

if mode == 0:
    for _, _, filelists in os.walk(viral_root):
        if len(filelists) == 0:
            print('Error: Empty root for viral genomes')
            exit(1)
    for _, _, filelists in os.walk(host_root):
        if len(filelists) == 0:
            print('Error: Empty root for host genomes')
            exit(1)
    if not os.path.exists(link_path):
        print('Error: Missing viral-host link file')
        exit(1)
    if result_root == '':
        print('Error: Not specified root for output')
        exit(1)
        
elif mode == 1:
    for _, _, filelists in os.walk(viral_root):
        if len(filelists) == 0:
            print('Error: Empty root for viral genomes')
            exit(1)
    if result_root == '':
        print('Error: Not specified root for output')
        exit(1)
        
elif mode == 2:
    for _, _, filelists in os.walk(viral_root):
        if len(filelists) == 0:
            print('Error: Empty root for viral KO sets')
            exit(1)
    for _, _, filelists in os.walk(host_root):
        if len(filelists) == 0:
            print('Error: Empty root for host KO sets')
            exit(1)
    if not os.path.exists(link_path):
        print('Error: Missing viral-host link file')
        exit(1)
    if result_root == '':
        print('Error: Not specified root for output')
        exit(1)
        
else:
    print('Error: invalid mode')
    exit(1)


print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Initialization successful')
print('-------------------------------------------------------------------------------------------------\n')


'''
=======================================================================
                Viral AMGs function interpretation
=======================================================================
'''

if mode==2:
    print('-------------------------------------------------------------------------------------------------')
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Interpretating function of viral AMGs')
    
    data = open('./database/modules.txt', 'r').readlines()
    dict_modules = {}
    for line in data:
        module_index = line.rstrip().split('\t')[0]
        module = line.rstrip().split('\t')[2]
        dict_modules[module_index] = module
    
    dict_module_definition = {}
    for line in data:
        module_index = line.rstrip().split('\t')[0]
        module_definition = line.rstrip().split('\t')[1]
        dict_module_definition[module_index] = module_definition
    
    output = []

    list_data = pd.read_csv(link_path, header=None)
    for i in range(len(list_data)):
        viral = list_data.iloc[i,0]
        host = list_data.iloc[i,1]

        viral_KO_path = viral_root + viral + '.txt'
        host_KO_path = host_root + host + '.txt'

        if os.path.exists(viral_KO_path) and os.path.exists(host_KO_path):
    
            viral_KO = []
            host_KO = []
            
            data = open(viral_KO_path, 'r').readlines()
            for line in data:
                viral_KO.append(line.rstrip().split('\t')[0])
                
            data = open(host_KO_path, 'r').readlines()
            for line in data:
                host_KO.append(line.rstrip().split('\t')[0])
            
            union_KO = list(set(viral_KO + host_KO))
                
            sorted_module = []
            for module_index, module in dict_modules.items():
                for item in viral_KO:
                    if item in module:
                        sorted_module.append(module_index)
            sorted_module = list(set(sorted_module))
            
            for module in sorted_module:
                try:
                    equation_module = dict_modules[module].strip()
                    viral_candidate = []
                    host_candidate = []
                    for item in viral_KO:
                        if item in equation_module:
                            viral_candidate.append(item)
                    for item in host_KO:
                        if item in equation_module:
                            host_candidate.append(item)
                            
                    if 'M' in equation_module:
                        continue
                        
                    if '-- ' in equation_module:
                        equation_module = equation_module.replace('-- ', '')
                    
                    if ' ' not in equation_module:
                        equation_module = '(' + equation_module + ')'
                        step_flag = [0]
                        step_host = []
                        step_viral = []
                        step_union = []
                        
                        if find_match_step(equation_module, viral_KO):
                            step_viral.append(0)
                        if find_match_step(equation_module, host_KO):
                            step_host.append(0)
                        if find_match_step(equation_module, union_KO):
                            step_union.append(0)
                        
                        if step_union==step_flag:
                            if len(step_viral)>0:
                                if set(step_viral).issubset(step_host):
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    if set(viral_candidate).issubset(host_candidate):
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                elif list(set(step_host + step_viral))==step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                if step_host!=step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if len(step_viral)>0:
                                if set(step_viral).issubset(step_host):
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    if set(viral_candidate).issubset(host_candidate):
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    min_index = min(step_union)
                                    max_index = max(step_union)
                                    continuous_flag = list(range(min_index, max_index+1))
                                    if continuous_flag == step_union:
                                        if list(set(step_host + step_viral))==step_union:
                                            step_host = list(map(str, step_host))
                                            step_viral = list(map(str, step_viral))
                                            step_union = list(map(str, step_union))
                                            step_flag = list(map(str, step_flag))
                                            output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 1)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                        else:
                                            # identify the inclusion relationship between union KO set and viral KO set and host KO set
                                            step_host = list(map(str, step_host))
                                            step_viral = list(map(str, step_viral))
                                            step_union = list(map(str, step_union))
                                            step_flag = list(map(str, step_flag))
                                            output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                if step_host!=step_union:
                                    
                                    min_index = min(step_union)
                                    max_index = max(step_union)
                                    continuous_flag = list(range(min_index, max_index+1))
                                    if continuous_flag == step_union:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        continue
                            
                        
                    step_host = []
                    step_viral = []
                    step_union = []
                    
                    steps = step_split(equation_module)
                    step_flag = list(range(0, len(steps)))
                    
                    for i in range(len(steps)):
                        if find_match_step(steps[i], viral_KO):
                            step_viral.append(i)
                        if find_match_step(steps[i], host_KO):
                            step_host.append(i)
                        if find_match_step(steps[i], union_KO):
                            step_union.append(i)
        
                    if step_union==step_flag:
                        if len(step_viral)>0:
                            if set(step_viral).issubset(step_host):
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                if set(viral_candidate).issubset(host_candidate):
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            elif list(sorted(set(step_host + step_viral)))==step_union:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if step_host!=step_union:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                    else:
                        if len(step_viral)>0:
                            if set(step_viral).issubset(step_host):
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                if set(viral_candidate).issubset(host_candidate):
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                min_index = min(step_union)
                                max_index = max(step_union)
                                continuous_flag = list(range(min_index, max_index+1))
                                if continuous_flag == step_union:
                                    if list(sorted(set(step_host + step_viral)))==step_union:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 1)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if step_host!=step_union:
                                # identify the step set corresponding to union KO set is continuous or not
                                min_index = min(step_union)
                                max_index = max(step_union)
                                continuous_flag = list(range(min_index, max_index+1))
                                if continuous_flag == step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (No supporteds)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                
                except Exception as e:
                    print('Find Error:')
                    traceback.print_exc()
                    print(viral)
                    print(host)
                    print(module)
                    print('\n')


    save_path = result_root + 'module_completeness_pattern.txt'
    f = open(save_path, 'w')
    for line in output:
        f.write(line)
    f.close()

    if len(output)>0:
        data = pd.read_csv(save_path, sep='\t', header=None)
        data.columns = ['Viral', 'Host', 'Module', 'Equation_module', 'Definition_module', 'Pattern', 'Pathway_blocks', 'Viral_supported_blocks', 'Host_supported_blocks', 'Union_supported_blocks']
        save_path = result_root + 'module_completeness_pattern.tsv'
        data.to_csv(save_path, sep='\t', index=False, header=True)
    else:
        data = pd.DataFrame(columns = ['Viral', 'Host', 'Module', 'Equation_module', 'Definition_module', 'Pattern', 'Pathway_blocks', 'Viral_supported_blocks', 'Host_supported_blocks', 'Union_supported_blocks'])
        save_path = result_root + 'module_completeness_pattern.tsv'
        data.to_csv(save_path, sep='\t', index=False, header=True)
     
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
    print('-------------------------------------------------------------------------------------------------\n')
    
    exit(1)

'''
=======================================================================
                Viral AMG candidates identification
=======================================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Identifying viral AMG candidates')


if not os.path.exists('./tmp/'):
    os.makedirs('./tmp/')

checkv_db = ''
for _, subdirs, _ in os.walk('./database/'):
    for subdir in subdirs:
        if subdir.startswith('checkv-db'):
            checkv_db = './database/' + subdir 

if checkv_db == '':
    print('Error: CheckV database missing or unreadable')
    exit(1)

output = []
for _, _, filelists in os.walk(viral_root):
    for file in filelists:
        path = viral_root + file
        records = list(SeqIO.parse(path, "fasta"))
        output.append(records[0])     
    
SeqIO.write(output, './tmp/viral_contigs.fa', "fasta")

    
viral_path = './tmp/viral_contigs.fa'
save_path = './tmp/protein/checkv/'

if not os.path.exists(save_path):
    os.makedirs(save_path)

command = 'checkv contamination ' + viral_path + ' ' + save_path + ' -t 16 -d ' + checkv_db
os.system(command)

gene_path = './tmp/protein/checkv/tmp/proteins.faa'
output_root = './tmp/diamond/viral/'

if not os.path.exists(output_root):
    os.makedirs(output_root)
    
output_path = output_root + 'viral.matches.txt'
command = 'diamond blastp -d ./db/reference -q ' + gene_path + ' -o ' + output_path + ' --sensitive'
os.system(command)

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')    
    
    
'''
====================================================================
         Viral AMG candidates purification preprocessing 
====================================================================
'''
    
print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Align viral AMG candidates with curated viral-benefit gene databases')

# 1. Align with NCBI Refseq Virus protein database and curated viral-benefit gene databases

gene_path = './tmp/protein/checkv/tmp/proteins.faa'
output_root = './tmp/diamond/viral/'
    
output_path = output_root + 'viral.matches.refseq.txt'
command = 'diamond blastp -d ./db/refseqvirus -q ' + gene_path + ' -o ' + output_path + ' --sensitive'
os.system(command)
    
output_path = output_root + 'viral.matches.cazyme.txt'
command = 'diamond blastp -d ./db/cazyme_protease -q ' + gene_path + ' -o ' + output_path + ' --sensitive'
os.system(command)
    
output_path = output_root + 'viral.matches.vogdb.txt'
command = 'diamond blastp -d ./db/vogdb -q ' + gene_path + ' -o ' + output_path + ' --sensitive'
os.system(command)

# 2. Preprocess checkV output

path = './tmp/protein/checkv/tmp/gene_features.tsv'
save_root = './tmp/viral_marker/checkv/'

if not os.path.exists(save_root):
    os.makedirs(save_root)

data = pd.read_csv(path, sep='\t')

seq = data.iloc[0,0]
start_index = 0
for i in range(1,len(data)):
    if data.iloc[i,0] !=seq:
        temp = data[start_index:i]
        save_path = save_root + seq + '.csv'
        temp.to_csv(save_path, header=True, index=False)
        start_index = i
        seq = data.iloc[i,0]
        
temp = data[start_index:]
save_path = save_root + seq + '.csv'
temp.to_csv(save_path, header=True, index=False)

# 3. Preprocess diamond output of NCBI RefSeq Virus protein database

dict_gene_length = {}
path = './database/virus.fasta'
records = list(SeqIO.parse(path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_gene_length:
        dict_gene_length[gene] = length

path = './tmp/diamond/viral/viral.matches.refseq.txt'
save_root = './tmp/viral_marker/refseq/'

if not os.path.exists(save_root):
    os.makedirs(save_root)
            
sorted_match = []
data = open(path, 'r').readlines()
seq = ('_').join(data[0].split('\t')[0].split('_')[:-1])

dict_query_gene = {}
gene_path = './tmp/protein/checkv/tmp/proteins.faa'
records = list(SeqIO.parse(gene_path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_query_gene:
        dict_query_gene[gene] = length 
             
for i in range(1,len(data)):
    if ('_').join(data[i].split('\t')[0].split('_')[:-1])==seq:
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
    else:
        if len(sorted_match)>0:
            save_path = save_root + seq + '_identity_50_coverage_50.txt'
            f = open(save_path, 'w')
            for line in sorted_match:
                f.write(line)
            f.close()  
            
        sorted_match = []
        seq = ('_').join(data[i].split('\t')[0].split('_')[:-1])  
        
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
            
if len(sorted_match)>0:
    save_path = save_root + seq + '_identity_50_coverage_50.txt'
        
    f = open(save_path, 'w')
    for line in sorted_match:
        f.write(line)
    f.close()

# 3. Preprocess diamond output of Cazyme and protease database  
 
dict_gene_length = {}
path = './database/cazyme_protease.fasta'
records = list(SeqIO.parse(path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_gene_length:
        dict_gene_length[gene] = length

path = './tmp/diamond/viral/viral.matches.cazyme.txt'
save_root = './tmp/viral_marker/cazyme/'  

if not os.path.exists(save_root):
    os.makedirs(save_root)
          
sorted_match = []
data = open(path, 'r').readlines()
seq = ('_').join(data[0].split('\t')[0].split('_')[:-1])

dict_query_gene = {}
gene_path = './tmp/protein/checkv/tmp/proteins.faa'
records = list(SeqIO.parse(gene_path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_query_gene:
        dict_query_gene[gene] = length 
             
for i in range(1,len(data)):
    if ('_').join(data[i].split('\t')[0].split('_')[:-1])==seq:
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
    else:
        if len(sorted_match)>0:
            save_path = save_root + seq + '_identity_50_coverage_50.txt'
            f = open(save_path, 'w')
            for line in sorted_match:
                f.write(line)
            f.close()  
            
        sorted_match = []
        seq = ('_').join(data[i].split('\t')[0].split('_')[:-1])  
        
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
            
if len(sorted_match)>0:
    save_path = save_root + seq + '_identity_50_coverage_50.txt'
        
    f = open(save_path, 'w')
    for line in sorted_match:
        f.write(line)
    f.close()

# 4: Proprocess diamond output of vogdb database
  
dict_gene_length = {}
path = './database/vogdb.faa'
records = list(SeqIO.parse(path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_gene_length:
        dict_gene_length[gene] = length

path = './tmp/diamond/viral/viral.matches.vogdb.txt'
save_root = './tmp/viral_marker/vogdb/'      

if not os.path.exists(save_root):
    os.makedirs(save_root)
      
sorted_match = []
data = open(path, 'r').readlines()
seq = ('_').join(data[0].split('\t')[0].split('_')[:-1])

dict_query_gene = {}
gene_path = './tmp/protein/checkv/tmp/proteins.faa'
records = list(SeqIO.parse(gene_path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_query_gene:
        dict_query_gene[gene] = length 
             
for i in range(1,len(data)):
    if ('_').join(data[i].split('\t')[0].split('_')[:-1])==seq:
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
    else:
        if len(sorted_match)>0:
            save_path = save_root + seq + '_identity_50_coverage_50.txt'
            f = open(save_path, 'w')
            for line in sorted_match:
                f.write(line)
            f.close()  
            
        sorted_match = []
        seq = ('_').join(data[i].split('\t')[0].split('_')[:-1])  
        
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
            
if len(sorted_match)>0:
    save_path = save_root + seq + '_identity_50_coverage_50.txt'
        
    f = open(save_path, 'w')
    for line in sorted_match:
        f.write(line)
    f.close()  
    
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')  

'''
========================================================================================================
        Remove viral AMG candidates with best hit against curated viral-benefit gene databases
========================================================================================================
'''

print('-------------------------------------------------------------------------------------------------')  
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Removing viral AMG candidates with best hit against curated viral-benefit gene databases')


# 1. sorting blast alignments according to coverage and identity
    
dict_gene_length = {}
path = './database/KEGG_sequences.faa'
records = list(SeqIO.parse(path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_gene_length:
        dict_gene_length[gene] = length

path = './tmp/diamond/viral/viral.matches.txt'
save_root = './tmp/KO/viral_identity/'

if not os.path.exists(save_root):
        os.makedirs(save_root)
        
sorted_match = []
data = open(path, 'r').readlines()
seq = ('_').join(data[0].split('\t')[0].split('_')[:-1])

dict_query_gene = {}
gene_path = './tmp/protein/checkv/tmp/proteins.faa'
records = list(SeqIO.parse(gene_path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_query_gene:
        dict_query_gene[gene] = length 
             
for i in range(1,len(data)):
    if ('_').join(data[i].split('\t')[0].split('_')[:-1])==seq:
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
    else:
        if len(sorted_match)>0:
            save_path = save_root + seq + '_identity_50_coverage_50.txt'
            f = open(save_path, 'w')
            for line in sorted_match:
                f.write(line)
            f.close()  
            
        sorted_match = []
        seq = ('_').join(data[i].split('\t')[0].split('_')[:-1])  
        
        if float(data[i].rstrip().split('\t')[2])<50:
            continue
        else:    
            query = data[i].rstrip().split('\t')[0]
            query_start = float(data[i].rstrip().split('\t')[6])
            query_end = float(data[i].rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = data[i].rstrip().split('\t')[1]
            match_start = float(data[i].rstrip().split('\t')[8])
            match_end = float(data[i].rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(data[i].rstrip() + '\t' + str(match_length) + '\n')
            
if len(sorted_match)>0:
    save_path = save_root + seq + '_identity_50_coverage_50.txt'
        
    f = open(save_path, 'w')
    for line in sorted_match:
        f.write(line)
    f.close()  

# 2: annotate all the sorted alignments with KO index
data = open('./database/sequences_KO.txt', 'r').readlines()
dict_KO = {}
for line in data:
    gene = line.rstrip().split('\t')[0]
    KO = line.rstrip().split('\t')[1]
    if gene not in dict_KO:
        dict_KO[gene] = KO

old_root = './tmp/KO/viral_identity/'
save_root = './tmp/KO/viral_annotation/'

if not os.path.exists(save_root):
        os.makedirs(save_root)
        
for _, _, filelists in os.walk(old_root):
    for file in filelists:
        old_path = old_root + file
        
        data = open(old_path, 'r').readlines()
        output = []
        for line in data:
            match = line.rstrip().split('\t')[1]
            KO = dict_KO[match]
            output.append(line.rstrip() + '\t' + KO + '\n')
        
        if len(output)==0:
            continue
            
        save_path = save_root + file.replace('.txt', '_annotation.txt')
        f = open(save_path, 'w')
        for line in output:
            f.write(line)
        f.close()


# 3: sorted the best hit
old_root = './tmp/KO/viral_annotation/'
save_root = './tmp/KO/viral_sorted/'

if not os.path.exists(save_root):
        os.makedirs(save_root)
        
for _, _, filelists in os.walk(old_root):
    for file in filelists:
        old_path = old_root + file        
        data = open(old_path, 'r').readlines()
        
        sorted_data = []
        dict_KO_count = {}
        dict_KO_diamond = {}
        match_query = ''
        for line in data:
            query = line.rstrip().split('\t')[0]
            if query == match_query:
                match_KO = line.rstrip().split('\t')[-1]
                if match_KO in dict_KO_count:
                    dict_KO_count[match_KO] += 1
                else:
                    dict_KO_count[match_KO] = 1
                    dict_KO_diamond[match_KO] = line
            else:
                if len(dict_KO_count)==0:
                    match_query = query
                    match_KO = line.rstrip().split('\t')[-1]
                    dict_KO_count[match_KO] = 1
                    dict_KO_diamond[match_KO] = line
                else:
                    ranking = sorted(dict_KO_count.items(), key=lambda x:x[1])
                    best_counts = 0
                    best_identity = 0
                    KO = ''
                    for key, value in ranking:
                        if value >= best_counts:
                            identity = float(dict_KO_diamond[key].split('\t')[2])
                            if identity > best_identity:
                                best_identity = identity
                                KO = key
                    sorted_data.append(dict_KO_diamond[KO])
                    
                    match_query = query
                    match_KO = line.rstrip().split('\t')[-1]
                    dict_KO_count = {}
                    dict_KO_diamond = {}
                    dict_KO_count[match_KO] = 1
                    dict_KO_diamond[match_KO] = line
                  
        ranking = sorted(dict_KO_count.items(), key=lambda x:x[1])
        best_counts = 0
        best_identity = 0
        KO = ''
        for key, value in ranking:
            if value >= best_counts:
                identity = float(dict_KO_diamond[key].split('\t')[2])
                if identity > best_identity:
                    best_identity = identity
                    KO = key
        sorted_data.append(dict_KO_diamond[KO])
                
        if len(sorted_data)==0:
            continue
                            
        save_path = save_root + file.replace('.txt', '_sorted.txt')        
        f = open(save_path, 'w')
        for line in sorted_data:
            f.write(line)
        f.close()
       
# 4: check the best hits from curated viral-benefit gene databases
old_root = './tmp/KO/viral_sorted/'
save_root = result_root + 'viral_KO/'
cazyme_root = './tmp/viral_marker/CAZYme/'
vogdb_root = './tmp/viral_marker/vogdb/'

if not os.path.exists(save_root):
        os.makedirs(save_root)
        
for _, _, filelists in os.walk(old_root):
    for file in filelists:
    
        dict_cazyme = {}
        dict_vogdb = {} 
        cazyme_path = cazyme_root + file.replace('_annotation_sorted.txt', '.txt')
        if os.path.exists(cazyme_path):
            data = open(cazyme_path, 'r').readlines()
            for line in data:
                protein = line.rstrip().split('\t')[0] 
                identity = float(line.rstrip().split('\t')[2])
                if protein not in dict_cazyme:
                    dict_cazyme[protein] = identity
                else:
                    dict_cazyme[protein] = max(identity, dict_cazyme[protein]) 
                
        vogdb_path = vogdb_root + file.replace('_annotation_sorted.txt', '.txt')
        if os.path.exists(vogdb_path):
            data = open(vogdb_path, 'r').readlines()
            for line in data:
                protein = line.rstrip().split('\t')[0] 
                identity = float(line.rstrip().split('\t')[2])
                if protein not in dict_cazyme:
                    dict_vogdb[protein] = identity
                else:
                    dict_vogdb[protein] = max(identity, dict_cazyme[protein]) 

        old_path = os.path.join(old_root, file)
        data = open(old_path, 'r').readlines()
        
        output = []  
        for line in data:
            protein = line.rstrip().split('\t')[0]
            identity = float(line.rstrip().split('\t')[2])
            if protein in dict_cazyme and identity<dict_cazyme[protein]:
                continue
            elif protein in dict_vogdb and identity<dict_vogdb[protein]:
                continue
            else:
                output.append(line)
                
        if len(output)==0:
            continue
        
        save_path = save_root + file.replace('.txt', '_besthit.txt')
        f = open(save_path, 'w')
        for line in output:
            f.write(line)
        f.close()

# 5: summarize the level C function category of each KO
dict_function = {}
function = ''
data = open('./database/metabolism.txt', 'r').readlines()
for line in data[3:]:
    if line.startswith('C'):
        function = line.split('C ')[1].split('[')[0].strip()
    elif line.startswith('D'):
        KO = 'K' + line.split('K')[1][:5]
        if KO not in dict_function:
            dict_function[KO] = [function]
        else:
            if function not in dict_function[KO]:
                dict_function[KO].append(function)

old_root = result_root + 'viral_KO/'
save_root = result_root + 'viral_KO_summary/'

if not os.path.exists(save_root):
        os.makedirs(save_root)
        
for _, _, filelists in os.walk(old_root):
    for file in filelists:
    
        old_path = old_root + file              
        data = open(old_path, 'r').readlines()

        count_KO = {}
        for line in data:
            
            KO = line.rstrip().split('\t')[-1]
            if KO not in count_KO:
                count_KO[KO] = 1
            else:
                count_KO[KO] += 1
        
        count_KO = sorted(count_KO.items(), key=lambda x:x[1], reverse=True)
        
        output = []    
        for key, value in count_KO:
            output.append(key + '\t' + str(value) + '\t' + ('\t').join(dict_function[key]) + '\n')
        
        if len(output)==0:
            continue
   
        save_path = save_root + file.replace('.txt', '_KO_summary.txt')
        f = open(save_path, 'w')
        for line in output:
            f.write(line)
        f.close()

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')  

'''
==================================================
        Purify viral AMG candidates 
==================================================
'''

if input_type==0:

    print('-------------------------------------------------------------------------------------------------')  
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Purifying identified viral AMG candidates')
    
    # 1: identify viral AMG candidates flanked by microbial hallmark gene on either side
    old_root = result_root + 'viral_KO/'
    save_root = result_root + 'viral_KO_no_microbial/'
    checkv_root = './tmp/viral_marker/checkv/'
    hallmark_root = './tmp/viral_marker/refseq/'
    
    if not os.path.exists(save_root):
                os.makedirs(save_root)
    
    for _, _, filelists in os.walk(old_root):
        for file in filelists:
            
            checkv_path = checkv_root + file.replace('_identity_50_coverage_50_annotation_sorted_besthit.txt', '.csv')
            hallmark_path = hallmark_root + file.replace('_identity_50_coverage_50_annotation_sorted_besthit.txt', '_identity_50_coverage_50.txt')
            old_path = old_root + file
            
            list_gene = pd.read_csv(checkv_path)
            if os.path.exists(hallmark_path):
                hallmarks = open(hallmark_path, 'r').readlines()
                for line in hallmarks:
                    index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                    list_gene.iloc[index, 6] = 1
                    
            data = open(old_path, 'r').readlines()
            for line in data:
                index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                list_gene.iloc[index, 6] = 0
            
            output = []
            for line in data:
                index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                left_microbial = False
                left_viral = False
                right_microbial = False
                right_viral = False
                
                for i in range(index, len(list_gene)):
                    if list_gene.iloc[i,6]==-1:
                        next_index = min(i+1, len(list_gene)-1)
                        if list_gene.iloc[next_index,6]==1:
                            right_viral = True
                        else:
                            right_microbial = True
                        break
                        
                    elif list_gene.iloc[i,6]==1:
                        right_viral = True
                        break
                    
                for i in range(index,-1,-1):
                    if list_gene.iloc[i,6]==-1:
                        next_index = max(i-1, 0)
                        if list_gene.iloc[next_index,6]==1:
                            left_viral = True
                        else:
                            left_microbial = True
                        break
                        
                    elif list_gene.iloc[i,6]==1:
                        left_viral = True
                        break
                        
                if left_microbial or right_microbial:
                    continue
                else:
                    output.append(line)
            
            if len(output)==0:
                continue
       
            save_path = save_root + file.replace('.txt', '_no_microbial.txt')
            f = open(save_path, 'w')
            for line in output:
                f.write(line)
            f.close()
            
    # 2: select viral AMG candidates flanked by viral hallmark gene on either side
    old_root = result_root + 'viral_KO/'
    save_root = result_root + 'viral_KO_category_1/'
    checkv_root = './tmp/viral_marker/checkv/'
    hallmark_root = './tmp/viral_marker/refseq/'
    
    if not os.path.exists(save_root):
        os.makedirs(save_root)
                
    for _, _, filelists in os.walk(old_root):
        for file in filelists:
            
            checkv_path = checkv_root + file.replace('_identity_50_coverage_50_annotation_sorted_besthit.txt', '.csv')
            hallmark_path = hallmark_root + file.replace('_identity_50_coverage_50_annotation_sorted_besthit.txt', '_identity_50_coverage_50.txt')
            old_path = old_root + file
            
            list_gene = pd.read_csv(checkv_path)
            if os.path.exists(hallmark_path):
                hallmarks = open(hallmark_path, 'r').readlines()
                for line in hallmarks:
                    index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                    list_gene.iloc[index, 6] = 1
                    
            data = open(old_path, 'r').readlines()
            for line in data:
                index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                list_gene.iloc[index, 6] = 0
            
            output = []
            for line in data:
                index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                left_microbial = False
                left_viral = False
                right_microbial = False
                right_viral = False
                
                for i in range(index, len(list_gene)):
                    if list_gene.iloc[i,6]==-1:
                        
                        next_index = min(i+1, len(list_gene)-1)
                        if list_gene.iloc[next_index,6]==1:
                            right_viral = True
                        else:
                            right_microbial = True
                        break
                       
                    elif list_gene.iloc[i,6]==1:
                        right_viral = True
                        break
                    
                for i in range(index,-1,-1):
                    if list_gene.iloc[i,6]==-1:
                        
                        next_index = max(i-1, 0)
                        if list_gene.iloc[next_index,6]==1:
                            left_viral = True
                        else:
                            left_microbial = True
                        break
                       
                    elif list_gene.iloc[i,6]==1:
                        left_viral = True
                        break
                        
                if left_microbial or right_microbial:
                    continue
                else:
                    if left_viral or right_viral:
                        output.append(line)
            
            if len(output)==0:
                continue
       
            save_path = save_root + file.replace('.txt', '_category_1.txt')
            f = open(save_path, 'w')
            for line in output:
                f.write(line)
            f.close()
            
    # 3: select viral AMG candidates flanked by viral hallmark gene on both sides
    old_root = result_root + 'viral_KO/'
    save_root = result_root + 'viral_KO_category_2/'
    checkv_root = './tmp/viral_marker/checkv/'
    hallmark_root = './tmp/viral_marker/refseq/'
    
    if not os.path.exists(save_root):
        os.makedirs(save_root)
                
    for _, _, filelists in os.walk(old_root):
        for file in filelists:
            
            checkv_path = checkv_root + file.replace('_identity_50_coverage_50_annotation_sorted_besthit.txt', '.csv')
            hallmark_path = hallmark_root + file.replace('_identity_50_coverage_50_annotation_sorted_besthit.txt', '_identity_50_coverage_50.txt')
            old_path = old_root + file
            
            list_gene = pd.read_csv(checkv_path)
            if os.path.exists(hallmark_path):
                hallmarks = open(hallmark_path, 'r').readlines()
                for line in hallmarks:
                    index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                    list_gene.iloc[index, 6] = 1
                    
            data = open(old_path, 'r').readlines()
            for line in data:
                index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                list_gene.iloc[index, 6] = 0
            
            output = []
            for line in data:
                index = int(line.rstrip().split('\t')[0].split('_')[-1]) - 1
                left_microbial = False
                left_viral = False
                right_microbial = False
                right_viral = False
                
                for i in range(index, len(list_gene)):
                    if list_gene.iloc[i,6]==-1:
                        
                        next_index = min(i+1, len(list_gene)-1)
                        if list_gene.iloc[next_index,6]==1:
                            right_viral = True
                        else:
                            right_microbial = True
                        break
                        
                    elif list_gene.iloc[i,6]==1:
                        right_viral = True
                        break
                    
                for i in range(index,-1,-1):
                    if list_gene.iloc[i,6]==-1:
                        
                        next_index = max(i-1, 0)
                        if list_gene.iloc[next_index,6]==1:
                            left_viral = True
                        else:
                            left_microbial = True
                        break
                        
                    elif list_gene.iloc[i,6]==1:
                        left_viral = True
                        break
                        
                if left_microbial or right_microbial:
                    continue
                else:
                    if left_viral and right_viral:
                        output.append(line)
            
            if len(output)==0:
                continue
       
            save_path = save_root + file.replace('.txt', '_category_2.txt')
            f = open(save_path, 'w')
            for line in output:
                f.write(line)
            f.close()
            
    # 4: summarize the level C function category of viral AMG candidates flanked by viral hallmark gene on either side
    
    dict_function = {}
    function = ''
    data = open('./database/metabolism.txt', 'r').readlines()
    for line in data[3:]:
        if line.startswith('C'):
            function = line.split('C ')[1].split('[')[0].strip()
        elif line.startswith('D'):
            KO = 'K' + line.split('K')[1][:5]
            if KO not in dict_function:
                dict_function[KO] = [function]
            else:
                if function not in dict_function[KO]:
                    dict_function[KO].append(function)
    
    old_root = result_root + 'viral_KO_category_1/'
    save_root = result_root + 'viral_KO_category_1_summary/'
    
    if not os.path.exists(save_root):
        os.makedirs(save_root)
    
    for _, _, filelists in os.walk(old_root):
        for file in filelists:
        
            old_path = old_root + file              
            data = open(old_path, 'r').readlines()
    
            count_KO = {}
            for line in data:
                
                KO = line.rstrip().split('\t')[-1]
                if KO not in count_KO:
                    count_KO[KO] = 1
                else:
                    count_KO[KO] += 1
            
            count_KO = sorted(count_KO.items(), key=lambda x:x[1], reverse=True)
            
            output = []    
            for key, value in count_KO:
                output.append(key + '\t' + str(value) + '\t' + ('\t').join(dict_function[key]) + '\n')
            
            if len(output)==0:
                continue
       
            save_path = save_root + file.replace('.txt', '_KO_summary.txt')
            f = open(save_path, 'w')
            for line in output:
                f.write(line)
            f.close()
    
    # 4: summarize the level C function category of viral AMG candidates flanked by viral hallmark gene on both sides
    
    old_root = result_root + 'viral_KO_category_2/'
    save_root = result_root + 'viral_KO_category_2_summary/'
    
    if not os.path.exists(save_root):
        os.makedirs(save_root)
    
    for _, _, filelists in os.walk(old_root):
        for file in filelists:
        
            old_path = old_root + file              
            data = open(old_path, 'r').readlines()
    
            count_KO = {}
            for line in data:
                
                KO = line.rstrip().split('\t')[-1]
                if KO not in count_KO:
                    count_KO[KO] = 1
                else:
                    count_KO[KO] += 1
            
            count_KO = sorted(count_KO.items(), key=lambda x:x[1], reverse=True)
            
            output = []    
            for key, value in count_KO:
                output.append(key + '\t' + str(value) + '\t' + ('\t').join(dict_function[key]) + '\n')
            
            if len(output)==0:
                continue
       
            save_path = save_root + file.replace('.txt', '_KO_summary.txt')
            f = open(save_path, 'w')
            for line in output:
                f.write(line)
            f.close()
        
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
    print('-------------------------------------------------------------------------------------------------')  


# Only perform viral AMG identification

if mode==1:
    exit(1)
    

'''
=======================================================================
            Romoving duplicated viral contigs in linked MAGs 
=======================================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Removing duplicated viral contigs in linked MAGs')

# 1: Copy host genomes for duplication

save_root = './tmp/MAG_with_removal/'
if not os.path.exists(save_root):
    os.makedirs(save_root)
    
for _, _, filelists in os.walk(host_root):
    for file in filelists:
        if file.endswith('.fa'):
            old_path = host_root + file
            save_path = save_root + file.replace('.fa','.fasta')
            shutil.copy(old_path, save_path)
        elif file.endswith('.fasta'):
            old_path = host_root + file
            save_path = save_root + file
            shutil.copy(old_path, save_path)
        else:
            print('Error: Host genome ' + file + ' unreadable')
            exit(1)

# 2: Copy viral genomes for duplication

save_root = './tmp/remove_duplicated/phage/'
if not os.path.exists(save_root):
    os.makedirs(save_root)
    
for _, _, filelists in os.walk(viral_root):
    for file in filelists:
        if file.endswith('.fa'):
            old_path = viral_root + file
            save_path = save_root + file
            shutil.copy(old_path, save_path)
        elif file.endswith('.fasta'):
            old_path = viral_root + file
            save_path = save_root + file.replace('.fasta','.fa')
            shutil.copy(old_path, save_path)
        else:
            print('Error: Viral genome ' + file + ' unreadable')
            exit(1)

# 3: Run blastn to find duplicated viral genomes in host genomes
   
phage_root = './tmp/remove_duplicated/phage/'
save_root = './tmp/remove_duplicated/blast/'
host_root = './tmp/MAG_with_removal/'

if not os.path.exists(save_root):
    os.makedirs(save_root)
    
data = pd.read_csv(link_path, header=None)
for i in range(len(data)):
    phage = data.iloc[i,0]
    host = data.iloc[i,1]
    
    phage_path = phage_root + phage + '.fa'
    host_path = host_root + host + '.fasta'
    if os.path.exists(phage_path) and os.path.exists(host_path):
        output_path = save_root + phage + '.txt'
        
        cmd = 'blastn -query ' + phage_path + ' -subject ' + host_path + ' -out ' + output_path + ' -max_hsps 1' + ' -outfmt ' + "'6 qseqid sseqid qlen slen pident length mismatch gaps qstart qend sstart send evalue bitscore'"
        subprocess.run(cmd, shell=True)
        
    else:
        print('Error: Viral genome ' + phage + ' and host genome ' + host + ' missing or unreadable')
        exit(1)

save_root = './tmp/remove_duplicated/blast/'
output = []
data = pd.read_csv(link_path, header=None)
for i in range(len(data)):
    phage = data.iloc[i,0]
    host = data.iloc[i,1]
            
    blast_path = save_root + phage + '.txt'
    if os.path.exists(blast_path):
        blast_data = open(blast_path, 'r').readlines()
        if len(blast_data)>0:
            blast_result = blast_data[0]
            if float(blast_result.split('\t')[4])==100.000 and blast_result.split('\t')[2]==blast_result.split('\t')[5]:
                output.append(host + '\t' + blast_result)
                
f = open('./tmp/remove_duplicated/blast_duplicated.txt', 'w')
for line in output:
    f.write(line)
f.close() 

shutil.rmtree('./tmp/remove_duplicated/phage/')

host_root = './tmp/MAG_with_removal/'
dict_remove = {}
dict_genome = {}

data = open('./tmp/remove_duplicated/blast_duplicated.txt', 'r').readlines()
for blast_result in data:
    phage_contig = blast_result.split('\t')[1]
        
    genome = blast_result.split('\t')[0]
    genome_index = blast_result.split('\t')[2]
    start = int(blast_result.split('\t')[11])
    end = int(blast_result.split('\t')[12])
    
    dict_genome[genome_index] = genome
    remove_index = [min(start, end)-1, max(start,end)]
    if genome_index in dict_remove:
        dict_remove[genome_index].append(remove_index)
    else:
        dict_remove[genome_index] = [remove_index]

for key, value in dict_remove.items():
    value.sort(reverse=True)
    host_path = host_root + dict_genome[key] + '.fasta'
    records = list(SeqIO.parse(host_path, "fasta"))
    for i in range(len(records)):
        if records[i].id == key:
            for item in value:
                start = item[0]
                end = item[1]
                if start==0 and end==len(records[i].seq):
                    del records[i]
                else:
                    if end== len(records[i].seq):
                        records[i].seq = records[i].seq[:start]
                    elif start==0:
                        records[i].seq = records[i].seq[end:]
                    else:
                        records[i].seq = records[i].seq[:start] + records[i].seq[end:]
            break
    SeqIO.write(records, host_path, "fasta")
    
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')

'''
=======================================================================
                Host AMG candidates identification
=======================================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Identifying host AMG candidates')

host_root = './tmp/MAG_with_removal/'
output_root = './tmp/protein/host/'
if not os.path.exists(output_root):
    os.makedirs(output_root)
    
for _, _, filelists in os.walk(host_root):
    for file in filelists:
        host_path = os.path.join(host_root, file)
        output_path = output_root + file.replace('.fasta', '.genes.faa')
        sentence = 'prodigal -p meta -i ' + host_path + ' -a ' + output_path + ' -q'
        os.system(sentence)

root = './tmp/MAG_with_removal/'
gene_root = './tmp/protein/host/'
output_root = './tmp/diamond/host/'

if not os.path.exists(output_root):
    os.makedirs(output_root)
    
for _, _, filelists in os.walk(root):
    for file in filelists:
        path = gene_root + file.replace('.fasta', '.genes.faa')
        output_path = output_root + file.replace('.fasta', '.matches.txt')
        if os.path.exists(output_path):
            continue
        sentence = 'diamond blastp -d ./db/reference -q ' + path + ' -o ' + output_path
        os.system(sentence)  

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')

'''
=======================================================================
            KO annotation on identified host AMG candidates
=======================================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Annotating host AMG candidates')

# 1: sorting blast alignments according to coverage and identity
dict_gene_length = {}
path = './database/KEGG_sequences.faa'
records = list(SeqIO.parse(path, "fasta"))
for i in range(len(records)):
    gene = records[i].id
    length = len(records[i].seq)
    if gene not in dict_gene_length:
        dict_gene_length[gene] = length

diamond_root = './tmp/diamond/host/'
gene_root = './tmp/protein/host/'
save_root = './tmp/KO/host_identity/'

if not os.path.exists(save_root):
    os.makedirs(save_root)
    
for _, _, filelists in os.walk(diamond_root):
    for file in filelists:
        
        # filter out alignments with identity>50 and coverage>50               
        sorted_match = []
        diamond_path = diamond_root + file
        data = open(diamond_path, 'r').readlines()
        if len(data)==0:
            continue
        
        dict_query_gene = {}
        gene_path = gene_root + file.replace('.matches.txt', '.genes.faa')
        records = list(SeqIO.parse(gene_path, "fasta"))
        for i in range(len(records)):
            gene = records[i].id
            length = len(records[i].seq)
            if gene not in dict_query_gene:
                dict_query_gene[gene] = length 
             
        for line in data:
            if float(line.rstrip().split('\t')[2])<=50:
                continue
                
            query = line.rstrip().split('\t')[0]
            query_start = float(line.rstrip().split('\t')[6])
            query_end = float(line.rstrip().split('\t')[7])
            query_length = dict_query_gene[query]
            
            match = line.rstrip().split('\t')[1]
            match_start = float(line.rstrip().split('\t')[8])
            match_end = float(line.rstrip().split('\t')[9])
            match_length = dict_gene_length[match]
            if abs(match_start - match_end) / match_length >= 0.5 and abs(query_start - query_end) / query_length >= 0.5:
                sorted_match.append(line.rstrip() + '\t' + str(match_length) + '\n')
        
        if len(sorted_match)==0:
            continue
            
        save_path = save_root + file.replace('.txt', '_identity_50_coverage_50.txt')
        f = open(save_path, 'w')
        for line in sorted_match:
            f.write(line)
        f.close() 


# 2: annotate all the sorted alignments with KO index 
data = open('./database/sequences_KO.txt', 'r').readlines()
dict_KO = {}
for line in data:
    gene = line.rstrip().split('\t')[0]
    KO = line.rstrip().split('\t')[1]
    if gene not in dict_KO:
        dict_KO[gene] = KO

old_root = './tmp/KO/host_identity/'
save_root = './tmp/KO/host_annotation/'

if not os.path.exists(save_root):
    os.makedirs(save_root)
    
for _, _, filelists in os.walk(old_root):
    for file in filelists:
        old_path = old_root + file
        
        data = open(old_path, 'r').readlines()
        output = []
        for line in data:
            match = line.rstrip().split('\t')[1]
            KO = dict_KO[match]
            output.append(line.rstrip() + '\t' + KO + '\n')
        
        if len(output)==0:
            continue
            
        save_path = save_root + file.replace('.txt', '_annotation.txt')
        f = open(save_path, 'w')
        for line in output:
            f.write(line)
        f.close()


# 3: select the best hit of alignments
old_root = './tmp/KO/host_annotation/'
save_root = result_root + 'host_KO/'

if not os.path.exists(save_root):
    os.makedirs(save_root)
    
for _, _, filelists in os.walk(old_root):
    for file in filelists:
        old_path = old_root + file        
        data = open(old_path, 'r').readlines()
        
        sorted_data = []
        dict_KO_count = {}
        dict_KO_diamond = {}
        match_query = ''
        for line in data:
            query = line.rstrip().split('\t')[0]
            if query == match_query:
                match_KO = line.rstrip().split('\t')[-1]
                if match_KO in dict_KO_count:
                    dict_KO_count[match_KO] += 1
                else:
                    dict_KO_count[match_KO] = 1
                    dict_KO_diamond[match_KO] = line
            else:
                if len(dict_KO_count)==0:
                    match_query = query
                    match_KO = line.rstrip().split('\t')[-1]
                    dict_KO_count[match_KO] = 1
                    dict_KO_diamond[match_KO] = line
                else:
                    ranking = sorted(dict_KO_count.items(), key=lambda x:x[1])
                    best_counts = 0
                    best_identity = 0
                    KO = ''
                    for key, value in ranking:
                        if value >= best_counts:
                            identity = float(dict_KO_diamond[key].split('\t')[2])
                            if identity > best_identity:
                                best_identity = identity
                                KO = key
                    sorted_data.append(dict_KO_diamond[KO])
                    
                    match_query = query
                    match_KO = line.rstrip().split('\t')[-1]
                    dict_KO_count = {}
                    dict_KO_diamond = {}
                    dict_KO_count[match_KO] = 1
                    dict_KO_diamond[match_KO] = line
                  
        ranking = sorted(dict_KO_count.items(), key=lambda x:x[1])
        best_counts = 0
        best_identity = 0
        KO = ''
        for key, value in ranking:
            if value >= best_counts:
                identity = float(dict_KO_diamond[key].split('\t')[2])
                if identity > best_identity:
                    best_identity = identity
                    KO = key
        sorted_data.append(dict_KO_diamond[KO])
                
        if len(sorted_data)==0:
            continue
                            
        save_path = save_root + file.replace('.txt', '_sorted.txt')        
        f = open(save_path, 'w')
        for line in sorted_data:
            f.write(line)
        f.close()


# 4: summarize the level C function category of each KO

dict_function = {}
function = ''
data = open('./database/metabolism.txt', 'r').readlines()
for line in data[3:]:
    if line.startswith('C'):
        function = line.split('C ')[1].split('[')[0].strip()
    elif line.startswith('D'):
        KO = 'K' + line.split('K')[1][:5]
        if KO not in dict_function:
            dict_function[KO] = [function]
        else:
            if function not in dict_function[KO]:
                dict_function[KO].append(function)
    
old_root = result_root + 'host_KO/'
save_root = result_root + 'host_KO_summary/'

if not os.path.exists(save_root):
    os.makedirs(save_root)
    
for _, _, filelists in os.walk(old_root):
    for file in filelists:
        old_path = old_root + file     
        data = open(old_path, 'r').readlines()
        
        count_KO = {}
        data = open(old_path, 'r').readlines()
        for line in data:
            KO = line.rstrip().split('\t')[-1]
            if KO not in count_KO:
                count_KO[KO] = 1
            else:
                count_KO[KO] += 1
        
        count_KO = sorted(count_KO.items(), key=lambda x:x[1], reverse=True)
        
        output = []    
        for key, value in count_KO:
            output.append(key + '\t' + str(value) + '\t' + ('\t').join(dict_function[key]) + '\n')
        
        if len(output)==0:
            continue
   
        save_path = save_root + file.replace('.txt', '_KO_summary.txt')
        f = open(save_path, 'w')
        for line in output:
            f.write(line)
        f.close()  

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')


'''
==============================================================
        Funtion interpretation of viral AMGs
==============================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Interpretating function of viral AMG')

data = open('./database/modules.txt', 'r').readlines()
dict_modules = {}
for line in data:
    module_index = line.rstrip().split('\t')[0]
    module = line.rstrip().split('\t')[2]
    dict_modules[module_index] = module

dict_module_definition = {}
for line in data:
    module_index = line.rstrip().split('\t')[0]
    module_definition = line.rstrip().split('\t')[1]
    dict_module_definition[module_index] = module_definition

output = []

if input_type==0:

    viral_root = result_root + 'viral_KO_category_1_summary/'
    host_root = result_root + 'host_KO_summary/'

    list_data = pd.read_csv(link_path, header=None)
    for i in range(len(list_data)):
        viral = list_data.iloc[i,0]
        host = list_data.iloc[i,1]

        viral_KO_path = viral_root + viral + '_identity_50_coverage_50_annotation_sorted_besthit_category_1_KO_summary.txt'
        host_KO_path = host_root + host + '.matches_identity_50_coverage_50_annotation_sorted_KO_summary.txt'

        if os.path.exists(viral_KO_path) and os.path.exists(host_KO_path):
    
            viral_KO = []
            host_KO = []
            
            data = open(viral_KO_path, 'r').readlines()
            for line in data:
                viral_KO.append(line.rstrip().split('\t')[0])
                
            data = open(host_KO_path, 'r').readlines()
            for line in data:
                host_KO.append(line.rstrip().split('\t')[0])
            
            union_KO = list(set(viral_KO + host_KO))
                
            sorted_module = []
            for module_index, module in dict_modules.items():
                for item in viral_KO:
                    if item in module:
                        sorted_module.append(module_index)
            sorted_module = list(set(sorted_module))
            
            for module in sorted_module:
                try:
                    equation_module = dict_modules[module].strip()
                    viral_candidate = []
                    host_candidate = []
                    for item in viral_KO:
                        if item in equation_module:
                            viral_candidate.append(item)
                    for item in host_KO:
                        if item in equation_module:
                            host_candidate.append(item)
                            
                    if 'M' in equation_module:
                        continue
                        
                    if '-- ' in equation_module:
                        equation_module = equation_module.replace('-- ', '')
                    
                    if ' ' not in equation_module:
                        equation_module = '(' + equation_module + ')'
                        step_flag = [0]
                        step_host = []
                        step_viral = []
                        step_union = []
                        
                        if find_match_step(equation_module, viral_KO):
                            step_viral.append(0)
                        if find_match_step(equation_module, host_KO):
                            step_host.append(0)
                        if find_match_step(equation_module, union_KO):
                            step_union.append(0)
                        
                        if step_union==step_flag:
                            if len(step_viral)>0:
                                if set(step_viral).issubset(step_host):
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    if set(viral_candidate).issubset(host_candidate):
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                elif list(set(step_host + step_viral))==step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                if step_host!=step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if len(step_viral)>0:
                                if set(step_viral).issubset(step_host):
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    if set(viral_candidate).issubset(host_candidate):
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    min_index = min(step_union)
                                    max_index = max(step_union)
                                    continuous_flag = list(range(min_index, max_index+1))
                                    if continuous_flag == step_union:
                                        if list(set(step_host + step_viral))==step_union:
                                            step_host = list(map(str, step_host))
                                            step_viral = list(map(str, step_viral))
                                            step_union = list(map(str, step_union))
                                            step_flag = list(map(str, step_flag))
                                            output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 1)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                        else:
                                            # identify the inclusion relationship between union KO set and viral KO set and host KO set
                                            step_host = list(map(str, step_host))
                                            step_viral = list(map(str, step_viral))
                                            step_union = list(map(str, step_union))
                                            step_flag = list(map(str, step_flag))
                                            output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                if step_host!=step_union:
                                    
                                    min_index = min(step_union)
                                    max_index = max(step_union)
                                    continuous_flag = list(range(min_index, max_index+1))
                                    if continuous_flag == step_union:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        continue
                            
                        
                    step_host = []
                    step_viral = []
                    step_union = []
                    
                    steps = step_split(equation_module)
                    step_flag = list(range(0, len(steps)))
                    
                    for i in range(len(steps)):
                        if find_match_step(steps[i], viral_KO):
                            step_viral.append(i)
                        if find_match_step(steps[i], host_KO):
                            step_host.append(i)
                        if find_match_step(steps[i], union_KO):
                            step_union.append(i)
        
                    if step_union==step_flag:
                        if len(step_viral)>0:
                            if set(step_viral).issubset(step_host):
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                if set(viral_candidate).issubset(host_candidate):
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            elif list(sorted(set(step_host + step_viral)))==step_union:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if step_host!=step_union:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                    else:
                        if len(step_viral)>0:
                            if set(step_viral).issubset(step_host):
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                if set(viral_candidate).issubset(host_candidate):
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                min_index = min(step_union)
                                max_index = max(step_union)
                                continuous_flag = list(range(min_index, max_index+1))
                                if continuous_flag == step_union:
                                    if list(sorted(set(step_host + step_viral)))==step_union:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 1)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if step_host!=step_union:
                                # identify the step set corresponding to union KO set is continuous or not
                                min_index = min(step_union)
                                max_index = max(step_union)
                                continuous_flag = list(range(min_index, max_index+1))
                                if continuous_flag == step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (No supporteds)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                
                except Exception as e:
                    print('Find Error:')
                    traceback.print_exc()
                    print(viral)
                    print(host)
                    print(module)
                    print('\n')
    
elif input_type==1:

    viral_root = result_root + 'viral_KO_summary/'
    host_root = result_root + 'host_KO_summary/'

    list_data = pd.read_csv(link_path, header=None)
    for i in range(len(list_data)):
        viral = list_data.iloc[i,0]
        host = list_data.iloc[i,1]

        viral_KO_path = viral_root + viral + '_identity_50_coverage_50_annotation_sorted_besthit_KO_summary.txt'
        host_KO_path = host_root + host + '.matches_identity_50_coverage_50_annotation_sorted_KO_summary.txt'

        if os.path.exists(viral_KO_path) and os.path.exists(host_KO_path):
    
            viral_KO = []
            host_KO = []
            
            data = open(viral_KO_path, 'r').readlines()
            for line in data:
                viral_KO.append(line.rstrip().split('\t')[0])
                
            data = open(host_KO_path, 'r').readlines()
            for line in data:
                host_KO.append(line.rstrip().split('\t')[0])
            
            union_KO = list(set(viral_KO + host_KO))
                
            sorted_module = []
            for module_index, module in dict_modules.items():
                for item in viral_KO:
                    if item in module:
                        sorted_module.append(module_index)
            sorted_module = list(set(sorted_module))
            
            for module in sorted_module:
                try:
                    equation_module = dict_modules[module].strip()
                    viral_candidate = []
                    host_candidate = []
                    for item in viral_KO:
                        if item in equation_module:
                            viral_candidate.append(item)
                    for item in host_KO:
                        if item in equation_module:
                            host_candidate.append(item)
                            
                    if 'M' in equation_module:
                        continue
                        
                    if '-- ' in equation_module:
                        equation_module = equation_module.replace('-- ', '')
                    
                    if ' ' not in equation_module:
                        equation_module = '(' + equation_module + ')'
                        step_flag = [0]
                        step_host = []
                        step_viral = []
                        step_union = []
                        
                        if find_match_step(equation_module, viral_KO):
                            step_viral.append(0)
                        if find_match_step(equation_module, host_KO):
                            step_host.append(0)
                        if find_match_step(equation_module, union_KO):
                            step_union.append(0)
                        
                        if step_union==step_flag:
                            if len(step_viral)>0:
                                if set(step_viral).issubset(step_host):
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    if set(viral_candidate).issubset(host_candidate):
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                elif list(set(step_host + step_viral))==step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                if step_host!=step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if len(step_viral)>0:
                                if set(step_viral).issubset(step_host):
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    if set(viral_candidate).issubset(host_candidate):
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    min_index = min(step_union)
                                    max_index = max(step_union)
                                    continuous_flag = list(range(min_index, max_index+1))
                                    if continuous_flag == step_union:
                                        if list(set(step_host + step_viral))==step_union:
                                            step_host = list(map(str, step_host))
                                            step_viral = list(map(str, step_viral))
                                            step_union = list(map(str, step_union))
                                            step_flag = list(map(str, step_flag))
                                            output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 1)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                        else:
                                            # identify the inclusion relationship between union KO set and viral KO set and host KO set
                                            step_host = list(map(str, step_host))
                                            step_viral = list(map(str, step_viral))
                                            step_union = list(map(str, step_union))
                                            step_flag = list(map(str, step_flag))
                                            output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                if step_host!=step_union:
                                    
                                    min_index = min(step_union)
                                    max_index = max(step_union)
                                    continuous_flag = list(range(min_index, max_index+1))
                                    if continuous_flag == step_union:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        continue
                            
                        
                    step_host = []
                    step_viral = []
                    step_union = []
                    
                    steps = step_split(equation_module)
                    step_flag = list(range(0, len(steps)))
                    
                    for i in range(len(steps)):
                        if find_match_step(steps[i], viral_KO):
                            step_viral.append(i)
                        if find_match_step(steps[i], host_KO):
                            step_host.append(i)
                        if find_match_step(steps[i], union_KO):
                            step_union.append(i)
        
                    if step_union==step_flag:
                        if len(step_viral)>0:
                            if set(step_viral).issubset(step_host):
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                if set(viral_candidate).issubset(host_candidate):
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            elif list(sorted(set(step_host + step_viral)))==step_union:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if step_host!=step_union:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (Complementary 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Complete (No supported blocks)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                    else:
                        if len(step_viral)>0:
                            if set(step_viral).issubset(step_host):
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                if set(viral_candidate).issubset(host_candidate):
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 1)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Alternative 2)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                min_index = min(step_union)
                                max_index = max(step_union)
                                continuous_flag = list(range(min_index, max_index+1))
                                if continuous_flag == step_union:
                                    if list(sorted(set(step_host + step_viral)))==step_union:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 1)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                    else:
                                        step_host = list(map(str, step_host))
                                        step_viral = list(map(str, step_viral))
                                        step_union = list(map(str, step_union))
                                        step_flag = list(map(str, step_flag))
                                        output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                        else:
                            if step_host!=step_union:
                                # identify the step set corresponding to union KO set is continuous or not
                                min_index = min(step_union)
                                max_index = max(step_union)
                                continuous_flag = list(range(min_index, max_index+1))
                                if continuous_flag == step_union:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Complementary 2)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                                else:
                                    step_host = list(map(str, step_host))
                                    step_viral = list(map(str, step_viral))
                                    step_union = list(map(str, step_union))
                                    step_flag = list(map(str, step_flag))
                                    output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (Supplementary)'+ '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                            else:
                                step_host = list(map(str, step_host))
                                step_viral = list(map(str, step_viral))
                                step_union = list(map(str, step_union))
                                step_flag = list(map(str, step_flag))
                                output.append(viral + '\t' + host + '\t' + module + '\t' + equation_module + '\t' + dict_module_definition[module] + '\t' + 'Incomplete (No supporteds)' + '\t' + (' ').join(step_flag) + '\t' + (' ').join(step_viral) + '\t' + (' ').join(step_host) + '\t' + (' ').join(step_union) + '\n')
                
                except Exception as e:
                    print('Error:')
                    traceback.print_exc()
                    print('Viral genome: ' + viral)
                    print('Host genome: ' + host)
                    print('Module: ' + module)
                    print('\n')
                    
else: 
    print('Error: invalid input type')
    exit(1)
    
save_path = result_root + 'module_completeness_pattern.txt'
f = open(save_path, 'w')
for line in output:
    f.write(line)
f.close()

if len(output)>0:
    data = pd.read_csv(save_path, sep='\t', header=None)
    data.columns = ['Viral', 'Host', 'Module', 'Equation_module', 'Definition_module', 'Pattern', 'Pathway_blocks', 'Viral_supported_blocks', 'Host_supported_blocks', 'Union_supported_blocks']
    save_path = result_root + 'module_completeness_pattern.tsv'
    data.to_csv(save_path, sep='\t', index=False, header=True)
else:
    data = pd.DataFrame(columns = ['Viral', 'Host', 'Module', 'Equation_module', 'Definition_module', 'Pattern', 'Pathway_blocks', 'Viral_supported_blocks', 'Host_supported_blocks', 'Union_supported_blocks'])
    save_path = result_root + 'module_completeness_pattern.tsv'
    data.to_csv(save_path, sep='\t', index=False, header=True)

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')