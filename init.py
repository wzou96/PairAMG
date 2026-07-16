import os
import time

'''
=======================================================================
                Initialize required database 
=======================================================================
'''

print('-------------------------------------------------------------------------------------------------')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Initializing Required database')

command = 'checkv download_database ./database/'
os.system(command)

checkv_db = ''
for _, subdirs, _ in os.walk('./database/'):
    for subdir in subdirs:
        if subdir.startswith('checkv-db'):
            checkv_db = './database/' + subdir 
                
command = 'export CHECKVDB=' + checkv_db
os.system(command)


diamond_root = './db/'
if not os.path.exists(diamond_root):
    os.makedirs(diamond_root)
    
sequence_path = './database/KEGG_sequences.faa'
db_path = './db/reference'
command = 'diamond makedb --in ' + sequence_path + ' -d ' + db_path
os.system(command)

sequence_path = './database/cazyme_protease.fasta'
db_path = './db/cazyme_protease'
command = 'diamond makedb --in ' + sequence_path + ' -d ' + db_path
os.system(command)

sequence_path = './database/vogdb.faa'
db_path = './db/vogdb'
command = 'diamond makedb --in ' + sequence_path + ' -d ' + db_path
os.system(command)

sequence_path = './database/virus.fasta'
db_path = './db/refseqvirus'
command = 'diamond makedb --in ' + sequence_path + ' -d ' + db_path
os.system(command)

print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish')
print('-------------------------------------------------------------------------------------------------\n')   