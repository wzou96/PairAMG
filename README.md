# PairAMG
PairAMG: host-context-aware interpretation on roles of viral auxiliary metabolic genes

## 1. Setup environment
*Note:* We suggest install all the required packages using Conda
* Clone the repository to the local.
```Linux
git clone https://github.com/wzou96/PairAMG
cd PairAMG
```
* Install and activate environment for PairAMG.
```Linux
conda env create -f requirements.yaml -n PairAMG
conda activate PairAMG
```
* Download the required database and unzip the file to folder database/.
https://drive.google.com/file/d/18KxGKXL0agmXpyO6zBtw_OyqdlCKrRWm/view?usp=drive_link
```Linux
unzip database.zip
```
* Download the CheckV database.
```Python
python init.py
```
## 2. Prepare the required input files
* Upload the viral sequences (.fa/ .fasta) and host genomes (.fa/ .fasta) to different foloders (e.g. viral/ and MAG/)  
  
*Note:* Each viral sequence should be uploaded seperately and renamed with sequence id
```
./viral/
k141_1.fa
k141_8.fa
k141_22.fa
```
```
./MAG/
bin.1.fasta
bin.20.fasta
bin.33.fasta
```
* If the user only wants to perform function interpretation, upload the viral KO sets (.txt) and host KO sets (.txt) to different folder (e.g. viral/ and MAG/)  

*Note:* Each KO set file should be uploaded seperately and renamed with sequence id. In the KO set file, **each KO should be separated into different rows**.

```
./viral/k141_1.txt
K00001
K00002
K00003
```
```
./viral/k141_2.txt
K00011
K00022
K00033
```

* Upload the predicted link file (.csv)  
  
*Note:* The link file **should not contain header** and the first column refers to the viral sequence id and the second column refers to the host genome file name

```
./link.csv
k141_1   bin.2
k141_2   bin.4
k141_5   bin.8
```
## 3. Quick start
* Download the sample data and unzip the file to foloder sample/
* ```Linux
unzip sample.zip
```
* Run PairAMG for end-to-end function interpretation of viral auxiliary metabolic genes
```Python
python [-viral VIRAL_ROOT] [-host HOST_ROOT] [-link LINK_PATH] [-d DATABASE] [-m MODE] [-t INPUT_TYPE] [-O OUTPUT_ROOT]
```
* Options: PairAMG supports different running modes for function interpretation of viral auxiliary genes
```
  --input_viral VIRAL_ROOT
                                 Root of the viral sequences or KO sets
  --input_host HOST_ROOT
                                 Root of the host genomes or KO sets
  --input_link LINK_PATH
                                 Path of the viral-host link file
  --d DATABASE
                                 Root of the required database (default: ./database/)
  --m MODE
                                 Mode of running PairAMG 0: End-to-end 1: Only viral AMG identification 2: Only viral AMG function interpretation (default: 0)
  --t INPUT_TYPE
                                 Type of input viral sequences 0: Require viral AMG candidates purification 1: Withput viral AMG candidates purification (default: 0)
  --o OUTPUT_ROOT
                                 Root of output files (default: ./result/)
```
## 4. Output explanation
```
  --host_KO/                          Identified AMG candidates from each host genome
  --host_KO_summary/                  Level C function category of AMG candidates identified from host genomes
  --viral_KO/                         Identified AMG candidates from each viral sequences
  --viral_KO_summary/                 Level C function category of AMG candidates identified from viral sequences
  --viral_KO_no_microbial/            Purified viral AMG candidates that not flanked by microbial marker genes on either or both sides
  --viral_KO_category_1/              Purified viral AMG candidates that flanked by viral marker genes on either or both sides
  --viral_KO_category_1_summary/      Level C function category of purified viral AMG candidates
  --viral_KO_category_2/              Purified viral AMG candidates that flanked by viral marker genes on both sides
  --viral_KO_category_2_summary/      Level C function category of purified viral AMG candidates
  --module_completeness_pattern.tsv   Identified completeness and patterns of metabolic pathways related to viral AMGs
  --module_completeness_pattern.txt   Identified completeness and patterns of metabolic pathways related to viral AMGs
```
## 5. Helps
* If you have any questions about the usage of PairAMG, please feel free to contact Wei Zou. (Email: [weizou-c@my.cityu.edu.hk](weizou-c@my.cityu.edu.hk))
