# PairAMG
PairAMG: host-context-aware interpretation on roles of viral auxiliary metabolic genes

## 1. Setup environment
Note: We suggest install all the required packages using Conda
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
Note: Each viral sequence should be uploaded seperately and renamed with sequence id
* Upload the predicted link file (.csv)  
Note: The link file should not contain header and the first column refers to the viral sequence id and the second column refers to the host genome file name
## 3. Quick start
* Run PairAMG for end-to-end function interpretation of viral auxiliary metabolic genes
```Python
python [-viral VIRAL_ROOT] [-host HOST_ROOT] [-link LINK_PATH] [-m MODE] [-t INPUT_TYPE]
```
* Options: PairAMG supports different running modes for function interpretation of viral auxiliary genes
```Linux
-viral root of viral genomes (mode 1 and 2) or viral KO lists (mode 3)
```
## 4. Helps
* If you have any questions about the usage of PairAMG, please feel free to contact Wei Zou. (Email: [weizou-c@my.cityu.edu.hk](weizou-c@my.cityu.edu.hk))
