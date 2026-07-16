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
* Download the CheckV database.
```Python
python init.py
```
