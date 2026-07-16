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
## 2. Quick start
  --input INPUT_FA
                        The name of your input file (FASTA format)
  --output OUTPUT_PTH
                        The path of the output directory
  --filename FILENAME
                        Custom name for output files (option)
  --database DATABASE
                        Model directory
  --len LEN
                        Predict only for sequences >= len bp (default: 500)
  --batch_size BATCH_SIZE
                        Batch size for prediction (default: 16)
  --threshold THRESHOLD
                        Threshold for prediction (default: 0.5)
  --force
                        Force overwrite of the output directory if it exists (option)
