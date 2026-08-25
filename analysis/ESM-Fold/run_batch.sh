#!/bin/bash
# set proxy and mirror (consistent with sub.sh)
export http_proxy=http://10.243.120.3:3128
export https_proxy=http://10.243.120.3:3128
export HF_ENDPOINT=https://hf-mirror.com

# load environment modules
module purge
module load cuda/13.0 miniforge3

# activate the conda environment
source activate esmfold

# run the batch prediction script
python run_batch.py
