#!/bin/bash
export http_proxy=http://10.243.120.3:3128
export https_proxy=http://10.243.120.3:3128
export HF_ENDPOINT=https://hf-mirror.com

module purge
module load cuda/12.8
module load nccl/2.26_cuda12.8_5090
source /data/home/scxj090/run/miniconda3/etc/profile.d/conda.sh
conda activate EVO2_new

# Show Python path
which python

# Run the batch processing script
# Process all data
python evo2_predict_advanced.py

# To process a subset or resume from a checkpoint, use the following commands:
# python batch_mutation.py --model_name evo2_7b --batch_size 32 --start_idx 0 --end_idx 1000 --output_excel mutation_results_part1.xlsx
# python batch_mutation.py --model_name evo2_7b --batch_size 32 --start_idx 1000 --end_idx 2000 --output_excel mutation_results_part2.xlsx
