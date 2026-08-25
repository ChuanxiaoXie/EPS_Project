#!/bin/bash
#SBATCH -J predict_muts
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --output=predict_%j.out

export http_proxy=http://10.243.120.3:3128
export https_proxy=http://10.243.120.3:3128
export HF_ENDPOINT=https://hf-mirror.com

module purge
module load cuda/12.8
module load nccl/2.26_cuda12.8_5090
conda activate SaProt

python saprot.py