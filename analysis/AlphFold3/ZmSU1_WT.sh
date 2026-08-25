#!/bin/bash
module load alphafold/3.0.0
python run_alphafold.py \
     --json_path=ZmSU1_WT.json \
     --model_dir=/data/public/datasets/alphafold3/3.0.0//model \
     --db_dir=/data/public/datasets/alphafold3/3.0.0/dataset \
      --pdb_database_path="/data/public/datasets/alphafold3/3.0.0/dataset/pdb_2022_09_28_mmcif_files.tar" \
     --output_dir=${SLURM_JOB_ID}_8aw3_output

