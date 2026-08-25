# select the models for ensemble prediction
DATASET=proteingym-benchmark
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --gnn_model_name k10_h512 k20_h512 k30_h512 k10_h768 k20_h768 k30_h768 k10_h1280 k20_h1280 k30_h1280 \
    --use_ensemble \
    --mutant_dataset_dir data/mutant_example/$DATASET \
    --result_dir result/$DATASET
