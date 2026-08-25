import argparse
import json
import warnings
import torch
import os
import sys
import yaml
import numpy as np
import pandas as pd
from torch import nn
from torch_geometric.loader import DataLoader
from numpy import nan
from typing import *
from tqdm import tqdm
from scipy.stats import spearmanr
from transformers import logging
from src.models import PLM_model, GNN_model
from src.data import build_mutant_dataset
from src.utils.utils import param_num

# set path
current_dir = os.getcwd()
sys.path.append(current_dir)
# ignore warning information
logging.set_verbosity_error()
warnings.filterwarnings("ignore")

amino_acids_type = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def label_row(rows, sequence, token_probs, offset_idx=1):
    s = []
    sep = ";"
    if ":" in rows:
        sep = ":"
    for row in rows.split(sep):
        if row.lower() == "wt":
            s.append(0)
            continue
        try:
            wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
        except Exception:
            print(f"Parsing error: row={row}, sequence={sequence}")
            return np.nan

        # Check if the index is valid
        if idx < 0 or idx >= len(sequence):
            print(f"Index out of range: idx={idx}, seq_len={len(sequence)}, row={row}")
            return np.nan

        # Check if the wild-type matches
        if sequence[idx] != wt:
            print(f"Wild-type mismatch: expected {wt}, found {sequence[idx]} at idx {idx}, row={row}")
            return np.nan

        try:
            wt_encoded = amino_acids_type.index(wt)
            mt_encoded = amino_acids_type.index(mt)
        except ValueError:
            print(f"Unknown amino acid: wt={wt}, mt={mt}, row={row}")
            return np.nan

        score = token_probs[idx, mt_encoded] - token_probs[idx, wt_encoded]
        s.append(score.item())

    return sum(s)


def predict(args, plm_model, gnn_model, loader):
    gnn_model.eval()
    softmax = nn.Softmax(dim=-1)
    result_dict = {"name": [], "count": [], args.score_name: []}

    with torch.no_grad():
        for data in loader:
            protein_name = data.protein_name[0]
            graph_data = plm_model(data)
            out, _ = gnn_model(graph_data)
            seq = "".join([amino_acids_type[i] for i in data.y])
            out = torch.log(softmax(out[:, :20]) + 1e-9)

            # Read the mutant data file
            mutant_file_tsv = os.path.join(args.mutant_dataset_dir, "DATASET", protein_name, f"{protein_name}.tsv")
            mutant_file_csv = os.path.join(args.mutant_dataset_dir, "DATASET", protein_name, f"{protein_name}.csv")
            if os.path.exists(mutant_file_tsv):
                mutant_df = pd.read_table(mutant_file_tsv)
            elif os.path.exists(mutant_file_csv):
                mutant_df = pd.read_csv(mutant_file_csv)
            else:
                raise ValueError(f"Invalid file: {mutant_file_tsv} or {mutant_file_csv}")

            # Compute the prediction score of the current model
            offset = 1
            mutant_df[args.score_name] = mutant_df[args.mutant_pos_col].apply(
                lambda x: label_row(x, seq, out.cpu().numpy(), offset)
            )

            # If all scores are NaN, skip this protein (do not save/update the file)
            if mutant_df[args.score_name].isna().all():
                print(f">>> {protein_name}: all mutations failed, skipping.")
                continue

            result_file = os.path.join(args.result_dir, protein_name + ".csv")

            # ========== Key change: append new column instead of overwriting the entire file ==========
            if os.path.exists(result_file):
                # File already exists -> read existing results and merge the new column
                existing_df = pd.read_csv(result_file)
                # Only add the new column if it does not already exist (avoid duplicates)
                if args.score_name not in existing_df.columns:
                    # Align by row index when merging (both DataFrames should have the same row order)
                    # Here we assume the mutation order is fixed, so assign directly
                    existing_df[args.score_name] = mutant_df[args.score_name].values
                else:
                    print(f"Warning: column {args.score_name} already exists in {result_file}, overwriting.")
                    existing_df[args.score_name] = mutant_df[args.score_name].values
                existing_df.to_csv(result_file, index=False)
                result = existing_df
            else:
                # File does not exist -> create a new file
                mutant_df.to_csv(result_file, index=False)
                result = mutant_df
            # ====================================================

            # Compute Spearman correlation coefficient (filter out NaN)
            valid = result[[args.mutant_score_col, args.score_name]].dropna()
            if len(valid) > 1:
                spearmanr_score = spearmanr(valid[args.mutant_score_col], valid[args.score_name]).correlation
                spearmanr_score = 0 if spearmanr_score is nan else spearmanr_score
            else:
                spearmanr_score = 0

            result_dict['count'].append(len(result))
            result_dict['name'].append(protein_name)
            result_dict[args.score_name].append(spearmanr_score)

            print(f">>> {protein_name}: {spearmanr_score}; mutant_num: {len(result)}")

    # Save the Spearman summary for all proteins
    if args.score_info is not None:
        if os.path.exists(args.score_info):
            total_result = pd.read_csv(args.score_info)
            total_result[args.score_name] = result_dict[args.score_name]
            total_result.to_csv(args.score_info, index=False)
        else:
            pd.DataFrame(result_dict).to_csv(args.score_info, index=False)

    print(f">>> {args.score_name} average spearmanr: {np.mean(result_dict[args.score_name])}\n")


def ensemble(args):
    print("----------------- Ensemble -----------------")
    result_files = os.listdir(args.result_dir)
    sp_scores = []
    for file in tqdm(result_files):
        result_file = os.path.join(args.result_dir, file)
        result_df = pd.read_csv(result_file)
        models_pred = [result_df[col].to_list() for col in result_df.columns if col.startswith("ProtSSN")]
        ensemble_pred = np.mean(models_pred, axis=0)
        result_df["ProtSSN_ensemble"] = ensemble_pred
        result_df.to_csv(result_file, index=False)
        sp_score = spearmanr(result_df[args.mutant_score_col], result_df["ProtSSN_ensemble"]).correlation
        sp_scores.append(sp_score)
    print(">>> Ensemble spearmanr: ", np.mean(sp_scores))


def prepare(args, dataset_name, k, h):
    # for build dataset
    args.mutant_name = f"{dataset_name}_k{k}"
    mutant_dataset = build_mutant_dataset(args)
    protein_names = mutant_dataset.protein_names
    print(f">>> Protein names: {protein_names}")
    mutant_loader = DataLoader(mutant_dataset, batch_size=1, shuffle=False)
    print(f">>> Number of proteins: {len(mutant_dataset)}")
    gnn_model = GNN_model(args)
    print(f">>> k{k}_h{h} {param_num(gnn_model)}")
    gnn_model_path = os.path.join(args.gnn_model_dir, f"protssn_k{k}_h{h}.pt")
    gnn_model.load_state_dict(torch.load(gnn_model_path))
    return args, mutant_loader, gnn_model


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gnn", type=str, default="egnn", help="gat, gcn, or egnn")
    parser.add_argument("--gnn_config", type=str, default="src/config/egnn.yaml", help="gnn config")
    parser.add_argument("--gnn_model_dir", type=str, default="model/", help="test model name")
    parser.add_argument("--gnn_model_name", type=str, default=None, nargs="+", help="test model name")

    parser.add_argument("--plm", type=str, default="facebook/esm2_t33_650M_UR50D", help="esm param number")
    parser.add_argument("--use_ensemble", action="store_true", help="use ensemble model")

    # dataset
    parser.add_argument("--mutant_dataset_dir", type=str, default="data/evaluation", help="mutation dataset")
    parser.add_argument("--mutant_name", type=str, default=None, help="name of mutation dataset")
    parser.add_argument("--mutant_pos_col", type=str, default="mutant", help="mutation column name")
    parser.add_argument("--mutant_score_col", type=str, default="score", help="the model output score column name")

    parser.add_argument("--score_info", type=str, default=None, help="the model output spearmanr score file")
    parser.add_argument("--result_dir", type=str, default="result/", help="the result output path")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = create_parser()
    args.gnn_config = yaml.load(open(args.gnn_config), Loader=yaml.FullLoader)[args.gnn]

    plm_model = PLM_model(args)
    args.plm_hidden_size = plm_model.model.config.hidden_size
    dataset_name = args.mutant_dataset_dir.split("/")[-1]
    os.makedirs(args.result_dir, exist_ok=True)

    for gnn in args.gnn_model_name:
        k, h = gnn.split("_")
        k, h = int(k[1:]), int(h[1:])
        print(f"--------------- ProtSSN k{k}_h{h} ---------------")
        assert k in [10, 20, 30], f"Invalid k: {k}"
        assert h in [512, 768, 1280], f"Invalid h: {h}"
        args.gnn_config["hidden_channels"] = h
        args.c_alpha_max_neighbors = k
        args.score_name = f"ProtSSN_k{k}_h{h}"
        args, mutant_loader, gnn_model = prepare(args, dataset_name, k, h)
        predict(args=args, plm_model=plm_model, gnn_model=gnn_model, loader=mutant_loader)
    if args.use_ensemble:
        ensemble(args)