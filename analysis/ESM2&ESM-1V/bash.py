#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import torch
import esm
from tqdm import tqdm


def parse_excel(args):
    df = pd.read_excel(args.input_excel, sheet_name=args.sheet_name, header=0)
    genes_data = []
    for _, row in df.iterrows():
        protein_id = row.iloc[0]
        if pd.isna(protein_id):
            continue
        protein_id = str(protein_id)

        seq = row.iloc[1]
        if pd.isna(seq):
            continue
        seq = str(seq).strip()

        muts = []
        # Modify here: start from index 2 up to the last column of this row
        for col in range(2, len(row)):
            mut = row.iloc[col]
            if pd.notna(mut) and str(mut).strip():
                muts.append(str(mut).strip())

        if muts:
            genes_data.append((protein_id, seq, muts))

    return genes_data

def validate_mutation(mut_str, sequence, offset_idx):
    """
    Parsing consistent with your source predict.py:
      wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    and verify sequence[idx] == wt (skip if mismatch)
    """
    mut_str = str(mut_str).strip()
    if len(mut_str) < 3:
        return None

    wt = mut_str[0]
    mt = mut_str[-1]
    try:
        idx = int(mut_str[1:-1]) - offset_idx
    except ValueError:
        return None

    if idx < 0 or idx >= len(sequence):
        return None
    if sequence[idx] != wt:
        return None

    return wt, idx, mt, mut_str


@torch.no_grad()
def score_gene_wt_marginals(model, alphabet, sequence, mutations, offset_idx):
    """
    Consistent with source predict.py:
      token_probs = log_softmax(model(batch_tokens)["logits"])
      score = token_probs[0, 1+idx, mt] - token_probs[0, 1+idx, wt]
    """
    batch_converter = alphabet.get_batch_converter()
    data = [("protein1", sequence)]
    _, _, batch_tokens = batch_converter(data)

    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()

    logits = model(batch_tokens)["logits"]  # [1, L, V]
    token_probs = torch.log_softmax(logits, dim=-1)

    valid_muts = []
    scores = []

    for mut in mutations:
        parsed = validate_mutation(mut, sequence, offset_idx)
        if parsed is None:
            continue
        wt, idx, mt, mut_str = parsed

        wt_idx = alphabet.get_idx(wt)
        mt_idx = alphabet.get_idx(mt)

        score = token_probs[0, 1 + idx, mt_idx].item() - token_probs[0, 1 + idx, wt_idx].item()
        valid_muts.append(mut_str)
        scores.append(score)

    return valid_muts, scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-excel", required=True)
    parser.add_argument("--sheet-name", default="Sheet1")
    parser.add_argument("--output-csv", default="all_mutations_predictions.csv")
    parser.add_argument("--model-location", nargs="+", required=True, help="List of ESM model names")
    parser.add_argument("--offset-idx", type=int, default=1, help="offset_idx consistent with source predict.py (1-based position)")
    parser.add_argument("--nogpu", action="store_true")
    args = parser.parse_args()

    if args.nogpu:
        torch.cuda.is_available = lambda: False

    # Consistent with your previous version: truncate at 1022 to avoid exceeding common token limit
    MAX_SEQ_LEN = 1022

    print(f"Reading {args.input_excel} (sheet={args.sheet_name})")
    genes_data = parse_excel(args)
    print(f"Total proteins: {len(genes_data)}")

    # key: ProteinID_MutString -> [scores from each model]
    all_mut_scores = {}

    for mi, model_name in enumerate(args.model_location, 1):
        print(f"\n--- Processing model {mi}/{len(args.model_location)}: {model_name} ---")
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
            print("Transferred model to GPU")

        skipped_proteins = 0

        for protein_id, seq, mutations in tqdm(genes_data, desc=f"Model {mi}"):
            # Truncate sequence (keep your previous logic)
            if len(seq) > MAX_SEQ_LEN:
                seq_trunc = seq[:MAX_SEQ_LEN]
            else:
                seq_trunc = seq

            # Filter mutations: keep only single-point mutations that remain valid after truncation and match WT
            filtered_muts = []
            for mut in mutations:
                if validate_mutation(mut, seq_trunc, args.offset_idx) is not None:
                    filtered_muts.append(mut)

            if not filtered_muts:
                skipped_proteins += 1
                continue

            valid_muts, scores = score_gene_wt_marginals(
                model, alphabet, seq_trunc, filtered_muts, offset_idx=args.offset_idx
            )
            if not valid_muts:
                skipped_proteins += 1
                continue

            for mut_str, score in zip(valid_muts, scores):
                key = f"{protein_id}_{mut_str}"
                all_mut_scores.setdefault(key, []).append(score)

        print(f"Model {mi} finished. Skipped {skipped_proteins} proteins in this model.")

        del model
        torch.cuda.empty_cache()

    # Average score output: keep only mutations that have scores from all models
    results = []
    for key, scores in all_mut_scores.items():
        if len(scores) == len(args.model_location):
            results.append({"ProteinID_Mut": key, "Score": float(np.mean(scores))})
        else:
            print(f"Warning: {key} only has {len(scores)}/{len(args.model_location)} scores, skipped.")

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output_csv, index=False)
    print(f"\nSaved {len(results)} mutations to {args.output_csv}")


if __name__ == "__main__":
    main()