import os
import sys
import glob
import re
import pandas as pd
from tqdm import tqdm

sys.path.append("/data/home/scxj090/run/SaProt-main")

from model.saprot.saprot_foldseek_mutation_model import SaprotFoldseekMutationModel
from utils.foldseek_util import get_struc_seq

# ========== Config ==========
ESMFOLD_PDB_ROOT = "/data/home/scxj090/run/esm/results"
EPS_MUT_XLSX = "/data/home/scxj090/run/SaProt-main/Sh2_MUT.xlsx"
OUTPUT_CSV = "/data/home/scxj090/run/SaProt-main/sh2_saprot_mut_scores.csv"

FOLDSEEK_BIN = "/data/home/scxj090/run/SaProt-main/bin/foldseek"
CHAIN_ID = "A"
PLDDT_MASK_FORCE_TRUE = True   # If PDB has no pLDDT, force False? Change to False?

MODEL_KWARGS = {
    "foldseek_path": None,
    "config_path": "/data/home/scxj090/run/SaProt-main/weights/PLMs/SaProt_1.3B_AFDB_OMG_NCBI",
    "load_pretrained": True,
}

DEVICE = "cuda"
# ======================================

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def has_plddt(pdb_file: str) -> bool:
    if pdb_file.endswith(".cif"):
        return True
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    b = line[60:66].strip()
                    try:
                        if float(b) != 0.0:
                            return True
                    except ValueError:
                        continue
        return False
    except Exception:
        return False

def find_pdb_for_gene(gene: str, pdb_root: str):
    gene = str(gene).strip()
    cand = (glob.glob(os.path.join(pdb_root, gene, "*.pdb")) +
            glob.glob(os.path.join(pdb_root, gene, "*.cif")))
    if cand:
        return cand[0]
    cand = (glob.glob(os.path.join(pdb_root, f"{gene}*.pdb")) +
            glob.glob(os.path.join(pdb_root, f"{gene}*.cif")))
    if cand:
        return cand[0]
    # Additional attempt: search directly within subdirectories
    cand = glob.glob(os.path.join(pdb_root, "*", f"{gene}*.pdb")) + \
           glob.glob(os.path.join(pdb_root, "*", f"{gene}*.cif"))
    if cand:
        return cand[0]
    return None

def clean_mutation_cell(x):
    if x is None:
        return None
    if isinstance(x, float) and pd.isna(x):
        return None
    s = str(x).strip()
    if s == "" or s.lower() in ["nan", "none"]:
        return None
    return s

def is_single_mut_format(m: str) -> bool:
    return re.match(r"^[A-Za-z]\d+[A-Za-z]$", m) is not None

def read_mutations_from_row(row, start_col_idx=1, end_col_idx=7):
    """Read columns B~H (indices 1 to 7)"""
    muts = []
    for idx in range(start_col_idx, min(end_col_idx + 1, len(row))):
        m = clean_mutation_cell(row[idx])
        if m is None:
            continue
        if m.startswith("#"):
            continue
        muts.append(m)
    return muts

def load_and_prepare_model():
    model = SaprotFoldseekMutationModel(**MODEL_KWARGS)
    model.eval()
    model.to(DEVICE)
    print("Model loaded.")
    return model

def run():
    model = load_and_prepare_model()
    df = pd.read_excel(EPS_MUT_XLSX)
    if df.shape[1] < 2:
        raise ValueError(f"EPS_MUT.xlsx has insufficient columns: {df.shape[1]}, at least 2 columns required")

    results = []
    failed_log = []  # Record failure information

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Predict genes"):
        gene = str(row.iloc[0]).strip()
        if gene.lower() in ["nan", "none", ""]:
            continue

        pdb_file = find_pdb_for_gene(gene, ESMFOLD_PDB_ROOT)
        if pdb_file is None:
            print(f"[Skipped] Gene {gene}: corresponding PDB not found")
            continue

        plddt_available = has_plddt(pdb_file)
        plddt_mask = (plddt_available and not PLDDT_MASK_FORCE_TRUE)

        try:
            parsed = get_struc_seq(
                FOLDSEEK_BIN,
                pdb_file,
                [CHAIN_ID],
                plddt_mask=plddt_mask,
            )
        except Exception as e:
            print(f"[Skipped] get_struc_seq failed gene={gene}, pdb={pdb_file}, err={e}")
            failed_log.append({"gene": gene, "error": str(e), "stage": "get_struc_seq"})
            continue

        if CHAIN_ID not in parsed:
            print(f"[Skipped] Chain {CHAIN_ID} not found in gene {gene}")
            continue

        seq, foldseek_seq, combined_seq = parsed[CHAIN_ID]

        # Read mutations: now columns B~H (indices 1~7)
        mutations = read_mutations_from_row(row, start_col_idx=1, end_col_idx=7)
        mutations = list(dict.fromkeys(mutations))  # Deduplicate
        if not mutations:
            continue

        for mut in mutations:
            # Check whether the mutation position exceeds the sequence length
            pos_match = re.search(r'\d+', mut)
            if pos_match:
                pos = int(pos_match.group())
                if pos > len(combined_seq):
                    msg = f"Mutation position {pos} exceeds sequence length {len(combined_seq)}"
                    print(f"[Failed] gene={gene} mut={mut} {msg}")
                    failed_log.append({"gene": gene, "mutation": mut, "error": msg, "stage": "predict"})
                    continue

            try:
                score = model.predict_mut(combined_seq, mut)
                results.append({
                    "gene": gene,
                    "chain": CHAIN_ID,
                    "mutation": mut,
                    "score": score,
                    "pdb_file": pdb_file,
                })
            except Exception as e:
                print(f"[Failed] gene={gene} mut={mut} err={e}")
                failed_log.append({"gene": gene, "mutation": mut, "error": str(e), "stage": "predict"})

    out_dir = os.path.dirname(OUTPUT_CSV)
    if out_dir:
        ensure_dir(out_dir)

    out_df = pd.DataFrame(results)
    out_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")
    print(f"Results saved to {OUTPUT_CSV}, total {len(results)} records")

    if failed_log:
        fail_df = pd.DataFrame(failed_log)
        fail_csv = OUTPUT_CSV.replace(".csv", "_failed.csv")
        fail_df.to_csv(fail_csv, index=False)
        print(f"Failed records saved to {fail_csv}")

if __name__ == "__main__":
    run()
