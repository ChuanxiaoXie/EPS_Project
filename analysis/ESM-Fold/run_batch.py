import os
import re
import torch
import esm
import pandas as pd
import biotite.structure.io as bsio
from tqdm import tqdm   # optional, for progress bar; comment out related lines if not installed

def sanitize_filename(name: str) -> str:
    """Replace illegal characters in the sample name with underscores; keep only letters, digits, and underscores"""
    return re.sub(r'[^a-zA-Z0-9_]', '_', name)

def main():
    # create the results folder
    os.makedirs("results", exist_ok=True)

    # read mutants.xlsx, assuming the first column is sample and the second is sequence
    df = pd.read_excel("CAD_mutation.xlsx", engine="openpyxl")
    # make sure column names are correct; modify if your actual column names differ
    df.columns = ["sample", "sequence"]   # force rename, or adjust according to your actual case

    # load the ESMFold model (loaded once globally)
    print("Loading ESMFold model...")
    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    model.set_chunk_size(128)   # adjust according to GPU memory
    print("Model loaded.")

    results = []   # store (sample, mean_plddt)

    # iterate over each row
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Predicting structures"):
        sample_raw = row["sample"]
        sequence = str(row["sequence"]).strip()   # ensure it is a string and strip leading/trailing whitespace

        # generate a safe file name
        sample_clean = sanitize_filename(sample_raw)
        pdb_path = os.path.join("results", f"{sample_clean}.pdb")

        # if the file already exists, you may skip it (we do not skip here, re-predict and overwrite)
        # to skip existing files, uncomment the next line
        # if os.path.exists(pdb_path): continue

        # predict
        with torch.no_grad():
            try:
                pdb_string = model.infer_pdb(sequence)
            except Exception as e:
                print(f"Error predicting {sample_raw}: {e}")
                continue

        # save the PDB file
        with open(pdb_path, "w") as f:
            f.write(pdb_string)

        # read the B-factor from the PDB file and compute the mean pLDDT
        try:
            struct = bsio.load_structure(pdb_path, extra_fields=["b_factor"])
            mean_plddt = struct.b_factor.mean()
        except Exception as e:
            print(f"Error reading {pdb_path}: {e}")
            mean_plddt = float('nan')

        results.append({"sample": sample_raw, "mean_plddt": mean_plddt})

        # optional: clear GPU cache every 100 predicted sequences
        if (idx + 1) % 100 == 0:
            torch.cuda.empty_cache()

    # save the summary CSV
    summary_df = pd.DataFrame(results)
    summary_df.to_csv("plddt_summary.csv", index=False)
    print("Done. Summary saved to plddt_summary.csv")

if __name__ == "__main__":
    main()
