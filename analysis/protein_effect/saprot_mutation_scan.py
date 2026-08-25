#!/usr/bin/env python3
"""Score all amino-acid substitutions at one structure position with SaProt."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--position", required=True, type=int)
    parser.add_argument("--protein-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--saprot-root", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--foldseek", required=True)
    parser.add_argument("--chain", default="A")
    parser.add_argument("--device", default="cpu")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    saprot_root = Path(args.saprot_root).resolve()
    if not saprot_root.is_dir():
        raise FileNotFoundError(saprot_root)
    sys.path.insert(0, str(saprot_root))
    from model.saprot.saprot_foldseek_mutation_model import SaprotFoldseekMutationModel
    from utils.foldseek_util import get_struc_seq

    structure_sequences = get_struc_seq(args.foldseek, args.structure, [args.chain], plddt_mask=False)
    if args.chain not in structure_sequences:
        raise ValueError(f"Chain not returned by Foldseek: {args.chain}")
    sequence, _foldseek_sequence, combined_sequence = structure_sequences[args.chain]
    if not 1 <= args.position <= len(sequence):
        raise ValueError(f"Position {args.position} is outside the sequence length {len(sequence)}")
    model = SaprotFoldseekMutationModel(
        foldseek_path=args.foldseek,
        config_path=args.model_path,
        load_pretrained=True,
    )
    model.eval()
    model.to(args.device)
    scores = model.predict_pos_mut(combined_sequence, args.position)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "position", "reference_amino_acid", "model_mutation_key", "score"])
        for mutation_key, score in sorted(scores.items(), key=lambda item: str(item[0])):
            writer.writerow([args.protein_id, args.position, sequence[args.position - 1], mutation_key, score])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
