#!/usr/bin/env python3
"""Plot a whole-genome Hi-C contact matrix from a Cooler multi-resolution file."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcool", required=True, help="Input .mcool file")
    parser.add_argument("--resolution", type=int, default=25000)
    parser.add_argument("--chromosomes", required=True, help="Comma-separated Cooler chromosome names")
    parser.add_argument("--output", required=True)
    parser.add_argument("--balance", default="", help="Cooler balance column; leave empty for raw counts")
    parser.add_argument("--cmap", default="Reds")
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    import cooler
    import cooltools
    import matplotlib.pyplot as plt
    import numpy as np

    chromosomes = [value.strip() for value in args.chromosomes.split(",") if value.strip()]
    if not chromosomes:
        raise ValueError("At least one chromosome is required")
    cooler_uri = f"{args.mcool}::/resolutions/{args.resolution}"
    contact_map = cooler.Cooler(cooler_uri)
    missing = [chromosome for chromosome in chromosomes if chromosome not in contact_map.chromnames]
    if missing:
        raise ValueError(f"Chromosomes absent from Cooler file: {', '.join(missing)}")
    balance: bool | str = args.balance if args.balance else False
    matrix = cooltools.matrix_whole_genome(contact_map, chrom_order=chromosomes, balance=balance)
    cumulative_lengths = np.cumsum([0, *[int(contact_map.chromsizes[chromosome]) for chromosome in chromosomes]])
    tick_positions = (cumulative_lengths[:-1] + cumulative_lengths[1:]) / 2

    figure, axis = plt.subplots(figsize=(14, 11), dpi=args.dpi)
    image = axis.imshow(matrix, cmap=args.cmap, vmin=args.vmin, vmax=args.vmax, aspect="auto")
    labels = [value if value.lower().startswith("chr") else f"Chr{value}" for value in chromosomes]
    axis.set_xticks(tick_positions, labels, rotation=90, fontsize=8)
    axis.set_yticks(tick_positions, labels, fontsize=8)
    axis.grid(linestyle="-", linewidth=0.2, color="#dddddd")
    color_bar = figure.colorbar(image, ax=axis, shrink=0.78)
    color_bar.set_label(f"{args.balance} normalized counts" if args.balance else "Contact counts")
    figure.tight_layout()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
