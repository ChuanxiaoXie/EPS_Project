#!/usr/bin/env python3
"""Render and optionally execute the manuscript five-track Circos configuration."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any


TRACK_FILES = {
    "genome_snp": "1_wai.txt",
    "gene_snp": "1_nei.txt",
    "gene_density": "2_genedensity.txt",
    "gc_count": "3_GC_wai.txt",
    "gc_to_at": "3_GC2AT_nei.txt",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-config", required=True, type=Path)
    parser.add_argument("--execute", action="store_true", help="Run Circos after rendering; rendering is the default.")
    return parser.parse_args()


def checked_path(value: object, label: str) -> str:
    text = str(value)
    if not text or any(character in text for character in "\r\n"):
        raise ValueError(f"{label} must be a non-empty single-line path")
    return text.rstrip("/\\")


def plot_range(config: dict[str, Any], name: str) -> tuple[float, float]:
    ranges = config.get("plot_ranges")
    if not isinstance(ranges, dict) or name not in ranges:
        raise ValueError(f"plot_ranges.{name} is required")
    values = ranges[name]
    if not isinstance(values, list) or len(values) != 2:
        raise ValueError(f"plot_ranges.{name} must contain [minimum, maximum]")
    minimum, maximum = float(values[0]), float(values[1])
    if maximum <= minimum:
        raise ValueError(f"plot_ranges.{name} maximum must exceed its minimum")
    return minimum, maximum


def render(config: dict[str, Any]) -> str:
    track_dir = checked_path(config.get("track_dir", ""), "track_dir")
    units = int(config.get("chromosomes_units", 5_000_000))
    if units <= 0:
        raise ValueError("chromosomes_units must be positive")
    ranges = {name: plot_range(config, name) for name in TRACK_FILES}

    def track(name: str) -> str:
        return f"{track_dir}/{TRACK_FILES[name]}"

    return f"""# Five-track chromosome plot used for the manuscript analysis.
karyotype = {track_dir}/karyotype.txt
chromosomes_units = {units}
chromosomes_display_default = yes

<colors>
genome_snp_blue = 27,136,175
gene_snp_blue = 76,188,204
</colors>

<image>
<<include etc/image.conf>>
</image>
<<include etc/ideogram.conf>>
<<include etc/colors_fonts_patterns.conf>>

<plots>
# Outer track: density of all selected SNPs.
<plot>
type = scatter
file = {track('genome_snp')}
r1 = 0.99r
r0 = 0.88r
min = {ranges['genome_snp'][0]:g}
max = {ranges['genome_snp'][1]:g}
color = genome_snp_blue
stroke_color = genome_snp_blue
glyph = circle
</plot>
<plot>
type = line
file = {track('genome_snp')}
r0 = 0.88r
r1 = 0.88r
color = black
thickness = 2
</plot>

# Inner part of track 1: density of selected SNPs overlapping genes.
<plot>
type = scatter
file = {track('gene_snp')}
r1 = 0.88r
r0 = 0.77r
orientation = in
min = {ranges['gene_snp'][0]:g}
max = {ranges['gene_snp'][1]:g}
color = gene_snp_blue
stroke_color = gene_snp_blue
glyph = circle
</plot>

# Track 2: gene-feature density.
<plot>
type = heatmap
file = {track('gene_density')}
r1 = 0.76r
r0 = 0.65r
min = {ranges['gene_density'][0]:g}
max = {ranges['gene_density'][1]:g}
color = greens-7-seq
</plot>
<plot>
type = line
file = {track('gene_density')}
r0 = 0.65r
r1 = 0.65r
color = black
thickness = 2
</plot>

# Outer part of track 3: reference G+C base count.
<plot>
type = histogram
file = {track('gc_count')}
r1 = 0.64r
r0 = 0.53r
min = {ranges['gc_count'][0]:g}
max = {ranges['gc_count'][1]:g}
color = rgb(164,222,209)
fill_color = rgb(164,222,209)
</plot>

# Inner part of track 3: G/C-to-A/T SNP density.
<plot>
type = histogram
file = {track('gc_to_at')}
r1 = 0.53r
r0 = 0.42r
orientation = in
min = {ranges['gc_to_at'][0]:g}
max = {ranges['gc_to_at'][1]:g}
color = rgb(105,180,219)
fill_color = rgb(105,180,219)
</plot>
<plot>
type = line
file = {track('gc_to_at')}
r0 = 0.53r
r1 = 0.53r
color = black
thickness = 2
</plot>
</plots>

<<include etc/ticks.conf>>
<<include etc/housekeeping.conf>>
data_out_of_range* = trim

<ideogram>
radius = 0.8r
label_size = 56
label_radius = 1.2r
</ideogram>
"""


def main() -> int:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError("Configuration root must be a JSON object")
    args.output_config.parent.mkdir(parents=True, exist_ok=True)
    args.output_config.write_text(render(config), encoding="utf-8")
    print(args.output_config)
    if args.execute:
        output_dir = Path(checked_path(config.get("output_dir", ""), "output_dir"))
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = str(config.get("output_file", "eps_chromosome_circos"))
        if not output_file or "/" in output_file or "\\" in output_file:
            raise ValueError("output_file must be a filename without directory separators")
        subprocess.run(
            [
                str(config.get("circos", "circos")),
                "-conf",
                str(args.output_config),
                "-png",
                "-svg",
                "-outputdir",
                str(output_dir),
                "-outputfile",
                output_file,
            ],
            check=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
