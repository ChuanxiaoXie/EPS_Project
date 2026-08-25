#!/usr/bin/env python3
"""Call adjacent plant/T-DNA segments on the same assembled contig from PAF."""
import argparse
import csv
from collections import defaultdict


def paf(path, min_aln, min_ident):
    out = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            aln = int(f[10]); matches = int(f[9])
            ident = matches / aln if aln else 0.0
            if aln < min_aln or ident < min_ident:
                continue
            out[f[0]].append({"qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]),
                "strand": f[4], "target": f[5], "ts": int(f[7]), "te": int(f[8]),
                "matches": matches, "aln": aln, "ident": ident, "mapq": int(f[11])})
    return out


def relation(a, b):
    if a["qe"] <= b["qs"]:
        return b["qs"] - a["qe"], "GENOME_THEN_TDNA", (a["qe"] + b["qs"]) // 2
    if b["qe"] <= a["qs"]:
        return a["qs"] - b["qe"], "TDNA_THEN_GENOME", (b["qe"] + a["qs"]) // 2
    ov = min(a["qe"], b["qe"]) - max(a["qs"], b["qs"])
    return -ov, "OVERLAP", (max(a["qs"], b["qs"]) + min(a["qe"], b["qe"])) // 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-paf", required=True); ap.add_argument("--tdna-paf", required=True)
    ap.add_argument("--homology-paf", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--sample", default="M1AA")
    ap.add_argument("--max-gap", type=int, default=300); ap.add_argument("--max-overlap", type=int, default=150)
    args = ap.parse_args()
    genome = paf(args.genome_paf, 80, 0.85); tdna = paf(args.tdna_paf, 50, 0.80)
    # Homology is a sensitive risk annotation, not a candidate filter. Use a
    # permissive threshold so divergent plant-homologous vector segments are
    # labelled rather than silently treated as junctions.
    homology = paf(args.homology_paf, 40, 0.60)
    hom_intervals = [(h["qs"], h["qe"], h["target"], h["ident"]) for hits in homology.values() for h in hits]
    rows=[]; seen=set()
    for contig in sorted(set(genome) & set(tdna)):
        for g in genome[contig]:
            for t in tdna[contig]:
                gap, order, bp = relation(g, t)
                if gap > args.max_gap or gap < -args.max_overlap: continue
                key=(contig, round(bp/20), g["target"], t["target"])
                if key in seen: continue
                seen.add(key); risks=[]; hom_detail=[]
                for hs,he,ht,hi in hom_intervals:
                    ov=min(t["te"],he)-max(t["ts"],hs)
                    if ov >= 20:
                        risks.append("TDNA_GENOME_HOMOLOGY"); hom_detail.append(f"{hs}-{he}:{ht}:{hi:.3f}:ov{ov}")
                if g["mapq"] < 20: risks.append("LOW_GENOME_MAPQ")
                if gap < -20: risks.append("LARGE_SEGMENT_OVERLAP")
                score=min(g["aln"],500)+min(t["aln"],500)+g["mapq"]+t["mapq"]
                rows.append({"candidate_id":"", "contig":contig, "contig_length":g["qlen"],
                    "contig_breakpoint_1based":bp+1, "segment_order":order, "query_gap_bp":gap,
                    "genome_target":g["target"], "genome_start_1based":g["ts"]+1, "genome_end_1based":g["te"],
                    "genome_strand":g["strand"], "genome_aln_bp":g["aln"], "genome_identity":f'{g["ident"]:.5f}',
                    "genome_mapq":g["mapq"], "tdna_target":t["target"], "tdna_start_1based":t["ts"]+1,
                    "tdna_end_1based":t["te"], "tdna_strand":t["strand"], "tdna_aln_bp":t["aln"],
                    "tdna_identity":f'{t["ident"]:.5f}', "tdna_mapq":t["mapq"],
                    "risk_labels":";".join(sorted(set(risks))) or "NONE",
                    "homology_detail":";".join(hom_detail) or ".", "candidate_score":score})
    rows.sort(key=lambda r:(-r["candidate_score"],r["contig"],r["contig_breakpoint_1based"]))
    for i,row in enumerate(rows,1): row["candidate_id"]=f"{args.sample}_J{i:04d}"
    fields=list(rows[0]) if rows else ["candidate_id","contig","contig_length","contig_breakpoint_1based",
        "segment_order","query_gap_bp","genome_target","genome_start_1based","genome_end_1based",
        "genome_strand","genome_aln_bp","genome_identity","genome_mapq","tdna_target","tdna_start_1based",
        "tdna_end_1based","tdna_strand","tdna_aln_bp","tdna_identity","tdna_mapq","risk_labels",
        "homology_detail","candidate_score"]
    with open(args.output,"w",newline="",encoding="utf-8") as out:
        w=csv.DictWriter(out,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__": main()
