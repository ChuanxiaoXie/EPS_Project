#!/usr/bin/env python3
"""Count distinct candidate-read QNAMEs supporting each contig breakpoint."""
import argparse,csv,re,subprocess
CIGAR=re.compile(r"(\d+)([MIDNSHP=X])")
def ref_span(pos,cigar):
    n=sum(int(x) for x,op in CIGAR.findall(cigar) if op in "MDN=X"); return pos,pos+max(0,n-1)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--samtools",required=True); ap.add_argument("--bam",required=True)
    ap.add_argument("--candidates",required=True); ap.add_argument("--output",required=True); args=ap.parse_args()
    with open(args.candidates,encoding="utf-8") as fh:
        reader=csv.DictReader(fh,delimiter="\t"); rows=list(reader); fields=reader.fieldnames or []
    extra=["mapped_qnames_100bp","spanning_qnames","endpoint_near_qnames","pair_qnames","raw_support_class","support_class","interpretation"]
    with open(args.output,"w",newline="",encoding="utf-8") as out:
        w=csv.DictWriter(out,fieldnames=fields+extra,delimiter="\t"); w.writeheader()
        for row in rows:
            bp=int(row["contig_breakpoint_1based"]); region=f'{row["contig"]}:{max(1,bp-100)}-{bp+100}'
            p=subprocess.run([args.samtools,"view",args.bam,region],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
            mapped=set(); spanning=set(); endpoint=set(); paired=set()
            for line in p.stdout.splitlines():
                f=line.split("\t")
                if len(f)<11 or f[5]=="*": continue
                qn,flag,pos,cigar=f[0],int(f[1]),int(f[3]),f[5]; start,end=ref_span(pos,cigar); mapped.add(qn)
                if start<=bp-5 and end>=bp+5: spanning.add(qn)
                if abs(start-bp)<=20 or abs(end-bp)<=20: endpoint.add(qn)
                if flag&1 and not flag&8: paired.add(qn)
            raw_support="HIGH" if len(spanning)>=3 else ("SUPPORTED" if len(spanning)>=1 or len(endpoint)>=2 else "WEAK")
            risks=set(filter(None,row.get("risk_labels","NONE").split(";"))); risks.discard("NONE")
            if len(mapped)>=10000: risks.add("EXTREME_LOCAL_DEPTH")
            row["risk_labels"]=";".join(sorted(risks)) or "NONE"
            non_specific=bool(risks & {"TDNA_GENOME_HOMOLOGY","LARGE_SEGMENT_OVERLAP","EXTREME_LOCAL_DEPTH"})
            interpretation=("Retained risk candidate: raw read support is non-specific in a plant-homologous, "
                "overlapping, or extreme-depth sequence context; not a confirmed insertion breakpoint."
                if non_specific else "Read support is compatible with the assembled breakpoint candidate.")
            row.update({"mapped_qnames_100bp":len(mapped),"spanning_qnames":len(spanning),
                "endpoint_near_qnames":len(endpoint),"pair_qnames":len(paired),
                "raw_support_class":raw_support,
                "support_class":"RISK_FLAGGED" if non_specific else raw_support,
                "interpretation":interpretation})
            w.writerow(row)
if __name__=="__main__": main()
