#!/usr/bin/env python3
import sys, os
from Bio import SeqIO

if len(sys.argv)!=2:
    sys.exit("Usage: extract_cluster.py <cluster_id>")

cid = sys.argv[1]
# load (or re‐use) your on‐disk FASTA index
idx = SeqIO.index_db("all_viral.idx","IMGVR_all_protein.faa","fasta")

outdir = "clusters/faa"
os.makedirs(outdir, exist_ok=True)
out_fa = f"{outdir}/Cluster_{cid}.faa"

with open("cluster_seq_map.csv") as map_fh, open(out_fa,"w") as of:
    for line in map_fh:
        cl, sid = line.strip().split(",",1)
        if cl!=cid: 
            continue
        of.write(idx[sid].format("fasta"))
print(f"[{cid}] wrote {out_fa}", file=sys.stderr)
