"""Median salmon mapping rate per quantified batch (flags any under the 20% discard rule).
Run on Hazel: python3 mapping_rule.py /rsstu/users/r/rrellan/sara/MaizeExpressD/results/Kremling2018"""
import glob, json, os, statistics, sys

R = sys.argv[1]
for d in sorted(x for x in os.listdir(R) if os.path.isdir(os.path.join(R, x)) and not x.startswith("_")):
    p = [json.load(open(f))["percent_mapped"]
         for f in glob.glob(os.path.join(R, d, "salmon/*/aux_info/meta_info.json"))]
    if p:
        med = statistics.median(p)
        flag = "DISCARD (<20%)" if med < 20 else ""
        print(f"{d:20s} n={len(p):3d}  median {med:5.1f}%  min {min(p):5.1f}%  {flag}")
