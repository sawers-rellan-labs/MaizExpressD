#!/usr/bin/env python3
"""Describe the structure of single-end 3' tag reads before choosing trimming settings.

For each FASTQ(.gz) file: read lengths, base composition at the 5' and 3' ends,
the most common 5' prefixes (a UMI + spacer design shows a fixed spacer, e.g. TATA),
poly(A) and poly(T) runs (QuantSeq FWD reads end in poly(A); REV reads start with
poly(T)), Illumina adapter content and exact duplicates.

Usage: python3 -I check_read_structure.py <out_prefix> <fastq.gz> [<fastq.gz> ...]
Writes <out_prefix>.summary.tsv and <out_prefix>.composition.tsv.
"""
import csv
import gzip
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

ADAPTERS = {"truseq": "AGATCGGAAGAGC", "nextera": "CTGTCTCTTATACACATCT"}
POLYA = re.compile(r"A{10,}")
POLYT_5P = re.compile(r"^.{0,12}T{10,}")
N_5P, N_3P = 20, 15


def reads(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            head = f.readline()
            if not head:
                return
            seq = f.readline().strip()
            f.readline()
            f.readline()
            yield seq


def describe(path):
    n = 0
    lengths = []
    comp5 = [Counter() for _ in range(N_5P)]
    comp3 = [Counter() for _ in range(N_3P)]
    prefix6, spacer = Counter(), Counter()
    polya_any = polya_end = polyt_start = 0
    polya_start = []
    adapter = Counter()
    adapter_pos = []
    seen = Counter()
    for s in reads(path):
        n += 1
        lengths.append(len(s))
        for i, b in enumerate(s[:N_5P]):
            comp5[i][b] += 1
        for i, b in enumerate(s[-N_3P:][::-1]):     # position 1 = last base
            comp3[i][b] += 1
        prefix6[s[:6]] += 1
        spacer[s[6:10]] += 1
        m = POLYA.search(s)
        if m:
            polya_any += 1
            polya_start.append(m.start())
            if m.end() >= len(s) - 2:
                polya_end += 1
        if POLYT_5P.match(s):
            polyt_start += 1
        for name, a in ADAPTERS.items():
            p = s.find(a)
            if p >= 0:
                adapter[name] += 1
                if name == "truseq":
                    adapter_pos.append(p)
        seen[s] += 1
    dup = sum(c - 1 for c in seen.values()) / n
    top_prefix, top_n = prefix6.most_common(1)[0]
    top_spacer, top_sp = spacer.most_common(1)[0]
    summary = {
        "file": Path(path).name.split(".")[0], "reads": n,
        "len_min": min(lengths), "len_median": statistics.median(lengths), "len_max": max(lengths),
        "pct_polyA10_any": round(100 * polya_any / n, 1),
        "pct_polyA10_at_end": round(100 * polya_end / n, 1),
        "polyA_start_median": statistics.median(polya_start) if polya_start else "",
        "pct_polyT10_at_start": round(100 * polyt_start / n, 2),
        "pct_truseq": round(100 * adapter["truseq"] / n, 2),
        "truseq_pos_median": statistics.median(adapter_pos) if adapter_pos else "",
        "pct_nextera": round(100 * adapter["nextera"] / n, 2),
        "top_prefix6": top_prefix, "pct_top_prefix6": round(100 * top_n / n, 2),
        "top_pos7_10": top_spacer, "pct_top_pos7_10": round(100 * top_sp / n, 2),
        "pct_exact_duplicates": round(100 * dup, 1),
    }
    comp = []
    for end, table in (("5p", comp5), ("3p", comp3)):
        for i, c in enumerate(table, 1):
            tot = sum(c.values())
            comp.append({"file": summary["file"], "end": end, "position": i,
                         **{b: round(100 * c[b] / tot, 1) for b in "ACGTN"}})
    return summary, comp


def main():
    out, files = sys.argv[1], sys.argv[2:]
    summaries, comps = [], []
    for f in files:
        s, c = describe(f)
        summaries.append(s)
        comps.extend(c)
        print("\t".join(f"{k}={v}" for k, v in s.items()))
    for name, rows in ((".summary.tsv", summaries), (".composition.tsv", comps)):
        with open(out + name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
