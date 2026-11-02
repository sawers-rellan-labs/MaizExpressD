#!/usr/bin/env python3
"""Check metadata/contrasts.tsv against metadata/samples.tsv (development).

Checks: columns; class and design vocabulary; sample ids exist, are included and are Kremling2018;
numerator and denominator are disjoint and n counts match; each side is one tissue, the contrast's;
panel contrasts hold the same pairs (line, or line and sampling group for day/night) on both sides,
leave out B73 and B73-10, and pair adult leaf within sampling group; B73 contrasts hold only B73;
no pilot L3Base, B73-10 or LMid library is used; *_qc counts match metadata/library_exclusions.tsv.
Also lists included libraries that no contrast uses, by tissue, for review.
Usage: python3 -I metadata/build/validate_contrasts.py
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_contrasts import COLS, NOT_PANEL, PILOT_FC, ROOT, STUDY, read_tsv  # noqa: E402

CLASSES = {"leaf_maturation", "day_night", "leaf_stage", "organ"}
DESIGNS = {"panel", "B73"}


def main():
    samples = {s["sample_id"]: s for s in read_tsv(ROOT / "metadata" / "samples.tsv")}
    excluded = {r["sample_id"] for r in read_tsv(ROOT / "metadata" / "library_exclusions.tsv")}
    rows = list(csv.reader(open(ROOT / "metadata" / "contrasts.tsv"), delimiter="\t"))
    problems = []
    if rows[0] != COLS:
        problems.append(f"header differs: {rows[0]}")
    contrasts = [dict(zip(rows[0], r)) for r in rows[1:]]
    if len({c["contrast_id"] for c in contrasts}) != len(contrasts):
        problems.append("duplicate contrast ids")
    used = set()
    for c in contrasts:
        cid = c["contrast_id"]
        if c["contrast_class"] not in CLASSES:
            problems.append(f"{cid}: class {c['contrast_class']}")
        if c["design"] not in DESIGNS or c["design"] not in c["tags"].split(";"):
            problems.append(f"{cid}: design {c['design']} (tags {c['tags']})")
        ns, ds = c["numerator_samples"].split(";"), c["denominator_samples"].split(";")
        if set(ns) & set(ds):
            problems.append(f"{cid}: numerator and denominator overlap")
        if len(ns) != int(c["n_numerator"]) or len(ds) != int(c["n_denominator"]):
            problems.append(f"{cid}: n counts do not match sample lists")
        if (c["replicated"] == "TRUE") != (len(ns) >= 2 and len(ds) >= 2):
            problems.append(f"{cid}: replicated flag does not match n")
        sides = {}
        for name, ids, tissue in (("num", ns, c["numerator_tissue"]), ("den", ds, c["denominator_tissue"])):
            ok = []
            for i in ids:
                s = samples.get(i)
                if s is None:
                    problems.append(f"{cid}: unknown sample {i}")
                    continue
                if s["include"] != "TRUE" or s["study_id"] != STUDY:
                    problems.append(f"{cid}: sample {i} not an included {STUDY} library")
                if s["tissue"] != tissue:
                    problems.append(f"{cid}: {i} is {s['tissue']}, not {tissue}")
                if (s["tissue"] == "L3Base" and s["flowcell"] == PILOT_FC) or s["genotype"] == "B73-10" or s["tissue"] == "LMid":
                    problems.append(f"{cid}: {i} is left out by design (pilot L3Base, B73-10 or LMid)")
                if c["design"] == "panel" and s["genotype"] in NOT_PANEL:
                    problems.append(f"{cid}: {i} is B73 in a panel contrast")
                if c["design"] == "B73" and s["genotype"] != "B73":
                    problems.append(f"{cid}: {i} is not B73")
                if c["sampling_group"] and s["sampling_group"] not in c["sampling_group"].split(";"):
                    problems.append(f"{cid}: {i} is from sampling group {s['sampling_group']!r}")
                used.add(i)
                ok.append(s)
            sides[name] = ok
        if c["design"] == "panel":
            key = (lambda s: f"{s['genotype']}|{s['sampling_group']}") if c["pair_key"] == "genotype|sampling_group" else (lambda s: s["genotype"])
            if c["contrast_class"] == "day_night" and c["pair_key"] != "genotype|sampling_group":
                problems.append(f"{cid}: day/night must pair within line and sampling group")
            kn, kd = {key(s) for s in sides["num"]}, {key(s) for s in sides["den"]}
            if kn != kd or len(kn) != int(c["n_pairs"]):
                problems.append(f"{cid}: sides hold different pairs ({len(kn)} vs {len(kd)}; n_pairs {c['n_pairs']})")
            nq = [s for s in sides["num"] if s["sample_id"] not in excluded]
            dq = [s for s in sides["den"] if s["sample_id"] not in excluded]
            if len({key(s) for s in nq} & {key(s) for s in dq}) != int(c["n_pairs_qc"]):
                problems.append(f"{cid}: n_pairs_qc does not match library_exclusions.tsv")
        for name in ("numerator", "denominator"):
            n_qc = sum(1 for s in sides[name[:3]] if s["sample_id"] not in excluded)
            if n_qc != int(c[f"n_{name}_qc"]):
                problems.append(f"{cid}: n_{name}_qc does not match library_exclusions.tsv")

    unused = Counter(s["tissue"] for i, s in samples.items() if s["include"] == "TRUE" and i not in used)
    print(f"{len(contrasts)} contrasts; {len(problems)} problems")
    for p in problems:
        print("  -", p)
    print("included libraries used by no contrast, by tissue:", ", ".join(f"{t} {n}" for t, n in sorted(unused.items())))
    return problems


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
