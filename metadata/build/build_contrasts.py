#!/usr/bin/env python3
"""Build metadata/contrasts.tsv for Kremling2018 (M4): developmental and organ contrasts between tissues
of the same line.

Contrasts (HANDOFF.md section 3; log2FC = log2(numerator / denominator)):
  leaf_maturation  L3Tip vs L3Base; L3Mid vs L3Base and L3Tip vs L3Mid (the 26 lines with all three
                   positions, all on one flowcell);
  day_night        LMAN vs LMAD, paired within line and sampling group (both groups together, and each
                   group alone as secondary contrasts);
  leaf_stage       LMAD vs L3Tip (adult leaf by day vs seedling leaf-3 tip);
  organ            GShoot vs GRoot.
Each has two versions:
  panel  the diversity panel, paired by line (DESeq2 ~ pair + tissue, as MaizExpressC's Zhou2022 panel):
         every line except B73 and B73-10 with libraries of both tissues; repeat libraries of a line enter
         as extra libraries of its pair; adult leaf pairs within line and sampling group;
  B73    the B73 libraries alone (replicated, ~ 0 + tissue), the reference line (M1 decision).
Left out by design (M1/M3 decisions): B73-10 (QC control only), the pilot L3Base libraries (flowcell
HHV55BGXX), LMid (display only); Kern has no developmental partner (it enters tissue specificity only).
Sample lists hold every included library of the design; library QC (metadata/library_exclusions.tsv) is
applied when DESeq2 runs, as in MaizExpressC, and the *_qc columns say what is left after it (panel: lines
with a library passing QC on both sides).
Usage: python3 -I metadata/build/build_contrasts.py [--write]
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STUDY = "Kremling2018"
PILOT_FC = "HHV55BGXX"
NOT_PANEL = {"B73", "B73-10"}
COLS = ["contrast_id", "study_id", "contrast_class", "design", "numerator_tissue", "denominator_tissue",
        "numerator_group", "denominator_group", "held_constant", "pair_key", "genotype", "sampling_group",
        "replicated", "n_numerator", "n_denominator", "n_pairs", "n_numerator_qc", "n_denominator_qc", "n_pairs_qc",
        "sign_note", "notes", "tags", "numerator_samples", "denominator_samples"]
LABEL = {"GRoot": "germinating root", "GShoot": "germinating shoot", "L3Base": "leaf 3 base", "L3Mid": "leaf 3 middle",
         "L3Tip": "leaf 3 tip", "LMAD": "adult leaf, day", "LMAN": "adult leaf, night"}
# class, numerator, denominator, sampling group ("" = any; "each" = per group), extra tags, sign note, notes
SPEC = [
    ("leaf_maturation", "L3Tip", "L3Base", "", "", "positive = higher at the tip of leaf 3 (mature, photosynthetic) than at its base (dividing, expanding cells)",
     "Seedling leaf 3, 2 cm from the tip vs 2 cm from the base."),
    ("leaf_maturation", "L3Mid", "L3Base", "", "subset_26_lines", "positive = higher in the middle of leaf 3 than at its base",
     "The 26 lines with base, middle and tip of leaf 3, all sequenced on one flowcell (HJ7FKBGXX)."),
    ("leaf_maturation", "L3Tip", "L3Mid", "", "subset_26_lines", "positive = higher at the tip of leaf 3 than in its middle",
     "The 26 lines with base, middle and tip of leaf 3, all sequenced on one flowcell (HJ7FKBGXX)."),
    ("day_night", "LMAN", "LMAD", "", "", "positive = higher at mid-night than at mid-day",
     "Mature mid-leaf after flowering, mid-night vs mid-day; paired within line and sampling group (8 or 26: dates matched "
     "to flowering time; germplasm differs between groups)."),
    ("day_night", "LMAN", "LMAD", "each", "secondary", "positive = higher at mid-night than at mid-day",
     "One sampling group only (sensitivity check of the combined contrast)."),
    ("leaf_stage", "LMAD", "L3Tip", "", "conditions_differ", "positive = higher in adult leaf (by day) than at the tip of seedling leaf 3",
     "Mature mid-leaf after flowering (by day) vs the tip of seedling leaf 3: stage and sampling date differ, and growth "
     "conditions may too (still NR), so read the result as adult vs seedling leaf in this experiment, not age alone."),
    ("organ", "GShoot", "GRoot", "", "", "positive = higher in the germinating shoot than in the germinating root",
     "Root and shoot of the same germinating seedlings."),
]


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build():
    samples = [s for s in read_tsv(ROOT / "metadata" / "samples.tsv") if s["include"] == "TRUE"]
    excluded = {r["sample_id"] for r in read_tsv(ROOT / "metadata" / "library_exclusions.tsv")}
    usable = [s for s in samples if not (s["tissue"] == "L3Base" and s["flowcell"] == PILOT_FC) and s["genotype"] != "B73-10"]
    rows = []

    def add(cls, a, b, design, group, tags, sign, notes, num, den, pair):
        lines = sorted({pair(s) for s in num} & {pair(s) for s in den}) if design == "panel" else []
        if design == "panel":
            num = [s for s in num if pair(s) in lines]
            den = [s for s in den if pair(s) in lines]
        nq = [s for s in num if s["sample_id"] not in excluded]
        dq = [s for s in den if s["sample_id"] not in excluded]
        lines_qc = sorted({pair(s) for s in nq} & {pair(s) for s in dq}) if design == "panel" else []
        who = "panel" if design == "panel" else "B73"
        cid = f"{STUDY}_{cls}_{a}_vs_{b}_{who}" + (f"_group{group}" if group else "")
        held = "line" + (" and sampling group" if cls == "day_night" else "") if design == "panel" else "genotype=B73"
        if group:
            held += f"; sampling_group={group}"
        same_fc = sum(1 for k in lines_qc if {s["flowcell"] for s in nq if pair(s) == k} & {s["flowcell"] for s in dq if pair(s) == k})
        t = [design] + [x for x in tags.split(";") if x]
        if design == "panel" and lines_qc and same_fc < 0.1 * len(lines_qc):
            t.append("flowcell_confounded")
        if design == "B73" and not {s["flowcell"] for s in nq} & {s["flowcell"] for s in dq}:
            t.append("flowcell_confounded")
        rows.append({
            "contrast_id": cid, "study_id": STUDY, "contrast_class": cls, "design": design,
            "numerator_tissue": a, "denominator_tissue": b,
            "numerator_group": f"{who} | {LABEL[a]}" + (f" | group {group}" if group else ""),
            "denominator_group": f"{who} | {LABEL[b]}" + (f" | group {group}" if group else ""),
            "held_constant": held, "pair_key": ("genotype|sampling_group" if cls == "day_night" else "genotype") if design == "panel" else "",
            "genotype": "B73" if design == "B73" else "", "sampling_group": group,
            "replicated": "TRUE" if len(num) >= 2 and len(den) >= 2 else "FALSE",
            "n_numerator": len(num), "n_denominator": len(den), "n_pairs": len(lines) if design == "panel" else "",
            "n_numerator_qc": len(nq), "n_denominator_qc": len(dq), "n_pairs_qc": len(lines_qc) if design == "panel" else "",
            "sign_note": f"log2FC = log2({a} / {b}); {sign}",
            "notes": notes + (f" {same_fc} of {len(lines_qc)} line pairs share a flowcell." if design == "panel" else ""),
            "tags": ";".join(t),
            "numerator_samples": ";".join(s["sample_id"] for s in num),
            "denominator_samples": ";".join(s["sample_id"] for s in den)})

    for cls, a, b, grp, tags, sign, notes in SPEC:
        groups = sorted({s["sampling_group"] for s in usable if s["tissue"] in (a, b)} - {""}, key=int) if grp == "each" else [""]
        for g in groups:
            side = lambda t, geno: [s for s in usable if s["tissue"] == t and geno(s["genotype"]) and (not g or s["sampling_group"] == g)]
            pair = (lambda s: f"{s['genotype']}|{s['sampling_group']}") if cls == "day_night" else (lambda s: s["genotype"])
            add(cls, a, b, "panel", g, tags, sign, notes, side(a, lambda x: x not in NOT_PANEL), side(b, lambda x: x not in NOT_PANEL), pair)
            if "subset_26_lines" in tags or "secondary" in tags:
                continue
            num, den = side(a, lambda x: x == "B73"), side(b, lambda x: x == "B73")
            if cls == "day_night":     # B73 by day was sampled in group 26 only
                gs = {s["sampling_group"] for s in num} & {s["sampling_group"] for s in den}
                num = [s for s in num if s["sampling_group"] in gs]; den = [s for s in den if s["sampling_group"] in gs]
                bnote = notes + f" B73 by day was sampled in group {', '.join(sorted(gs))} only."
            else:
                bnote = notes
            add(cls, a, b, "B73", "", tags, sign, bnote, num, den, None)
            if cls == "day_night":
                rows[-1]["sampling_group"] = ";".join(sorted(gs))
                rows[-1]["held_constant"] += f"; sampling_group={rows[-1]['sampling_group']}"
    return rows


def main(write=False):
    rows = build()
    print(f"{len(rows)} contrasts")
    for r in rows:
        print(f"  {r['contrast_id']:55s} {r['n_numerator']:>4}/{r['n_denominator']:<4} libraries"
              + (f", {r['n_pairs']} pairs ({r['n_pairs_qc']} after QC)" if r["design"] == "panel" else
                 f" ({r['n_numerator_qc']}/{r['n_denominator_qc']} after QC)") + f"  [{r['tags']}]")
    if write:
        with open(ROOT / "metadata" / "contrasts.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print("wrote metadata/contrasts.tsv")


if __name__ == "__main__":
    main(write="--write" in sys.argv[1:])
