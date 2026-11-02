#!/usr/bin/env python3
"""Build metadata/runs.tsv, metadata/samples.tsv and metadata/lines.tsv for Kremling2018.

All 1,960 runs of BioProject PRJNA383416 (3' RNA-seq of the maize 282 panel in seven
tissues; Kremling et al. 2018) share one BioSample, so the design is read from the
library names: [plate_][submission_]<tissue>_<line>_<barcode>, e.g.
KAKRNA23_10346932_L3Base_A556_ACAGAT or 10343264_LMAN26_CI21E_ATGAAC. The flowcell
comes from the submitted file name (run_alias).

Inputs (metadata/build/): ena_runs/PRJNA383416.tsv (ENA run report), tissues.tsv
(tissue definitions, hand-curated), line_aliases.tsv (library-name spellings that
match a different name in Lei et al. 2025), label_fixes.tsv (tissue labels that the
expression contradicts, with the label used instead; the sample_id keeps the label in
the name, since the quantified tables are keyed by it), source/Lei2025_TableS1_lines.tsv
(the 298 lines with subpopulation, from extract_lei_tables.R).

Usage: python3 metadata/build/build_metadata.py [--write]
Without --write it only validates and prints the summary.
"""
import csv
import hashlib
import json
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "metadata" / "build"
CACHE = BUILD / "cache"
STUDY = "Kremling2018"
BIOPROJECT = "PRJNA383416"

RUN_COLS = ["run_accession", "study_id", "also_used_by", "bioproject", "biosample", "experiment",
            "sample_title", "library_name", "organism", "library_strategy", "library_source",
            "library_layout", "read_length", "instrument", "n_spots", "n_bases", "fastq_urls",
            "fastq_md5", "biosample_attributes", "include", "exclude_reason"]
SAMPLE_COLS = ["sample_id", "study_id", "run_accessions", "library_name", "genotype", "genotype_sra",
               "genotype_lei", "genotype_class", "subpopulation", "origin", "tissue", "tissue_lei",
               "tissue_label", "tissue_group", "tissue_detail", "po_term", "stage", "stage_class",
               "stage_order", "sampling_time", "sampling_group", "growth_system", "seq_submission",
               "plate", "flowcell", "barcode", "replicate", "include", "curation_note"]
LINE_COLS = ["genotype", "genotype_sra_spellings", "genotype_lei", "lei_match", "genotype_class",
             "subpopulation", "origin", "pedigree", "n_libraries", "tissues_with_data",
             "adult_leaf_group", "note"]
VOCAB = {
    "genotype_class": {"B73", "other inbred", "landrace", "NR"},
    "subpopulation": {"ss", "nss", "ts", "sweet", "popcorn", "mixed", "landrace", "NR"},
    "tissue_group": {"Root", "Shoot", "Leaf", "Kernel"},
    "stage_class": {"germination", "seedling", "adult", "kernel", "NR"},
    "include": {"TRUE", "FALSE"},
}
SKIP_TAGS = {"INSDC secondary accession", "INSDC center name", "INSDC first public",
             "INSDC last update", "INSDC status", "ENA-CHECKLIST", "ENA-FIRST-PUBLIC",
             "ENA-LAST-UPDATE", "SRA accession", "Submitter Id", "organism", "BioSampleModel"}


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def fetch(url):
    """GET with an on-disk cache (metadata/build/cache, not committed)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / hashlib.sha1(url.encode()).hexdigest()
    if path.exists():
        return path.read_text(errors="replace")
    req = urllib.request.Request(url, headers={"User-Agent": "MaizExpressD-metadata/0.1"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                text = r.read().decode("utf-8", errors="replace")
            path.write_text(text)
            return text
        except Exception:
            time.sleep(2 ** attempt)
    sys.exit(f"could not fetch {url}")


def biosample_json(acc):
    """The BioSample's alias, title and attributes as JSON (as in MaizExpressN/C runs.tsv)."""
    s = ET.fromstring(fetch("https://www.ebi.ac.uk/ena/browser/api/xml/" + acc)).find("SAMPLE")
    out = {"alias": s.get("alias", ""), "title": s.findtext("TITLE", "")}
    for a in s.iter("SAMPLE_ATTRIBUTE"):
        tag = a.findtext("TAG", "")
        if tag not in SKIP_TAGS:
            out[tag] = a.findtext("VALUE", "")
    return json.dumps(out, ensure_ascii=False)


def norm(name):
    """Comparison key for line names: lower case, letters and digits only."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def lei_short(name):
    """Lei's Table S1 name without the 282set_ prefix and Goodman-Buckler suffix."""
    return re.sub(r"Goodman-Buckler$", "", re.sub(r"^282set_", "", name))


def parse_library(name, label_to_tissue):
    """Split a library name into plate, submission, tissue label, line and barcode."""
    toks = name.split("_")
    hits = [i for i, t in enumerate(toks) if t in label_to_tissue]
    if len(hits) != 1 or hits[0] + 2 > len(toks):
        sys.exit(f"cannot parse library name {name}")
    i = hits[0]
    pre = toks[:i]
    plate = [p for p in pre if re.fullmatch(r"KAKRNA\d+", p)]
    sub = [p for p in pre if p.isdigit()]
    if len(plate) > 1 or len(sub) > 1 or len(plate) + len(sub) != len(pre):
        sys.exit(f"unexpected prefix in library name {name}")
    return {"plate": plate[0] if plate else "", "seq_submission": sub[0] if sub else "",
            "label": toks[i], "genotype_sra": "_".join(toks[i + 1:-1]), "barcode": toks[-1]}


def flowcell(run_alias):
    """Flowcell ID from the submitted file name, e.g. 5999_8601_26956_HLTVYBGXX_..."""
    toks = run_alias.split("_")
    return toks[3] if len(toks) > 3 and re.fullmatch(r"H[A-Z0-9]{8}", toks[3]) else ""


def main(write=False):
    ena = read_tsv(BUILD / "ena_runs" / f"{BIOPROJECT}.tsv")
    tissues = {t["tissue"]: t for t in read_tsv(BUILD / "tissues.tsv")}
    label_to_tissue = {lab: t for t, row in tissues.items() for lab in row["sra_labels"].split(";")}
    lei = {r["inbred_lei"]: r for r in read_tsv(BUILD / "source" / "Lei2025_TableS1_lines.tsv")}
    lei_by_norm = {norm(lei_short(n)): n for n in lei}
    aliases = {a["genotype_sra"]: a for a in read_tsv(BUILD / "line_aliases.tsv")}
    fixes = {f["run_accession"]: f for f in read_tsv(BUILD / "label_fixes.tsv")}
    problems = []

    # Parse every run.
    runs = []
    for r in sorted(ena, key=lambda r: r["run_accession"]):
        p = parse_library(r["library_name"], label_to_tissue)
        p["name_tissue"] = label_to_tissue[p["label"]]
        f = fixes.pop(r["run_accession"], None)
        if f:
            if f["label_in_name"] != p["label"] or f["label_used"] not in label_to_tissue:
                problems.append(f"label_fixes.tsv: {r['run_accession']} does not match its library name")
            p["fix"] = f
            p["label"] = f["label_used"]
        p.update(run=r, tissue=label_to_tissue[p["label"]], flowcell=flowcell(r["run_alias"]))
        m = re.fullmatch(r"LMA[DN](\d+)", p["label"])
        p["sampling_group"] = m.group(1) if m else ""
        runs.append(p)
    problems += [f"label_fixes.tsv: {acc} is not a run of {BIOPROJECT}" for acc in fixes]

    # Match each library-name spelling to Lei's Table S1.
    match = {}
    for g in sorted({p["genotype_sra"] for p in runs}):
        if g in aliases:
            a = aliases[g]
            if a["genotype_lei"] not in lei:
                problems.append(f"line_aliases.tsv: {a['genotype_lei']} is not in Lei Table S1")
            match[g] = (a["genotype_lei"], a["kind"])
        elif norm(g) in lei_by_norm:
            match[g] = (lei_by_norm[norm(g)], "direct")
        else:
            match[g] = ("", "none")

    # One canonical genotype per Lei line: the spelling with most libraries (ties: alphabetical).
    libs_per_spelling = Counter(p["genotype_sra"] for p in runs)
    spellings = defaultdict(list)
    for g, (ln, _) in match.items():
        spellings[ln or g].append(g)
    canon = {}
    for key, gs in spellings.items():
        best = sorted(gs, key=lambda g: (-libs_per_spelling[g], g))[0]
        for g in gs:
            canon[g] = best

    def genotype_class(g):
        if g == "B73":
            return "B73"
        if re.match(r"MR\d\d", g):
            return "landrace"
        if g == "B73-10":
            return "NR"
        return "other inbred"

    def subpopulation(lei_name):
        s = lei[lei_name]["subpopulation"] if lei_name else ""
        return s if s in VOCAB["subpopulation"] - {"NR"} else "NR"

    # Samples: one library per run; replicates numbered within genotype x tissue by run accession.
    rep = Counter()
    samples = []
    for p in runs:
        r, t = p["run"], tissues[p["tissue"]]
        g = canon[p["genotype_sra"]]
        lei_name, kind = match[p["genotype_sra"]]
        rep[(g, p["name_tissue"])] += 1
        notes = [f"library {r['library_name']}"]
        if "fix" in p:
            notes.append(f"labelled {p['fix']['label_in_name']} in the library name; {p['fix']['label_used']} by expression "
                         f"(see label_fixes.tsv)")
        if kind in ("spelling", "alias", "probable"):
            notes.append(f"Lei name {lei_name} ({kind} match)")
        if not lei_name:
            notes.append("not in Lei et al. 2025 Table S1")
        if p["tissue"] == "LMid":
            notes.append("tissue unclear (LMid, not L3Mid); see needs_review.tsv")
        if g == "B73-10":
            notes.append("genotype label B73-10; see needs_review.tsv")
        samples.append({
            "sample_id": f"{STUDY}_{p['name_tissue']}_{g}_R{rep[(g, p['name_tissue'])]}",
            "study_id": STUDY, "run_accessions": r["run_accession"], "library_name": r["library_name"],
            "genotype": g, "genotype_sra": p["genotype_sra"], "genotype_lei": lei_name,
            "genotype_class": genotype_class(g), "subpopulation": subpopulation(lei_name),
            "origin": lei[lei_name]["origin"] if lei_name else "",
            "tissue": p["tissue"], "tissue_lei": t["lei_code"], "tissue_label": t["tissue_label"],
            "tissue_group": t["tissue_group"], "tissue_detail": t["tissue_detail"], "po_term": t["po_term"],
            "stage": t["stage"], "stage_class": t["stage_class"], "stage_order": t["stage_order"],
            "sampling_time": t["sampling_time"], "sampling_group": p["sampling_group"],
            "growth_system": t["growth_system"], "seq_submission": p["seq_submission"], "plate": p["plate"],
            "flowcell": p["flowcell"], "barcode": p["barcode"], "replicate": rep[(g, p["name_tissue"])],
            "include": "TRUE", "curation_note": "; ".join(notes)})

    # Runs table, in the MaizExpressN/C layout so the pipeline scripts read it unchanged.
    attrs = biosample_json(ena[0]["sample_accession"])
    run_rows = []
    for p in runs:
        r = p["run"]
        n_spots, n_bases = int(r["read_count"]), int(r["base_count"])
        run_rows.append({
            "run_accession": r["run_accession"], "study_id": STUDY, "also_used_by": "",
            "bioproject": r["study_accession"], "biosample": r["sample_accession"],
            "experiment": r["experiment_accession"], "sample_title": r["sample_title"],
            "library_name": r["library_name"], "organism": r["scientific_name"],
            "library_strategy": r["library_strategy"], "library_source": r["library_source"],
            "library_layout": r["library_layout"], "read_length": round(n_bases / n_spots) if n_spots else "",
            "instrument": r["instrument_model"], "n_spots": n_spots, "n_bases": n_bases,
            "fastq_urls": ";".join("https://" + u for u in r["fastq_ftp"].split(";") if u),
            "fastq_md5": r["fastq_md5"], "biosample_attributes": attrs, "include": "TRUE", "exclude_reason": ""})

    # Lines table: every genotype in the library names, then Lei lines with no libraries.
    by_line = defaultdict(list)
    for s in samples:
        by_line[s["genotype"]].append(s)
    order = [t for t in sorted(tissues, key=lambda t: float(tissues[t]["display_order"] or 99))]
    lines = []
    for g in sorted(by_line, key=str.lower):
        ss = by_line[g]
        lei_name = ss[0]["genotype_lei"]
        kinds = {match[s["genotype_sra"]][1] for s in ss}
        per_t = Counter(s["tissue"] for s in ss)
        groups = sorted({s["sampling_group"] for s in ss if s["sampling_group"]}, key=int)
        lines.append({
            "genotype": g, "genotype_sra_spellings": ";".join(sorted({s["genotype_sra"] for s in ss})),
            "genotype_lei": lei_name, "lei_match": ";".join(sorted(kinds)),
            "genotype_class": ss[0]["genotype_class"], "subpopulation": ss[0]["subpopulation"],
            "origin": ss[0]["origin"], "pedigree": lei[lei_name]["pedigree"] if lei_name else "",
            "n_libraries": len(ss), "tissues_with_data": ";".join(f"{t}:{per_t[t]}" for t in order if per_t[t]),
            "adult_leaf_group": ";".join(groups), "note": ""})
    matched_lei = {s["genotype_lei"] for s in samples if s["genotype_lei"]}
    for n in sorted(set(lei) - matched_lei):
        lines.append({"genotype": "", "genotype_sra_spellings": "", "genotype_lei": n, "lei_match": "none",
                      "genotype_class": "", "subpopulation": subpopulation(n), "origin": lei[n]["origin"],
                      "pedigree": lei[n]["pedigree"], "n_libraries": 0, "tissues_with_data": "",
                      "adult_leaf_group": "", "note": "in Lei et al. 2025 Table S1; no library name matches"})

    # Checks.
    if len(run_rows) != len({r["run_accession"] for r in run_rows}):
        problems.append("duplicate run accessions")
    if len(samples) != len({s["sample_id"] for s in samples}):
        problems.append("duplicate sample IDs")
    for s in samples:
        for col, allowed in VOCAB.items():
            if s[col] not in allowed:
                problems.append(f"{s['sample_id']}: {col} = {s[col]!r}")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", s["sample_id"]):
            problems.append(f"sample_id with unexpected characters: {s['sample_id']}")
        if not s["flowcell"]:
            problems.append(f"{s['sample_id']}: no flowcell in run_alias")
    for r in run_rows:
        if len(r["fastq_urls"].split(";")) != len(r["fastq_md5"].split(";")) or not r["fastq_md5"]:
            problems.append(f"{r['run_accession']}: FASTQ URLs and MD5s do not match")

    # Summary.
    print(f"{len(run_rows)} runs, {len(samples)} libraries, {len(by_line)} genotypes "
          f"({sum(1 for l in lines if l['genotype'] and l['genotype_lei'])} matched to Lei Table S1; "
          f"{len(set(lei) - matched_lei)} Lei lines without libraries)")
    print("libraries and genotypes per tissue:")
    for t in order:
        ss = [s for s in samples if s["tissue"] == t]
        print(f"  {t:7s} {tissues[t]['lei_code'] or '-':3s} {len(ss):4d} libraries, "
              f"{len({s['genotype'] for s in ss}):3d} genotypes, "
              f"{sum(1 for s in ss if s['genotype_lei']):4d} libraries of Lei lines")
    print("genotype classes:", dict(Counter(l["genotype_class"] for l in lines if l["genotype"])))
    print("subpopulations:", dict(Counter(l["subpopulation"] for l in lines if l["genotype"])))
    print(f"total reads {sum(r['n_spots'] for r in run_rows) / 1e9:.2f} G; "
          f"FASTQ {sum(int(b) for r in ena for b in r['fastq_bytes'].split(';') if b) / 1e9:.0f} GB")
    if problems:
        print(f"\n{len(problems)} problem(s):")
        print("\n".join("  " + p for p in problems[:50]))
        sys.exit(1)
    if write:
        write_tsv(ROOT / "metadata" / "runs.tsv", RUN_COLS, run_rows)
        write_tsv(ROOT / "metadata" / "samples.tsv", SAMPLE_COLS, samples)
        write_tsv(ROOT / "metadata" / "lines.tsv", LINE_COLS, lines)
        print("wrote metadata/runs.tsv, samples.tsv, lines.tsv")


if __name__ == "__main__":
    main(write="--write" in sys.argv[1:])
