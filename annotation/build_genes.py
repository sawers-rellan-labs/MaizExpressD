#!/usr/bin/env python3
"""M5: build the gene annotation table (HANDOFF section 5, annotation/genes.parquet).

Inputs in annotation/source/ (downloaded 19 Oct 2026):
  MaizeGDB  https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/
    Zm00001eb.1.fulldata.txt.gz         symbol, full name, gene product, GO terms, location (no header)
    Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz                    strand (fallback)
    Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.genemodel_locus.txt.gz     NCBI GeneID per gene
    Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt.gz   pre-release Zm00001e IDs
  MaizeGDB  https://download.maizegdb.org/Pan-genes/B73_gene_xref/
    B73v5_to_B73v4.tsv, B73v5_to_B73v3.tsv                          v4 and v3 IDs
  Ensembl Plants BioMart (plants_mart, zmays_eg_gene)
    ensembl_zmays_genes.tsv                                         strand, biotype, name, description
    ensembl_zmays_athaliana_orthologs.tsv                           Arabidopsis orthologs and % identity
  MaizeGDB  Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_entap_results.tsv.gz  predicted descriptions (EnTAP:
    EggNOG description, else an informative similarity hit), canonical transcript only

Writes annotation/genes.tsv and annotation/go_terms.tsv; annotation/build_parquet.R turns them
into Parquet. Columns (genes): gene_id, v4_ids, v3_ids, symbol, full_name, synonyms, description,
description_source (MaizeGDB, Ensembl, EnTAP-EggNOG, EnTAP-similarity), chrom, start, end, strand,
biotype, go_terms, ath_best_hit, ath_symbol, ath_orthologs, quantified.
Multi-valued fields are comma-separated. ath_best_hit is the Ensembl ortholog with the highest mean
of the two % identities; ath_orthologs lists every ortholog as AGI or AGI:symbol, for search.
Usage: python3 annotation/build_genes.py
"""
import csv
import glob
import gzip
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "annotation" / "source"
NONE = {"", "-", "NULL"}

def opn(name):
    p = SRC / name
    return gzip.open(p, "rt", errors="replace") if name.endswith(".gz") else open(p, errors="replace")

genes = {}
go_names = {}

def gene(gid):
    if gid not in genes:
        genes[gid] = {"gene_id": gid, "v4_ids": "", "v3_ids": "", "symbol": "", "full_name": "", "synonyms": [],
                      "description": "", "description_source": "", "chrom": "", "start": "", "end": "", "strand": "", "biotype": "",
                      "go_terms": [], "ath_best_hit": "", "ath_symbol": "", "ath_orthologs": [], "quantified": False}
    return genes[gid]

# MaizeGDB gene table
for line in opn("Zm00001eb.1.fulldata.txt.gz"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 14 or not f[1].startswith("Zm00001eb"):
        continue
    g = gene(f[1])
    g["chrom"], g["start"], g["end"] = f[7], f[8], f[9]
    if f[10] not in NONE:
        g["symbol"] = f[10]
    if f[11] not in NONE:
        g["full_name"] = f[11]
    if f[12] not in NONE:
        g["description"] = ", ".join(dict.fromkeys(p.strip() for p in f[12].split(",") if p.strip()))
        g["description_source"] = "MaizeGDB"
    for term in re.findall(r"(GO:\d{7})=([^,]*(?:,(?!GO:)[^,]*)*)", f[13]):
        g["go_terms"].append(term[0])
        go_names.setdefault(term[0], term[1].strip())

# v4 / v3 IDs
for col, fname in (("v4_ids", "B73v5_to_B73v4.tsv"), ("v3_ids", "B73v5_to_B73v3.tsv")):
    for line in opn(fname):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0].startswith("Zm00001eb") and f[1] not in NONE:
            gene(f[0])[col] = ",".join(dict.fromkeys(x for x in f[1].split(",") if x))

# Ensembl gene features
ens_name = {}
with opn("ensembl_zmays_genes.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or not f[0].startswith("Zm00001eb"):
            continue
        g = gene(f[0])
        if not g["chrom"]:
            g["chrom"], g["start"], g["end"] = "chr" + f[1] if f[1].isdigit() else f[1], f[2], f[3]
        g["strand"] = {"1": "+", "-1": "-"}.get(f[4], "")
        g["biotype"] = f[5]
        if f[6]:
            ens_name[f[0]] = f[6]
        if not g["description"] and f[7]:
            g["description"] = re.sub(r"\s*\[Source:.*?\]$", "", f[7])
            g["description_source"] = "Ensembl"

# predicted descriptions from EnTAP for genes still without one (canonical transcript)
canonical = {line.strip() for line in opn("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts.gz")}
with opn("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_entap_results.tsv.gz") as fh:
    head = next(fh).rstrip("\n").split("\t")
    i_desc, i_inf, i_egg = head.index("Description"), head.index("Informative"), head.index("EggNOG Description")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] not in canonical:
            continue
        g = genes.get(f[0].split("_")[0])
        if g is None or g["description"]:
            continue
        egg = f[i_egg] if len(f) > i_egg else ""
        if egg and egg != "NaN":
            g["description"], g["description_source"] = egg, "EnTAP-EggNOG"
        elif len(f) > i_inf and f[i_inf].lower().startswith("yes") and f[i_desc] not in ("", "NaN"):
            d = re.sub(r"^\S+\s+", "", f[i_desc])          # drop the subject accession
            d = re.sub(r"\s*\[[^\]]*\]$", "", d)          # drop the [species] suffix
            g["description"], g["description_source"] = d, "EnTAP-similarity"

# strand fallback from the GFF3
for line in opn("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if len(f) > 8 and f[2] == "gene":
        m = re.search(r"ID=(?:gene:)?(Zm00001eb\d+)", f[8])
        if m and m.group(1) in genes and not genes[m.group(1)]["strand"]:
            genes[m.group(1)]["strand"] = f[6]

# synonyms: Ensembl name if different, NCBI GeneID, pre-release Zm00001e ID
for gid, nm in ens_name.items():
    g = genes[gid]
    if nm and nm != g["symbol"]:
        g["synonyms"].append(nm)
with opn("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.genemodel_locus.txt.gz") as fh:
    next(fh)
    for line in fh:
        f = line.split()
        if len(f) >= 2 and f[0] in genes and f[1].isdigit():
            genes[f[0]]["synonyms"].append("NCBI:" + f[1])
for line in opn("Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt.gz"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[8] in genes and f[3].startswith("Zm00001e"):
        genes[f[8]]["synonyms"].append(f[3])

# Arabidopsis orthologs
best = {}
with opn("ensembl_zmays_athaliana_orthologs.tsv") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or not f[1] or f[0] not in genes:
            continue
        agi, nm = f[1], f[2]
        label = agi if (not nm or nm == agi) else f"{agi}:{nm}"
        g = genes[f[0]]
        if label not in g["ath_orthologs"]:
            g["ath_orthologs"].append(label)
        score = (float(f[4] or 0) + float(f[5] or 0)) / 2
        if f[0] not in best or score > best[f[0]][0]:
            best[f[0]] = (score, agi, "" if nm == agi else nm)
for gid, (score, agi, nm) in best.items():
    genes[gid]["ath_best_hit"], genes[gid]["ath_symbol"] = agi, nm

# genes quantified by salmon (the gene set of the nf-core gene tables)
tabs = sorted(glob.glob(str(ROOT / "results" / "*" / "salmon.merged.gene_counts.tsv")))
quant = set()
if tabs:
    with open(tabs[0]) as fh:
        next(fh)
        quant = {line.split("\t", 1)[0] for line in fh}
for gid in quant:
    gene(gid)["quantified"] = True

cols = ["gene_id", "v4_ids", "v3_ids", "symbol", "full_name", "synonyms", "description", "description_source", "chrom", "start", "end",
        "strand", "biotype", "go_terms", "ath_best_hit", "ath_symbol", "ath_orthologs", "quantified"]
with open(ROOT / "annotation" / "genes.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for gid in sorted(genes):
        g = genes[gid]
        row = []
        for c in cols:
            v = g[c]
            if isinstance(v, list):
                v = ",".join(dict.fromkeys(v))
            elif isinstance(v, bool):
                v = "TRUE" if v else "FALSE"
            row.append(v)
        w.writerow(row)
with open(ROOT / "annotation" / "go_terms.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["go_id", "go_name"])
    w.writerows(sorted(go_names.items()))

n = len(genes)
cnt = lambda c: sum(1 for g in genes.values() if g[c])
print(f"{n} genes ({len(quant)} quantified, {sum(1 for q in quant if q not in genes or not genes[q]['chrom'])} of them without a location)")
import collections
print("  description sources:", dict(collections.Counter(g["description_source"] or "none" for g in genes.values())))
for c in ("symbol", "description", "v4_ids", "v3_ids", "strand", "go_terms", "ath_best_hit", "synonyms"):
    print(f"  {c:13s} {cnt(c):6d}  ({cnt(c) / n:.0%})")
print(f"  {len(go_names)} GO terms")
