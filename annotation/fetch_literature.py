#!/usr/bin/env python3
"""Gene-specific literature for the gene card, from NCBI gene2pubmed.

gene2pubmed links NCBI Gene records to PubMed papers. For maize (taxid 4577) the NCBI GeneID maps
to B73 v5 gene models through NCBI's own cross-references in Zea_mays.gene_info (dbXrefs
"Ensembl:Zm00001eb..."), or else through its v4 locus tag (ZEAMMB73_Zm00001d...) and the v4-to-v5
table in annotation/genes.tsv.
Genome-scale papers link thousands of genes, so a paper counts as gene-specific when it is linked
to at most MAX_GENES maize v5 genes.

Downloads, cached in annotation/source/ncbi/ (--refresh fetches them again):
  gene2pubmed_maize.tsv.gz   maize rows of https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz
                             (the full file, ~290 MB, is streamed and not kept)
  Zea_mays.gene_info.gz      https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/Zea_mays.gene_info.gz
  pubmed_summaries.tsv.gz    pmid, year, first author, title, journal (E-utilities esummary)
  release.json               download date
Writes annotation/literature.tsv: gene_id, pmid, year, first_author, title, journal, n_genes
(gene-specific papers only). annotation/build_parquet.R turns it into Parquet.
Usage: python3 annotation/fetch_literature.py [--refresh]
"""
import csv
import datetime
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "annotation" / "source" / "ncbi"
G2P = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GENE_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/Zea_mays.gene_info.gz"
MAX_GENES = 10


def maize_rows(refresh):
    out = SRC / "gene2pubmed_maize.tsv.gz"
    if out.exists() and not refresh:
        return out
    n = 0
    tmp = out.with_suffix(".part")
    with urllib.request.urlopen(G2P, timeout=600) as r, gzip.open(r, "rt") as src, gzip.open(tmp, "wt") as dst:
        for line in src:
            if line.startswith("4577\t") or line.startswith("#"):
                dst.write(line)
                n += 1
    tmp.rename(out)
    print(f"gene2pubmed: {n:,} maize rows kept")
    return out


def summaries(pmids, refresh):
    path = SRC / "pubmed_summaries.tsv.gz"
    have = {}
    if path.exists() and not refresh:
        with gzip.open(path, "rt") as f:
            have = {r[0]: r for r in csv.reader(f, delimiter="\t")}
    todo = sorted(set(pmids) - set(have))
    for i in range(0, len(todo), 200):
        batch = todo[i:i + 200]
        data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "json",
                                       "tool": "maizexpressn"}).encode()
        for attempt in range(5):
            try:
                with urllib.request.urlopen(urllib.request.Request(ESUMMARY, data=data), timeout=120) as r:
                    res = json.load(r)["result"]
                break
            except Exception:
                time.sleep(2 + 3 * attempt)
        else:
            raise SystemExit(f"esummary failed for batch starting {batch[0]}")
        for p in batch:
            d = res.get(p, {})
            authors = d.get("authors") or []
            have[p] = [p, (d.get("pubdate") or "")[:4], authors[0]["name"] if authors else "",
                       (d.get("title") or "").strip(), d.get("source") or ""]
        time.sleep(0.4)                      # NCBI: at most 3 requests per second without an API key
    with gzip.open(path, "wt", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows(have[p] for p in sorted(have))
    if todo:
        print(f"pubmed: {len(todo):,} summaries fetched ({len(have):,} cached)")
    return have


def main():
    refresh = "--refresh" in sys.argv
    SRC.mkdir(parents=True, exist_ok=True)
    info = SRC / "Zea_mays.gene_info.gz"
    if refresh or not info.exists():
        urllib.request.urlretrieve(GENE_INFO, info)
    v4to5 = defaultdict(set)
    with open(ROOT / "annotation" / "genes.tsv") as f:
        for g in csv.DictReader(f, delimiter="\t"):
            for v4 in filter(None, g["v4_ids"].split(",")):
                v4to5[v4].add(g["gene_id"])
    locus = defaultdict(set)                 # NCBI GeneID -> v5 genes
    with gzip.open(info, "rt") as f:
        for r in csv.reader(f, delimiter="\t"):
            if r[0].startswith("#"):
                continue
            v5 = {x.split(":", 1)[1] for x in r[5].split("|") if x.startswith("Ensembl:Zm00001eb")}
            if not v5 and r[3].startswith("ZEAMMB73_"):
                v5 = v4to5.get(r[3][len("ZEAMMB73_"):], set())
            locus[r[1]] |= v5
    links = defaultdict(set)                  # pmid -> v5 genes
    with gzip.open(maize_rows(refresh), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            tax, gid, pmid = line.rstrip("\n").split("\t")[:3]
            for g in locus.get(gid, ()):
                links[pmid].add(g)
    specific = {p: gs for p, gs in links.items() if len(gs) <= MAX_GENES}
    summ = summaries(specific, refresh)
    out = sorted(([g, p] + summ[p][1:] + [len(gs)] for p, gs in specific.items() for g in gs),
                 key=lambda r: (r[0], -int(r[2] or 0), r[1]))
    with open(ROOT / "annotation" / "literature.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "pmid", "year", "first_author", "title", "journal", "n_genes"])
        w.writerows(out)
    (SRC / "release.json").write_text(json.dumps({"downloaded": datetime.date.today().isoformat(), "source": G2P,
                                                  "max_genes": MAX_GENES}, indent=1) + "\n")
    print(f"literature: {len(links):,} papers linked to maize v5 genes; {len(specific):,} gene-specific "
          f"(<= {MAX_GENES} genes) on {len({r[0] for r in out}):,} genes")


if __name__ == "__main__":
    main()
