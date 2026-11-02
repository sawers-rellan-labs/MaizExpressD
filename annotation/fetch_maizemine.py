#!/usr/bin/env python3
"""Gene-card extras from MaizeMine, MaizeGDB's InterMine data warehouse: metabolic pathways,
GWAS hits near the gene, protein-changing variants, and pan-gene presence across the B73 and NAM
founder genomes.

MaizeGDB's own API (alpha.maizegdb.org/api/v1) sits behind a Cloudflare bot check that refuses
scripted requests (22 Oct 2026), so this uses MaizeMine's web services, which answer them:
  https://maizemine.rnet.missouri.edu/maizemine/service/query/results   (InterMine path queries)

Downloads, cached in annotation/source/maizemine/ (--refresh fetches them again):
  pathways.tsv.gz         gene, pathway id, pathway name (MaizeGDB E2P2 pathways data set)
  gwas.tsv.gz             GWAS hits (MaizeMine "QTL": GWAS Atlas and Wallace 2014 GWAS): id, trait,
                          source, chromosome, position, p-value, R2, PubMed id
  variants_high.tsv.gz    stop-gained, frameshift, splice acceptor/donor, start-lost and stop-lost
                          consequences (Ensembl variants): gene, variant id, type
  variants_missense.tsv.gz  missense consequences: gene, variant id, SIFT prediction
  pangene_members.tsv.gz  pan-gene group membership: group id, gene id, gene source (genome)
  atlas.tsv.gz            B73 expression atlas (MaizeMine Expression, mean FPKM per sample): gene, FPKM, sample
  atlas_samples.tsv.gz    atlas sample annotation: name, organ group, tissue, growth stage, qTeller id
  release.json            MaizeMine release, download date and the pan-genome genome order
Writes (B73 v5 genes only):
  annotation/pathways.tsv   gene_id, pathway_id, pathway_name
  annotation/gwas_hits.tsv  gene_id, hit_id, trait, source, chrom, position, distance_bp (0 = inside
                            the gene), p_value, r2, pubmed; hits within WINDOW bp of the gene
  annotation/variants.tsv   gene_id and the number of distinct variants per class: stop_gained,
                            frameshift, splice_site, start_lost, stop_lost, missense,
                            missense_deleterious (SIFT "deleterious", low-confidence calls included)
  annotation/pangenes.tsv   gene_id, pangene_id, genomes_present, n_genomes, missing (genomes without
                            a member), b73_genes_in_group
  annotation/atlas_samples.tsv  sample, organ_group, tissue, stage, order (organ group, then tissue and stage)
  annotation/atlas.tsv.gz   gene_id, sample, fpkm (B73 samples only)
annotation/build_parquet.R turns them into Parquet.
Usage: python3 annotation/fetch_maizemine.py [--refresh]
"""
import bisect
import csv
import datetime
import gzip
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "annotation" / "source" / "maizemine"
SERVICE = "https://maizemine.rnet.missouri.edu/maizemine/service"
V5 = re.compile(r"^Zm00001eb\d{6}$")
WINDOW = 10_000          # GWAS hits within 10 kb of the gene (maize LD mostly decays within a few kb)
ORGANS = ["Root", "Leaf", "Shoot", "Reproductive", "Seed"]
HIGH = ["STOP_GAINED", "FRAMESHIFT_VARIANT", "SPLICE_ACCEPTOR_VARIANT", "SPLICE_DONOR_VARIANT",
        "START_LOST", "STOP_LOST"]

QUERIES = {
    "pathways": """<query model="genomic" view="Gene.primaryIdentifier Gene.pathways.identifier Gene.pathways.name">
      <constraint path="Gene.source" op="=" value="B73 Zm00001eb.1"/></query>""",
    "gwas": """<query model="genomic" view="QTL.primaryIdentifier QTL.name QTL.source QTL.chromosome.primaryIdentifier
      QTL.chromosomeLocation.start QTL.pValue QTL.rSquared QTL.publications.pubMedId"
      ><join path="QTL.publications" style="OUTER"/></query>""",
    "variants_high": """<query model="genomic" view="Consequence.transcript.gene.primaryIdentifier Consequence.snpId
      Consequence.consequenceTypes.name">
      <constraint path="Consequence.consequenceTypes.name" op="ONE OF">%s</constraint></query>"""
      % "".join(f"<value>{t}</value>" for t in HIGH),
    "variants_missense": """<query model="genomic" view="Consequence.transcript.gene.primaryIdentifier Consequence.snpId
      Consequence.siftQualitativePrediction">
      <constraint path="Consequence.consequenceTypes.name" op="=" value="MISSENSE_VARIANT"/></query>""",
    "pangene_members": """<query model="genomic" view="PangeneGroup.primaryIdentifier
      PangeneGroup.genes.primaryIdentifier PangeneGroup.genes.source"></query>""",
    "atlas": """<query model="genomic" view="Expression.gene.primaryIdentifier Expression.FPKM Expression.sample.name">
      <constraint path="Expression.gene.source" op="=" value="B73 Zm00001eb.1"/></query>""",
    "atlas_samples": """<query model="genomic" view="Sample.name Sample.organGroup Sample.tissueDescription Sample.growthStage
      Sample.qTellerId"></query>""",
}


def fetch(name, refresh):
    out = SRC / f"{name}.tsv.gz"
    if out.exists() and not refresh:
        return out
    data = urllib.parse.urlencode({"query": " ".join(QUERIES[name].split()), "format": "tab"}).encode()
    req = urllib.request.Request(SERVICE + "/query/results", data=data)
    tmp = out.with_suffix(".part")
    n = 0
    with urllib.request.urlopen(req, timeout=3600) as r, gzip.open(tmp, "wb") as f:
        for line in r:
            f.write(line)
            n += 1
    tmp.rename(out)
    print(f"{name}: downloaded {n:,} rows")
    return out


def rows(path):
    with gzip.open(path, "rt") as f:
        yield from csv.reader(f, delimiter="\t")


def write(name, header, data):
    with open(ROOT / "annotation" / name, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(data)


def clean(x):
    return "" if x in ('""', "None") else x


def main():
    refresh = "--refresh" in sys.argv
    SRC.mkdir(parents=True, exist_ok=True)
    paths = {k: fetch(k, refresh) for k in QUERIES}
    rel = SRC / "release.json"
    if refresh or not rel.exists():
        with urllib.request.urlopen(SERVICE + "/version/release", timeout=60) as r:
            release = r.read().decode().strip().strip('"')
        rel.write_text(json.dumps({"maizemine_release": release, "downloaded": datetime.date.today().isoformat(),
                                   "service": SERVICE}, indent=1) + "\n")

    # pathways
    pw = sorted({(g, pid, name) for g, pid, name in rows(paths["pathways"]) if V5.match(g) and pid})
    write("pathways.tsv", ["gene_id", "pathway_id", "pathway_name"], pw)
    print(f"pathways: {len(pw):,} gene-pathway links, {len({p[0] for p in pw}):,} genes, {len({p[1] for p in pw}):,} pathways")

    # GWAS hits near each gene (gene coordinates from annotation/genes.tsv, B73 v5)
    hits = defaultdict(dict)                        # chrom -> hit id -> record (one per hit, first PubMed id)
    for hid, trait, source, chrom, pos, p, r2, pmid in rows(paths["gwas"]):
        if not pos.isdigit():
            continue
        h = hits[chrom].setdefault(hid, [int(pos), hid, trait, source, clean(p), clean(r2), clean(pmid)])
        if not h[6] and clean(pmid):
            h[6] = clean(pmid)
    by_chrom = {c: sorted(v.values()) for c, v in hits.items()}
    pos_idx = {c: [h[0] for h in v] for c, v in by_chrom.items()}
    gh = []
    with open(ROOT / "annotation" / "genes.tsv") as f:
        for g in csv.DictReader(f, delimiter="\t"):
            c = g["chrom"]
            if c not in by_chrom or not g["start"]:
                continue
            s, e = int(g["start"]), int(g["end"])
            lo = bisect.bisect_left(pos_idx[c], s - WINDOW)
            hi = bisect.bisect_right(pos_idx[c], e + WINDOW)
            for pos, hid, trait, source, p, r2, pmid in by_chrom[c][lo:hi]:
                d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
                gh.append([g["gene_id"], hid, trait, source, c, pos, d, p, r2, pmid])
    gh.sort(key=lambda r: (r[0], r[6], r[2]))
    write("gwas_hits.tsv", ["gene_id", "hit_id", "trait", "source", "chrom", "position", "distance_bp", "p_value", "r2", "pubmed"], gh)
    print(f"gwas: {sum(len(v) for v in by_chrom.values()):,} hits; {len(gh):,} hit-gene pairs within {WINDOW:,} bp "
          f"on {len({r[0] for r in gh}):,} genes")

    # protein-changing variants: distinct variant ids per gene and class
    cls = {"STOP_GAINED": "stop_gained", "FRAMESHIFT_VARIANT": "frameshift", "SPLICE_ACCEPTOR_VARIANT": "splice_site",
           "SPLICE_DONOR_VARIANT": "splice_site", "START_LOST": "start_lost", "STOP_LOST": "stop_lost"}
    seen = defaultdict(set)
    for g, snp, t in rows(paths["variants_high"]):
        if V5.match(g) and t in cls:
            seen[(g, cls[t])].add(snp)
    sift = defaultdict(int)
    for g, snp, pred in rows(paths["variants_missense"]):
        if V5.match(g):
            seen[(g, "missense")].add(snp)
            if "deleterious" in pred.lower():
                seen[(g, "missense_deleterious")].add(snp)
            sift[clean(pred) or "none"] += 1
    cols = ["stop_gained", "frameshift", "splice_site", "start_lost", "stop_lost", "missense", "missense_deleterious"]
    genes = sorted({g for g, _ in seen})
    write("variants.tsv", ["gene_id"] + cols, [[g] + [len(seen.get((g, c), ())) for c in cols] for g in genes])
    print(f"variants: {len(genes):,} genes; SIFT calls on missense rows: {dict(sift)}")

    # pan-gene groups: which genomes have a member
    members = defaultdict(list)
    genomes = set()
    for grp, gene, source in rows(paths["pangene_members"]):
        genome = source.split()[0] if source else ""
        if genome:
            genomes.add(genome)
            members[grp].append((gene, genome))
    order = ["B73"] + sorted(genomes - {"B73"})
    out = []
    for grp, mem in members.items():
        present = {gm for _, gm in mem}
        b73 = [g for g, gm in mem if gm == "B73" and V5.match(g)]
        for g in b73:
            out.append([g, grp, len(present), len(order), ",".join(x for x in order if x not in present), len(b73)])
    out.sort()
    write("pangenes.tsv", ["gene_id", "pangene_id", "genomes_present", "n_genomes", "missing", "b73_genes_in_group"], out)
    info = json.loads(rel.read_text())
    if info.get("genomes") != order:
        info["genomes"] = order
        rel.write_text(json.dumps(info, indent=1) + "\n")
    print(f"pan-genes: {len(members):,} groups over {len(order)} genomes; {len({r[0] for r in out}):,} B73 v5 genes")

    # B73 expression atlas: one mean FPKM per gene and sample; samples ordered by organ group, tissue, stage
    ann = {r[0]: [clean(x) for x in r[1:]] for r in rows(paths["atlas_samples"])}
    vals = [(g, smp, f) for g, f, smp in rows(paths["atlas"]) if V5.match(g) and smp.startswith("B73")]
    used = sorted({smp for _, smp, _ in vals},
                  key=lambda x: (ORGANS.index(ann.get(x, [""])[0]) if ann.get(x, [""])[0] in ORGANS else 99,
                                 ann.get(x, ["", ""])[1].lower(), ann.get(x, ["", "", ""])[2], x))
    write("atlas_samples.tsv", ["sample", "organ_group", "tissue", "stage", "order"],
          [[x] + (ann.get(x, ["", "", "", ""])[:3]) + [i + 1] for i, x in enumerate(used)])
    with gzip.open(ROOT / "annotation" / "atlas.tsv.gz", "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "sample", "fpkm"])
        w.writerows((g, smp, f"{float(fp):.3g}" if fp not in ("", '""') else "") for g, smp, fp in sorted(vals))
    print(f"atlas: {len(vals):,} values, {len({v[0] for v in vals}):,} genes x {len(used)} B73 samples")


if __name__ == "__main__":
    main()
