#!/usr/bin/env python3
"""Write per-batch download manifests and nf-core/rnaseq samplesheets for Kremling2018.

The 1,960 libraries run as eight batches by tissue (L3Mid and LMid share one), so each
download fits on scratch (~40-85 GB) and each Nextflow run stays a manageable size. From
metadata/runs.tsv and metadata/samples.tsv (include = TRUE):

  pipeline/manifests/Kremling2018_<batch>.downloads.tsv   run, url, md5, dest            (ENA route)
  pipeline/manifests/Kremling2018_<batch>.sra.tsv         run, url, out1, out2, reads_per_file (AWS route)
  pipeline/manifests/Kremling2018_<batch>.samplesheet.csv sample,fastq_1,fastq_2,strandedness

Destination paths are on Hazel scratch (FASTQ_ROOT). Strandedness is 'forward': salmon found
stranded-forward reads in every test run (pipeline/read_structure/).

The AWS route (download_sra.sh, convert_sra.sbatch) fetches each run's .sra object from NCBI's
open-data bucket and converts it to the same FASTQ path, so the samplesheet serves both routes.
reads_per_file is ENA's read_count (n_spots): one read per spot in these single-end runs.

Usage: python3 -I pipeline/make_manifests.py
"""
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "pipeline" / "manifests"
STUDY = "Kremling2018"
FASTQ_ROOT = f"/share/maize/rrellan/MaizExpressD/fastq/{STUDY}"
BATCH = {"GRoot": "GRoot", "GShoot": "GShoot", "L3Base": "L3Base", "L3Tip": "L3Tip", "LMAD": "LMAD",
         "LMAN": "LMAN", "Kern": "Kern", "L3Mid": "LeafMid", "LMid": "LeafMid"}


def main():
    runs = {r["run_accession"]: r for r in csv.DictReader(
        open(ROOT / "metadata" / "runs.tsv"), delimiter="\t") if r["include"] == "TRUE"}
    samples = [s for s in csv.DictReader(open(ROOT / "metadata" / "samples.tsv"), delimiter="\t")
               if s["include"] == "TRUE"]
    OUTDIR.mkdir(parents=True, exist_ok=True)
    by_batch = {}
    for s in samples:
        by_batch.setdefault(BATCH[s["tissue"]], []).append(s)
    n_runs = Counter()
    for batch, ss in sorted(by_batch.items()):
        unit = f"{STUDY}_{batch}"
        dl, sra, sheet, nbytes = [], [], [], 0
        for s in sorted(ss, key=lambda s: s["sample_id"]):
            for acc in s["run_accessions"].split(";"):
                r = runs[acc]
                urls, md5s = r["fastq_urls"].split(";"), r["fastq_md5"].split(";")
                if r["library_layout"] != "SINGLE" or len(urls) != 1 or len(md5s) != 1:
                    sys.exit(f"{acc}: expected one single-end FASTQ")
                dest = f"{FASTQ_ROOT}/{batch}/{urls[0].rsplit('/', 1)[1]}"
                dl.append([acc, urls[0], md5s[0], dest])
                sra.append([acc, f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}", dest, "", r["n_spots"]])
                sheet.append([s["sample_id"], dest, "", "forward"])
                nbytes += int(r["n_bases"])
                n_runs[acc] += 1
        with open(OUTDIR / f"{unit}.downloads.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["run", "url", "md5", "dest"])
            w.writerows(dl)
        with open(OUTDIR / f"{unit}.sra.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["run", "url", "out1", "out2", "reads_per_file"])
            w.writerows(sra)
        with open(OUTDIR / f"{unit}.samplesheet.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
            w.writerows(sheet)
        print(f"{unit}: {len(sheet)} libraries, {nbytes / 1e9:.1f} Gbases", file=sys.stderr)
    if len(n_runs) != len(runs) or max(n_runs.values()) != 1:
        sys.exit("every included run must appear in exactly one batch")


if __name__ == "__main__":
    main()
