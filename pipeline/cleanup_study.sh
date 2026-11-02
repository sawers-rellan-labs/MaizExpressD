#!/usr/bin/env bash
# Delete one batch's FASTQ and Nextflow work files from scratch after its results are
# confirmed. Refuses to run unless the salmon gene table exists. The run list and MD5s stay
# in pipeline/manifests and in the download and conversion logs (copied to the results
# folder), so everything can be re-downloaded. For the AWS route, the (by then empty) .sra
# folder goes too.
#   bash cleanup_study.sh <unit>        e.g. Kremling2018_GRoot
set -euo pipefail

UNIT="${1:?unit, e.g. Kremling2018_GRoot}"
STUDY="${UNIT%%_*}"
BATCH="${UNIT#*_}"
SCRATCH=/share/maize/rrellan/MaizExpressD
FASTQ="$SCRATCH/fastq/$STUDY/$BATCH"
OUT=/rsstu/users/r/rrellan/sara/MaizeExpressD/results/$STUDY/$BATCH
GENES="$OUT/salmon/salmon.merged.gene_counts.tsv"

[[ -s "$GENES" ]] || { echo "no $GENES: not cleaning $UNIT"; exit 1; }
cp -n "$FASTQ/download.log" "$OUT/download.log" 2>/dev/null || true
cp -n "$FASTQ/convert.log" "$OUT/convert.log" 2>/dev/null || true
cp -n "$SCRATCH/sra/$UNIT/download.log" "$OUT/sra_download.log" 2>/dev/null || true
du -sh "$FASTQ" "$SCRATCH/work/$UNIT" "$SCRATCH/sra/$UNIT" 2>/dev/null || true
read -r -p "Delete FASTQ and work files for $UNIT? [y/N] " ok
[[ "$ok" == "y" ]] || { echo "kept"; exit 0; }
rm -rf "$FASTQ" "$SCRATCH/work/$UNIT" "$SCRATCH/sra/$UNIT"
echo "removed scratch files for $UNIT"
