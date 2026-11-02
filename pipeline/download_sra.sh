#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/share/maize/rrellan/MaizExpressD/logs/sra_%x_%j.out
# Download one batch's .sra files from NCBI's AWS open-data bucket (xfer node). From MaizExpressC.
#   sbatch --job-name=sra_<unit> download_sra.sh <unit> [parallel, default 8]
# The speed test (speed_test.out, 29 Oct 2026) gave ~250 MB/s with 8 connections, against
# 13.5 MB/s at best from ENA. Each file is resumed until its size equals the bucket's
# Content-Length; integrity is then checked by vdb-validate in convert_sra.sbatch (SRA files
# carry their own checksums). Runs missing from the bucket (6 of 1,960 return 404, 29 Oct
# 2026) are fetched from ENA instead, MD5-checked, straight to their FASTQ path (download_fastq.sh
# with those rows of <unit>.downloads.tsv); convert_sra.sbatch then skips them.
# Re-submitting is safe: complete files are skipped, and runs whose FASTQ already exist are not
# fetched. Writes <sra_dir>/download.log.
set -uo pipefail
UNIT="${1:?unit, e.g. Kremling2018_GRoot}"
JOBS="${2:-8}"
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline
SRA_DIR=/share/maize/rrellan/MaizExpressD/sra/$UNIT
mkdir -p "$SRA_DIR"

fetch() {
  local run="$1" url="$2" dest="$SRA_DIR/$1.sra" hdr want got try
  hdr=$(curl -sI --max-time 60 "$url" | tr -d '\r')
  if [[ "$(awk 'NR == 1 {print $2}' <<< "$hdr")" == 404 ]]; then
    echo -e "$(date +%F_%T)\t$run\tNOT_ON_AWS" >> "$SRA_DIR/download.log"; return 0
  fi
  want=$(awk 'tolower($1) == "content-length:" {print $2}' <<< "$hdr")
  [[ -z "$want" ]] && { echo -e "$(date +%F_%T)\t$run\tNO_SIZE" >> "$SRA_DIR/download.log"; return 1; }
  if [[ -s "$dest" && "$(stat -c %s "$dest")" == "$want" ]]; then
    echo -e "$(date +%F_%T)\t$run\tSKIP" >> "$SRA_DIR/download.log"; return 0
  fi
  for try in $(seq 1 30); do
    curl -fsSL --connect-timeout 60 --speed-limit 100000 --speed-time 120 -C - -o "$dest.part" "$url"
    got=$(stat -c %s "$dest.part" 2>/dev/null || echo 0)
    [[ "$got" == "$want" ]] && { mv "$dest.part" "$dest"; echo -e "$(date +%F_%T)\t$run\tOK\t$want" >> "$SRA_DIR/download.log"; return 0; }
    (( got > want )) && rm -f "$dest.part"    # overshoot: start again
    sleep $(( try < 6 ? try * 10 : 60 ))
  done
  echo -e "$(date +%F_%T)\t$run\tINCOMPLETE\t$got/$want" >> "$SRA_DIR/download.log"; return 1
}
export -f fetch; export SRA_DIR

status=0
# read splits on runs of tabs, so an empty out2 (single-end runs) is passed as "-"
tail -n +2 "$PIPE/manifests/$UNIT.sra.tsv" | awk -F'\t' '{print $1, $2, $3, ($4 == "" ? "-" : $4)}' |
while read -r run url out1 out2; do
  out2=${out2#-}
  # runs already converted need no .sra
  [[ -s "$out1" && ( -z "$out2" || -s "$out2" ) ]] && continue
  printf '%s\t%s\n' "$run" "$url"
done | xargs -P "$JOBS" -L 1 bash -c 'fetch "$0" "$1"' || status=1

# runs not on AWS: ENA FASTQ instead (header line plus their rows of the ENA manifest)
missing=$(awk -F'\t' '$3 == "NOT_ON_AWS" {print $2}' "$SRA_DIR/download.log" 2>/dev/null | sort -u)
if [[ -n "$missing" ]]; then
  awk -F'\t' 'NR == FNR {m[$1]; next} FNR == 1 || $1 in m' <(echo "$missing") \
    "$PIPE/manifests/$UNIT.downloads.tsv" > "$SRA_DIR/ena_fallback.tsv"
  echo "not on AWS, fetching from ENA: $(echo $missing)"
  bash "$PIPE/download_fastq.sh" "$SRA_DIR/ena_fallback.tsv" 6 || status=1
fi
echo "$(find "$SRA_DIR" -name '*.sra' | wc -l) .sra files ($(du -sh "$SRA_DIR" | cut -f1)) for $(( $(wc -l < "$PIPE/manifests/$UNIT.sra.tsv") - 1 )) runs; failures: $(grep -c -P "\tNO_SIZE|\tINCOMPLETE" "$SRA_DIR/download.log")"
exit $status
