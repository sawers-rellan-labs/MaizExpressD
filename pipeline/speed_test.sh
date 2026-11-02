#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=speedtest
#SBATCH --output=/share/maize/rrellan/MaizExpressD/logs/speedtest_%j.out
# Compare download routes from the xfer node for Kremling2018 (one-off check, 29 Oct 2026).
# ENA FASTQ over HTTPS and FTP, and NCBI's SRA objects on AWS (sra-pub-run-odp), with 1 to 24
# parallel connections. Each connection fetches a different GRoot run to /dev/null, capped
# at 60 s; throughput = bytes received / wall time of the batch. Then one run is fetched from
# both sources and converted, to check the SRA object holds the same reads and qualities.
#   sbatch speed_test.sh
set -uo pipefail
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline
TMP=/share/maize/rrellan/MaizExpressD/tmp/speedtest
SIF=/share/maize/rrellan/MaizeExpressC/singularity/sra-tools_3.2.1.sif
mkdir -p "$TMP"
mapfile -t RUNS < <(tail -n +2 "$PIPE/manifests/Kremling2018_GRoot.downloads.tsv" | cut -f1 | head -24)
mapfile -t ENA < <(tail -n +2 "$PIPE/manifests/Kremling2018_GRoot.downloads.tsv" | cut -f2 | head -24)

batch() {   # batch <label> <n> <url-maker>
  local label="$1" n="$2" mk="$3" i t0 t1 bytes=0
  t0=$(date +%s.%N)
  for ((i = 0; i < n; i++)); do
    curl -s --max-time 60 -o /dev/null -w '%{size_download}\n' "$($mk "$i")" > "$TMP/$label.$i" &
  done
  wait
  t1=$(date +%s.%N)
  for ((i = 0; i < n; i++)); do bytes=$((bytes + $(cat "$TMP/$label.$i"))); done
  awk -v l="$label" -v n="$n" -v b="$bytes" -v t="$(echo "$t1 - $t0" | bc)" \
    'BEGIN {printf "%-10s %3d conn  %7.1f MB in %5.1f s  %6.1f MB/s total  %5.1f MB/s per conn\n", l, n, b/1e6, t, b/1e6/t, b/1e6/t/n}'
  rm -f "$TMP/$label".*
}
ena_https() { echo "${ENA[$1]}"; }
ena_ftp()   { echo "${ENA[$1]/https:\/\//ftp://}"; }
aws()       { echo "https://sra-pub-run-odp.s3.amazonaws.com/sra/${RUNS[$1]}/${RUNS[$1]}"; }

echo "== $(date '+%F %T') on $(hostname); other jobs: $(squeue -h -u "$USER" -t R -o %j | tr '\n' ' ')"
for n in 1 6 12 24; do batch ena_https "$n" ena_https; done
for n in 1 6; do batch ena_ftp "$n" ena_ftp; done
for n in 1 4 8; do batch aws "$n" aws; done

echo "== same run from both sources: ${RUNS[0]}"
cd "$TMP"
curl -s -o ena.fastq.gz "${ENA[0]}"
curl -s -o "${RUNS[0]}.sra" "$(aws 0)"
source /etc/profile.d/modules.sh; module load apptainer >/dev/null 2>&1
t0=$(date +%s)
singularity exec -B /share -B /gpfs_common "$SIF" fasterq-dump --split-files -e 1 -t "$TMP" -O "$TMP" "${RUNS[0]}.sra" > /dev/null 2>&1
echo "fasterq-dump on 1 CPU: $(( $(date +%s) - t0 )) s"
ls -la ena.fastq.gz "${RUNS[0]}.sra" "${RUNS[0]}"*.fastq 2>&1
echo "records: ENA $(( $(zcat ena.fastq.gz | wc -l) / 4 ))  SRA $(( $(cat "${RUNS[0]}"*.fastq | wc -l) / 4 ))"
# compare sequences and qualities (lines 2 and 4 of each record), ignoring read names
cmp <(zcat ena.fastq.gz | awk 'NR % 4 == 2 || NR % 4 == 0') \
    <(cat "${RUNS[0]}"*.fastq | awk 'NR % 4 == 2 || NR % 4 == 0') \
  && echo "sequences and qualities identical" || echo "sequences or qualities differ"
echo "ENA first record:"; zcat ena.fastq.gz | head -4
echo "SRA first record:"; head -4 "${RUNS[0]}"*.fastq
echo "distinct quality characters: ENA $(zcat ena.fastq.gz | awk 'NR % 4 == 0' | head -20000 | fold -w1 | sort -u | tr -d '\n')  SRA $(cat "${RUNS[0]}"*.fastq | awk 'NR % 4 == 0' | head -20000 | fold -w1 | sort -u | tr -d '\n')"
cd / && rm -rf "$TMP"
echo "== done $(date '+%F %T')"
