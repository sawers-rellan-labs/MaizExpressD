#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/share/maize/rrellan/MaizExpressD/logs/download_%x_%j.out
# Download one batch's FASTQ on Hazel's data-transfer node (xfer partition).
#   sbatch --job-name=dl_<unit> download_study.sbatch <unit> [connections, default 6]
# Re-submitting is safe: verified files are skipped, partial files resume.
# QoS xfer allows 4 CPUs / 24 GB per user in total: 1 CPU per download (curl is I/O-bound)
# leaves room for two downloads plus a Nextflow head job.
set -uo pipefail
UNIT="${1:?unit, e.g. Kremling2018_GRoot}"
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline
bash "$PIPE/download_fastq.sh" "$PIPE/manifests/$UNIT.downloads.tsv" "${2:-6}" || exit 1
