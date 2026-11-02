#!/usr/bin/env bash
# Queue the AWS SRA route for one or more batches: sra_<unit> downloads the .sra files on the
# xfer node (download_sra.sh), then cv_<unit> converts them to FASTQ on a compute node
# (convert_sra.sbatch) once the download has succeeded. Run queue_studies.sh afterwards: each
# head job then waits for its cv_<unit> job.
#   bash queue_sra.sh <unit> [<unit> ...]     e.g. Kremling2018_GRoot Kremling2018_GShoot
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p /share/maize/rrellan/MaizExpressD/logs      # Slurm needs the log folder to exist
for u in "$@"; do
  [[ -s "manifests/$u.sra.tsv" ]] || { echo "$u: no manifests/$u.sra.tsv"; exit 1; }
  sra=$(sbatch --parsable --job-name="sra_$u" download_sra.sh "$u")
  cv=$(sbatch --parsable --job-name="cv_$u" --dependency="afterok:$sra" --kill-on-invalid-dep=yes \
         convert_sra.sbatch "$u")
  echo -e "$u\tsra_$u $sra\tcv_$u $cv"
done
