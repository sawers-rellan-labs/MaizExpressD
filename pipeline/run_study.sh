#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=7G
#SBATCH --output=/share/maize/rrellan/MaizExpressD/logs/%x_%j.out
# Nextflow head job for one batch of Kremling2018 (e.g. Kremling2018_GRoot): salmon
# quantification with nf-core/rnaseq. Runs on the xfer node (the only nodes with internet:
# pipeline code, nf-core configs and container pulls); each pipeline step is submitted as its
# own Slurm job on the compute partition (-profile hazel), which needs no internet.
#   sbatch --job-name=<unit> run_study.sh <unit>
# Needs: FASTQ downloaded (download_study.sbatch). The reference comes from MaizExpressN.
set -euo pipefail

UNIT="${1:?unit, e.g. Kremling2018_GRoot}"
STUDY="${UNIT%%_*}"
BATCH="${UNIT#*_}"
SHARE=/rsstu/users/r/rrellan/sara/MaizeExpressD
PIPE="$SHARE/pipeline"
SCRATCH=/share/maize/rrellan/MaizExpressD
OUT="$SHARE/results/$STUDY/$BATCH"
REF_SAVED="$SCRATCH/reference/nfcore_genome"   # transcript FASTA, GTF, salmon index

mkdir -p "$SCRATCH/logs" "$SCRATCH/work/$UNIT" "$SCRATCH/singularity"
export NXF_HOME="$SCRATCH/.nextflow"            # keep Nextflow's cache out of $HOME
export NXF_SINGULARITY_CACHEDIR="$SCRATCH/singularity"
# peak RSS was 5.8 GB at -Xmx6g; at 7G two head jobs and two downloads fit xfer's 24 GB per user
export NXF_OPTS='-Xms1g -Xmx4g'
export APPTAINER_TMPDIR="$SCRATCH/tmp" SINGULARITY_TMPDIR="$SCRATCH/tmp"   # xfer /tmp is only 16 GB
export APPTAINER_CACHEDIR="$SCRATCH/singularity/apptainer_cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

source /etc/profile.d/modules.sh                # batch shells don't get 'module' by default
module purge
module load nextflow/25                         # 25.10.3 launcher; NXF_VER fetches the pinned release
module load singularity                         # the head job pulls any image missing from the BRC library
export NXF_VER=25.10.8                          # nf-core/rnaseq 3.27.0 needs >= 25.10.4
NXF_VERSION=$(nextflow -version 2>/dev/null | grep -oE 'version [0-9.]+' | awk '{print $2}')
if printf '25.10.4\n%s\n' "$NXF_VERSION" | sort -V -C; then RNASEQ=3.27.0; else RNASEQ=3.24.0; fi
echo "unit=$UNIT nextflow=$NXF_VERSION nf-core/rnaseq=$RNASEQ"

# stop before Nextflow starts if any FASTQ in the sample sheet is missing (e.g. purged)
missing=$(tail -n +2 "$PIPE/manifests/$UNIT.samplesheet.csv" | cut -d, -f2,3 | tr ',' '\n' |
          while read -r f; do [[ -z "$f" || -s "$f" ]] || echo "$f"; done)
if [[ -n "$missing" ]]; then
  echo "missing FASTQ, download first:"; echo "$missing"; exit 1
fi

# the B73 v5 reference bundle (genome, gene models, filtered GTF, transcript FASTA, salmon index)
# is MaizExpressN's, same as for MaizExpressC; restore anything scratch purged (30-day rule)
REF_BACKUP=/rsstu/users/r/rrellan/sara/MaizeExpressN/reference
rsync -a --ignore-existing "$REF_BACKUP/" "$SCRATCH/reference/"
[[ -d "$REF_SAVED/index/salmon" ]] || { echo "no salmon index in $REF_SAVED"; exit 1; }
REF_ARGS=(--salmon_index "$REF_SAVED/index/salmon"
          --gtf "$(ls "$REF_SAVED"/*.gtf | head -1)"     # the filtered GTF, as in N and C
          --transcript_fasta "$REF_SAVED/genome.transcripts.fa")

# study settings (3' tag trimming, no length correction) for every batch
STUDY_CONF=()
[[ -f "$PIPE/conf/$STUDY.config" ]] && STUDY_CONF=(-c "$PIPE/conf/$STUDY.config")

# each batch gets its own launch folder, so -resume only ever resumes that batch's run
mkdir -p "$SCRATCH/launch/$UNIT" && cd "$SCRATCH/launch/$UNIT"
nextflow run nf-core/rnaseq -r "$RNASEQ" \
  -profile hazel \
  -c "$PIPE/conf/maizexpressd.config" \
  ${STUDY_CONF[@]+"${STUDY_CONF[@]}"} \
  --input "$PIPE/manifests/$UNIT.samplesheet.csv" \
  --outdir "$OUT" \
  -work-dir "$SCRATCH/work/$UNIT" \
  -resume \
  "${REF_ARGS[@]}"

echo "done: $UNIT -> $OUT (FASTQ kept; remove with cleanup_study.sh after checking results)"
