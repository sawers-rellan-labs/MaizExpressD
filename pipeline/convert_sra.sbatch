#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/share/maize/rrellan/MaizExpressD/logs/convert_%x_%j.out
# Convert one batch's .sra files to gzipped FASTQ on a compute node (no internet needed).
# From MaizExpressC.
#   sbatch --job-name=cv_<unit> --dependency=afterok:<sra job> convert_sra.sbatch <unit>
# Per run: vdb-validate (SRA checksums), fasterq-dump with ENA-style read names
# (@<run>.<spot> <spot>/1, bare '+' line), record count must equal reads_per_file in
# manifests/<unit>.sra.tsv, qualities must not be flattened (SRA Lite), pigz, then the .sra
# is removed. The FASTQ go to the samplesheet paths, so run_study.sh is the same for both
# routes. 8 runs at a time, 2 threads each (runs are small: ~16 s of fasterq-dump on 1 CPU).
# Re-running skips runs whose FASTQ exist. Log: <fastq_dir>/convert.log.
set -uo pipefail
UNIT="${1:?unit, e.g. Kremling2018_GRoot}"
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline
SCR=/share/maize/rrellan/MaizExpressD
SIF=/share/maize/rrellan/MaizeExpressC/singularity/sra-tools_3.2.1.sif
export TMPDIR=$SCR/tmp/convert_$UNIT HOME=$SCR/tmp     # set before the module: it reads TMPDIR
mkdir -p "$TMPDIR"
source /etc/profile.d/modules.sh
module load apptainer >/dev/null 2>&1

one() {
  local run="$1" out1="$2" out2="${3#-}" n="$4"          # "-" = no second file
  local sra="$SCR/sra/$UNIT/$run.sra" dir; dir=$(dirname "$out1"); mkdir -p "$dir"
  local log="$dir/convert.log" t="$TMPDIR/$run"
  if [[ -s "$out1" && ( -z "$out2" || -s "$out2" ) ]]; then echo -e "$(date +%F_%T)\t$run\tSKIP" >> "$log"; return 0; fi
  [[ -s "$sra" ]] || { echo -e "$(date +%F_%T)\t$run\tMISSING_SRA" >> "$log"; return 1; }
  singularity exec -B /share -B /gpfs_common "$SIF" vdb-validate "$sra" > "$t.validate" 2>&1 \
    || { echo -e "$(date +%F_%T)\t$run\tVALIDATE_FAILED" >> "$log"; return 1; }
  mkdir -p "$t"
  singularity exec -B /share -B /gpfs_common "$SIF" fasterq-dump --split-files -e 2 -t "$t" -O "$t" \
      --seq-defline '@$ac.$si $si/$ri' --qual-defline '+' "$sra" > "$t.dump" 2>&1 \
    || { echo -e "$(date +%F_%T)\t$run\tDUMP_FAILED" >> "$log"; rm -rf "$t"; return 1; }
  local f c q ok=1
  for f in "$t/${run}_1.fastq" ${out2:+"$t/${run}_2.fastq"}; do
    [[ -s "$f" ]] || { ok=0; echo -e "$(date +%F_%T)\t$run\tNO_OUTPUT\t$(basename "$f"); got $(ls "$t" | tr '\n' ' ')" >> "$log"; continue; }
    c=$(( $(wc -l < "$f") / 4 ))
    [[ "$c" == "$n" ]] || { ok=0; echo -e "$(date +%F_%T)\t$run\tCOUNT_MISMATCH\t$(basename "$f") $c != $n" >> "$log"; }
    q=$(awk 'NR % 4 == 0' "$f" | head -2000 | fold -w1 | sort -u | wc -l)
    (( q >= 3 )) || { ok=0; echo -e "$(date +%F_%T)\t$run\tFLAT_QUALITIES\t$(basename "$f") $q distinct" >> "$log"; }
  done
  if [[ $ok -eq 1 ]]; then
    pigz -p 2 -c "$t/${run}_1.fastq" > "$out1.part" && mv "$out1.part" "$out1"
    [[ -n "$out2" ]] && pigz -p 2 -c "$t/${run}_2.fastq" > "$out2.part" && mv "$out2.part" "$out2"
    rm -rf "$t" "$t.validate" "$t.dump" "$sra"
    echo -e "$(date +%F_%T)\t$run\tOK\t$n reads per file" >> "$log"; return 0
  fi
  rm -rf "$t"; return 1
}
export -f one; export SCR UNIT SIF TMPDIR

t0=$(date +%s)
# xargs splits on whitespace, so an empty out2 (single-end runs) is passed as "-"
tail -n +2 "$PIPE/manifests/$UNIT.sra.tsv" | awk -F'\t' '{print $1, $3, ($4 == "" ? "-" : $4), $5}' |
  xargs -P 8 -L 1 bash -c 'one "$0" "$1" "$2" "$3"'
st=$?
log="$(dirname "$(awk -F'\t' 'NR==2{print $3}' "$PIPE/manifests/$UNIT.sra.tsv")")/convert.log"
echo "converted: $(awk -F'\t' '$3 == "OK" || $3 == "SKIP" {print $2}' "$log" | sort -u | wc -l) of $(( $(wc -l < "$PIPE/manifests/$UNIT.sra.tsv") - 1 )) runs in $(( ($(date +%s) - t0) / 60 )) min"
awk -F'\t' '$3 != "OK" && $3 != "SKIP"' "$log" | tail -20     # problems, if any
exit $st
