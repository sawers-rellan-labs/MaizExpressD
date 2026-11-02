#!/usr/bin/env bash
# Queue quantification (run_study.sh head jobs) for several batches as parallel chains.
# A batch waits for its download job (dl_<unit>, ENA) or conversion job (cv_<unit>, AWS route,
# queued by queue_sra.sh) if one is still queued or running,
# and for the previous batch in its chain, so at most <chains> head jobs hold xfer slots at once.
# If a download or conversion fails, its head job is cancelled (--kill-on-invalid-dep) and
# the chain moves on; requeue that batch after re-downloading it.
#   bash queue_studies.sh <chains> <unit> [<unit> ...]     e.g. 1 Kremling2018_GRoot Kremling2018_GShoot
set -euo pipefail
CHAINS="${1:?number of parallel chains}"; shift
cd "$(dirname "$0")"

declare -a last
i=0
for s in "$@"; do
  dl=$(squeue -h -u "$USER" -n "dl_$s,cv_$s" -o %i | paste -sd:)
  if [[ -z "$dl" ]]; then
    # no download pending: every FASTQ in the sample sheet must be there
    n_missing=$(tail -n +2 "manifests/$s.samplesheet.csv" | cut -d, -f2,3 | tr ',' '\n' |
                while read -r f; do [[ -z "$f" || -s "$f" ]] || echo "$f"; done | wc -l)
    if [[ $n_missing -gt 0 ]]; then
      echo -e "$s\tSKIPPED\t$n_missing FASTQ missing and no download or conversion queued"; continue
    fi
  fi
  c=$((i % CHAINS)); i=$((i + 1))
  dep=()
  [[ -n "$dl" ]] && dep+=("afterok:$dl")
  [[ -n "${last[$c]:-}" ]] && dep+=("afterany:${last[$c]}")
  opt=(--kill-on-invalid-dep=yes)
  [[ ${#dep[@]} -gt 0 ]] && opt+=(--dependency="$(IFS=,; echo "${dep[*]}")")
  jid=$(sbatch --parsable --job-name="$s" "${opt[@]}" run_study.sh "$s")
  last[$c]=$jid
  echo -e "$s\t$jid\tchain $c\t${dep[*]:-no dependency}"
done
