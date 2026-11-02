#!/usr/bin/env bash
# Download the FASTQ files of one study from a manifest made by make_manifests.py.
#   bash download_fastq.sh manifests/<study>.downloads.tsv [parallel_jobs]
# Columns: run, url, md5, dest. Resumes partial files, checks every MD5, skips
# files already verified, and writes <dest_dir>/download.log. Re-run to retry
# failures. Exit status is non-zero if any file is missing or fails its MD5.
set -uo pipefail

MANIFEST="${1:?manifest .downloads.tsv}"
JOBS="${2:-6}"

fetch_one() {
  local run="$1" url="$2" md5="$3" dest="$4"
  local dir; dir=$(dirname "$dest")
  mkdir -p "$dir"
  if [[ -s "$dest" && -n "$md5" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
    echo -e "$(date +%F_%T)\t$run\tSKIP\t$(basename "$dest")" >> "$dir/download.log"
    return 0
  fi
  # No curl --retry: with -C - it can resume from a stale offset and duplicate bytes
  # (seen with curl 7.76 on Hazel). Each attempt below resumes from the real file size;
  # a stalled transfer (<50 KB/s for 3 min) is dropped and resumed. CNCB often closes
  # transfers after a few hundred MB, so a file may need many resumes: an unfinished
  # .part is kept for the next round. Only a file that curl finished and that fails
  # its MD5 is deleted and downloaded from scratch (3 rounds).
  local round try got done_ok
  for round in 1 2 3; do
    done_ok=0
    for try in $(seq 1 40); do
      curl -fsSL --connect-timeout 60 --speed-limit 50000 --speed-time 180 \
           -C - -o "$dest.part" "$url" && { done_ok=1; break; }
      sleep $(( try < 8 ? try * 15 : 120 ))
    done
    if [[ $done_ok -eq 0 ]]; then
      echo -e "$(date +%F_%T)\t$run\tINCOMPLETE(round $round)\t$(basename "$dest")\t$(stat -c %s "$dest.part" 2>/dev/null || echo 0) bytes" >> "$dir/download.log"
      continue
    fi
    got=$(md5sum "$dest.part" 2>/dev/null | cut -d' ' -f1)
    if [[ -n "$got" && ( -z "$md5" || "$got" == "$md5" ) ]]; then
      mv "$dest.part" "$dest"
      echo -e "$(date +%F_%T)\t$run\tOK\t$(basename "$dest")\t$got" >> "$dir/download.log"
      return 0
    fi
    echo -e "$(date +%F_%T)\t$run\tMD5_MISMATCH(round $round)\t$(basename "$dest")" >> "$dir/download.log"
    rm -f "$dest.part"
  done
  echo -e "$(date +%F_%T)\t$run\tFAILED\t$(basename "$dest")" >> "$dir/download.log"
  return 1
}
export -f fetch_one

tail -n +2 "$MANIFEST" | tr -d '\r' | \
  xargs -P "$JOBS" -L 1 bash -c 'fetch_one "$0" "$1" "$2" "$3"'

# summary
missing=0
while IFS=$'\t' read -r run url md5 dest; do
  [[ -s "$dest" ]] || { echo "missing: $dest"; missing=$((missing + 1)); }
done < <(tail -n +2 "$MANIFEST")
n=$(($(wc -l < "$MANIFEST") - 1))
echo "$(basename "$MANIFEST"): $((n - missing))/$n files present and verified"
[[ $missing -eq 0 ]]
