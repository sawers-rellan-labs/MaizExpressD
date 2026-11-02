"""Collect per-library salmon and fastp QC for every quantified batch into one table.
Run on Hazel:  python3 library_qc.py /rsstu/users/r/rrellan/sara/MaizeExpressD/results/Kremling2018 > library_qc.tsv
Columns: study_id (the batch here), sample_id, reads_processed, reads_mapped, percent_mapped, library_type,
compatible_ratio (as in MaizExpressC), then reads_raw (before fastp), adapter_trimmed (reads cut by the poly(A)
adapter), duplication, q30_raw, mean_length (after fastp), decoy_fragments and percent_decoy (salmon reads that
fit the genome better than any transcript, as % of reads_processed)."""
import glob
import json
import os
import sys

R = sys.argv[1]
print("\t".join(["study_id", "sample_id", "reads_processed", "reads_mapped", "percent_mapped", "library_type", "compatible_ratio",
                 "reads_raw", "adapter_trimmed", "duplication", "q30_raw", "mean_length", "decoy_fragments", "percent_decoy"]))
for study in sorted(os.listdir(R)):
    if study.startswith("_"):
        continue
    for mi in sorted(glob.glob(os.path.join(R, study, "*", "*", "aux_info", "meta_info.json"))):
        d = os.path.dirname(os.path.dirname(mi))
        if os.path.basename(os.path.dirname(d)) not in ("salmon", "star_salmon"):
            continue
        m = json.load(open(mi))
        lf = os.path.join(d, "lib_format_counts.json")
        l = json.load(open(lf)) if os.path.exists(lf) else {}
        fp = os.path.join(R, study, "fastp", os.path.basename(d) + ".fastp.json")
        f = json.load(open(fp)) if os.path.exists(fp) else {}
        before = f.get("summary", {}).get("before_filtering", {})
        after = f.get("summary", {}).get("after_filtering", {})
        n, decoy = m.get("num_processed", 0), m.get("num_decoy_fragments", 0)
        print("\t".join(map(str, [study, os.path.basename(d), m.get("num_processed", ""), m.get("num_mapped", ""),
                                  round(m.get("percent_mapped", 0), 2), l.get("expected_format", ""),
                                  round(l.get("compatible_fragment_ratio", 0), 3),
                                  before.get("total_reads", ""), f.get("adapter_cutting", {}).get("adapter_trimmed_reads", ""),
                                  round(f.get("duplication", {}).get("rate", 0), 4), round(before.get("q30_rate", 0), 4),
                                  after.get("read1_mean_length", ""), decoy, round(100 * decoy / n, 2) if n else ""])))
