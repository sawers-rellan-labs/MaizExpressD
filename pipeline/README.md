# pipeline: download and quantification on Hazel (M2)

Same route as MaizExpressN and MaizExpressC: Slurm, **nf-core/rnaseq** 3.27.0 in salmon pseudo-alignment
mode with the nf-core `hazel` profile, plus `conf/maizexpressd.config` and the study settings in
`conf/Kremling2018.config`. FASTQ come either from ENA with our MD5-checked downloader or, much faster,
as .sra files from NCBI's AWS bucket converted on compute nodes (MaizExpressC's AWS route; see below). The B73 v5 reference bundle,
including the decoy-aware salmon index, is MaizExpressN's (restored by `run_study.sh` from
`sara/MaizeExpressN/reference`), so no index is built here.

Kremling2018 has 1,960 libraries (532 GB of FASTQ), so it runs as eight batches by tissue. Each batch is a
*unit* named `Kremling2018_<batch>`:

| Unit | Libraries | FASTQ |
|---|---|---|
| Kremling2018_GRoot | 291 | 64.6 GB |
| Kremling2018_GShoot | 295 | 74.7 GB |
| Kremling2018_L3Base | 302 | 84.3 GB |
| Kremling2018_L3Tip | 295 | 86.6 GB |
| Kremling2018_LMAD | 210 | 60.4 GB |
| Kremling2018_LMAN | 276 | 78.5 GB |
| Kremling2018_Kern | 254 | 73.5 GB |
| Kremling2018_LeafMid (L3Mid + LMid) | 37 | 9.2 GB |

## Locations (approved at the M1 checkpoint, 29 Oct 2026)

| What | Where |
|---|---|
| Pipeline files (copied from the project) | `/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline` (Mac: `/Volumes/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline`) |
| Final outputs (salmon, gene tables, MultiQC) | `/rsstu/users/r/rrellan/sara/MaizeExpressD/results/Kremling2018/<batch>` |
| FASTQ (temporary), Nextflow work, containers, reference, logs | `/share/maize/rrellan/MaizExpressD/{fastq,work,singularity,reference,logs}` |

## Files

| File | Purpose |
|---|---|
| `check_read_structure.py` | Read lengths, 5′/3′ base composition, poly(A)/poly(T), adapters, duplicates for a few FASTQ heads; results of the 29 Oct 2026 check in `read_structure/` |
| `make_manifests.py` | (Mac) Writes `manifests/<unit>.downloads.tsv` (ENA), `<unit>.sra.tsv` (AWS) and `<unit>.samplesheet.csv` from `metadata/runs.tsv` and `samples.tsv` |
| `download_fastq.sh` | ENA route: downloads one unit's FASTQ: resumable, parallel, MD5-checked, re-runnable; log in the FASTQ folder |
| `download_study.sbatch` | ENA route: runs `download_fastq.sh` as a job on the `xfer` node (`dl_<unit>`) |
| `download_sra.sh` | AWS route: downloads one unit's .sra files on the `xfer` node (`sra_<unit>`), 8 connections, size-checked; runs missing from AWS come from ENA instead |
| `convert_sra.sbatch` | AWS route: on a compute node (`cv_<unit>`), vdb-validate, fasterq-dump with ENA-style read names, read count and quality checks, pigz to the samplesheet paths; removes each .sra |
| `queue_sra.sh` | AWS route: queues `sra_<unit>` then `cv_<unit>` for one or more units |
| `speed_test.sh` | One-off comparison of ENA (HTTPS, FTP) and AWS download speeds; log in `speed_test.out` |
| `conf/maizexpressd.config` | nf-core/rnaseq parameters: reference, salmon only, fastp, resources (small salmon jobs for 3′ tags) |
| `conf/Kremling2018.config` | 3′ tag settings: poly(A) cut as an adapter, poly(G)/poly(X) trimming, no length correction; no fq lint or FastQC (see below) |
| `run_study.sh` | Nextflow head job for one unit (`sbatch --job-name=<unit> run_study.sh <unit>`) |
| `queue_studies.sh` | Queues several head jobs as chains, each waiting for its download (`dl_<unit>`) or conversion (`cv_<unit>`) job |
| `cleanup_study.sh` | Deletes a unit's FASTQ and work files after its gene table exists (asks first) |
| `library_qc.py`, `mapping_rule.py` | Per-library salmon and fastp QC table (`metadata/library_qc.tsv`, read by `metadata/build/library_exclusions.R`) and per-batch median mapping |

## Settings chosen from the read-structure check (29 Oct 2026)

First 200,000 reads of 11 runs covering every tissue, 11 flowcells, read lengths 85/86/90 and the pilot
batches (`read_structure/`):

- Reads are untrimmed, fixed length, with no UMI or spacer; no poly(T) starts. Salmon finds **stranded
  forward** (SF) reads, so the samplesheets say `forward`.
- 25–42% of reads run into a **poly(A) tail** (median start nt 44–51) followed by a degraded TruSeq adapter
  and G runs. fastp is given an 18-nt poly(A) as the adapter, which cuts the tail and everything after it,
  then trims 3′ poly(G) and other homopolymers.
- The first ~12 nt show random-primer composition bias, but removing them lowers salmon mapping by about 3
  points, so they are kept (MaizExpressC's Xue2021 removed them).
- Mapping of the test reads (transcript-only index, % of input reads): 35–53% untrimmed, 37–58% with the
  Xue2021 settings, **42–68% with these settings**. The decoy-aware index on Hazel may give slightly less.
- Only 3–4% of trimmed tags are shorter than salmon's k = 31, so the shared index (k = 31) is kept.
- Extending transcripts 500 bp downstream (for tags beyond annotated 3′ ends) added 5–7 points of mapping,
  but much of the gain went to a few genes on unplaced scaffolds and single loci, and some well-expressed
  genes lost most of their reads. Not used: the shared MaizExpressN annotation is kept.

## Download routes (speed test, 29 Oct 2026)

From the `xfer` node (`speed_test.sh`, log in `speed_test.out`), ENA gave ~0.7 MB/s per connection (13.5 MB/s
with 24 connections; FTP slower), while NCBI's open-data bucket on AWS (`sra-pub-run-odp`) gave 56 MB/s on one
connection and ~250 MB/s on eight. The .sra files total 403 GB (75% of the ENA FASTQ), so the whole set downloads
in well under an hour from AWS against ~11 h at ENA's best. One run converted with fasterq-dump has the same
reads and qualities as ENA's FASTQ (all 4.5 M records compared); `convert_sra.sbatch` also writes ENA's read
names (`@<run>.<spot> <spot>/1`, bare `+`), so the files match ENA's content too. ENA's MD5s cannot be checked on
converted files: instead vdb-validate checks the SRA object's own checksums and the read count must equal ENA's.

6 of the 1,960 runs are not in the AWS bucket (404: SRR5911565 GRoot, SRR5911272 GShoot, SRR5910075 and
SRR5912038 Kern, SRR5911633 L3Tip, SRR5911302 LMAN); `download_sra.sh` fetches those few from ENA.

LeafMid ran on the ENA route; the seven tissue units use the AWS route.

## Fewer Slurm jobs per library (29 Oct 2026)

Every Slurm job on Hazel runs a task epilog (the "Resource usage summary" at the end of each log) that calls
`scontrol show job` 13 times as the user, and the controller allows each user 6 requests per second. With
~6 short jobs per library our account hit that limit, which capped all our pipelines together at ~12 finished
jobs a minute and stalled Nextflow's `sbatch` calls for minutes. `conf/Kremling2018.config` therefore skips
fq lint and FastQC (`skip_linting`, `skip_fastqc`): each library now runs only fastp and salmon. fastp still
reports read quality, adapter content and duplication; the FASTQ were checked on download; library QC uses
salmon's output. LeafMid ran before this change and has lint and FastQC outputs too.

## Running (on a Hazel login node)

```bash
# Mac: copy pipeline/ to sara first (rsync -a pipeline/ /Volumes/rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline/)
cd /rsstu/users/r/rrellan/sara/MaizeExpressD/pipeline
# AWS route: download (xfer) -> conversion (compute) -> head job, per unit
bash queue_sra.sh Kremling2018_GRoot Kremling2018_GShoot
bash queue_studies.sh 2 Kremling2018_GRoot Kremling2018_GShoot     # each head job waits for its cv_<unit>
# ENA route (used for LeafMid)
mkdir -p /share/maize/rrellan/MaizExpressD/logs           # Slurm needs the log folder to exist
sbatch --job-name=dl_Kremling2018_LeafMid download_study.sbatch Kremling2018_LeafMid
bash queue_studies.sh 1 Kremling2018_LeafMid                # head job waits for the download
```

The small LeafMid unit (37 libraries, 9 GB) ran first as the check (29 Oct 2026, 34 min): 33–54% of reads
adapter-trimmed by fastp (median 39%), salmon `SF`, 47–62% of the trimmed reads mapped (1.3–4.3 M mapped reads).
Clean up each unit's FASTQ once its gene table is in `results/` (`bash cleanup_study.sh <unit>`); the full
set of FASTQ would be ~530 GB on scratch, which has room (100 TB free), so cleanup need not wait on the queue.
