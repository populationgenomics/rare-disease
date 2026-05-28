# GVCF Truncation Check

Utility scripts that scan GCS buckets for `.tbi` index files and verify they contain the expected number of contigs
(up to chrX). Files missing expected contigs are reported as truncated.

Requires `gcloud` authentication and read access to the target buckets as well as bcftools:

```bash
gcloud auth application-default login
```

## Scripts

### `check_truncated_gvcfs_from_dirs.sh`

Scans bucket directories listed in a paths file, finds all `.tbi` files, and checks each one.
Processes 16 files in parallel.
(This will make your laptop cry so feel free to tone it down if you don't need to do this fast)

Can be run locally or in a cloud shell environment.

```bash
./check_truncated_gvcfs_from_dirs.sh bucket_paths.txt > results.txt
```

### `check_gvcfs_from_files.sh`

Checks specific `.tbi` file paths listed one per line. Also reports the creation date and last contig for each file.

```bash
./check_gvcfs_from_files.sh file_paths.txt > results.txt
```

### `summarize_report.sh`

Reads the output of either check script and writes a summary
(total, OK, truncated counts and truncated file list) to `summary_report.txt`.

```bash
./summarize_report.sh results.txt
```
