# Handoff: WU347825 `prolfqua_qc` failure — fix ready, needs validation

**Date:** 2026-06-26
**Status:** Fix written + regression test added in `/scratch/wolski/prolfquapp`.
NOT validated — this host (fgcz-r-038) has no R packages and the docker image is a
different arch (`Rscript: cannot execute binary file`). Validate on Mac or via image rebuild.

## Symptom
`A386_DIANN/WU347825` failed at the `prolfqua_qc` Snakemake rule (exit 1).
- `.snakemake/log/2026-06-25T181636.707691.snakemake.log` → `Error in rule prolfqua_qc`.
- All `work/logs/*.log` are 0 bytes; `Rqc.1.log` was deleted by Snakemake on failure.
- Reproduced manually with `prolfquapp-docker --image prolfqua/prolfquapp:2.2.6 -- prolfqua_qc.sh ...`
  (binary: `/scratch/wolski/diann-runner-wolski/.venv/bin/prolfquapp-docker`).

Real error captured from the rerun:
```
INFO  Files data: out-DIANN_quantB/WU347825_report-first-pass.site_report.parquet; out-DIANN_quantB/WU347825_report.parquet; out-DIANN_quantB/WU347825_report.site_report.parquet
ERROR the condition has length > 1
```

## Root cause
`get_DIANN_files()` in `R/preprocess_DIANN.R` greps DIA-NN reports with the **unanchored**
pattern `report\.parquet$`. The new DIA-NN run emits PTM site-report files
(`*.site_report.parquet`) which also END in `report.parquet`, so 3 files matched instead of 1.
The length-3 vector flows into `read_diann_report()`, whose `if (grepl("\\.parquet$", path))`
errors under R >= 4.2 ("the condition has length > 1").

## Fix (already applied)
`R/preprocess_DIANN.R`, inside `get_DIANN_files()`, right after the initial `grep`:
```r
# DIA-NN PTM "site_report.parquet" files also end in "report.parquet" and
# would otherwise be matched above; drop them so only the main report remains.
diann.path <- diann.path[!grepl("site_report\\.(parquet|tsv)$", diann.path)]
```
(`report-first-pass.parquet` does NOT match `report\.parquet$`, so it was never a problem.)

Verified the corrected logic against the real dir on this host (plain R, no pkgs needed):
resolves the 3 matches down to exactly `WU347825_report.parquet`.

## Regression test (already added)
`tests/testthat/test-preprocess_DIANN-native2x.R` →
`test_that("get_DIANN_files ignores DIA-NN PTM site_report files", ...)`.
Creates `WU1_report.parquet` + the two `*.site_report.parquet` variants and asserts
`get_DIANN_files()` returns exactly the main report.

## To validate (Mac or rebuild)
```bash
cd /scratch/wolski/prolfquapp      # or your Mac clone with the same diff
Rscript -e 'devtools::test(filter = "preprocess_DIANN")'   # expect new test PASS
# full suite:
Rscript -e 'devtools::test()'
```
Then rebuild the image and re-run the real step:
```bash
docker build -t prolfqua/prolfquapp:2.2.6 .   # bakes in the fixed R source
cd /scratch/A386_DIANN/WU347825/work
PATH=/scratch/wolski/diann-runner-wolski/.venv/bin:$PATH \
  prolfquapp-docker --runtime docker --image prolfqua/prolfquapp:2.2.6 -- \
  prolfqua_qc.sh --indir out-DIANN_quantB -s DIANN --dataset dataset.csv \
  --project 37485 --order 37485 --workunit 347825 --outdir qc_result --flat_outdir
# expect exit 0 and a populated qc_result/
```

## Notes / open question
- There's also a benign warning: `there are more than one column for sample: Relative Path, File`
  (annotation has both a `Relative Path` and a `File` column). Not the cause of the failure, but
  worth a look later if QC annotation behaves oddly.
- Only one file changed: `R/preprocess_DIANN.R` (+ the test). `git diff` in the repo shows it.
