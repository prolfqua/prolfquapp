# Handoff: WU347825 `prolfqua_qc` failure — TWO fixes

**Date:** 2026-06-26
**Status:** TWO separate bugs hit in sequence at the `prolfqua_qc` step.

- **Bug 1 — site_report.parquet (R code).** FIXED in `R/preprocess_DIANN.R` + test.
  **CONFIRMED WORKING**: rebuilt as `prolfqua/prolfquapp:2.2.7`, rerun now shows
  `Files data: out-DIANN_quantB/WU347825_report.parquet` (single file) — no more
  "the condition has length > 1".
- **Bug 2 — arrow zstd codec (Docker image).** Exposed once Bug 1 was past. The DIA-NN 2.x
  `report.parquet` is zstd-compressed but the image's `arrow` was a minimal build:
  `NotImplemented: Support for codec 'zstd' not built`. FIXED in `Dockerfile` (install
  arrow with `LIBARROW_MINIMAL=false` + zstd, with a build-time `stopifnot` capability
  check). NEEDS A REBUILD (2.2.7 still has minimal arrow) → next tag (2.2.8).

This host (fgcz-r-038) has no R packages and the image is a different arch, so the image
must be rebuilt (Mac / CI) to validate Bug 2.

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

## Bug 1 — root cause
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

## Bug 2 — root cause
DIA-NN 2.x writes `report.parquet` with **zstd** compression. The image's `arrow` R package
was a *minimal* build (no compression codecs) because arrow was pulled in implicitly as a
prolfquapp dependency with no feature flags. Reading the parquet fails with:
```
NotImplemented: Support for codec 'zstd' not built
```

## Bug 2 — fix (already applied, NEEDS REBUILD)
`Dockerfile`, in the `build` stage right after `mkdir /opt/r-libs-site` (before pak / dep
resolution):
```dockerfile
ENV LIBARROW_MINIMAL=false
ENV ARROW_WITH_ZSTD=ON
RUN R -e 'options(warn=2); install.packages("arrow", repos = "https://stat.ethz.ch/CRAN/"); stopifnot(arrow::arrow_info()$capabilities[["zstd"]])'
```
Installed first so `pak::local_install_deps(upgrade = FALSE)` keeps this full build. The
`stopifnot` fails the build if zstd is still missing.

## To validate (Mac or rebuild)
```bash
cd /scratch/wolski/prolfquapp      # or your Mac clone with the same diff
Rscript -e 'devtools::test(filter = "preprocess_DIANN")'   # expect new test PASS (Bug 1)
Rscript -e 'devtools::test()'                              # full suite
```
Then rebuild the image (bakes in BOTH fixes) and re-run the real step:
```bash
docker build -t prolfqua/prolfquapp:2.2.8 .   # build fails fast if arrow lacks zstd
cd /scratch/A386_DIANN/WU347825/work
PATH=/scratch/wolski/diann-runner-wolski/.venv/bin:$PATH \
  prolfquapp-docker --runtime docker --image prolfqua/prolfquapp:2.2.8 -- \
  prolfqua_qc.sh --indir out-DIANN_quantB -s DIANN --dataset dataset.csv \
  --project 37485 --order 37485 --workunit 347825 --outdir qc_result --flat_outdir
# expect exit 0 and a populated qc_result/
```
NOTE: bump the image tag the Snakefile requests (currently `prolfqua/prolfquapp:2.2.7`)
to whatever you build, or the runner keeps using the old minimal-arrow image.

## Notes / open question
- Benign warning: `there are more than one column for sample: Relative Path, File`
  (annotation has both a `Relative Path` and a `File` column). Not the cause of either failure;
  worth a later look only if QC annotation misbehaves.
- Files changed: `R/preprocess_DIANN.R` + its test (Bug 1), and `Dockerfile` (Bug 2).
