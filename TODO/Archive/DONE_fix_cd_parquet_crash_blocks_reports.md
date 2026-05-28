# Fix CD Parquet crash blocking reports

## Problem

The CompoundDiscoverer DEA command can finish the classic R Markdown reports and
then exit on Windows with status `3221225477` before creating
`SummarizedExperiment.rds`, the Quarto report, and `index.html`.

The output listing shows the run stops before `lfqdata_normalized.parquet` and
`lfqdata.yaml`. In `inst/application/CMD_DEA_CD.R`, the next statement after the
classic reports is `arrow::write_parquet(...)`, followed by YAML,
`SummarizedExperiment`, Quarto, and index generation. Exit status `3221225477`
is a native Windows access violation, which `tryCatch()` cannot catch inside the
same R process.

## Plan

1. [x] Move the required report chain (`SummarizedExperiment.rds`, Quarto HTML,
   `index.html`) before the optional Parquet export in `CMD_DEA_CD.R`.
2. [x] Run the Arrow Parquet write in a separate `Rscript` process so a native
   Arrow crash cannot terminate the main DEA command.
3. [x] Keep YAML export in the main process because it is pure R and cheap.
4. [x] Log Parquet subprocess output and status so Windows failures are visible
   in the run log without hiding the root failing operation.
5. [x] Add focused tests for the isolated Parquet helper and run the relevant CD
   tests.

## Verification

- `Rscript -e "devtools::test(filter = '^CMD_DEA_CD$')"`
- `R CMD INSTALL .`
- `NOT_CRAN=true Rscript -e "testthat::test_file('tests/testthat/test-CMD_DEA_CD.R')"`
- `git diff --check`
