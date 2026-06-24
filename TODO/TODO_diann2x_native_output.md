# TODO: Read native DIA-NN 2.x output (parquet / `Run` column)

## Goal

Make the prolfquapp DIA-NN reader consume **default DIA-NN 2.x output** directly, so
upstream pipelines no longer have to rename `Run` → `File.Name` or down-convert the
parquet to a prolfqua-specific TSV.

This is the prolfquapp counterpart to the `diann_runner` change in
`gstore/diann_runner/TODO/TODO_drop_prolfqua_tsv_conversion.md`. The two must be
coordinated.

## Background (why this exists)

`gstore/diann_runner` currently runs DIA-NN 2.x, whose main report is a **parquet** with:

- a bare `Run` column (e.g. `20260623_010_C42222_S1172811_Plate_7001_H1` — basename, no
  path, no extension),
- a `Run.Index` column,
- **no `File.Name` column** (DIA-NN 1.x had `File.Name`; 2.x dropped it from the main report).

prolfquapp was written for DIA-NN 1.x and depends on `File.Name`. To bridge this, the
runner converts the parquet to TSV and renames `Run` → `File.Name`
(`diann_runner/src/diann_runner/snakemake_helpers.py:525`). That renamed TSV is the only
reason prolfquapp works today — and it breaks other consumers (pmultiqc needs native
`Run`). Teaching prolfquapp the native schema lets the runner drop the shim.

## Verified Current State (in this repo)

All in `prolfquapp/R/preprocess_DIANN.R`:

- **`diann_read_output()` (lines 40-44)** derives `raw.file` from `File.Name`:
  ```r
  report2$raw.file <- gsub("^x|\\.d\\.zip$|\\.d$|\\.raw$|\\.mzML$", "",
                           basename(gsub("\\\\", "/", report2$File.Name)))
  ```
  → hard dependency on `File.Name`.
- **`preprocess_DIANN()` (line 162)** and **`dataset_template_diann()` (line 253)** read
  via `readr::read_tsv(quant_data)` → TSV only, no parquet path.
- **`get_DIANN_files()` (lines 93-98)** discovers the report by grepping
  `report\\.tsv$|diann-output\\.tsv` → never finds a parquet.
- **`arrow` is already a dependency** (`prolfquapp/DESCRIPTION:26`) and is already used to
  *write* parquet (`R/cmd_helpers.R:451`, `R/parquet_helpers.R`). Reading parquet adds no
  new dependency.

Other columns the reader uses are native DIA-NN names and unaffected: `Protein.Group`,
`Stripped.Sequence`, `Precursor.Quantity`, `Precursor.Normalised`, `PEP`,
`Lib.PG.Q.Value`, `PG.Q.Value`, `Protein.Names`, `PG.Quantity`/`PG.MaxLFQ`.

## Proposed Changes

### 1. `diann_read_output()` — derive `raw.file` from `Run`, fall back to `File.Name`

```r
run_col <- if ("Run" %in% names(report2)) "Run"
           else if ("File.Name" %in% names(report2)) "File.Name"
           else stop("DIA-NN report has neither 'Run' nor 'File.Name'")
report2$raw.file <- gsub("^x|\\.d\\.zip$|\\.d$|\\.raw$|\\.mzML$", "",
                         basename(gsub("\\\\", "/", report2[[run_col]])))
```

Safe because the DIA-NN 2.x `Run` value is already the bare basename: `basename()`,
the backslash `gsub`, and the extension strip are all no-ops on it, so `raw.file` is
**identical** to the value produced today from the renamed `File.Name`. (For DIA-NN 1.x
`File.Name`, which is a full path with extension, the existing logic still applies via
the fallback.)

### 2. `get_DIANN_files()` — discover the parquet (prefer it)

Extend the grep to also match `report\\.parquet$` and prefer the parquet when both are
present:

```r
diann.path <- grep("report\\.parquet$|report\\.tsv$|diann-output\\.tsv",
                   dir(path = path, recursive = TRUE, full.names = TRUE), value = TRUE)
# prefer parquet if available
```

### 3. `preprocess_DIANN()` / `dataset_template_diann()` — read parquet or TSV

Replace the bare `readr::read_tsv(quant_data)` with a small dispatch:

```r
read_diann <- function(path) {
  if (grepl("\\.parquet$", path)) arrow::read_parquet(path) else readr::read_tsv(path)
}
```

`arrow::read_parquet` returns a tibble/data.frame; the rest of the pipeline is unchanged.

## Tests

- Add a small DIA-NN 2.x **parquet** fixture (native `Run`/`Run.Index`, no `File.Name`)
  under `prolfquapp/tests/...` (or `inst/application/DIANN/...`).
- `diann_read_output()` produces the same `raw.file` values from the parquet (`Run`) as
  from a legacy `File.Name` TSV.
- `get_DIANN_files()` discovers and prefers the parquet.
- `preprocess_DIANN()` builds an `LFQData` from the parquet fixture.
- Keep a legacy `File.Name` TSV test to confirm the fallback path still works.
- See existing `prolfquapp/tests/testthat/test-preprocess_DIANN-empty.R` for the
  empty-report path.

## Coordination

- After this lands and is released (image `prolfqua/prolfquapp:<new>`), `diann_runner`
  can drop `convert_parquet_to_tsv` and point `prolfqua_qc` at the parquet
  (`gstore/diann_runner/TODO/TODO_drop_prolfqua_tsv_conversion.md`).
- The DIA-NN-runner `diann-qc` plotter already reads parquet natively, so only prolfquapp
  blocks the cutover.
- Bump the prolfquapp version / Docker image tag so the runner can pin the version that
  understands native output.

## Acceptance Criteria

- prolfquapp reads a default DIA-NN 2.x report (parquet with `Run`, no `File.Name`)
  without any upstream renaming.
- `raw.file` values are identical to the current renamed-TSV path (regression-tested).
- The legacy DIA-NN 1.x `File.Name` TSV path still works (fallback).
- No new package dependency (uses the existing `arrow`).
