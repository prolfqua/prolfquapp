# TODO: option for static (vs dynamic dated) output folder naming

Created: 2026-06-15

## Problem

QC (and DEA) output is written into a **dynamically named** folder built by
[`zipdir_name()`](../R/R6_AppConfiguration.R) (`R/R6_AppConfiguration.R:103-145`):

```
<prefix>_<YYYYMMDD>_PI<project>_O<order>_WU<workunit>_<transform>
# e.g. QC_20260509_PI37485_O37485_WU345497_none
```

This subdir is created **even when the caller passes an explicit `--outdir`**.
So `prolfqua_qc.sh --outdir qc_result …` does not produce
`qc_result/proteinAbundances.html`; it produces
`qc_result/QC_20260509_PI37485_O37485_WU345497_none/proteinAbundances.html`.

The folder name embeds `Sys.Date()` — the **run date on the compute node**, not
the submit date — so the path is non-deterministic and **cannot be reconstructed
by a caller ahead of time** (a queued job can run on a different day than it was
submitted).

## Why it matters (downstream breakage)

Consumers that need a **static** path break:

- **SUSHI** — `sushi/master/lib/DIANNApp.rb#next_dataset` registers
  `'Protein Abundances [Link]' => qc_result/proteinAbundances.html` and
  `'Sample Sizes [Link]' => qc_result/QC_sampleSizeEstimation.html`. A static
  `[Link]` can't point into a dated, dynamic subdir.
- **diann_runner** — the Snakefile QC step calls
  `prolfqua_qc.sh … --outdir qc_result` and expects the htmls directly under
  `qc_result/`.
- **legacy FGCZ Makefile** — works around this with
  `find qc_result -name proteinAbundances.html`, which only papers over it.

Observed inconsistency in real output (same pipeline, different prolfquapp runs):

| run | layout |
|-----|--------|
| 2026‑02‑05 | **flat**: `qc_result/proteinAbundances.html` |
| 2026‑05‑09 | **nested**: `qc_result/QC_20260509_PI37485_O37485_WU345497_none/proteinAbundances.html` |

The nested layout has **no** top-level html at all (not even `index.html`).

## Root cause

- `R/R6_AppConfiguration.R::zipdir_name()` — builds the dated `<prefix>_…` name.
- `R/ProjectStructure.R` — `outpath` / `qc_dir` (default `"qc_results"`) /
  `qc_path()` (`file.path(self$outpath, self$qc_dir)`).
- QC entry: `inst/application/CMD_QUANT_QC.R` (via
  `inst/application/bin/prolfqua_qc.sh`), which forwards `--outdir`.

## Proposed fix (opt-in, non-breaking)

Add a flag that makes QC/DEA write outputs **directly into `--outdir`** with **no
dynamic subdir** — e.g. `--flat-outdir` (boolean) or `--outdir-mode static|dynamic`.

- **Default = current (dynamic)** behavior, so multi-run namespacing in shared
  output dirs is preserved for existing users.
- When set, skip the `zipdir_name()` subdir: write `proteinAbundances.html`,
  `QC_sampleSizeEstimation.html`, `index.html`, … directly under `--outdir`.
- Wire the flag through `CMD_QUANT_QC.R` → the QC generator / `ProjectStructure`
  so `qc_path()` returns `--outdir` itself when static.

Then diann_runner calls `prolfqua_qc.sh --outdir qc_result --flat-outdir …` and
the static `qc_result/proteinAbundances.html` contract holds for SUSHI + AppRunner.

## Acceptance criteria

- `prolfqua_qc.sh --outdir qc_result --flat-outdir …` produces, directly:
  - `qc_result/proteinAbundances.html`
  - `qc_result/QC_sampleSizeEstimation.html`
  - `qc_result/index.html`
- Without the flag, output layout is unchanged (dated subdir as today).
- Same option should apply to the DEA entry for consistency (DEA also uses
  `zipdir_name()`), even if only QC is needed first.

## Implementation (done 2026-06-15)

Implemented as opt-in `--flat_outdir` boolean (snake_case to match existing
flags like `--pattern_decoys`; maps cleanly to `opt$flat_outdir`).

Single root-cause switch in `ProlfquAppConfig$get_zipdir()`: when
`flat_outdir` is TRUE it returns `path` directly, skipping the dated subdir.
Everything downstream (QC `output_dir`, DEA `ZIPDIR` / `get_result_dir()` /
`get_input_dir()`) derives from `get_zipdir()`, so no call sites needed
touching.

Changes:

- `R/R6_AppConfiguration.R` — new `flat_outdir = FALSE` field + `initialize`
  param; `get_zipdir()` honors it; `list_to_R6_app_config()` round-trips it.
- `R/utils.R` — `sync_opt_config()` applies `opt$flat_outdir` (DEA path).
- `R/cmd_helpers.R` — `run_qc_preprocess()` gained a `flat_outdir` param.
- `inst/application/CMD_QUANT_QC.R` — `--flat_outdir` option, passed through.
- `inst/application/CMD_DEA_V2.R` — `--flat_outdir` option; `dir.create`
  calls made idempotent (`showWarnings = FALSE, recursive = TRUE`).
- `tests/testthat/test-flat_outdir.R` — new unit tests (config switch,
  list round-trip, `sync_opt_config` override, `run_qc_preprocess`).

Shell scripts unchanged: they forward `"$@"`, so the new flag flows through.

Default behavior (no flag) is unchanged: dated subdir as today.

## Notes

- Caller-side reconstruction of the dynamic name (e.g. in DIANNApp.rb) was
  considered and rejected: the name embeds the run-time `Sys.Date()`, which the
  caller cannot know at submit time. Fixing it here (the producer) is the correct
  layer.
- If the dated subdir is intentional for namespacing concurrent runs into one
  output dir, that's exactly why the static mode must be **opt-in**.
