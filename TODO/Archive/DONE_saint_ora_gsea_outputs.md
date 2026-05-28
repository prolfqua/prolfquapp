# TODO: Enable ORA and RNK Outputs for SAINT Reports

## Problem

`DEAReportGenerator$write_DEA()` currently skips ORA and GSEA/rnk generation when `deanalyse$default_model == "saint"`. The SAINT workflow needs these files for downstream enrichment tools.

## Plan

1. Inspect the SAINT contrast table columns produced by `ContrastsSAINTexpress` in the current prolfquapp output.
2. Update `DEAReportGenerator` so ORA foreground/background files and rank files can be written for SAINT.
3. Keep regular DEA behavior unchanged.
4. For SAINT rank files, use SAINT-compatible ranking columns rather than the DEA-only `statistic` column. Prefer `log2_EFCs`, then `SaintScore`, then `BFDR` if needed.
5. Add focused tests for SAINT ORA/RNK generation and regular DEA rank behavior.
6. Regenerate documentation if roxygen comments change, run targeted tests, reinstall via package-local `make install`, and rerun one SAINT example to verify files appear.
