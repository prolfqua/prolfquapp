# Preserve child-count assays in SummarizedExperiment reports

## Problem

The tabbed SummarizedExperiment report reconstructs `LFQData` from `rawData` and
`transformedData` assays. `make_SummarizedExperiment()` currently does not store
the configured `nr_children` matrix, so `se_report_lfqdata()` falls back to one
child per protein/sample. As a result, report sections requiring at least two
child measurements show an empty placeholder even when peptide-count data were
available in the original `LFQData`.

## Plan

1. Store the raw LFQData child-count wide matrix in the generated
   `SummarizedExperiment` as a dedicated assay.
2. Reconstruct the configured `nr_children` column from that assay in
   `se_report_lfqdata()` before creating the `LFQData` object.
3. Keep the fallback to one child only for older SE files that do not contain
   the child-count assay.
4. Use experiment-level peptide support from rowData-derived feature annotation
   for the report section that describes features supported by at least two
   peptides in the experiment.
5. Keep the per-sample child-count assay as a separate round-trip feature for
   reconstructed `LFQData`.
6. Add tests that the generated SE includes the assay and that reconstructed
   report data exposes experiment-level peptide support.
7. Reinstall `prolfquapp` and rerun the SAINT reports.
