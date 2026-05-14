# TODO: SAINT as a prolfquapp model path

Implement the intermediate package graph:

```text
prolfquapp -> prolfquasaint -> prolfqua
prolfquapp -> prolfqua
```

Planned changes:

- Remove `prolfquapp` from `prolfquasaint` runtime imports.
- Keep SAINT core execution and result adapter code in `prolfquasaint`.
- Move parser/report orchestration responsibility to `prolfquapp`.
- Add `model = "saint"` support to the existing `prolfquapp` DEA command path.
- Store SAINT input/output artifacts on the `DEAnalyse` object.
- Make report/XLSX output write SAINT sheets and skip DEA-only ORA/GSEA exports.
- Add focused tests for dependency direction, SAINT dispatch, and SAINT report outputs.

Follow-up:

- Adapt the checked-in `prolfquasaint` WU337670 example to the current
  `prolfquapp` YAML format with `processing_options$model: saint`.
- Update the root Snakemake package install graph to match
  `prolfquapp -> prolfquasaint -> prolfqua` and keep
  `prolfquabenchmark` last in the R package list.

Volcano/report follow-up:

- Root cause: the legacy DEA report hard-codes volcano colors for DEA model
  names. SAINT produces `modelName = "ContrastSaint"`, so the fallback color is
  assigned to a palette entry whose name remains `NA`. Plotly then treats
  `ContrastSaint` as outside the color scale and renders an empty-looking
  volcano.
- Fix the report palette construction in the source vignette so fallback colors
  are named with the actual model levels.
- Keep the fix model-agnostic; do not special-case only `ContrastSaint`.
- Add a focused test for the palette helper or extracted report logic so unknown
  model names retain named colors.

Tabbed/SAINT report follow-up:

- Root cause of tiny headings in the legacy R Markdown report is the shared
  `inst/templates/fgcz.css` rule `h1 { font-size: 8px; }`. Replace it with
  sensible heading sizes that preserve large report titles and readable section
  headings.
- For `model = "saint"`, skip the legacy differential-expression introduction
  because the linear-model explanation is misleading for SAINTexpress.
- Add the existing SummarizedExperiment-backed Quarto tabbed report to the
  normal `CMD_DEA_V2.R` output path, mirroring the existing CompoundDiscoverer
  command behavior: write `SummarizedExperiment.rds`, render
  `DE_WU*_quarto.html` when Quarto is available, and include it in `index.html`.
- Fix `render_quarto_se_report()` to normalize the `SummarizedExperiment.rds`
  path before entering the temporary Quarto render directory, so command-line
  runs with relative output paths work.
- Fix index link generation for Quarto output after `render_quarto_se_report()`
  returns an absolute file path while the ZIP/result directory is relative.
- Keep the R Markdown report as the legacy output for now; add the tabbed report
  as an additional artifact rather than replacing report generation in this
  phase.

SummarizedExperiment contrast padding follow-up:

- Root cause: `DEAReportGenerator$make_SummarizedExperiment()` stores each contrast table aligned to every assay row. Proteins without a result become all-`NA` rows in the nested rowData table.
- The Quarto SE extraction helper `.se_report_contrast_table()` currently returns those padded all-`NA` rows, so plots facet by a real contrast plus an empty `NA` contrast.
- Fix `.se_report_contrast_table()` to drop padded contrast rows where both the feature identifier and contrast label are missing after reconstructing each nested contrast table.
- Add a focused test covering a nested rowData contrast table with padded all-`NA` rows.
- Rerender the WU345302 SAINT example and verify the report has no empty `NA` contrast facet.

R coding skill compliance cleanup:

- Remove `library()` attachment calls from the SE Quarto template and use `requireNamespace()` checks plus explicit package namespaces.
- Replace new broad boolean test expectations for padded contrast rows with specific count expectations.
- Re-run targeted SE Quarto report tests and reinstall/rerender if template behavior changed.
