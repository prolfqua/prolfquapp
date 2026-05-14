# TODO: Superseded - Add Contrast Agreement UpSet Tabs to Quarto Results

Status: superseded by `TODO_restructure_quarto_dea_report.md`, which places contrast agreement under
`Differential Expression` rather than under the current top-level `Results` tab.

## Context

The Quarto DEA report source is `inst/templates/quarto/Grp2Analysis_V2_SE.qmd`.

The current `Results` tab contains two DT tables:

- Number of samples per condition
- Significant proteins by contrast and direction

The report already includes an UpSet-style significant protein overlap plot later under `Differential Expression >
Significant Protein Overlap`, but the user expects the contrast-agreement visualization directly on the Results page,
next to the summary tables. Because the table summary is split into up/down regulated proteins, the overlap view should
also support all significant proteins, up-regulated proteins, and down-regulated proteins.

## Plan

1. Update only the tracked Quarto source template `inst/templates/quarto/Grp2Analysis_V2_SE.qmd`.
2. Restructure the `Results` section into a nested tabset:
   - `Tables`: keep the existing sample-count and significant-summary DT tables.
   - `Contrast Overlap`: add UpSet plots for significant proteins shared across contrasts.
3. Build the UpSet input directly from existing `significant_contrasts`, using `protein_Id`, `contrast`, and `diff`.
   Keep the implementation local to the QMD template and do not add public package API.
4. Render three overlap views where possible:
   - all significant proteins
   - up-regulated proteins (`diff > 0`)
   - down-regulated proteins (`diff < 0`)
   Each view should show a compact explanatory placeholder when fewer than two non-empty contrast sets are available.
5. Remove the existing duplicate `Differential Expression > Significant Protein Overlap` section or replace it with a
   pointer-style sentence if needed, so the report does not show the same overlap plot twice.
6. Verify with the existing Quarto report test path if practical:
   - `Rscript -e "testthat::test_file('tests/testthat/test-quarto-se-report.R')"`
   If full rendering is blocked by local dependencies or Quarto runtime, at minimum run a focused R check for the helper
   data transformation and report the limitation.
