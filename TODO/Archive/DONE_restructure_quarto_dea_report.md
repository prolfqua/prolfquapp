# TODO: Restructure SummarizedExperiment Quarto DEA Report

Status: implemented in `inst/templates/quarto/Grp2Analysis_V2_SE.qmd`.

## Context

The active Quarto source is `inst/templates/quarto/Grp2Analysis_V2_SE.qmd`. It is rendered by
`prolfquapp:::render_quarto_se_report()` from a `SummarizedExperiment` written by
`DEAReportGenerator$make_SummarizedExperiment()`.

There is currently only one source QMD in `inst/templates/quarto/`. The `DiffExpQC_SE.html` file in the same directory is
generated/static and should not be treated as source unless a matching QMD is added later.

The current top-level tabs are:

1. `Settings`
2. `Results`
3. `Protein Identification`
4. `Quality Control`
5. `Differential Expression`
6. `Tables`
7. `Session Info`

Vocabulary problem:

- The report is used for proteins, peptides, and metabolites. User-facing section names should therefore avoid
  protein-only terms where the section is conceptually generic.
- Use `feature` for the measured row entity in generic sections.
- Use `differential abundance` or `differential analysis` rather than `differential expression` in generic report
  navigation. `Differential expression` can remain only in method text if a proteomics-specific context explicitly
  requires it.

Problems with the current structure:

- `Results` appears before the reader has seen identification and quality-control context.
- The sample grouping table is under `Results`, but it is experimental design/settings context.
- Differential-abundance result summaries are split across `Results`, `Differential Expression`, and `Tables`.
- The significant-feature UpSet plot is currently under `Differential Expression > Significant Protein Overlap`, while
  the related increased/decreased summary table is under `Results`.
- `Protein Identification` contains both feature-count summaries and sample annotation; sample annotation belongs closer to
  settings/design.
- The current layout gives equal top-level weight to `Tables`, although this report currently contains one full result
  table. The navigation label should be singular and specific.

## Design Goal

Make the report follow the analysis workflow and reader expectations:

1. What was run and on which samples?
2. What was identified and how complete is the data?
3. Are abundance distributions and sample structure acceptable?
4. What differential-abundance results were found?
5. Where is the full machine-readable result table?

## Proposed Top-Level Structure

### 1. Settings

Purpose: all run metadata and experimental design inputs.

Keep:

- B-Fabric links
- FDR threshold
- Difference threshold
- Model formula

Move here from current `Results` / `Protein Identification`:

- `tbl-sample-counts` as `Sample counts per condition`
- `tbl-sample-annotation` as `Sample annotation`

Add if available from metadata:

- Contrast definitions from `meta$contrasts` / `meta$formula`, shown as a compact table.

Suggested subtabs:

- `Run`
- `Design`
- `Contrasts`

### 2. Feature Detection

Purpose: what was measured before any statistical testing.

Keep:

- `fig-protein-counts`
- `fig-protein-counts-two-peptides`

Move here from current `Quality Control > Missing Values`:

- `fig-protein-detection-overlap`, because it describes feature detection overlap across experimental groups, not
  differential expression.

Suggested subtabs:

- `Counts`
- `Overlap`

Notes:

- Preserve compact placeholder behavior for unavailable two-peptide plots.
- If contaminant/decoy summary is available in future `SummarizedExperiment` metadata, this section is the correct place
  to add it, but keep the top-level tab label generic.

### 3. Quality Control

Purpose: evaluate whether the abundance data and sample structure are suitable for differential-abundance analysis.

Keep, but reorganize:

- `Missing Values`
  - `fig-missing-heatmap`
  - `fig-missingness-per-group`
- `Abundance`
  - `fig-abundance-density`
- `Variance`
  - `fig-variance-violin`
  - `tbl-variance-summary`
- `Sample Structure`
  - `fig-pca`
  - `fig-correlation`
  - `fig-abundance-heatmap`

Notes:

- Keep this before differential-abundance results so users see data quality before interpreting hits.
- Consider making placeholder figure heights conditional in this QMD too, matching the Rmd compact-placeholder fix.

### 4. Differential Abundance

Purpose: all result summaries, diagnostic plots, significant-feature views, and contrast agreement.

Suggested subtabs:

- `Summary`
  - significant features by contrast and direction table, moved from current `Results`
  - optional count table including total tested, significant, not significant by contrast
- `Contrast Agreement`
  - UpSet plot for all significant features
  - UpSet plot for increased-abundance significant features (`diff > 0`)
  - UpSet plot for decreased-abundance significant features (`diff < 0`)
  - compact placeholder when fewer than two non-empty contrast sets are available
- `Volcano`
  - current volcano plot
- `MA Plot`
  - current MA and ranked MA plot
- `Significant Features`
  - current significant-feature heatmap
  - current significant-feature table or a filtered table view of significant rows
- `Distributions`
  - current fold-change and p-value histograms from `fig-fc-pvalue`

Important implementation detail:

- Build UpSet inputs locally from `significant_contrasts` using `protein_Id`, `contrast`, and `diff`.
- Do not add public package API.
- Keep set names as contrast names.
- Use `unique(protein_Id)` within each contrast/direction to avoid duplicate row effects.
- For direction-specific sets, use:
  - increased: `diff > 0`
  - decreased: `diff < 0`

### 5. Result Table

Purpose: full data table and export-oriented view.

Keep:

- Full `contrast_table` as `Differential abundance results`

Add or move if useful:

- Significant-only filtered table
- Sample annotation mirror only if users strongly expect all tables in one place; otherwise keep sample annotation in
  `Settings > Design` and avoid duplication.

Suggested subtabs:

- none by default; this is a single top-level table
- add subtabs only if additional full-detail tables are added later, for example `All Features` and `Significant Features`
- `Samples` only if duplication is explicitly desired

### 6. Session Info

Keep as final tab.

## Implementation Plan

1. Update only `inst/templates/quarto/Grp2Analysis_V2_SE.qmd`.
2. Do not edit generated HTML files under `inst/templates/quarto/*.html`.
3. Add local helper functions in the setup chunk:
   - `significant_sets_by_direction(significant_contrasts, direction = c("all", "increased", "decreased"))`
   - `draw_upset_or_empty(sets, empty_label)`
   These helpers stay private to the QMD.
4. Move existing chunks without changing their computational behavior:
   - Move sample counts and sample annotation into `Settings`.
   - Move feature detection overlap into `Feature Detection`.
   - Move significant summary table into `Differential Abundance > Summary`.
   - Move or replace the existing significant overlap section with the new `Contrast Agreement` section.
5. Keep chunk labels stable where possible to reduce generated asset churn. Rename labels only when their meaning changes.
6. Update captions to match the new section context:
   - Detection overlap: feature detection across groups.
   - Contrast agreement: significant differential-abundance calls shared across contrasts.
7. Render locally through the existing test pathway if possible:
   - `Rscript -e "testthat::test_file('tests/testthat/test-quarto-se-report.R')"`
8. If Quarto rendering is slow or blocked, run a focused R validation of:
   - helper functions on synthetic `significant_contrasts`
   - no duplicate protein IDs inside UpSet sets
   - empty-placeholder behavior for zero/one contrast

## Self-Review Pass 1

This structure fixes the main semantic issue: experimental design moves out of `Results`, identification comes before
statistical results, and all differential-abundance summaries and overlap plots live together. It also keeps the
machine-readable full result table available without letting it dominate the report narrative.

Risk: a top-level `Results` tab disappears. That is acceptable because the report title is already differential
abundance/analysis, and the meaningful results are now under `Differential Abundance`. If a top-level result landing page is
preferred, it should be a short `Summary` tab after `Quality Control`, not a mixed table dump before identification.

## Self-Review Pass 2

The plan avoids public API expansion and confines changes to the QMD source. It also avoids changing the
`SummarizedExperiment` schema, which is important because the report already reconstructs enough data from rowData,
colData, assays, and metadata. The only helper logic needed is local set construction for UpSetR.

The singular `Result Table` label is better than `Tables` because the current report has one full result table. Summary
tables should live with their corresponding narrative sections rather than inside the full result-table view.
