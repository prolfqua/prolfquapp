# Replace pheatmap with ComplexHeatmap

## Goal

Replace the `pheatmap`-based heatmaps in the prolfqua ecosystem with
`ComplexHeatmap`, using the root plotting implementation rather than report-level
wrappers.

Target visual behavior:

- Abundance heatmaps use a green-black-red palette.
- Missing abundance values are displayed in light gray.
- Missingness heatmaps keep a clear observed/missing encoding.
- Existing report-facing methods stay minimal and stable unless a signature
  change is explicitly needed.

Also review the PCA implementation for correctness issues and propose targeted
improvements.

## Current Findings

- The root heatmap implementation is in `prolfqua/R/tidyMS_plotting.R`. **Four**
  functions call `pheatmap::pheatmap` directly: `plot_heatmap_cor()`,
  `plot_heatmap()`, `plot_na_heatmap()`, and `plot_raster()`.
- `prolfqua/R/tidyMS_plotting.R` also defines and exports an S3 method
  `print.pheatmap()` (registered in `NAMESPACE` via
  `S3method(print,pheatmap)`). It calls `grid::grid.newpage()` +
  `grid::grid.draw(x$gtable)` so heatmaps render correctly when knitr
  auto-prints them in a chunk. This is the actual rendering path the reports
  rely on and must be removed/replaced when `pheatmap` is dropped.
- `prolfqua/R/LFQDataPlotter.R` exposes the report-facing methods:
  `heatmap()`, `heatmap_cor()`, `na_heatmap()`, and `raster()`. Its roxygen
  examples (`stopifnot(class(...) == "pheatmap")`) and `@return pheatmap` tags
  must be updated.
- `prolfquapp` consumes these only via the R6 methods (`pl$heatmap()`,
  `pl$na_heatmap()`) and renders by knitr auto-print — there is **no**
  `inherits(plot, "pheatmap")` handling to remove, and **no** downstream package
  (`prolfquapp`, `prophosqua`, `prolfquappPTMreaders`, `prolfquasaint`) declares
  or calls `pheatmap`. So the change is confined to `prolfqua`.
- `prolfqua` currently declares `pheatmap` in `DESCRIPTION` (Imports) and imports
  it via roxygen in `R/A_dataset_docu.R`.
- Existing tests in `prolfqua/tests/testthat/test-plotting_functions.R` assert
  `pheatmap` classes (`expect_s3_class(p, "pheatmap")`) **and** inspect the
  `pheatmap` gtable directly (`p$gtable$grobs`, filtered by
  `inherits(grob, "text")`). ComplexHeatmap objects have no `$gtable`, so the
  label-inspection tests need a different access path, not just a class change.
- `plot_pca()` currently drops all features with any missing value, transposes
  the matrix, runs `stats::prcomp()` (default `center = TRUE`, `scale. = FALSE`),
  then joins scores to sample annotation.

## Implementation Plan

### 1. Confirm Heatmap Scope

- Search all ecosystem packages for direct `pheatmap::pheatmap`,
  `class(...) == "pheatmap"`, and report-specific pheatmap handling.
- Separate true source files from generated outputs such as rendered HTML,
  `doc/`, and copied report snapshots.
- Treat `prolfqua` plotting helpers as the root cause location for shared
  heatmap behavior.

### 2. Add ComplexHeatmap Dependencies

- Add `ComplexHeatmap` and `circlize` to the correct field in
  `prolfqua/DESCRIPTION`.
- Remove `pheatmap` only after all source usages and tests are migrated.
- Update roxygen imports in `prolfqua/R/A_dataset_docu.R`.
- Regenerate `NAMESPACE` with `make document`; do not edit `NAMESPACE`
  manually.

### 3. Create Small Internal Heatmap Helpers

Add internal helpers in `prolfqua`, scoped to plotting only:

- Build column annotations from sample metadata and `factor_keys`.
- Apply sample label suffix truncation consistently.
- Build row label truncation consistently.
- Render or return `ComplexHeatmap::Heatmap` objects in a way that works in
  reports, PDFs, and tests.
- Define shared colors:
  - Abundance values: green-black-red, centered around zero for row-scaled
    heatmaps.
  - Missing values in abundance heatmaps: light gray via `na_col`.
  - Missingness heatmap: observed and missing states with missing displayed
    distinctly; prefer light gray for missing if that is the desired final
    convention.

Keep these helpers internal unless a public API change is explicitly needed.

### 4. Migrate Source Heatmap Functions

Update these source functions in `prolfqua/R/tidyMS_plotting.R`:

- `plot_heatmap_cor()`
  - Compute correlations as before.
  - Preserve the current row/column ordering behavior unless a correctness issue
    is found.
  - Return a ComplexHeatmap object.

- `plot_heatmap()`
  - Preserve row scaling and NA-row filtering semantics. Note the current code
    z-scores twice — once via `t(scale(t(matrix)))` (used for NA filtering and
    `hclust` ordering) and again via pheatmap's `scale = "row"`. Re-z-scoring is
    idempotent, so the displayed values are a single row z-score. ComplexHeatmap
    has no `scale=` argument, so scale the matrix once yourself and pass it with
    `cluster_rows = FALSE` (rows are already ordered by `hclust`) to keep the
    same displayed values and ordering.
  - Use green-black-red colors for scaled abundance values, centered at zero via
    `circlize::colorRamp2()`.
  - Use light gray for missing values (`na_col`).
  - Preserve row and sample label truncation behavior.

- `plot_na_heatmap()`
  - Preserve the row filtering, row limit, and clustering intent.
  - Use an explicit observed/missing encoding and update the report prose if the
    colors change from white/black.

- `plot_raster()`
  - This function calls `pheatmap::pheatmap` directly and is exposed via
    `LFQDataPlotter$raster()`, so it must be migrated to satisfy the "no source
    code depends on pheatmap" acceptance criterion (see resolved scope note
    below). Preserve the mean/var ordering and the `not_na` ordering; it uses
    `cluster_rows = FALSE` and `cluster_cols = FALSE` (no clustering).

Remove the exported `print.pheatmap()` S3 method (and its roxygen
`@method print pheatmap` / `@export`), then regenerate `NAMESPACE` with
`make document` so `S3method(print,pheatmap)` is dropped. ComplexHeatmap objects
auto-print via their own `show`/`draw` mechanism, so no replacement S3 method is
needed for knitr chunks.

Update `LFQDataPlotter` return documentation (`@return pheatmap` on `heatmap()`
and `heatmap_cor()`) and the roxygen examples (the
`stopifnot(class(...) == "pheatmap")` lines, including the `raster()` example)
from `pheatmap` to the new object class.

### 5. Update Rendering Utilities and Reports

- Update `.render_plot_to_device()` (in `prolfqua/R/utilities.R`, used by
  `LFQDataPlotter$write_pdf()`/`write_boxplots()`) if `ComplexHeatmap` objects
  require `ComplexHeatmap::draw()` rather than the current `print()`. Printing a
  `Heatmap` object normally dispatches to its `show` method and draws, so verify
  with a real `write_pdf()` smoke test before changing it.
- No `inherits(plot, "pheatmap")` handling exists in `prolfquapp` templates, so
  there is nothing to replace there. The reports render heatmaps purely by knitr
  auto-print (e.g. `ph` / `plot_or_empty(nah, ...)` in
  `Grp2Analysis_V2_R6.Rmd`), which keeps working for ComplexHeatmap objects.
- Do not edit generated HTML output or the copied report snapshots under
  `inst/application/.../DEA_*/` and `ptm-pipeline/test_data/`; fix only the
  source template `vignettes/Grp2Analysis_V2_R6.Rmd` (and its metabo variant).
- Update prose in `prolfquapp` reports where it states old heatmap colors — the
  missingness heatmap caption currently says "White: Protein is observed, black:
  Protein is not observed" (`Grp2Analysis_V2_R6.Rmd`, ref:naHeat).

### 6. Update Tests

In `prolfqua/tests/testthat/test-plotting_functions.R`:

- Replace the four `expect_s3_class(p, "pheatmap")` assertions (one each in the
  `plot_raster` truncation test, the correlation, abundance, and raster suffix
  tests) with the ComplexHeatmap class (an S4 `Heatmap` object before drawing,
  `HeatmapList` after `draw()`) — use `expect_s4_class()` / `methods::is()`.
- The label tests currently read `p$gtable$grobs` and filter
  `inherits(grob, "text")`. ComplexHeatmap objects have no `$gtable`, so this
  inspection must be rewritten — read the labels off the heatmap object
  (e.g. row/column name slots) or capture grobs after `draw()` via
  `grid::grid.grab()`. This is the heaviest part of the test migration; budget
  for it rather than treating it as a class-name swap.
- Add a test that abundance heatmaps set missing values to light gray (`na_col`).
- Add a test that the abundance color mapping spans green, black, and red.
- Add or adapt smoke tests for writing heatmaps to PDF through
  `LFQDataPlotter$write_pdf()`.

Run focused validation first:

```bash
cd prolfqua
Rscript -e "testthat::test_file('tests/testthat/test-plotting_functions.R')"
```

Then run package-level validation:

```bash
cd prolfqua
make test
```

After installing `prolfqua`, validate downstream report usage:

```bash
cd ../prolfquapp
make test
```

### 7. Review PCA Implementation

**Status: implemented (2026-06-18).** Focused fixes applied to `plot_pca()` in
`prolfqua/R/tidyMS_plotting.R` and regression tests added to
`prolfqua/tests/testthat/test-plotting_functions.R` (45 plotting tests pass via
`devtools::load_all()`). See the per-item status and "implemented" notes below.

Review `plot_pca()` in `prolfqua/R/tidyMS_plotting.R`. Status against the code:

- **Real risk:** when all/nearly all features have an NA, `na.omit()` yields a
  0-row matrix and `prcomp()` errors. Needs a preflight check.
- **Partly handled:** the PC-count case is guarded by
  `if (max(PC) > (ncol(xx) - 1))`, but it `return(NULL)`s silently and the caller
  `pca_plotly()` then feeds `NULL` to `ggplotly()` and errors. The few-samples
  case still needs an explicit guard with an actionable message.
- **Already correct:** `prcomp(ff)` already uses `center = TRUE`,
  `scale. = FALSE` (the defaults), which is what's wanted for transformed LFQ
  data. Making these explicit is cosmetic, not a behavior change.
- **Conditional:** removing zero-variance features only matters if scaling is
  enabled (`scale. = TRUE`); with the current `scale. = FALSE` they are harmless,
  so only address this if scaling is added.
- **Valid defensive fix:** `inner_join(annotation, xx)` has no `by`, so it joins
  on common columns (currently `sample_name`) and emits a dplyr "Joining, by"
  message. Specify `by = sample_name` explicitly.
- **Already handled — no change:** `factor_keys[2]` is `NA` when only one factor
  key is present, and both the `geom_point(aes(shape = ...))` branch and the
  `scale_shape_manual()` branch are guarded by `if (!is.na(sh))`. A regression
  test is still worthwhile, but no fix is required here.
- **Valid:** duplicated sample names become duplicated row names and cause a
  cartesian expansion in the join; error early.
- **Valid (verified):** the imputation message points to `impute_with_zcomp()`,
  which exists (`R/LFQDataImp.R`). Wording can be made more actionable.

PCA follow-up changes (all done unless marked deferred):

- **Done.** Preflight checks now `stop()` with actionable messages for the
  empty-matrix (0 complete features after NA filtering) and insufficient-samples
  (`ncol < max(PC) + 1`) cases instead of returning a bare `NULL` that breaks
  `pca_plotly()`. The post-`prcomp()` PC-count guard also `stop()`s now.
- **Done.** Scores join annotation with explicit `by = sample_name` (no more
  dplyr "Joining, by" message).
- **Done.** Errors early on duplicated sample names (before the transpose/join).
- **Done.** `prcomp()` made explicit (`center = TRUE, scale. = FALSE`) to
  document intent; PCA output stays `ggplot`, no public API change.
- **Deferred (unchanged from analysis above):** zero-variance feature removal —
  only matters under `scale. = TRUE`; harmless with the current `scale. = FALSE`.
- **Done.** Added focused tests: missing-heavy data (error), too-few-samples
  (error), duplicate sample names (error), and a single-factor / explicit-join
  regression guard for the `factor_keys[2]` NA path.

## Acceptance Criteria

- No source code in `prolfqua` directly depends on `pheatmap` (downstream
  packages already have no `pheatmap` dependency). `pheatmap` is removed from
  `DESCRIPTION`, the roxygen `@importFrom pheatmap pheatmap`, and the
  `print.pheatmap()` S3 method; `NAMESPACE` is regenerated.
- `prolfqua` heatmap functions return/draw ComplexHeatmap objects.
- `plot_raster()` calls `pheatmap::pheatmap` directly, so it **is** in scope and
  must be migrated (resolves the earlier "out of scope unless confirmed" note —
  the dependency is confirmed).
- Abundance heatmaps use green-black-red colors and show missing values in light
  gray.
- Downstream `prolfquapp` reports continue to render the heatmaps with no
  special pheatmap handling (they already use knitr auto-print).
- PCA review findings are either implemented as focused fixes or recorded with
  a clear reason for deferring.
- Focused plotting tests and package tests pass for `prolfqua`; downstream
  `prolfquapp` tests pass after installing the updated `prolfqua`.
