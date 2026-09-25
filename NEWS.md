# prolfquapp 2.10.5

- For the `lm_impute` default model, the DEA SummarizedExperiment and AnnData carry an `imputedData` assay/layer, `transformedData` with every missing cell filled by `prolfqua::impute_from_model()`, and an `imputation` rowData block with `n_observed`, `n_imputed` and `route` per feature. Decoys are not modelled and stay NA in both. `schema_version` is 2.1.0.
- The DEA SummarizedExperiment contrast frames carry the feature key columns again, so every rowData frame can be joined on the keys. Since 2.10.0 they were stripped with the annotation columns, which broke exploreDE. The contrast frames are now the model's contrasts without the annotation join; the xlsx export is unchanged.
- On the protein path, the DEA SummarizedExperiment and AnnData carry an `ibaq` assay/layer with the IBAQ values also written to `IBAQ_<workunit>.xlsx`; features without an IBAQ value are NA.
- `DEAResultReader` keys every table by the artifact's `feature_keys`: `lfq_raw` and `lfq_transformed` carry the real key columns (`protein_Id` and `site` for a site DEA), joined from the `annotation` frame on the feature id, instead of the flattened row name under `protein_Id`, and `subject_id` is the feature keys. With the keys back in the contrast frames, the tabbed report on a site DEA drew an empty significant-feature heatmap and UpSet sets of proteins; both are per site again.
- The AnnData stores every nested rowData frame (contrasts, statistics, imputation) as a `varm` data frame indexed by the feature ids, with its key columns, instead of numeric matrices plus `uns$prolfquapp$varm_columns`/`varm_annotations`; the order of every metadata list is kept in `uns$prolfquapp$uns_list_order`. Needs the anndataR fork `wolski/anndataR` (scverse/anndataR#525), which writes the data-frame index so Python `anndata` can read the file. AnnData files written by 2.10.0-2.10.4 are not read.
- `DEAResultReader` requires the artifact's `analysis_configuration`; the heuristic configuration for SummarizedExperiments without one is removed.
- `DEAResultReader` gains `annotation`, `samples`, `lfq_imputed` (from `imputedData`) and `imputation`, so downstream packages such as prophosqua read a DEA artifact only through it. `write_summarized_experiment_h5ad()` is exported.

# prolfquapp 2.10.4

- The tabbed DEA report shows the volcano plot again for SAINT analyses. It
  looked the panel up under the fixed name `FDR`, while the plotter names each
  panel after the backend's own column, so a SAINT run (`BFDR`) got "No plot is
  available for this data set" where the volcano belongs. Models whose FDR
  column is called `FDR` were unaffected.

- Align the report-template dependency with protsea so downstream packages can resolve their combined dependencies in CI.

- Read and atomically write complete MuData containers from R, preserving independent AnnData modalities, shared sample annotations, and container metadata for downstream PTM analysis.

- The per-feature peptide count is now called `nrPeptides` everywhere, in the
  annotation, the XLSX sheet and the AnnData artifact. Simulated and MSstats
  analyses used to call it `nr_peptides`, and two readers carried a line
  copying one spelling to the other, so a tool reading the artifact had to know
  which reader had run. Both aliases are gone. The `nr_peptides` *option*
  (minimum peptides per protein) is unchanged.

# prolfquapp 2.10.3

- New exported `write_h5ad_atomic()`: writes an AnnData to a temporary file
  beside the destination, reads it back, hands the restored object to a
  caller-supplied validator, and only then moves it into place. A crash or a
  failed validation leaves no partial `.h5ad` behind. The DEA artifact writer
  uses it, and downstream packages writing their own h5ad can share the same
  protocol instead of copying it.

- The Docker image now carries the prolfquasaint build that stamps
  `estimate_type` on SAINT contrast tables, so a SAINT analysis writes an
  artifact whose per-row estimate provenance downstream enrichment tools can
  read. Without it, asking to exclude imputed estimates failed on a SAINT
  result.

# prolfquapp 2.10.2

- The Docker image now refreshes the unpinned `prolfquasaint` and `saintexpress` dependency chain, so SAINT GSEA rank generation uses the current `score = NULL` contract instead of the dependencies embedded in the preceding image.

# prolfquapp 2.10.1

- GSEA rank files (`.rnk`) now carry the model's test statistic rather than the
  signed `-log10(p)` derived from it. The statistic is what the modelling
  backend actually reports, and it distinguishes features that share a p-value.
  SAINTexpress is unchanged: it reports no p-value, and its `SaintScore` is a
  bounded probability that carries no direction, so its rank stays `log2_EFCs`.
  Rank files written by earlier versions are not comparable to these.

# prolfquapp 2.10.0

- The `SummarizedExperiment` and the AnnData now record `identifier_key`, the
  annotation column holding the identifier enrichment tools are given (STRING,
  ORA). Which column that is depends on the reader that produced the analysis,
  so a consumer of the artifact no longer has to guess at column names.
- A metadata table (the contrast definitions, the model formula) is now written
  into the AnnData column by column rather than as a dataframe group: anndataR
  encodes a one-row table's columns as HDF5 scalars, which `anndata` in Python
  refuses to read, so every single-contrast analysis produced an `.h5ad` that
  Python could not open. Reading the file back in R is unchanged.
- New `SIM_PEPTIDE` input, the peptide-level twin of `SIM`, so the nested
  modelling facades (`lmer_nested`, `ropeca_nested`, `firth_nested`,
  `limpa_nested`, `binomial_nested`) can be exercised on simulated data.

- The example analysis now carries a fold-change threshold its own data can
  clear (0.4 rather than 1, against simulated effects that top out near 0.8), so
  the example and vignette renders of the tabbed report show populated
  significance tables and contrast-agreement plots instead of empty panels.
  Real analyses are unaffected: the threshold comes from the run's
  configuration.

- The tabbed differential-expression report now renders from `AnnData.h5ad`
  rather than from `SummarizedExperiment.rds`, so every analysis run that
  produces a report has demonstrated that its AnnData can be read back. Both
  files are still written, and the report reads either one.
- Fixed silent metadata loss when writing AnnData: a data frame stored in the
  metadata (the contrast definitions and the model formula) was flattened to a
  list of columns, and a list without names -- which AnnData cannot represent at
  all -- was written as an empty group. Data frames now survive the round-trip,
  unnamed lists are refused with an error naming the offending entry, and the
  order of layers, contrast tables and their columns is preserved. This also
  repairs `LFQData_from_anndata()`, which silently dropped every extra
  abundance layer when reading an `.h5ad` back, because the layer names it
  needed had been written as an unnamed list.

- New `DEAResultReader` reads a differential-expression result artifact — a
  `SummarizedExperiment`, an AnnData, or a path to an `.rds` or `.h5ad` file —
  back into the familiar prolfqua objects: the raw and transformed `LFQData`,
  the contrast table with the column names the backend produced, and a
  `ContrastsTable` that already knows the backend's column roles, so
  `significant()` and the contrast plots resolve columns by role. The tabset
  report is built on it and no longer reimplements any of that.
- Both differential-expression reports now decide what to show from the
  backend's column roles rather than from which backend ran: the significance
  filter, the volcano and score panels, and the difference-test introduction all
  follow the recorded roles. Results from a backend without a p-value (SAINT)
  are presented correctly without the report naming it.
- AnnData files written by a differential-expression run can now be read back
  into a `SummarizedExperiment`, so the `.h5ad` is a full-fidelity artifact
  rather than an export: layers, sample and feature annotation, every contrast
  table, and the analysis metadata all survive the round-trip.
- The `SummarizedExperiment` now stores feature annotation once, in its own
  `rowData` frame named `annotation`, and the per-contrast frames carry results
  only. In the AnnData that annotation frame is `var` and each contrast frame is
  its own `varm` entry under the name prolfquapp gave it. Its metadata states
  the artifact type and schema version (`2.0.0`), and carries the analysis
  configuration and provenance under the names `analysis_configuration` and
  `provenance`; readers of artifacts written by earlier versions must be
  updated.

# prolfquapp 2.9.1

- The `SummarizedExperiment` written by a differential-expression run now
  records the contrast column roles in its metadata as
  `contrast_configuration`, so downstream tools can find the contrast, effect,
  score and FDR columns by role instead of guessing their names. This makes
  results from backends with their own naming — SAINTexpress uses `Bait`,
  `log2_EFCs`, `SaintScore` and `BFDR` — readable by the same generic code that
  handles the standard schema.

- The `SummarizedExperiment` written by a differential-expression run now
  carries the processing options in its metadata, so reports rendered from it
  can state the settings the analysis actually used.
- The tabset report reads its settings straight from the `SummarizedExperiment`
  metadata and no longer falls back when they are absent, so it requires a
  `SummarizedExperiment` written by this version or later. Its example render
  now builds an example on the fly instead of loading a shipped snapshot, and
  the 4.4 MB `inst/extdata/3106962.rds` fixture it used has been dropped.
- The processing-parameter table is now `ProcessingOptions$parameters_table()`,
  shared by both differential-expression reports instead of being built inline
  in one of them. A cleared decoy pattern survives the round-trip through
  `SummarizedExperiment` metadata and is reported as "not identified".
- The differential-expression report now opens with a table of every processing
  parameter it was run with (peptide filter, aggregation, normalization,
  thresholds, contaminant and decoy handling). The contaminant and decoy rows
  describe what the pipeline actually does rather than the unused `remove_cont`
  / `remove_decoys` flags, and the peptide-filter and significance-threshold
  sections point at the table instead of restating the values.
- The redundant "at least two peptides" feature-count panel is gone from both
  differential-expression reports; it ignored the configured minimum peptides
  per protein and, at any minimum of two or more, simply duplicated the panel
  above it. The remaining per-sample count plot now says which protein matrix
  it is drawn from.

# prolfquapp 2.9.0

- Differential-expression runs now write `AnnData.h5ad` alongside
  `SummarizedExperiment.rds`, preserving sample and feature axes, abundance
  layers, feature annotations, contrast statistics, and analysis provenance for
  typed R/Python downstream workflows.
- `SummarizedExperiment` contrast tables now retain the authoritative feature
  names for empty decoy rows instead of synthetic `NA` row names.
- Protein annotation now keys on the level the analysis reports on: a site-level
  annotation stays one row per protein and site instead of being collapsed to
  one row per protein, so per-site columns such as the sequence window reach the
  result tables. Protein- and peptide-level analyses are unaffected.
- Full Docker releases now refresh their GitHub package dependencies once per
  release version, preventing cached older SAINT packages from replacing the
  constant-control-variance fix.

# prolfquapp 2.8.0

- Quarto reports now vendor and install the synchronized FGCZ template from
  `fgczQuartoTemplate@dbeb852`, including the current report DPI and
  View Source/per-figure code controls.

# prolfquapp 2.7.1

- The Docker image now installs the SAINTexpress constant-control-variance fix,
  allowing SAINT intensity analyses with constant complete control profiles to
  finish instead of failing with `NA/NaN/Inf in 'y'`.

# prolfquapp 2.7.0

- Command-line failures now retain and log their originating R call stack
  instead of reporting that no traceback is available.
- Volcano, MA, and score plots now consistently display observed estimates in
  black, LOD-imputed estimates in green, and group-mean fallback estimates in
  blue.
- Repeated-measures annotation processing now evaluates its subject-column
  and repeated-design checks as a single scalar condition.
- Shared Quarto reports now include the current FGCZ responsive figure-grid
  styling and opt-in full-width layout.
- Protein annotations are now joined to peptide and PTM results by protein ID,
  so every quantified feature retains its protein metadata without losing or
  multiplying result rows.
- The FGCZ Quarto dependency now follows its canonical
  `fgcz/fgczQuartoTemplate` upstream and camel-cased package name, including
  the latest shared report assets and tab/download controls.
- The experimental-design survey now links to NIST's maintained guidance on
  blocking factors instead of an obsolete external URL.
- Package builds now install Quarto visual abstracts as dedicated runtime
  assets instead of vignette output, eliminating an R CMD check NOTE while
  preserving runtime-rendered report overviews.
- Quarto reports now vendor the dynamic horizontal FGCZ toolbar: it stays below
  the visible banner, pins to the top-right while scrolling, and expands its
  Find/Download text labels on hover or keyboard focus.

# prolfquapp 2.6.1

- The Docker image now installs the `prolfqua` Firth-model performance fix, so
  `firth_nested` analyses with very large peptide effects complete instead of
  spending days on unused coefficient profile-likelihood intervals.

# prolfquapp 2.6.0

- The full Docker image now includes the fixed-size Quarto report figures and
  interactive widgets introduced in 2.5.1.

# prolfquapp 2.5.1

- Quarto report figures now use fixed display dimensions based on the historic
  six-inch figures, so QC and differential-expression diagnostics remain a
  readable size and fit vertically on ultrawide screens.

# prolfquapp 2.5.0

- Docker releases now derive their build mode from the semantic version: `X.Y.0`
  tags run the full multi-architecture build and runtime checks, while later
  `X.Y.Z` patch tags update the packages declared in `Remotes` and rebuild
  prolfquapp on the matching `X.Y.0` image. An urgent patch release fails if
  that full base image is unavailable or if the Docker environment or package
  dependencies changed since `X.Y.0`.
- Docker builds now use Posit's Ubuntu Noble R base image and native Posit
  Package Manager binaries on AMD64 and ARM64, avoiding lengthy source
  compilation while retaining Arrow's zstd support.
- Docker images are now published to the GitHub Container Registry
  (`ghcr.io/prolfqua/prolfquapp`) instead of Docker Hub. Images from earlier
  versions remain available on Docker Hub; new releases are published to GHCR.

# prolfquapp 2.4.1

- Quarto analysis reports now render standard static and interactive figures at a
  centred two-thirds of the content width. Intentional multi-column figure
  layouts, data tables, and compact Overview visual abstracts are unchanged.
- All five Quarto reports now open with a compact, report-specific visual Overview: three summary cards show the number of samples, experimental groups, and quantified proteins (or the analysed feature type), followed by the visual abstract. The reports finish with a two-subtab Session Info area: Report provenance records the B-Fabric/input context, creator, timestamp, software/model, and package version, while R session info contains only `sessionInfo()`. The visual abstracts are packaged with the reports, so the same layout is retained in runtime-rendered HTML output.
- Quarto visual abstracts are now copied as individual vignette assets, so `make build-vignettes` no longer fails when `devtools` stages report assets into `doc/`.
- Package builds now exclude Quarto's transient `vignettes/.quarto` freeze cache, avoiding nonportable tar-path warnings and cache files in source tarballs.
- Package and vignette builds now synchronize the four FGCZ Quarto assets from `fgczquartotemplate`, so generated reports use the current shared toolbar and styling.
- Quarto vignette extraction now uses safe defaults for report-local conditional metadata, preventing spurious missing-object errors while creating the companion `.R` sources.
- Differential Expression Analysis Quality Control report (`DiffExpQC_R6_tabset.qmd`) now presents its report/analysis
  metadata (Workunit, Order, Project, generated-by, timestamp, software, model, package version) once — as a table in a
  final Session Info tab, alongside `sessionInfo()` — instead of duplicating it in a top-of-page callout.
- Quarto reports now ship the updated right-aligned **Find / Download** toolbar asset, positioned below the FGCZ banner
  around one-quarter of the viewport height from the top; plot ZIP downloads can include Order/Workunit metadata and a
  current timestamp.
- The documentation website now uses Quarto repository source links instead of the broken altdoc `code-links: true`
  sidebar entry, so the source link points to GitHub instead of `undefined`.
- Interactive Plotly subplots now fade non-hovered keyed traces, so abundance-density curves from prolfqua become easier
  to inspect sample by sample.
- Quarto reports now render reliably under `R CMD check` and in fresh installs. The reports are rendered via `fgczquartotemplate::fgcz_render()`, which stages the FGCZ template assets (`_metadata.yml`, `fgcz.scss`, `fgcz_header_quarto.html`, `fgcz-plot-finder.html`) next to the report from the installed `fgczquartotemplate` package. Previously the reports used the `fgczquartotemplate-html` Quarto *extension* and relied on the `_extensions/` tree being shipped into the installed `doc/`, but the `vignettes/.install_extras` rule never actually shipped it, so rendering failed when no `_extensions/` directory was present (e.g. under `R CMD check`). The reports are now plain `format: html` documents styled by a directory-level `_metadata.yml`, with the Find/Download toolbar wired via `include-after-body: fgcz-plot-finder.html`; `fgczquartotemplate` was added to `Imports`.
- The explanatory info callouts in the Quarto reports (e.g. "Why look at the fold-change and p-value distributions?" in the tabbed DEA report, and the "About this report" note in the QC & sample-size report) are now collapsed by default, so the reports open with a cleaner overview and readers expand a note only when they want it.
- Restructured the "Protein Signal Intensities within Groups" report (`QC_ProteinAbundances_tabset.qmd`): the Protein abundances tab now leads with the iBAQ signal figure, followed by the table, with the column descriptions and the "What is the iBAQ signal?" explanation moved into collapsed info boxes; Order/Workunit metadata moved from the top of the page into a new final "Session Info" tab (alongside `sessionInfo()`). Also fixed the column-description list, where the `nrMeasured_<GroupName>` / `meanAbundance_<GroupName>` / `signal_percent_<GroupName>` names had their `<GroupName>` suffix silently dropped (it was parsed as an HTML tag); the group-suffix placeholders now render.
- Documentation website Articles menu cleaned up: the entries are now ordered Differential Expression Analysis, its Quality Control report, the tabbed Differential Expression report, Protein Signal Intensities within Groups, Quality Control & Sample Size Estimation, and the Auxiliary meeting-agenda article last. The `Grp2Analysis_V2_SE_tabset` report was retitled from "Differential Abundance Analysis" to "Differential Expression Analysis (Tabbed Report)" so it is clearly the tabbed variant of the main DEA report. The stub `prolfquapp.Rmd` introduction vignette (an incomplete skeleton duplicating the README) was removed.
- Tabbed DEA report (`Grp2Analysis_V2_SE_tabset.qmd`): rewrote the cryptic UpSet-plot captions (feature-detection overlap between groups; significant / increased / decreased features shared between contrasts) so each states that it is an UpSet plot and explains what the intersection bars, dot matrix, and set-size bars represent, and clarified the MA-plot and significant-feature heatmap captions.
- The default contrast model is now `lm_impute` instead of `lm_missing`. `lm_impute` refits proteins whose per-protein linear model failed or was singular by imputing at the limit of detection with borrowed variance (flagging rescued rows as `lod_imputed`), whereas the deprecated `lm_missing` substituted group means without a model fit. This changes default DEA results for proteins that could not be fit directly, removes the prolfqua deprecation warning emitted on every default run, and applies to `make_DEA_config_R6()`, `run_make_yaml()`, the `prolfqua_yaml.sh --model` default, and the CompoundDiscoverer DEA entry point. `lm_missing` remains available as an explicit `model =` / `--model` choice.
- The documentation website is now built with [altdoc](https://altdoc.etiennebacher.com/) (Quarto Website backend) instead of pkgdown. pkgdown's `tweak_tabsets` step crashes on the Quarto `panel-tabset`s used by the tabbed report vignettes; altdoc renders the vignettes natively through Quarto, so the tabset reports appear on the site with their tabs intact. Removed the obsolete `make quarto`/`render-quarto` preview targets from the Makefile.
# prolfquapp 2.4.0

- The DEA results `index.html` (`write_index_html`) is now rendered as an FGCZ Quarto entry page with separate deliverable tables for HTML reports, Excel workbooks, ORA input gene lists, and GSEA rank files, including captions, file sizes, report descriptions, and Excel-content descriptions.
- Report HTML output filenames now match their Quarto source rather than the workunit-based `DE_`/`QC_` scheme: a DEA run writes `Grp2Analysis_V2_R6.html`, `Grp2Analysis_V2_SE_tabset.html`, `DiffExpQC_R6_tabset.html`, and `QCandSSE_tabset.html`; a QC run writes `QC_ProteinAbundances_tabset.html` and `QCandSSE_tabset.html` — all inside the per-workunit `Results_WU_<workunit>/` (or QC output) folder. The index page continues to label links with descriptive report titles, not filenames.
- Renamed the Quarto report vignettes to drop the now-redundant `_quarto` suffix (every report is Quarto): the primary DEA report is `Grp2Analysis_V2_R6.qmd`, and the tabbed reports carry a `_tabset` suffix — `DiffExpQC_R6_tabset.qmd`, `QCandSSE_tabset.qmd`, and `QC_ProteinAbundances_tabset.qmd` (`Grp2Analysis_V2_SE_tabset.qmd` already followed this convention). The report HTML filenames produced by the DEA and QC pipelines are unchanged.
- Reviewed every figure and table caption across the five Quarto reports against the FGCZ searchable-caption rule and rewrote the 16 that were vague or inaccurate into specific scientific labels (naming the measured quantity, what points/bars represent, the axes/encoding, grouping, and transformations). This corrects two DEA figures that were labelled "Venn diagram" but actually draw UpSet plots, disambiguates the two previously identical per-sample protein-count captions (proteins with ≥1 vs ≥2 peptides), and fixes the protein-abundance figure caption (its x-axis is the abundance-rank percentile, not the signal contribution).
- SE tabset DEA report (`Grp2Analysis_V2_SE_tabset.qmd`): the fold-change / p-value figure now carries a specific, searchable scientific caption (replacing the vague "Fold-change and p-value summaries.") and is preceded by a callout explaining why both the fold-change distribution and the p-value distribution are inspected as model diagnostics (centred-near-zero fold-changes; approximately uniform p-values under the null).
- Removed the retired R Markdown report sources now that the pipelines render Quarto only: `Grp2Analysis_V2_R6.Rmd`, `DiffExpQC_R6.Rmd`, `QC_ProteinAbundances.Rmd`, and `QCandSSE.Rmd`. The legacy `DEAReportGenerator$render_DEA()` method and the `render`/`markdown`/`markdown_qc` arguments of `write_DEA_all()` (which drove the R Markdown rendering) were removed, and the exported helper `copy_DEA_R6_Files()` (which copied the R Markdown templates into the run's input folder) was removed. The two non-report R Markdown vignettes (`prolfquapp.Rmd`, `Auxiliary_ExDesignSurvey.Rmd`) are unaffected.
- The DEA and QC command-line pipelines now render **Quarto reports only**; the R Markdown reports are no longer rendered. The DEA run produces `DE_<workunit>.html` (the R6 Quarto report, primary), `DE_<workunit>_tabset.html` (the SummarizedExperiment tabset overview), `QC_<workunit>.html` (the differential-expression QC report, now a tabbed Quarto report), and `QC_sampleSizeEstimation.html` (the sample-size estimation report, now also produced during a DEA run). The QC run produces `proteinAbundances.html` and `QC_sampleSizeEstimation.html` from their Quarto reports. Report rendering is centralized in one place and each report renders independently, so a single report failure warns without aborting the run or dropping the others. The `.Rmd` report sources remain in the package but are no longer used by the pipelines.
- New Quarto report `DiffExpQC_R6_quarto.qmd`: a tabbed port of the differential-expression QC report (`DiffExpQC_R6.Rmd`), with Missing Values, Variance, and Differential Expression tabs, styled with the FGCZ Quarto extension. It builds from a serialized `DEAnalyse` `.rds` (falling back to `example_deanalyse()`).
- `QCandSSE_quarto.qmd`: the Sample Size Calculation section now splits into one sub-tab per tested log2 fold-change (0.59, 1, 2), each showing the sample-size bar chart and a smaller per-effect-size table, instead of stacked sub-figures and one wide table.
- `QCandSSE_quarto.qmd`: when the data contain no missing values (so the missing-value heatmap is empty) the report now shows an explicit "no missing values" placeholder in place of the heatmap, instead of dropping the figure and leaving a broken "Figure ??" cross-reference.
- `QC_ProteinAbundances_quarto.qmd` is now a two-tab report with no table of contents: the first tab (shown by default) holds the interactive protein table and the protein signal-contribution plot (the crosstalk table stays interactive), and the per-sample protein-count barplot moves to a second tab.
- QC runs now also produce the Quarto sample-size report (`QC_sampleSizeEstimation_quarto.html`) alongside the existing R Markdown `QC_sampleSizeEstimation.html`. It is rendered from a serialized copy of the QC data and receives the B-fabric project, order, and workunit identifiers, so its Workunit/Project/Order header is populated (falling back to "n/a" when an identifier is unset). Rendering is skipped with a warning, without failing the QC run, when the Quarto CLI or the installed report sources are unavailable.
- `example_deanalyse()` now sets example B-fabric identifiers (project/order/workunit), so the differential-expression Quarto report renders with a populated project context instead of blank fields when built from the bundled example.
- `QCandSSE_quarto.qmd` layout tweaks: the per-sample feature counts and the sample-overlap UpSet plot are now shown side by side, as are the missing-value histogram and the missingness heatmap. The sample-size figure is split into one sub-panel per log2 fold-change (0.59, 1, 2) instead of a single cramped facet grid, giving each effect size full height.
- The `QCandSSE_quarto.qmd` report is now laid out as a tabbed report (Quarto `panel-tabset`), matching `Grp2Analysis_V2_SE_tabset.qmd`. Top-level tabs are Introduction, Quality Control, Sample Size Calculation, Sample Mapping, and Session Info; the Quality Control tab groups Feature Detection, Missing Values, Abundance Distributions, Variance, and Sample Structure into sub-tabs, and the Variance sub-tab splits the raw-intensity coefficient of variation and the transformed-scale standard deviation into their own sub-tabs. All figures and tables keep their cross-reference numbers across tabs.
- New Quarto ports of the QC vignettes, `QC_ProteinAbundances_quarto.qmd` and `QCandSSE_quarto.qmd`, styled with the vendored FGCZ Quarto extension (left table of contents, Find/Save toolbar) to match the other Quarto reports. They build from example data by default (`example_qc_generator()` and the `data_ionstar` example) and take file parameters (`pap_file`, `qc_data_file`) for real inputs. The existing `QC_ProteinAbundances.Rmd` / `QCandSSE.Rmd` are unchanged and the QC pipeline still renders them.
- Restructured the `QCandSSE_quarto.qmd` report for clearer scientific storytelling: consistent section hierarchy (Quality Control → Feature Detection / Missing Values / Abundance Distributions / Variance / Sample Structure, then Sample Size Calculation), an introduction roadmap and a visible caveat callout, every figure and table now referenced in the text, and an explicit coefficient-of-variation → standard-deviation → sample-size narrative. The sample correlation heatmap and the transformed-intensity overview heatmap are now shown side by side. Coefficient-of-variation and standard-deviation summary tables are rounded (2 and 3 decimals), sample-size counts are shown as integers, and a Session Info section with a provenance line was added.
- Both Quarto reports now live in `vignettes/` as the single source of truth. The SummarizedExperiment tabbed report (`Grp2Analysis_V2_SE_tabset.qmd`) moved out of `inst/templates/quarto/` into `vignettes/` and adopts the vendored FGCZ Quarto extension styling (with the Find/Save toolbar and a left table of contents), matching `Grp2Analysis_V2_R6_quarto.qmd`. Both reports are shipped into the installed package's `doc/` directory, and the DEA CLI renders them from there: each run now also writes `DEAnalyse.rds` next to `SummarizedExperiment.rds` and produces the Grp2 differential-expression Quarto report (`*_quarto_dea.html`) alongside the SE report. Quarto rendering requires the package to be installed with vignettes built (the default for `make install`); when the Quarto CLI or the built sources are absent, report rendering is skipped with a warning instead of failing the run.
- DEA report (`Grp2Analysis_V2_R6.Rmd` and its Quarto port `Grp2Analysis_V2_R6_quarto.qmd`): the interactive PCA and volcano plots now use a responsive width with a fixed height instead of a fixed pixel width, fixing the tall/narrow rendering; the abundance-density subplot is given an explicit height so it no longer collapses.
- New Quarto vignette `Grp2Analysis_V2_R6_quarto.qmd`: a Quarto port of the differential-expression report, styled with the vendored FGCZ Quarto template extension and built through the package vignette machinery. It loads the analysis object from a serialized `.rds` via the `deanalyse_file` parameter (falling back to a generated example when none is supplied), places the table of contents in a wider right-hand margin, and enables the FGCZ plot-finder toolbar (the floating "Find" and "Save"/download buttons). The existing R Markdown report is unchanged.
- CLI logging: the command-line entry points now route stray `message()` and `warning()` output (from prolfqua core, dplyr joins, and prolfquapp itself) through the logger, so those lines appear with the standard `INFO`/`WARN [timestamp]` layout and are captured in the run log file instead of printing untagged. readr's column-specification chatter (`Rows:/Columns:`, `Use spec()`) is silenced. New exported helper `route_messages_to_logger()`.
# prolfquapp 2.3.3

- Protein abundance QC report (`QC_ProteinAbundances.Rmd`): contaminant and decoy proteins are now drawn on top of the regular points instead of being hidden behind them, and are rendered with a separate, higher opacity so they stand out. `plot_abundance_vs_percent()` gains a `highlight_alpha` argument (default `1`, so existing callers are unchanged); the report sets it to `0.8`.
- DiffExpQC report (`DiffExpQC_R6.Rmd`): reverted the four combined figures from interactive `plotly` subplots back to static `gridExtra::grid.arrange()` panels.
