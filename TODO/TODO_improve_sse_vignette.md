# Improve the QCandSSE Quarto vignette (QCandSSE_quarto.qmd)

## Goal and scope

This is a **planning document only** — no code is written yet. The single target file is
`/Users/wolski/projects/prolfqua_fml/prolfquapp/vignettes/QCandSSE_quarto.qmd`, the FGCZ Quarto QC + sample-size vignette (`format: fgczquartotemplate-html`). The work has three drivers from the user: (1) place the transformed-intensity overview heatmap next to the sample-correlation heatmap (or pair the correlation heatmap with the pairwise scatter matrix); (2) fix the ~7-significant-digit rounding in `tbl-cv-quantiles` and audit the sibling tables for the same defect; (3) improve the scientific storytelling and section consistency.

Hard constraints every change must respect:

- The vignette stays **self-contained** and must keep building from the `data_ionstar` fallback in the `setup` chunk. `data_ionstar` is an untransformed data set, so `show_text` is `TRUE` in the default build — this is the branch we must test.
- Every `prolfquapp::plotly_ggplot_subplot` chunk (`fig-plot-distributions`, `fig-intensity-distribution`, `fig-sd-violin`, `fig-sd-ecdf`) **must keep `#| fig-dpi: 96`** or the plotly subplots inflate to ~2000px. Do not add `fig-dpi: 96` to the static `print()` heatmap chunks — they should keep inheriting the extension's 300 dpi.
- Cross-references use `@fig-*` / `@tbl-*`. The FGCZ extension sets a global `out.width` of 40% and `fig-dpi` 300 for static figures. **Quarto renders an unresolved cross-reference as visible `?@tbl-...` text (a warning, not a build failure)** — so a broken reference degrades the default build's output even though it does not crash it.
- Do not break rendering or the cross-reference graph. Reordering chunks is safe (refs resolve by label, not position) **only if** every reordered chunk stays after `transformIntensities`, which defines `dataTransformed` / `plotter_transformed`, and after the chunk that defines `stats` / `st`.

Scope is deliberately kept at the report call site (prose, chunk options, chunk order, kable formatting). Anything that would add package API (a shared `example_qc_sse()` helper, factoring `empty_report_plot` into `prolfquapp`, renaming the `project_conf` param across the CLI) is out of scope for this pass and listed under Open questions for explicit sign-off. Per the repo Changelog Discipline, a `NEWS.md` bullet is required once the user-visible report changes land.

## Current structure

As-is section / figure / table map (line numbers approximate, chunk labels exact):

- **Front matter** — `author: "WEW@FGCZ.ETHZ.CH"`; params `qc_data_file`, `project_conf`, `plot_density`, `plot_sd_vs_mean`, `target_type`.
- **`setup`** (`include=FALSE`) — loads data or `data_ionstar` fallback, sets theme, builds `strings` translation list.
- **`# Introduction`** — B-fabric bullets (Workunit/Project/Order, blank in the vignette build); the "no automated conclusions, discuss with your bioinformatician" caveat is **hidden inside an HTML comment** (line 95); two-sentence purpose statement, no roadmap.
- **`# Quality Control: Identifications and Quantifications`** (one mega-H1 covering ~70% of the report):
  - `typicalObservations` + prose on identification counts.
  - `hierarchyCounts` — totals table via `prolfqua::table_facade` (line 128), **no `#| label`**; also runs `library(rlang)` / `library(prolfqua)` mid-document (lines 123-124).
  - `fig-hierarchy-barplot` — per-sample counts (orphan: never `@`-referenced).
  - `fig-upset-missing` — detection overlap (orphan).
  - `fig-missing-histogram` (orphan), `fig-missingness-heatmap` (orphan; `preparemissingnessHeatmap` uses `include=FALSE` header arg).
  - `checktransformation` → `show_text`; `transformIntensities` builds `dataTransformed` / `plotter_transformed`.
  - `variability` (`#| eval: !expr show_text`, `#| results: asis`) — **emits the `## Variability of raw intensities` heading and its intro via `cat()`**, forward-referencing `@fig-plot-distributions` and `@fig-intensity-distribution`.
  - `fig-plot-distributions` (plotly, `fig-dpi:96`, `eval: show_text`) — CV densities; defines `stats`.
  - `computeCVQuantiles` — **has both `include=FALSE` header arg and a `#| eval: !expr show_text` hash-pipe option** (fragile mix, lines 268-269); reads `stats`, builds `cv_quantiles_res`.
  - `tbl-cv-quantiles` — `knitr::kable(cv_quantiles_res, digits = getOption("digits"))` → **~7 sig figs (Table 1, user's flagged defect)**.
  - **`## Variability of transformed intensities`** (static markdown; one dense intro paragraph referencing five targets):
    - `fig-intensity-distribution` (plotly, `fig-dpi:96`, `eval: show_text`).
    - `preparecorrelationHeat` → `fig-correlation-heat` (static, square sample×sample; reads `width`, line 317).
    - `fig-pairs-smooth` (static, dense 12×12 matrix; orphan).
    - `fig-sd-violin` (plotly, `fig-dpi:96`; defines `st`), `fig-sd-ecdf` (plotly, `fig-dpi:96`).
    - `fig-sd-vs-mean` (`eval: !expr plot_sd_vs_mean`; orphan).
    - `computeSDQuantiles` → `tbl-sd-quantiles` — `kable(..., digits = getOption("digits"))` (~7 sig figs).
    - `prepareoverviewHeat` → `fig-overview-heat` (static, tall target×sample; reads `width`, line 402; **orphan; stranded at the very end after the SD material — this is the "Figure 11" the user wants moved**).
- **`# Sample Size Calculation`** — intro re-references `@fig-sd-violin`; effect-size paragraph; `fig-sample-size`; `tbl-sample-size` (`kable(..., digits = getOption("digits"))`, `sdtrimmed` fractional + N counts; guarded by `if (exists("tmp") && !is.null(tmp))`, lines 452-454); the power/confidence/fold-change **definitions arrive after** the figure and table (lines 459-464).
- **`# Appendix`** — `tbl-sample-mapping` (character, `table_facade`, no `#| tbl-cap`), `tbl-hierarchy-counts-sample` (integer counts, `table_facade`; duplicates the barplot content).
- **`resetTheme`** (`include=FALSE`). **No Session Info section, no closing provenance line** (unlike both sibling DEA reports).

The heading tree is lopsided: identification-count and missingness content have no subheadings, the only H2s are the two "Variability" siblings, and those two are authored by two different mechanisms (`cat()` vs static markdown).

**Two verified traps that any reorder/de-orphan plan must handle (see the relevant sections below):**

1. **`table_facade` tables cannot be cross-referenced as-is.** `prolfqua::table_facade(df, caption)` is literally `knitr::kable(df, digits = getOption("digits"), caption = caption)` (verified in `prolfqua/R/utilities.R:370`). The caption is baked into the kable HTML rather than supplied via Quarto's `#| tbl-cap`, so Quarto assigns the table no number and any `@tbl-...` reference renders as `?@tbl-...`. This is exactly why every **working** table cross-reference in the file (`tbl-cv-quantiles`, `tbl-sd-quantiles`, `tbl-sample-size`) uses bare `knitr::kable()` + `#| tbl-cap:`, while all three `table_facade` tables are currently un-referenced. Adding `#| label` alone is **not** sufficient.

2. **The `width` variable is shared and clobbered.** `preparecorrelationHeat` sets `width <- max(8, nrSamples * 0.15)` (line 317) and `prepareoverviewHeat` sets `width = max(10, nrSamples * 0.15)` (line 402) — the same name. `fig-correlation-heat` reads `#| fig-width: !expr width` **and** `#| fig-height: !expr width` (lines 323-324, square); `fig-overview-heat` reads `#| fig-width: !expr width` (line 408). Both prep chunks are `include=FALSE`; to place the two figure cells adjacent inside a layout div, both preps must run before the div, so `prepareoverviewHeat` overwrites `width` before `fig-correlation-heat` evaluates its size expressions — the correlation heatmap would render at the overview's width.

## Proposed structure

Adopt the SE-tabset family's QC vocabulary as **static** H2/H3 subsections (not tabsets — see Open questions; the cross-reference prose and linear vignette read-through favour a flat scroll). Exact target outline (heading text and levels):

```
# Introduction
  (roadmap paragraph + visible callout with the "no conclusions" caveat + forward ref to @tbl-sample-mapping)

# Quality Control
## Feature Detection
    tbl-hierarchy-counts (totals, converted to kable + tbl-cap) + fig-hierarchy-barplot + fig-upset-missing
## Missing Values
    fig-missing-histogram + fig-missingness-heatmap
## Abundance Distributions
    fig-intensity-distribution   (raw-vs-transformed density; keep eval: show_text)
## Variance
### Coefficient of variation (raw intensities)
    fig-plot-distributions + tbl-cv-quantiles   (bodies gated eval: show_text)
### Standard deviation (transformed intensities)
    fig-sd-violin + fig-sd-ecdf + fig-sd-vs-mean + tbl-sd-quantiles
## Sample Structure
    fig-correlation-heat + fig-overview-heat  (side by side)  ...  fig-pairs-smooth (full width, below)

# Sample Size Calculation
  (collapsible "Key concepts" callout with power/significance/fold-change BEFORE the results)
  fig-sample-size + tbl-sample-size

# Appendix
    tbl-sample-mapping (primary) + tbl-hierarchy-counts-sample (or drop — see checklist)

# Session Info
    sessionInfo() + closing italic provenance line
```

Rationale and rules:

- Rename the mega-H1 to **`# Quality Control`**; the identification/quantification content now lives under the descriptive `## Feature Detection`, so the H1 no longer mis-scopes its content.
- The two variability siblings become **`### Coefficient of variation (raw intensities)`** and **`### Standard deviation (transformed intensities)`** under one static **`## Variance`** H2. This is the central structural fix: both headings become static markdown; only the *body* chunks keep `#| eval: !expr show_text`. This removes the `cat()`-emitted heading and its string-embedded `@`-refs.
- The `## Sample Structure` subsection groups the three sample-similarity figures consecutively so the requested side-by-side pairing (`fig-correlation-heat` + `fig-overview-heat`) becomes natural. `prepareoverviewHeat` and `fig-overview-heat` move up to sit immediately after `fig-correlation-heat`. Both are static `print()` chunks, so this cannot error — **but it is not free of the `width` clobber described above; the fix is mandatory (see Layout changes).** This subsection ends the QC section, immediately before Sample Size Calculation.
- **Pre-transformed case (`show_text == FALSE`, e.g. SE-derived input) — reference gating.** Making the CV and Abundance-Distributions lead-ins static must not leave dangling refs on this path. Any prose that references a `show_text`-gated exhibit (`@fig-plot-distributions`, `@tbl-cv-quantiles`, and `@fig-intensity-distribution`) must itself be gated to the same condition — emit the referencing sentence from a small `#| eval: !expr show_text` `#| results: asis` `cat()` block — with a **parallel static fallback sentence** for the `show_text == FALSE` case, for **both** the `### Coefficient of variation (raw intensities)` subsection **and** `## Abundance Distributions` (which, under the new outline, holds only `fig-intensity-distribution` and would otherwise be a heading + dangling `@fig-intensity-distribution` ref with no figure). Keep the H2/H3 headings themselves static and generic; do not interpolate `strings[["targets"]]` into headings via asis. Suggested fallback text: "The data were already transformed on import, so raw-intensity CVs are not reported." and "The data were already transformed on import; only the transformed distribution is shown." *(If the team prefers, an acceptable alternative is to declare `show_text == FALSE` explicitly out of scope and guarantee only the `TRUE` branch — but the default plan is to gate the refs as above.)*

## Storytelling improvements

The report already contains the right latent arc — identifications → reproducibility of detection → intensity distributions → variability (CV then SD) → variance estimate → sample size — but never states it and partially scrambles it. The prose fixes below are all low-risk (prose, `@`-refs, chunk reordering) and touch no plotting code.

**Introduction roadmap.** After the purpose sentence, add a short roadmap paragraph naming the two parts and pre-connecting them: a QC part (how many features were identified, how reproducibly detected, how their distributions compare, how variable they are) and a sample-size part that *consumes* the variability estimated in QC to size the main experiment. End it with a forward reference to the sample-mapping key: "Sample names used throughout are defined in @tbl-sample-mapping." **This forward-ref only resolves once `tbl-sample-mapping` is converted from `table_facade` to `knitr::kable` + `#| tbl-cap` (see Formatting & rounding); otherwise it renders as `?@tbl-sample-mapping`.**

**Surface the hidden caveat.** Promote the HTML-commented "no automated conclusions / discuss with your bioinformatician" note (line 95) to a visible Quarto callout right after the roadmap, matching the SE report's visible-callout convention:

```
::: {.callout-note}
This report is generated automatically and contains no experiment-specific conclusions.
Please discuss it with your project bioinformatician or statistician.
:::
```

**De-orphan every figure.** Seven figures and both appendix tables currently render with no `@`-reference. Add one motivating sentence (and an `@`-ref) for each in the nearest prose, so each numbered exhibit has a narrative anchor and "what good looks like" cue:

- `@fig-hierarchy-barplot`, `@fig-upset-missing` — cite from the Feature Detection intro. Also convert the `hierarchyCounts` totals table from `table_facade` to `knitr::kable` + `#| label: tbl-hierarchy-counts` + `#| tbl-cap:` (see below), then add one sentence linking totals (`@tbl-hierarchy-counts`), the per-sample barplot (`@fig-hierarchy-barplot`), and the appendix per-sample counts (`@tbl-hierarchy-counts-sample`).
- `@fig-missing-histogram`, `@fig-missingness-heatmap` — cite from the Missing Values intro (LOD/reproducibility framing already exists nearby).
- `@fig-pairs-smooth` — "Pairwise scatter plots (@fig-pairs-smooth) let us inspect sample-to-sample agreement in detail; tight diagonal clouds indicate reproducible samples, curvature or outliers flag deviating samples."
- `@fig-overview-heat` — "@fig-overview-heat gives an overview of all transformed intensities; consistent columns and replicate clustering indicate a homogeneous data set."
- `@fig-sd-vs-mean` — gate the reference behind the same `plot_sd_vs_mean` condition as the chunk: "When shown, @fig-sd-vs-mean confirms that vsn has stabilized the variance — SD should be roughly flat across mean intensity."

**Split the dense intro paragraph.** The single paragraph at the top of the old transformed-variability section (line 284) forward-references five targets that, after regrouping, live in three subsections. Split it so each new H2/H3 opens with a 1–2 sentence lead-in that references only its own figures: Abundance Distributions gets the vsn/normalization sentence + `@fig-intensity-distribution` (gated, per the pre-transformed rule); the SD H3 gets the "cannot report CV on the transformed scale, use SD" sentence + `@fig-sd-violin` / `@fig-sd-ecdf` / `@tbl-sd-quantiles`; Sample Structure gets the correlation/clustering sentence + `@fig-correlation-heat` / `@fig-pairs-smooth` / `@fig-overview-heat`. This also makes prose order match render order (currently the intro promises SD material next but the page shows the correlation heatmap and pairs matrix first).

**The CV→SD→sample-size bridge (highest-value narrative fix).** The pivotal transition is compressed into one clause. Expand it into a bridge paragraph in the Variance section: the coefficient of variation is only interpretable on a ratio scale with a meaningful zero; after variance-stabilizing normalization the intensities are on an approximately log2 scale where `sd/mean` is no longer meaningful, so we quantify variability directly as the standard deviation on the transformed scale; because that scale is roughly log2, this SD approximates the relative (fold-change) error on the original scale — exactly the input a fold-change-based power calculation requires. Then open Sample Size Calculation by naming the hand-off explicitly: "We take the SD at the 50th and 75th percentile from @tbl-sd-quantiles as the assumed within-group standard deviation and feed it into a two-sample t-test power calculation." This stitches the two halves into one argument.

**Move the definitions before the results.** The power (1−β), confidence (1−α), and fold-change definitions (lines 459-464) currently arrive as a postscript after `fig-sample-size` / `tbl-sample-size`. Move them above the results, ideally in a collapsible callout so they prime interpretation without interrupting expert readers:

```
::: {.callout-note collapse="true"}
## Key concepts: power, significance, effect size
... the three existing definition paragraphs ...
:::
```

**Generic closing synthesis.** The QC section currently ends abruptly. Add a brief *generic* interpretation guide (compatible with "no automated conclusions") at the end of the QC section: reproducible QC data show comparable intensity distributions across samples, high sample-to-sample correlations, replicate clustering in the overview heatmap, and missingness confined to low-intensity features near the LOD; deviations should be discussed before proceeding. Frame it as criteria, not a verdict.

## Layout changes (user request 1)

The three sample-structure figures have conflicting aspect ratios: `fig-correlation-heat` is square (sample×sample), `fig-overview-heat` is tall (target×sample), and `fig-pairs-smooth` is a dense 12×12 matrix that is illegible below full width. All three are static `print()` outputs (no plotly), so a column layout carries no `fig-dpi:96` concern.

**Recommended arrangement:** pair **`fig-correlation-heat` next to `fig-overview-heat`** (the user's first option), and keep `fig-pairs-smooth` on its own full-width line. Reject the correlation+pairs pairing (the user's alternative): halving the 12×12 pairs matrix destroys its legibility.

**Mandatory pre-step — fix the shared `width` clobber.** Before wrapping the two chunks side by side, rename the shared variable so each figure keeps its own size. In `preparecorrelationHeat` set `width_cor <- max(8, nrSamples * 0.15)` and in `prepareoverviewHeat` set `width_overview <- max(10, nrSamples * 0.15)`; update `fig-correlation-heat` to `#| fig-width: !expr width_cor` / `#| fig-height: !expr width_cor` and `fig-overview-heat` to `#| fig-width: !expr width_overview`. (Alternatively hard-code the two figure cells' sizes.) Without this, `prepareoverviewHeat` overwrites `width` before `fig-correlation-heat` reads it and the correlation heatmap renders at the overview's width.

**Exact Quarto syntax.** After moving `prepareoverviewHeat` + `fig-overview-heat` up next to `fig-correlation-heat`, wrap the two figure chunks in an unequal-column layout div so each keeps its own label and caption (preserving `@fig-correlation-heat` / `@fig-overview-heat` cross-refs). Give the tall overview more width and the square correlation less, and add per-chunk `#| out-width: "100%"` so each fills its column instead of shrinking to the extension's global 40%:

```
::: {layout="[[40,60]]"}
```{r}
#| label: fig-correlation-heat
#| out-width: "100%"
#| fig-cap: ...
print(chmap)
```
```{r}
#| label: fig-overview-heat
#| out-width: "100%"
#| fig-cap: ...
print(hm)
```
:::
```

Reduce `fig-overview-heat`'s `fig-height` from 10 toward ~7 so it does not tower over the square correlation heatmap in the 60% column. Keep `fig-pairs-smooth` full width with `#| out-width: "100%"` on its own line below the div.

**Render risk (medium — must be verified once):** the interaction between a `layout` div holding two computational figure cells and the extension's global 40% `out.width` is the one thing to confirm on a trial render — check that both heatmaps size sensibly, that cross-refs still resolve, and that the page body does **not** scroll horizontally.

**Family-consistent fallback (see Open questions):** if the side-by-side render proves fragile, group the three sample-structure figures in a `::: {.panel-tabset}` with tabs Correlation / Pairwise scatter / Overview heatmap, mirroring the SE report's "Sample Structure" tabset — each figure then renders at its natural size and no cross-ref is lost.

**Independent sizing fix:** regardless of the side-by-side layout, the global 40% width crushes the dense exhibits. Add `#| out-width: "100%"` to `fig-pairs-smooth`, `fig-overview-heat`, and `fig-missingness-heatmap` so the information-dense matrices/heatmaps use full text width while small single-panel plots keep the 40% default.

## Formatting & rounding (user request 2)

Three tables use `knitr::kable(..., digits = getOption("digits"))`, which prints ~7 significant figures. **This is a pre-existing defect, not a Quarto-port regression:** `prolfqua::table_facade` also defaults to `digits = getOption("digits")` (verified in `prolfqua/R/utilities.R:370`), so the original `.Rmd`'s CV/SD tables were equally unrounded. The family standard is to round, so fix each of the three:

- **`tbl-cv-quantiles` (Table 1, the flagged defect):** CV is a percentage (`sd/mean*100`, hence `xlim(0,150)`). Round to **2 dp** to match the sibling `Grp2Analysis_V2_R6_quarto.qmd` `tbl-cv`.
- **`tbl-sd-quantiles`:** transformed-scale SDs (~0.1–0.5). Round to **3 dp** to match the SE report's `tbl-variance-summary`.
- **`tbl-sample-size`:** mixed — `sdtrimmed` is fractional, the per-delta N columns are sample counts and must read as integers. A single `digits=` cannot do both.

**Exact recipe.** Prefer the family idiom (rounding the frame before `kable` so it is robust to integer key columns and column reordering):

```r
# tbl-cv-quantiles
cv_quantiles_res |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 2))) |>
  knitr::kable()

# tbl-sd-quantiles
sd_quantile_res2 |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3))) |>
  knitr::kable()

# tbl-sample-size: KEEP the existing exists()/!is.null() guard; round sdtrimmed, keep N integer
if (exists("tmp") && !is.null(tmp)) {
  tmp |>
    dplyr::mutate(sdtrimmed = round(sdtrimmed, 3)) |>
    dplyr::mutate(dplyr::across(where(\(x) is.numeric(x) && !identical(x, sdtrimmed)), ceiling)) |>
    knitr::kable()
}
```

Preserve the `if (exists("tmp") && !is.null(tmp))` guard around `tbl-sample-size` (lines 452-454): `tmp` is created only in the else-branch of `fig-sample-size`, so an unguarded pipe errors when `tmp` is absent. Never leave `digits = getOption("digits")`. Avoid a positional `digits=` vector for `tbl-sample-size` — it is fragile if `relevant_factor_keys()` changes column order. `ceiling()` the N columns so a fractional N from `power_t_test_quantiles` still reads as an integer sample count. *(Simplify the `ceiling()` selector to the concrete N column names if that reads more clearly at implementation time.)*

**Convert the cross-referenced `table_facade` tables to `knitr::kable` + `#| tbl-cap`.** Any table that must acquire a resolvable `@tbl-` reference has to leave `table_facade` behind (its caption is baked into the kable HTML, so Quarto never numbers it). Convert:

- the `hierarchyCounts` totals table → `knitr::kable(...)` with `#| label: tbl-hierarchy-counts` and the caption moved into `#| tbl-cap:` (matching the family idiom in `Grp2Analysis_V2_R6_quarto.qmd`);
- `tbl-sample-mapping` → `knitr::kable(...)` with `#| tbl-cap:` (so the Introduction's `@tbl-sample-mapping` forward-ref resolves).

`tbl-hierarchy-counts-sample` needs the same conversion only if it acquires a `@`-reference; if it stays un-referenced it can remain `table_facade`.

**Appendix tables need no numeric rounding.** `tbl-sample-mapping` renders character raw-file/sample names and `tbl-hierarchy-counts-sample` renders integer counts, so the default digits has no visible effect on values — the only appendix change is the `table_facade`→`kable` conversion needed for cross-referencing, not a rounding fix. Do not "fix" what is not broken.

**Units and caption clarity (small, do with the rounding).** CV is a percentage but neither the axis nor the caption says so: set `x = "CV [%]"` in both `stats$density()` panels of `fig-plot-distributions` (currently `x = "CV"`, lines 250/254), and append "(in %)" to the `tbl-cv-quantiles` caption. Align the two quantile-table captions in parallel phrasing and note the scale, e.g. "CV quantiles (in %) at the 50th–90th percentiles" and "Standard-deviation quantiles (transformed / vsn scale) at the 50th–90th percentiles."

## Family consistency & robustness

**Front matter.** Change `author: "WEW@FGCZ.ETHZ.CH"` to a team string consistent with the family — `"Functional Genomics Center Zurich"` (as in `Grp2Analysis_V2_R6_quarto.qmd`). `date: today` and the title already match.

**Empty B-fabric fields.** In the default vignette build `local_project_conf` is `list()`, so the Introduction bullets render as bare "Workunit:/Project:/Order:". Guard with `%||%`, e.g. `` `r local_project_conf$workunit_Id %||% "n/a"` ``, or omit the bullets when the fields are NULL. **Dependency:** the Introduction inline code (lines 91-93) executes *before* `hierarchyCounts` (line 123), which is the only place `library(rlang)` is currently loaded, and base-R `%||%` exists only in R ≥ 4.4. So this guard requires the "move `library(rlang)` into `setup`" item to land first (or an explicit R ≥ 4.4 assumption). Sequence the two checklist items accordingly.

**Normalize chunk options to hash-pipe.** The file mixes old-style header args with hash-pipe options. The dangerous case is `computeCVQuantiles`, which carries **both** `include=FALSE` (header) and `#| eval: !expr show_text` (hash-pipe, lines 268-269): if the hash-pipe `eval` is not honoured, the chunk runs when `show_text` is FALSE and references `stats`, which is only created under the same gate in `fig-plot-distributions` — an undefined-object abort. Convert `computeCVQuantiles`, `preparemissingnessHeatmap`, `preparecorrelationHeat`, `prepareoverviewHeat`, `transformIntensities`, `checktransformation`, and `hierarchyCounts` to hash-pipe-only options (`#| label:`, `#| include: false`, `#| eval:`), matching the family, and guard the CV chunk with the same gate as its figure (or `exists("stats")`).

**Library loading.** Move `library(rlang)` / `library(prolfqua)` out of the `hierarchyCounts` body (lines 123-124) into `setup` (the family loads in setup; SE validates a `required_packages` vector). This is also the prerequisite for the `%||%` guard above.

**Session Info + provenance.** Add a trailing `# Session Info` H1 with a `sessionInfo()` chunk (SE style) and a closing italic provenance line naming `QCandSSE_quarto.qmd` and `packageVersion("prolfquapp")` (Grp2 R6 style). This closes the only footer gap versus the DEA reports.

**fig-dpi discipline (leave as-is; document only).** Every plotly-capable chunk correctly carries `#| fig-dpi: 96`; static heatmaps correctly omit it. Do not change this. Optionally add a one-line comment near the first plotly chunk explaining why `fig-dpi: 96` is mandatory for `plotly_ggplot_subplot`, so a future editor does not strip it. (Minor wrinkle: `fig-plot-distributions` and `fig-sd-violin` fall back to static gridExtra when `plot_density` is FALSE, where 96 dpi is slightly low-res but harmless; `plot_density` defaults TRUE, so no action needed.)

## Prioritized action checklist

### P0 — user-flagged defects and the structural root cause

- [ ] Round `tbl-cv-quantiles` to 2 dp via `dplyr::across(where(is.numeric), \(x) round(x, 2))` before `kable` (drop `digits = getOption("digits")`). *(user request 2, Table 1)*
- [ ] Rename the shared `width` variable — `width_cor` in `preparecorrelationHeat`, `width_overview` in `prepareoverviewHeat` — and update the `!expr` in `fig-correlation-heat` (fig-width + fig-height) and `fig-overview-heat` (fig-width). **Prerequisite for the reorder.** *(user request 1)*
- [ ] Move `prepareoverviewHeat` + `fig-overview-heat` up next to `fig-correlation-heat`, wrap the two chunks in `::: {layout="[[40,60]]"}` with per-chunk `#| out-width: "100%"`, and reduce `fig-overview-heat` `fig-height` ~10→~7; keep `fig-pairs-smooth` full width. *(user request 1)*
- [ ] Rename the mega-H1 to `# Quality Control` and introduce the static H2/H3 outline (Feature Detection / Missing Values / Abundance Distributions / Variance{CV raw, SD transformed} / Sample Structure). *(user request 3)*
- [ ] Make both variability headings static markdown under `## Variance`; remove the `cat()`-emitted heading in `variability`; gate only the body chunks (`fig-plot-distributions`, `computeCVQuantiles`, `tbl-cv-quantiles`) with `#| eval: !expr show_text`. Gate the CV and Abundance-Distributions lead-in sentences (which reference `@fig-plot-distributions`, `@tbl-cv-quantiles`, `@fig-intensity-distribution`) behind the same `show_text` condition, with a static fallback sentence each for `show_text == FALSE` so no `?@...` dangling refs appear on the pre-transformed path.
- [ ] Add an Introduction roadmap paragraph + a visible `::: {.callout-note}` for the "no automated conclusions" caveat (promote it out of the HTML comment) + forward ref to `@tbl-sample-mapping` (which requires the `tbl-sample-mapping` kable conversion below).
- [ ] Add the CV→SD→sample-size bridge paragraph in `## Variance` and the explicit hand-off sentence opening `# Sample Size Calculation`.

### P1 — consistency, robustness, and remaining storytelling

- [ ] Round `tbl-sd-quantiles` to 3 dp; round `tbl-sample-size` (`sdtrimmed` 3 dp, N columns `ceiling()`/integer) using the `across`/`mutate` idiom **inside the existing `exists("tmp") && !is.null(tmp)` guard**.
- [ ] Convert the `hierarchyCounts` totals table and `tbl-sample-mapping` from `prolfqua::table_facade(df, caption)` to `knitr::kable(df)` + `#| tbl-cap:` (+ `#| label: tbl-hierarchy-counts` on the totals table) so their `@tbl-` cross-references resolve; keep the rounding idiom on the same pipe where numeric. `#| label` alone does **not** make a `table_facade` table referenceable.
- [ ] Split the dense transformed-variability intro paragraph so each new subsection opens with a lead-in referencing only its own figures; make prose order match render order.
- [ ] De-orphan `fig-hierarchy-barplot`, `fig-upset-missing`, `fig-missing-histogram`, `fig-missingness-heatmap`, `fig-pairs-smooth`, `fig-overview-heat` with `@`-refs + one motivating sentence each; gate the `@fig-sd-vs-mean` ref behind `plot_sd_vs_mean`.
- [ ] Add a sentence linking `@tbl-hierarchy-counts`, `@fig-hierarchy-barplot`, and `@tbl-hierarchy-counts-sample`.
- [ ] Add `#| out-width: "100%"` to `fig-pairs-smooth`, `fig-overview-heat`, `fig-missingness-heatmap`.
- [ ] Normalize `computeCVQuantiles` (and the other header-arg chunks) to hash-pipe-only options; guard the CV chunk with the same gate as `fig-plot-distributions` (or `exists("stats")`).
- [ ] Set `x = "CV [%]"` in both panels of `fig-plot-distributions` and note "(in %)" / scale in the `tbl-cv-quantiles` and `tbl-sd-quantiles` captions (parallel phrasing).
- [ ] Move `library(rlang)` / `library(prolfqua)` from `hierarchyCounts` into `setup` (prerequisite for the `%||%` guard in P2).

### P2 — polish

- [ ] Change front-matter `author` to `"Functional Genomics Center Zurich"`.
- [ ] Guard empty B-fabric bullets in the Introduction with `%||% "n/a"` (only after the `library(rlang)`→`setup` move, or under an R ≥ 4.4 assumption).
- [ ] Move the power/confidence/fold-change definitions above `fig-sample-size`/`tbl-sample-size`, ideally in a `::: {.callout-note collapse="true"}`.
- [ ] Add a generic "what good QC looks like" closing paragraph at the end of `# Quality Control`.
- [ ] Add a `# Session Info` H1 (`sessionInfo()`) + closing provenance line.
- [ ] Retitle `# Appendix` (e.g. "Appendix: Sample Mapping"), put `tbl-sample-mapping` first, and decide whether `tbl-hierarchy-counts-sample` stays (as the numeric backing for `@fig-hierarchy-barplot`) or is dropped.
- [ ] Drop the empty `#| fig-alt: ""` on `fig-sample-size` (or add real alt text for family consistency).
- [ ] Optionally add a one-line comment explaining the mandatory `#| fig-dpi: 96` near the first plotly chunk.

## Open questions / decisions

- **Tabsets vs linear flow.** The SE report is fully tabsetted; this report is a top-to-bottom narrative with heavy cross-reference prose that reads better as a flat scroll. **Recommended default:** keep a flat, static H2/H3 structure and do not tabset the report. The one place a tabset could earn its keep is the three Sample Structure figures if the side-by-side layout proves fragile.
- **Side-by-side vs panel-tabset for Sample Structure.** **Recommended default:** the `layout="[[40,60]]"` side-by-side for `fig-correlation-heat` + `fig-overview-heat` (direct answer to the user), with the panel-tabset as the fallback if the trial render misbehaves with the extension's global 40% width. Decide after one trial render.
- **Rounding precision (2 vs 3 dp).** CV to 2 dp (Grp2 R6) and SD to 3 dp (SE) mixes two family precedents. **Recommended default:** keep CV at 2 dp and SD at 3 dp — they follow the closest sibling for each quantity; revisit only if the user wants a single uniform precision.
- **`show_text == FALSE` support.** **Recommended default:** gate the referencing prose (as in P0) so the pre-transformed path stays ref-clean. If the team decides the pre-transformed path is not a supported use of this vignette, the alternative is to state that explicitly and guarantee only the `show_text == TRUE` branch — but then the dense-intro split and de-orphaning must avoid referencing gated exhibits in static prose.
- **Keep `tbl-hierarchy-counts-sample`?** It duplicates `fig-hierarchy-barplot`. **Recommended default:** keep it in the appendix but frame it explicitly as the numeric backing for the barplot; drop only if the user wants a leaner report.
- **Error-handling parity (needs explicit sign-off, not a default action).** The family DEA reports degrade failing plots to a placeholder (`empty_report_plot` + `plot_or_empty`, or the SE `safe_plot`/`show_plot` set) and wrap the vsn transform / heatmap builders in `tryCatch`. Adopting this here is defensible family parity but (a) is beyond the three explicit user asks and (b) sits against the standing instruction to avoid `tryCatch`/skip-condition bandaids unless explicitly requested. **Recommended default:** do not add it in this pass; raise it as a separate, explicitly approved change.
- **Package-API items (deferred, need explicit sign-off).** Factoring `empty_report_plot`/`plot_or_empty`/`safe_plot` into `prolfquapp`; adding a `prolfquapp::example_qc_sse()` fallback to replace the hand-rolled tempdir RDS bundle; standardizing the B-fabric identifier param name (`project_conf` vs `project_info` vs `project_spec`) across the CLI callers (`CMD_QUANT_QC.R` and render callers); unifying the table renderer on `DT::datatable`/`dt_table` across the family. **Recommended default:** leave all four out of this pass and raise them separately, since each touches package API or CLI wiring and the repo rules require explicit approval. Add a `NEWS.md` bullet when the in-scope report changes land.

## Follow-up

Once this pass lands, `QC_ProteinAbundances_quarto.qmd` (the sibling QC report) may receive the analogous treatment — rounding audit, de-orphaning, section-consistency, and layout conventions — as a separate follow-up.
