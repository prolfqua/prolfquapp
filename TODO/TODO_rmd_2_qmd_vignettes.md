# Migrate `Grp2Analysis_V2_R6.Rmd` to Quarto (`.qmd`)

> Establish a repeatable migration pattern for params-driven R Markdown reports by producing a Quarto version of
> `Grp2Analysis_V2_R6.Rmd` — the largest and most reference-heavy report in the package, and therefore the one that
> proves the pattern.

## Goal and scope

Add `prolfquapp/vignettes/Grp2Analysis_V2_R6.qmd` as the single source of truth for the Quarto report path, built through
the standard package vignette machinery. The work is **additive**: the existing Rmd and its rendering code are left
untouched, and the QMD is not wired into `CMD_DEA_V2.R` or any `inst/application/` script in this pass.

The two paths differ only in how the analysis object arrives:

- The **Rmd** keeps passing a live R object through `rmarkdown::render(params = list(deanalyse = ...))` and resolves it
  with `dea <- eval(params$deanalyse)`. This is intentional and stays as-is.
- The **QMD** takes a file parameter `deanalyse_file` pointing at a serialized `.rds`, and resolves it with
  `dea <- readRDS(deanalyse_file)`. External callers pass a path, never a live object.

Everything downstream of the loaded `dea` object — structure, prose, citations, figure/table intent, analysis
semantics — stays as analogous to the Rmd as Quarto syntax allows.

If a genuine bug turns up in the existing Rmd rendering path, report it and wait for explicit permission before changing
any Rmd rendering code.

### Verified against the codebase

These assumptions were checked before writing the plan, so the implementation can rely on them:

- `prolfquapp::example_deanalyse()` exists and is exported (`R/example_deanalyse.R`; signature `example_deanalyse(Nprot = 100)`).
- The `fgczquartotemplate-html` **Quarto extension** exists at `fgczquartotemplate/_extensions/fgczquartotemplate/` and
  contains exactly `_extension.yml`, `fgcz.scss`, `fgcz_header_quarto.html`, and `fgcz-plot-finder.html`. The extension
  already sets `embed-resources`/`self-contained`, `code-fold`, `execute: {echo, warning, message} = false`, and a
  `crossref` block — matching the Rmd's `self_contained: true` and global `echo = FALSE`.
- The installed `quarto` R package registers the vignette engine `quarto::format` for `*.qmd` (confirmed on this
  machine, quarto R 1.5.1). The extension declares `quarto-required: ">=1.4.0"`; the local Quarto CLI is 1.8.25.

## Environment prerequisites

Building a Quarto vignette needs the **Quarto CLI on `PATH`** (>= 1.4.0), not just the `quarto` R package. The `quarto`
R package only provides the R-side engine; it shells out to the binary.

Behavior when the binary is absent, which matters for CI:

- Under `R CMD check` / `R CMD build`: the engine emits an informational message and produces an **empty** vignette —
  the check does **not** fail. This means a quarto-less CI will *silently* ship an empty `Grp2Analysis_V2_R6.html`.
- Under a direct `devtools::build_vignettes()`: the engine **aborts** with an error.

Action: ensure the ecosystem/CI environment installs Quarto so the vignette is actually rendered rather than silently
emptied. This is an ecosystem-Makefile / CI concern, not a `DESCRIPTION` one.

## Design

### Input contract

The report body works from a single loaded object:

```r
dea <- readRDS(deanalyse_file)
```

The invariant across both scenarios below is: report code consumes a loaded `dea`, and callers hand it a file path.

### Scenario 1 — package documentation vignette

When no `deanalyse_file` is supplied (the vignette-build default), the QMD creates the example object once, caches it to
an `.rds` under `tempdir()`, and then loads it through the same `readRDS()` path as any other render. Writing to
`tempdir()` (not the package/working directory) keeps `R CMD check` clean.

```r
deanalyse_file <- params$deanalyse_file
if (is.null(deanalyse_file) || !nzchar(deanalyse_file)) {
  deanalyse_file <- file.path(tempdir(), "prolfquapp-example-deanalyse.rds")
  if (!file.exists(deanalyse_file)) {
    saveRDS(prolfquapp::example_deanalyse(), deanalyse_file)
  }
}
dea <- readRDS(deanalyse_file)
```

Note on determinism: the Rmd regenerates a fresh example on every render, so its output is not byte-stable run-to-run
either. That is fine for a vignette — output stability is **not** an invariant here. If we ever want the example report
to be reproducible (e.g. for snapshot tests), the clean fix is to add a `set.seed()` inside `example_deanalyse()` so the
generated object is constant; do that in the package, not in the report.

### Scenario 2 — application / report rendering (deferred)

The application scripts already build a completed `DEAnalyse` object. When QMD rendering is later integrated, that path
should: save the real `DEAnalyse` object beside the report outputs, pass its path as `deanalyse_file`, render the same
`Grp2Analysis_V2_R6.qmd` (in place or as a staged copy), and keep the Rmd path available as a selectable legacy option.

This is explicitly out of scope for the first pass — see [Follow-up](#follow-up-application-integration).

### FGCZ Quarto template deployment — extension style

`fgczquartotemplate` supports two deployment styles: an R-helper style (`fgcz_render()` stages `_metadata.yml`,
`fgcz.scss`, `fgcz_header_quarto.html`, `fgcz-plot-finder.html` next to the QMD, then calls `quarto::quarto_render()`),
and a Quarto **extension** style (`format: fgczquartotemplate-html`, resolved from `_extensions/fgczquartotemplate/`).

Use the **extension style** for the vignette. The Quarto vignette engine calls Quarto directly and cannot run an
external `fgcz_render()` first to stage assets, so the styling must be self-contained in a vendored extension. The
target layout:

```text
prolfquapp/vignettes/
  Grp2Analysis_V2_R6.qmd
  bibliography.bib               # already present; reused
  _extensions/
    fgczquartotemplate/
      _extension.yml
      fgcz.scss
      fgcz_header_quarto.html
      fgcz-plot-finder.html
```

These four files are **vendored** copies of the `fgczquartotemplate` extension, treated as generated assets, not hand-
edited in `prolfquapp`. (`_metadata.yml` belongs to the R-helper style and is deliberately *not* part of the extension
directory.)

### Front matter

```yaml
---
title: "Differential Expression Analysis."
format: fgczquartotemplate-html
params:
  deanalyse_file: null
vignette: >
  %\VignetteIndexEntry{Differential Expression Analysis}
  %\VignetteEngine{quarto::format}
  %\VignetteEncoding{UTF-8}
bibliography: bibliography.bib
---
```

`quarto::format` is the correct engine: it honors the declared `format:` (here `fgczquartotemplate-html`). The
alternative `quarto::html` would override the format and defeat the FGCZ styling, so do not use it. The Rmd's LaTeX
`fancyhdr`/`header-includes` block has no effect on HTML output and is simply dropped.

### Cross-reference and caption conversion — the hard part

This is where most of the effort and risk lives; "convert the references" is not a one-liner. The Rmd relies on four
bookdown-specific mechanisms that Quarto does **not** support 1:1, and each fails *silently* (dangling `?@ref`
placeholders) rather than erroring:

1. **`(ref:label)` text references — 10 of them.** Pure bookdown; **no Quarto equivalent**. Each caption must move
   inline into a `#| fig-cap:` chunk option. Several captions embed inline R (e.g. the volcano-plot caption references
   `cfg$processing_options$FDR_threshold`/`diff_threshold`); express those with `#| fig-cap: !expr <r-expression>` or
   build the caption string on a preceding line.
2. **Chunk labels must be renamed to Quarto's reserved prefixes.** A referenceable figure needs a `fig-` label, a table
   needs a `tbl-`. Rename via `#| label: fig-...` / `#| label: tbl-...` (hyphenated; Quarto dislikes camelCase) and
   reference with `@fig-...` / `@tbl-...`.
3. **`knitr::kable(caption = ...)` does not yield a `tbl-` cross-reference in Quarto.** Set the caption via the
   `#| tbl-cap:` chunk option (paired with a `#| label: tbl-...`), not the `kable()` argument.
4. **Drop the literal "Figure"/"Table"/"Plot" word before a reference.** In Quarto `@fig-x` already renders as
   "Figure 1"; writing "Figure @fig-x" produces "Figure Figure 1". The Rmd writes "Figure \@ref(fig:x)",
   "Table \@ref(tab:x)", and once "Plot \@ref(fig:pca)" — all become just `@fig-x` / `@tbl-x`.

**Known risk — DT datatables referenced as figures.** Two results tables are `DT::datatable` htmlwidgets referenced as
`Figure \@ref(fig:tableAllProt)` and `Figure \@ref(fig:SigPrey)`. Quarto's cross-reference machinery does not reliably
number/link htmlwidgets as `fig-`. Expect these to need either dropping the cross-reference (refer to "the table below")
or a wrapping workaround; treat as an explicit test case, not an assumed success.

**Conditionally-evaluated chunks.** `fig:SigPrey`, `fig:sigroteins`, and `fig:vennDiagramSig` run under
`eval = showSignificant`, and their reference text is emitted from inside `cat(...)`. Verify Quarto resolves a
cross-reference whose target chunk did not evaluate (dangling ref) and whose reference is produced by `cat()` into the
knitted markdown.

**Minor porting fix.** The `nrPerSample`/`nrPerSample2` chunks use `fig.with=10` — a typo for `fig.width` that knitr
silently ignores. Use the correct `#| fig-width: 10` in the QMD.

#### Figure / table inventory (bookdown → Quarto)

Full conversion checklist. Line numbers reference `vignettes/Grp2Analysis_V2_R6.Rmd`.

| Bookdown label | Quarto label | Kind | Chunk / ref lines | Notes |
| --- | --- | --- | --- | --- |
| `tab:samples` | `tbl-samples` | kable | ~L114 / L108 | `#| tbl-cap:` + label |
| `tab:annotation` | `tbl-annotation` | kable | ~L134 / L108 | `#| tbl-cap:` + label |
| `tab:CVtable` | `tbl-cv` | kable | ~L308 / L301 | `#| tbl-cap:` + label |
| `tab:contrtable` | `tbl-contrasts` | kable | ~L350 / L345 | `#| tbl-cap:` + label |
| `tab:nrsignificant` | `tbl-nr-significant` | kable | ~L458 / L437 | `#| tbl-cap:` + label |
| `fig:nrPerSample` | `fig-nr-per-sample` | ggplot | L157 / L155 | fix `fig.with`→`fig-width` |
| `fig:nrPerSample2` | `fig-nr-per-sample-2` | ggplot | L173 / L171 | fix `fig.with`→`fig-width` |
| `fig:naHeat` | `fig-na-heat` | plot | L213 / L207 | `(ref:naHeat)` L211 → `fig-cap` |
| `fig:vennProteins` | `fig-venn-proteins` | plot | L224 / L220 | `(ref:vennProteins)` L222 |
| `fig:normalized` | `fig-normalized` | plot | L260 / L233 | `(ref:normalized)` L258 |
| `fig:SDViolin` | `fig-sd-violin` | plot | L294 / L290 | `(ref:SDViolin)` L292 |
| `fig:heatmap` | `fig-heatmap` | plot | L327 / L320 | `(ref:heatmap)` L323 |
| `fig:pca` | `fig-pca` | plot | L335 / L331 | ref says "Plot"; `(ref:pca)` L333 |
| `fig:volcanoplot` | `fig-volcano` | plot | L418 / L377 | `(ref:volcanoplot)` L415 has inline R |
| `fig:sigroteins` | `fig-sig-proteins` | plot | L527 / L518 | `eval=showSignificant`; ref via `cat()` |
| `fig:vennDiagramSig` | `fig-venn-sig` | plot | L534 / — | `eval=showSignificant`; `(ref)` L531 |
| `fig:tableAllProt` | `fig-table-all` | **DT** | L394 / L358,L363 | htmlwidget crossref risk |
| `fig:SigPrey` | `fig-sig-prey` | **DT** | L481 / L475 | htmlwidget crossref risk; `eval=showSignificant` |

The report's first setup chunk keeps `cfg <- dea$prolfq_app_config`, `is_saint_report <- identical(dea$default_model,
"saint")`, and the `empty_report_plot`/`plot_or_empty` helpers unchanged — only the `dea <- eval(params$deanalyse)` line
is replaced by the `deanalyse_file` resolution above.

## Package configuration changes (`DESCRIPTION`)

- Add `quarto` to **`Suggests:`**.
- Change `VignetteBuilder: knitr` to `VignetteBuilder: knitr, quarto` — **additive**. Keep `knitr`; the six existing
  Rmd vignettes still build through it. Dropping `knitr` would break them.
- `fgczquartotemplate` is already in `Imports:` and `Remotes:` — leave it. Because the extension is vendored, it is not
  strictly required at vignette-build time, but it remains needed for the optional sync check below and is already a
  dependency.

## Implementation plan

- [ ] Vendor the `fgczquartotemplate` extension into `prolfquapp/vignettes/_extensions/fgczquartotemplate/` (the four
      files listed above), copied from the local `fgczquartotemplate` package.
- [ ] Add `vignettes/Grp2Analysis_V2_R6.qmd`, starting from the Rmd structure and prose.
- [ ] Replace the object-resolution line with the `deanalyse_file` setup handling (Scenario 1 fallback).
- [ ] Convert front matter to the Quarto/extension form (format, params, vignette engine `quarto::format`).
- [ ] Perform the cross-reference/caption conversion using the inventory table:
  - [ ] chunk option syntax (`#|` cell options, `fig-width`, etc.);
  - [ ] rename chunk labels to `fig-*` / `tbl-*`;
  - [ ] inline all 10 `(ref:...)` captions into `#| fig-cap:` (use `!expr` where inline R is needed);
  - [ ] convert `\@ref(fig:...)`/`\@ref(tab:...)` to `@fig-*`/`@tbl-*` and drop the literal "Figure"/"Table"/"Plot";
  - [ ] `knitr::kable` captions moved to `#| tbl-cap:`;
  - [ ] decide and test the DT-datatable references (`fig-table-all`, `fig-sig-prey`).
- [ ] Preserve report prose and analysis chunks otherwise unchanged.
- [ ] Point the QMD at `vignettes/bibliography.bib` (already present).
- [ ] Update `DESCRIPTION` (`Suggests: quarto`; `VignetteBuilder: knitr, quarto`).
- [ ] Confirm `R CMD build` keeps `vignettes/_extensions/fgczquartotemplate/` in the tarball and `.Rbuildignore` does
      not exclude it (inspect the built tarball for `_extension.yml`).
- [ ] Verify the existing Rmd vignette path still builds unchanged.

## Validation

Direct render smoke test from the package source tree (fast iteration; requires the Quarto CLI):

```bash
quarto render vignettes/Grp2Analysis_V2_R6.qmd
```

Full package vignette build — the real acceptance gate:

```bash
Rscript -e "devtools::build_vignettes()"
```

Explicit `deanalyse_file` path (Scenario-2 rehearsal), rendering from a real `.rds`:

```bash
Rscript -e 'saveRDS(prolfquapp::example_deanalyse(), "/tmp/dea.rds")'
quarto render vignettes/Grp2Analysis_V2_R6.qmd -P deanalyse_file:/tmp/dea.rds
```

Checks to eyeball in the rendered HTML: every figure/table is numbered and every `@ref` resolves (search the output for
`?@` to catch dangling references), the FGCZ styling/header is applied, and the two DT tables render and are referenced
sensibly.

## Follow-up: application integration

Deferred to a later task — do **not** touch `CMD_DEA_V2.R`, `CMD_DEA_CD.R`, or `write_dea_run_outputs()` in this pass.
The later design must let callers select the Rmd (legacy) or QMD (new) report path. Recommended placement: write the
real `DEAnalyse.rds` into the results directory **beside `SummarizedExperiment.rds`**, so it is archived alongside the
other RDS/parquet outputs and passed as `deanalyse_file`.

## Open questions

- **Sync guard for the vendored extension.** A test should keep `vignettes/_extensions/fgczquartotemplate/` from
  drifting. Reuse the existing `fgczquartotemplate/data-raw/sync_assets.R` pattern: compare the three shared assets
  (`fgcz.scss`, `fgcz_header_quarto.html`, `fgcz-plot-finder.html`) against the installed copies at
  `system.file("quarto/<file>", package = "fgczquartotemplate")`. Caveat: `_extension.yml` is the extension manifest and
  has no `inst/quarto/` counterpart, so it cannot be covered by that comparison — decide whether to check it separately
  or accept it as unguarded.
- **Plot-finder toolbar.** Recommend keeping it **disabled** for vignette builds. The `_extension.yml` already makes it
  opt-in (`include-after-body`), and a search/download toolbar is application-report UX, not documentation UX. Enable it
  only on the Scenario-2 application path if wanted.
