# Integrate `Grp2Analysis_V2_R6_quarto.qmd` into the CMD_DEA application

> Follow-up to `TODO/TODO_rmd_2_qmd_vignettes.md` ("Scenario 2 — application / report rendering", deferred there).
> Goal: render the Quarto DEA report from the CLI pipeline (`CMD_DEA_V2.R` / `CMD_DEA_CD.R`) so a run produces the
> Quarto HTML alongside the existing outputs, with the FGCZ styling and the Find/Save toolbar.

## Status: IMPLEMENTED (Path A)

Done. Both Quarto reports live in `vignettes/` in the extension style; the SE tabbed report was moved from
`inst/templates/quarto/`. `vignettes/.install_extras` ships them (and the `_extensions/` tree) into installed `doc/`,
and `R/quarto_report_helpers.R` renders both from there (`render_quarto_se_report`, `render_quarto_dea_report`) via a
direct `quarto::quarto_render()` — no `fgcz_render`. The CLI (`CMD_DEA_V2.R`, `CMD_DEA_CD.R`, and
`write_dea_run_outputs()`) writes `DEAnalyse.rds` and renders `*_quarto_dea.html` beside the SE report.

**The `make build-vignettes` / devtools caveat is resolved.** The `.install_extras` pattern `vignettes/_extensions$`
ships the extension tree correctly under base-R `R CMD build` (it matches the full-path directory entry, so
`file.copy(recursive = TRUE)` preserves `doc/_extensions/fgczquartotemplate/...`) while staying invisible to
`devtools::build_vignettes()` (which matches only top-level basenames, so it never attempts the directory copy that
previously failed with `ENOTSUP`). No ecosystem-Makefile change was needed. A bare `_extensions` pattern ships under
base R but breaks devtools; a `_extensions/` pattern is devtools-safe but flattens/loses the tree under base R — only the
anchored full-path form satisfies both.

The rest of this document is the original design analysis, retained for context.

## The constraint that drives the design

The obvious "copy the qmd and its `_extensions/` folder into the render dir and call `quarto::quarto_render()`" plan
does **not** work at application runtime, because **there is nowhere to copy the extension from**:

- In `fgczquartotemplate`, the Quarto extension lives at the **repo root** (`_extensions/fgczquartotemplate/`), and it is
  **excluded from the built/installed package** by `.Rbuildignore` (`^_extensions$`). At runtime
  `system.file("_extensions", package = "fgczquartotemplate")` returns `""`.
- What the installed `fgczquartotemplate` *does* ship is `inst/quarto/` — the **R-helper style** assets
  (`_metadata.yml`, `fgcz.scss`, `fgcz_header_quarto.html`, `fgcz-plot-finder.html`, `template.qmd`), reachable via
  `fgczquartotemplate::fgcz_quarto_dir()`. This is what `fgcz_render()` stages.
- `prolfquapp` has no `_extensions/` under `inst/`. The vignette's copy lives in `vignettes/_extensions/`, and by
  default `vignettes/` is not installed into the library.

### Escape hatch: `.install_extras` ships the vignette's extension into `doc/` (verified)

The last point is fixable. R's `vignettes/.install_extras` file (a list of regexps matched against paths under
`vignettes/`) tells `R CMD build`/`R CMD INSTALL` to copy matching files into the installed package's `doc/` directory.
The file already exists here (it lists `bibliography.bib`); Path A adds an `_extensions` line.
**Verified on this machine:** with `_extensions` in `vignettes/.install_extras`, `R CMD build` copies the whole
nested tree into the tarball at `inst/doc/_extensions/fgczquartotemplate/` (all four files) alongside
`inst/doc/Grp2Analysis_V2_R6_quarto.qmd`. After install these are reachable at runtime via
`system.file("doc/_extensions/fgczquartotemplate", package = "prolfquapp")` and
`system.file("doc/Grp2Analysis_V2_R6_quarto.qmd", package = "prolfquapp")`.

That means the **extension style is available to the application after all** — the app copies the qmd and its
`_extensions/` out of `prolfquapp`'s own installed `doc/` (not from `fgczquartotemplate`) into a temp render dir and
renders directly. This unlocks a genuine single source of truth: the vignette file *is* the application report.

**Caveat (must be resolved with this approach):** `devtools::build_vignettes()` — which `make build-vignettes` calls —
chokes on a *directory* extra. It matches `.install_extras` against top-level names non-recursively, matches the bare
`_extensions` directory, and then `fs` errors copying a directory (`ENOTSUP`); a file-only pattern (`_extensions/`) makes
devtools match nothing. Base-R `R CMD build`/`INSTALL` (the real install path, and what `make install` uses via the
tarball) handles it correctly. So adopting `.install_extras` requires either pointing `make build-vignettes` at a
base-R build step instead of `devtools::build_vignettes()`, or copying the extra via base R. Do **not** commit
`.install_extras` until that is settled, or it breaks `make build-vignettes`.

So there are two workable application paths, detailed under Decision below: the **R-helper style** via
`fgcz_render(buttons = TRUE)` (what the existing SE report uses, no `.install_extras` needed), or the **extension style**
via `.install_extras` + a direct `quarto_render` (true single source of truth, with the caveat above).

## What already exists (reuse it)

`CMD_DEA_V2.R:259` and `CMD_DEA_CD.R:321` already call `prolfquapp:::render_quarto_se_report()`
(`R/quarto_report_helpers.R`), driven from `write_dea_run_outputs()` in `R/cmd_helpers.R:484-504`:

```r
se_file <- file.path(reporter$resultdir, "SummarizedExperiment.rds")
SE <- reporter$make_SummarizedExperiment()
saveRDS(SE, file = se_file)
outdir$quarto_file <- render_quarto_se_report(
  se_file = se_file,
  output_dir = reporter$resultdir,
  output_file = paste0(reporter$fname, "_quarto.html")
)
```

`render_quarto_se_report()` copies a template `.qmd` into a temp dir, `setwd()`s there, and calls
`fgczquartotemplate::fgcz_render(template, buttons = TRUE, execute_params = list(se_file = ...))`, then copies the HTML
to the results dir. Its template (`inst/templates/quarto/Grp2Analysis_V2_SE_tabset.qmd`) is **R-helper style**: no
`format:` line, styling supplied by the staged `_metadata.yml`.

The DEA Quarto report is the same pattern with a different serialized input.

## Decision — two viable paths

**Path A (recommended) — extension style, single source of truth, via `.install_extras`.**
Ship the vignette's `_extensions/` and the qmd into installed `doc/` (the verified `.install_extras` hack above). The app
copies `system.file("doc/Grp2Analysis_V2_R6_quarto.qmd", ...)` and `system.file("doc/_extensions", ...)` into a temp
render dir and calls `quarto::quarto_render(input, execute_params = list(deanalyse_file = <path>))` directly. The
Find/Save toolbar and the `toc-location`/`fig-dpi` tweaks are already in the vignette front matter, so **no `fgcz_render`
and no `buttons` argument** — the one file is both vignette and app report. Requires resolving the
`make build-vignettes`/devtools caveat.

**Path B (fallback) — R-helper style via `fgcz_render(buttons = TRUE)`.**
Maintain a separate app template `inst/templates/quarto/Grp2Analysis_V2_R6_dea.qmd` with R-helper front matter (no
`format:`), rendered exactly like `render_quarto_se_report()` but with `execute_params = list(deanalyse_file = ...)`. No
`.install_extras`, no devtools caveat; assets come from the installed `fgczquartotemplate/inst/quarto/`. Here
`buttons = TRUE` is the correct and only way to get the toolbar (it merges `include-after-body` pointing at the staged
`fgcz-plot-finder.html`). Downside: two front-matter variants of the same body → the single-source concern below.

### Single-source-of-truth (only relevant to Path B)

The report **body** (~660 lines) should have one source. Under Path A this is automatic. Under Path B, only the front
matter differs, so keep the body shared:

1. **Shared included body.** Extract the body into `_grp2_body.qmd` (no front matter); two thin wrappers
   `{{< include _grp2_body.qmd >}}` — the vignette (extension front matter) and the app template (R-helper front matter).
   Keep the vignette and `inst/` copies of the body in sync via a check (see Open questions).
2. **Two full copies** of the qmd. Simplest to wire, worst for drift. Only if 1 is rejected.

Under Path B, the app-facing front matter must reproduce the vignette tweaks: `toc-location: left`, the wider sidebar,
and `#| fig-dpi: 96` on the three plotly chunks (`fig-normalized`, `fig-pca`, `fig-volcano`) — the last is essential, or
the global 300 dpi inflates the plotly widgets to ~2100px. The toolbar comes from `buttons = TRUE`, so the app front
matter should **not** hard-code `include-after-body`.

## Data plumbing

The Rmd path passes a live object (`rmarkdown::render(params = list(deanalyse = self$deanalyse))`,
`R6_DEAReportGenerator.R:366`). The QMD consumes a file (`dea <- readRDS(deanalyse_file)`), so:

- [ ] Serialize the completed `DEAnalyse` object to `.rds` in the results dir, beside `SummarizedExperiment.rds`
      (e.g. `DEAnalyse.rds`). Note this is the **DEAnalyse R6**, not the SummarizedExperiment — a new save step. In
      `write_dea_run_outputs()` the object is in scope as `deanalyse` (`cmd_helpers.R:449`).
- [ ] (Path A only) Add an `_extensions` line to the existing `vignettes/.install_extras` (it already lists
      `bibliography.bib`), and resolve the `make build-vignettes` caveat (point that target at a base-R build or copy
      extras via base R rather than `devtools::build_vignettes()`).
- [ ] Add `render_quarto_dea_report()` in `R/quarto_report_helpers.R`. Path A: copy the qmd + `_extensions/` out of
      `system.file("doc", package = "prolfquapp")` into a temp dir and `quarto::quarto_render(execute_params =
      list(deanalyse_file = <path>))` (no `buttons`). Path B: near-copy of `render_quarto_se_report()` rendering the app
      template with `execute_params = list(deanalyse_file = <path>)` and `buttons = TRUE`.
- [ ] Guard on the Quarto CLI with the existing `Sys.which("quarto")` warn/skip already in `render_quarto_se_report()`.
- [ ] Wire the call into `write_dea_run_outputs()` (`cmd_helpers.R:484-504`), inside the same `tryCatch` that already
      guards the SE/Quarto export, so a missing Quarto CLI or render error degrades gracefully.

## Report selection (Rmd legacy vs QMD)

- [ ] Let callers choose the DEA report path — Rmd (legacy, `render_DEA`) or QMD (new). A config/CLI flag threaded from
      `CMD_DEA_V2.R` / `CMD_DEA_CD.R`, defaulting to keep current behavior until the QMD path is validated. Keep the Rmd
      render code untouched (the migration TODO's invariant).

## Out of scope / do not touch in a first pass

- The vignette's report body itself — Path A reuses it verbatim; Path B only adds a wrapper/copy.
- The legacy Rmd rendering path (`render_DEA`, `Grp2Analysis_V2_R6.Rmd`).

## Open questions

- **`make build-vignettes` under `.install_extras` (Path A).** Decide how to keep that target working given
  `devtools::build_vignettes()` fails on the directory extra: switch it to a base-R build/`pkgbuild` step, or copy the
  extra via base R. This is an ecosystem-Makefile coordination item, not just a prolfquapp one.
- **Body sync guard (Path B).** If the shared-include option is used, the body must be identical in `vignettes/` and
  `inst/`. Add a test comparing the two, or generate one from the other at build time.
- **Does `make_SummarizedExperiment()`'s data suffice instead?** The SE report already exists. Confirm the DEA report
  genuinely needs the full `DEAnalyse` object (it uses `$prolfq_app_config`, `$lfq_data_raw`, `$lfq_data`,
  `$annotated_contrasts`, `$contrast_results`, `$default_model`, etc.) and cannot be reconstructed from the SE — it
  cannot, so a separate `DEAnalyse.rds` is required.
- **Reproducibility.** For snapshot/output stability, add `set.seed()` inside `example_deanalyse()` (package-side), as
  noted in the migration TODO — only relevant if we snapshot the app output.

## Validation

- [ ] Run `CMD_DEA_V2.R` end-to-end on a test dataset; confirm `DEAnalyse.rds` and the `*_quarto.html` (DEA) appear in
      the results dir with FGCZ styling, the Find/Save toolbar, resolved cross-references, and sane plotly heights.
- [ ] Confirm graceful skip (logged warning, no run failure) when the Quarto CLI is absent.
- [ ] Confirm the legacy Rmd report still renders when the selector picks it.
- [ ] `make check-fast` / package tests still pass.
