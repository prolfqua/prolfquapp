# TODO: prolfquapp/Dockerfile — vignette + Quarto build-time checks

## Scope (final, narrowed)

Working on **`prolfquapp/Dockerfile`** — the deployment image
`docker.io/prolfqua/prolfquapp:latest` consumed by
[`integration_test/prolfquapp_docker.sh`](https://prolfqua.github.io/integration_test/prolfquapp_docker.sh).

The image is **almost fine** as-is. The one thing to add: **build-time
checks that prove every vignette can actually be re-rendered with
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
given the packages baked into the image**, plus the same end-to-end
check for the Quarto SE template. The user-described verification path
is “go into the image, run
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
on every vignette — this is the fast smoke test that proves vignette
completeness.” Step 1 below executes exactly that path at build time so
a broken image never ships.

Out of scope for this round:

- Quarto stays (qmd template under `inst/templates/quarto/` is rendered
  at analysis time).
- TinyTeX stays.
- No final-stage / source-baking changes.
- No slimming of the `build` stage `install.packages` line.
- Root `prolfqua_fml/Dockerfile` edits made earlier this session are
  unrelated to this TODO and untouched.

Folded in (was parked, now in scope):

- **Drop `prolfquadata` from the
  [`pak::pkg_install`](https://pak.r-lib.org/reference/pkg_install.html)
  line.** Verified unused: not in `DESCRIPTION`, not referenced from
  `R/`, `vignettes/`, or `tests/`. Only hits are
  [`inst/samples/maxquant_txt/tinydata.R:3`](https://prolfqua.github.io/prolfquapp/inst/samples/maxquant_txt/tinydata.R)
  (a sample helper script — not run by build/install/tests/vignettes)
  and
  [`inst/poster/prolfqua.bib:799`](https://prolfqua.github.io/prolfquapp/inst/poster/prolfqua.bib)
  (a citation entry, text only). Safe to remove from the build-stage
  [`pak::pkg_install`](https://pak.r-lib.org/reference/pkg_install.html)
  call.

## Current `prolfquapp/Dockerfile` (the bits that matter)

Multi-stage build:

1.  **`base`** — `r-base:4.5.2` + `pandoc`, `gdebi`, Quarto.
2.  **`build`** — dev libs, `R_LIBS_USER=/opt/r-libs-site`, then
    upstream deps via
    [`pak::pkg_install`](https://pak.r-lib.org/reference/pkg_install.html),
    then `COPY ./DESCRIPTION` +
    [`pak::local_install_deps`](https://pak.r-lib.org/reference/local_install_deps.html),
    then `COPY . /opt/prolfqua`,
    `install.packages(knitr, rmarkdown, DT, gridExtra, KernSmooth, plotly)`,
    then `R CMD build` (renders Rmd vignettes) + `R CMD INSTALL`.
    **Existing lines 39–40 do `stopifnot(file.exists(... .Rmd ...))`
    only — they check the Rmd *source* survived install, not that any
    HTML was rendered, and not that re-rendering works.** That is what
    we are replacing.
3.  **Final stage** — TinyTeX + `COPY --from=build /opt/r-libs-site`,
    `PATH=/opt/r-libs-site/prolfquapp/application/bin:...`,
    `ENTRYPOINT ["/bin/bash"]`. Untouched in this round.

## Vignettes that must render

All six Rmd files under
[`vignettes/`](https://prolfqua.github.io/prolfquapp/vignettes/) have
proper `%\VignetteIndexEntry` + `%\VignetteEngine{knitr::rmarkdown}`
declarations and are rendered by `R CMD build`:

- `Auxiliary_ExDesignSurvey.Rmd` (standalone)
- `DiffExpQC_R6.Rmd` (uses `params:` with defaults)
- `Grp2Analysis_V2_R6.Rmd` (uses `params:` with defaults)
- `QC_ProteinAbundances.Rmd` (uses `params:` with defaults)
- `QCandSSE.Rmd` (uses `params:` with defaults)
- `prolfquapp.Rmd` (standalone)

After `R CMD INSTALL`, each lands at
`system.file("doc", "<name>.Rmd"/".html", package = "prolfquapp")`. The
Rmd source survives in the installed `doc/`, so
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
against it works **without** needing the original `/opt/prolfqua` source
tree — i.e. the check is robust across the build → final stage
transition.

The **Quarto template** at
[`inst/templates/quarto/Grp2Analysis_V2_SE.qmd`](https://prolfqua.github.io/prolfquapp/inst/templates/quarto/Grp2Analysis_V2_SE.qmd)
is rendered at *analysis* time (not by `R CMD build`) via
`prolfquapp:::render_quarto_se_report()`
([R/quarto_report_helpers.R:3](https://prolfqua.github.io/prolfquapp/R/quarto_report_helpers.R)).
It needs: the `quarto` CLI on `PATH`, the support files in
`inst/templates/quarto/` (`_fgcz-report.yml`, `fgcz_header_quarto.html`,
the template), and the R packages listed in its setup chunk: `DT`,
`dplyr`, `ggplot2`, `grid`, `gridExtra`, `prolfqua`, `prolfquapp`,
`tibble`, `UpSetR`. The bundled fixture
[`inst/extdata/3106962.rds`](https://prolfqua.github.io/prolfquapp/inst/extdata/3106962.rds)
is the default `se_file` and is the build-time render target.

## Plan

Two new `RUN` steps in the **`build` stage**, replacing the existing
[`Dockerfile`](https://prolfqua.github.io/prolfquapp/Dockerfile) lines
39–40, plus a one-line edit on the
[`pak::pkg_install`](https://pak.r-lib.org/reference/pkg_install.html)
call to drop `prolfquadata`. The new steps sit after `R CMD INSTALL` and
before the `data.table` smoke check (which we leave alone).

### Step 0 — Drop unused `prolfquadata` from upstream deps

In [`Dockerfile`](https://prolfqua.github.io/prolfquapp/Dockerfile) line
31, remove `"git::https://gitlab.bfabric.org/wolski/prolfquadata.git"`
from the `pak::pkg_install(c(...))` vector. The line becomes:

``` dockerfile
RUN R -e 'options(warn=2); pak::pkg_install(c("any::seqinr", "any::prozor", "any::logger", "any::lubridate", "github::fgcz/prolfqua", "github::prolfqua/prolfquasaint"))'
```

Justification in the “Folded in” bullet above.

### Step 1 — Re-render every Rmd vignette with `rmarkdown::render()`

This is the “completeness” / fast-smoke test the user described: for
every `.Rmd` in the installed `doc/`, run
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
to a tmp output and assert the HTML appears. Failure surfaces the exact
vignette + the underlying rmarkdown error.

``` dockerfile
RUN Rscript -e ' \
  doc <- system.file("doc", package = "prolfquapp", mustWork = TRUE); \
  rmds <- list.files(doc, pattern = "\\.Rmd$", full.names = TRUE); \
  if (length(rmds) == 0) stop("No Rmd vignettes installed under doc/"); \
  out_dir <- tempfile("vignette_check_"); dir.create(out_dir, recursive = TRUE); \
  for (rmd in rmds) { \
    cat("Rendering: ", basename(rmd), "\n", sep = ""); \
    out <- rmarkdown::render(rmd, output_dir = out_dir, quiet = TRUE, envir = new.env()); \
    stopifnot(file.exists(out)); \
  } \
  cat("All ", length(rmds), " Rmd vignettes rendered via rmarkdown::render().\n", sep = "") \
  '
```

Notes: - `envir = new.env()` keeps each render isolated. - The four
`params:`-driven vignettes have default values in their YAML headers;
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
uses those defaults when no `params` argument is passed. - `output_dir`
is a tempdir away from the installed `doc/`, so we don’t perturb the
image state.

### Step 2 — Render the Quarto SE template end-to-end

Render the Quarto template against the bundled `3106962.rds` to verify
the full Quarto toolchain + every R package the template needs is
present and works end-to-end. `render_quarto_se_report()` already
handles `Sys.which("quarto")`, support-file copying, `quarto render`,
and an expected-HTML existence check internally; this step just calls it
and re-asserts the output.

``` dockerfile
RUN Rscript -e ' \
  out <- tempfile("qmd_check_"); dir.create(out, recursive = TRUE); \
  se  <- system.file("extdata", "3106962.rds", package = "prolfquapp", mustWork = TRUE); \
  prolfquapp:::render_quarto_se_report(se_file = se, output_dir = out); \
  stopifnot(file.exists(file.path(out, "Grp2Analysis_V2_SE.html"))); \
  cat("Quarto SE report rendered OK ->", file.path(out, "Grp2Analysis_V2_SE.html"), "\n") \
  '
```

Note: `render_quarto_se_report()` issues
`logger::log_warn(...); return(NULL)` if `Sys.which("quarto")` is empty.
The post-call `stopifnot(file.exists(...))` still fails in that case, so
a missing or broken Quarto CLI is caught.

## Files touched

- [`prolfquapp/Dockerfile`](https://prolfqua.github.io/prolfquapp/Dockerfile)
  - Line 31: remove `prolfquadata` from the
    [`pak::pkg_install`](https://pak.r-lib.org/reference/pkg_install.html)
    vector.
  - Lines 39–40: replace with the two new `RUN` blocks (Steps 1 and 2).

No Makefile change. No README change. No final-stage change.

## Verification

1.  **Build succeeds and prints the new check banners.**

    ``` bash
    docker build -t prolfqua/prolfquapp:dev -f prolfquapp/Dockerfile prolfquapp/
    ```

    Build log must contain:

    - `Rendering: <each>.Rmd` lines +
      `All 6 Rmd vignettes rendered via rmarkdown::render().`
    - `Quarto SE report rendered OK -> ...`.

2.  **CLI contract intact** (PATH + ENTRYPOINT for
    `prolfquapp_docker.sh`):

    ``` bash
    docker run --rm prolfqua/prolfquapp:dev prolfqua_dea.sh --help     # exit 0
    ```

3.  **Interactive fast smoke test** (the user-described path — Step 1
    reproduced by hand):

    ``` bash
    docker run --rm -it --entrypoint bash prolfqua/prolfquapp:dev
    # inside the container:
    Rscript -e ' \
      doc <- system.file("doc", package="prolfquapp"); \
      for (f in list.files(doc, pattern="\\.Rmd$", full.names=TRUE)) { \
        cat(basename(f), "\n"); \
        rmarkdown::render(f, output_dir=tempdir(), quiet=TRUE, envir=new.env()); \
      }'
    ```

    Plus the Quarto template:

    ``` bash
    Rscript -e 'prolfquapp:::render_quarto_se_report( \
      se_file    = system.file("extdata","3106962.rds", package="prolfquapp"), \
      output_dir = "/tmp/qmd_out")'
    ```

4.  **Real DEA against fixture data** (the user marked this “YES
    absolutely essential” — *not* covered by the build-time checks, must
    be run separately via `prolfquapp_docker.sh`):

    ``` bash
    cd integration_test
    prolfquapp_docker.sh \
      --image-repo prolfqua/prolfquapp \
      --image-version dev \
      -- prolfqua_dea.sh <args for one of the fixture datasets>
    ```

    See `integration_test/CLAUDE.md` for the exact fixture invocation
    and the existing `save-references-docker` / `test-dea-regression`
    flow.

Failure modes the build-time checks must catch:

- Any Rmd vignette that
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
  cannot complete → step 1 fails with the rmarkdown error and the
  offending vignette name.
- Quarto CLI missing or broken → step 2 fails on the post-call
  `file.exists`.
- Any R package required by an Rmd or the Quarto template setup chunk
  that isn’t in the image → corresponding step fails inside the render
  call.
