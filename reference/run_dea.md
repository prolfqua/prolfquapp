# Run differential expression analysis pipeline

Reads annotation, preprocesses quantification data, then runs the full
DEA pipeline: aggregation, transformation, model fitting, and contrast
computation. Returns the
[`DEAnalyse`](https://prolfqua.github.io/prolfquapp/reference/DEAnalyse.md)
object and supporting data needed for report generation.

## Usage

``` r
run_dea(indir, dataset, software, config)
```

## Arguments

- indir:

  directory containing quantification output files

- dataset:

  path to annotation CSV/TSV/XLSX

- software:

  software key as returned by
  [`get_procfuncs`](https://prolfqua.github.io/prolfquapp/reference/get_procfuncs.md)
  (e.g. "prolfquapp.DIANN")

- config:

  a
  [`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md)
  object (already merged with CLI overrides via
  [`sync_opt_config`](https://prolfqua.github.io/prolfquapp/reference/sync_opt_config.md))

## Value

list with `deanalyse` (DEAnalyse), `xd` (preprocessed data),
`annotation` (parsed annotation including contrasts), and `files`
(discovered file paths)
