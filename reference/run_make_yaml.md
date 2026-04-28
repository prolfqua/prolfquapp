# Generate a DEA configuration list

Creates a
[`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md)
object, converts it to a plain list, and reorders fields so that
verbose/internal sections appear at the bottom of the YAML output.

## Usage

``` r
run_make_yaml(
  project = "",
  order = "",
  workunit = "",
  norm = "vsn",
  model = "lm_missing",
  outdir = NULL
)
```

## Arguments

- project:

  project ID

- order:

  order ID

- workunit:

  workunit ID

- norm:

  normalization method (e.g. "vsn", "none", "robscale")

- model:

  contrast facade method (see `names(prolfqua::FACADE_REGISTRY)`)

- outdir:

  optional output directory; if it exists, stored in the config so
  downstream scripts know where to write results

## Value

named list suitable for
[`yaml::write_yaml()`](https://yaml.r-lib.org/reference/write_yaml.html)

## Examples

``` r
cfg <- run_make_yaml(project = "p100", workunit = "WU123")
cfg$project_spec$workunit_Id
#> [1] "WU123"
```
