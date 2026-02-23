# Generate zip directory name based on project information and date

Generate zip directory name based on project information and date

## Usage

``` r
zipdir_name(
  prefix = "DEA",
  project_id = "",
  order_id = "",
  workunit_id = "",
  transform = "vsn",
  date = Sys.Date()
)
```

## Arguments

- prefix:

  Analysis prefix (e.g., "DEA" or "QC")

- project_id:

  Project identifier (optional)

- order_id:

  Order identifier (optional)

- workunit_id:

  Workunit identifier (optional)

- transform:

  Data transformation method (e.g., "vsn", "quantile", "robscale")

- date:

  Date to use for naming (defaults to current date)

## Value

Generated zip directory name

## Examples

``` r
zipdir_name("DEA", "12345", "67890", "11111", "vsn")
#> [1] "DEA_20260223_PI12345_O67890_WU11111_vsn"
zipdir_name("QC", transform = "quantile")
#> [1] "QC_20260223_quantile"
zipdir_name("DEA", workunit_id = "99999", transform = "robscale")
#> [1] "DEA_20260223_WU99999_robscale"
```
