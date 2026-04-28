# Generate contrasts for a single-factor design

Reads an annotation file, validates the control level exists, and adds a
CONTROL column (C for control, T for treatment).

## Usage

``` r
run_contrasts_single(annotation_file, control, group = NULL)
```

## Arguments

- annotation_file:

  path to annotation CSV/TSV/XLSX

- control:

  reference level name (e.g. "WT")

- group:

  group column name, or NULL to auto-detect

## Value

data.frame with CONTROL column added

## Examples

``` r
csv <- system.file("application/contrasts/scenario1_single_factor.csv",
  package = "prolfquapp")
result <- run_contrasts_single(csv, control = "WT")
#> Rows: 9 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): file, name, group
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
table(result$group, result$CONTROL)
#>     
#>      C T
#>   KO 0 3
#>   OE 0 3
#>   WT 3 0
```
