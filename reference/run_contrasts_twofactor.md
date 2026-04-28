# Generate contrasts for a two-factor design

Reads an annotation file and adds ContrastName/Contrast columns using
[`annotation_add_contrasts`](https://wolski.github.io/prolfqua/reference/annotation_add_contrasts.html).

## Usage

``` r
run_contrasts_twofactor(annotation_file, f1, f2, interactions = TRUE)
```

## Arguments

- annotation_file:

  path to annotation CSV/TSV/XLSX

- f1:

  primary factor column name

- f2:

  secondary factor column name

- interactions:

  logical; include interaction contrasts? Default TRUE.

## Value

data.frame with ContrastName and Contrast columns added

## Examples

``` r
csv <- system.file("application/contrasts/scenario2_two_factor.csv",
  package = "prolfquapp")
result <- run_contrasts_twofactor(csv, f1 = "treatment", f2 = "time")
#> Rows: 12 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (4): file, name, treatment, time
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
unique(result[!is.na(result$ContrastName), c("ContrastName", "Contrast")])
#> # A tibble: 4 × 2
#>   ContrastName                           Contrast                               
#>   <chr>                                  <chr>                                  
#> 1 MINOCA_vs_MI                           ( (G_MINOCA_T0 + G_MINOCA_T150)/2 - (G…
#> 2 MINOCA_vs_MI_at_T0                     G_MINOCA_T0 - G_MI_T0                  
#> 3 MINOCA_vs_MI_at_T150                   G_MINOCA_T150 - G_MI_T150              
#> 4 interaction_MINOCA_vs_MI_at_T150_vs_T0 (G_MINOCA_T150 - G_MI_T150) - (G_MINOC…
```
