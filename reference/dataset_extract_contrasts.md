# extect contrasts from dataset

extect contrasts from dataset

## Usage

``` r
dataset_extract_contrasts(annot, GRP2)
```

## Arguments

- annot:

  annotation data frame

- GRP2:

  configuration list

## Examples

``` r
file <- system.file("application/dataset_csv/dataset_25.csv", package = "prolfquapp")
res <- readr::read_csv(file)
#> Rows: 143 Columns: 14
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (11): Relative.Path, Name, Description, PatientID, SpinalSegment, BMsect...
#> dbl  (2): RunID, Biopsie
#> lgl  (1): Dot
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
GRP2 <- make_DEA_config()
#> Warning: DEPRECATED
GRP2 <- dataset_extract_contrasts(res,GRP2)
#> Warning: DEPRECATED
stopifnot(length(GRP2$pop$Contrasts) == 0)

file <- system.file("application/dataset_csv/dataset_26.csv", package = "prolfquapp")
res <- readr::read_csv(file)
#> Rows: 18 Columns: 8
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (8): Relative Path, Name, Animal, Treatment, SubjectID, Group, ContrastN...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
GRP2 <- dataset_extract_contrasts(res,GRP2)
#> Warning: DEPRECATED
stopifnot(length(GRP2$pop$Contrasts) == 3)
```
