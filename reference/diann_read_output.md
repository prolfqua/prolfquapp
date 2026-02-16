# read DiaNN diann-output.tsv file

filter for 2 peptides per protein, and for Q.Value \< 0.01 (default)

## Usage

``` r
diann_read_output(data, Lib.PG.Q.Value = 0.01, PG.Q.Value = 0.05)
```

## Arguments

- data:

  data frame of DIA-NN report

- Lib.PG.Q.Value:

  library protein group q-value threshold

- PG.Q.Value:

  protein group q-value threshold

## Examples

``` r
if (FALSE) { # \dontrun{
xx <- readr::read_tsv("WU292720_report.tsv")
report2 <- prolfquapp::diann_read_output(xx)
nrow(report2)
} # }
```
