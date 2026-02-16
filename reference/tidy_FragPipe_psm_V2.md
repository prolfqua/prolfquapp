# read psm.tsv produced by FragPipe and convert into long format

read psm.tsv produced by FragPipe and convert into long format

## Usage

``` r
tidy_FragPipe_psm_V2(
  psm_files,
  purity_threshold = 0.5,
  PeptideProphetProb = 0.9,
  abundance_threshold = 0,
  quan_column_prefix = "^Intensity",
  aggregate = TRUE
)
```

## Arguments

- psm_files:

  path(s) to psm.tsv file(s)

- purity_threshold:

  purity threshold default = 0.5

- PeptideProphetProb:

  default 0.9

- abundance_threshold:

  minimum abundance threshold

- quan_column_prefix:

  regex prefix for quantitative columns

- aggregate:

  aggregate spectra to psm level
