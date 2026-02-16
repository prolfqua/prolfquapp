# read psm.tsv produced by FragPipe and convert into long format

read psm.tsv produced by FragPipe and convert into long format

## Usage

``` r
tidy_FragPipe_psm(
  psm_files,
  purity_threshold = 0.5,
  PeptideProphetProb = 0.9,
  abundance_threshold = 0,
  column_before_quants = c("Quan Usage", "Mapped Proteins"),
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

- column_before_quants:

  describes the last column before the quantitative values (this is not
  consistent with in different versions of FP, default "Quan Usage"

- aggregate:

  aggregate spectra to psm level
