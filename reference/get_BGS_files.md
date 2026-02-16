# get BGS and fasta file location in folder

get BGS and fasta file location in folder

get BGS and fasta file location in folder

## Usage

``` r
get_BGS_files(
  path,
  bgs_pattern = "*BGS Factory Report \\(Normal\\).tsv|_Report.tsv"
)

get_BGS_files(
  path,
  bgs_pattern = "*BGS Factory Report \\(Normal\\).tsv|_Report.tsv"
)
```

## Arguments

- path:

  path to data directory

- bgs_pattern:

  glob pattern for BGS report files

## Value

list with paths to data and fasta

list with paths to data and fasta

## Examples

``` r
if (FALSE) { # \dontrun{
x <- get_DIANN_files("inst/application/DIANN/2517219/")
} # }
if (FALSE) { # \dontrun{
x <- get_DIANN_files("inst/application/DIANN/2517219/")
} # }
```
