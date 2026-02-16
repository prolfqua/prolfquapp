# Compute number of tryptic peptides

Compute number of tryptic peptides

## Usage

``` r
nr_tryptic_peptides(sequence, min_length = 6, max_length = 30)
```

## Arguments

- sequence:

  amino acid sequence

- min_length:

  minimum peptide length

- max_length:

  maximum peptide length

## Examples

``` r
# example code

sequence <- "MKGLPRAKSHGSTGWGKRKRNKPK"
nr_tryptic_peptides(sequence, min_length=5)
#> [1] 1
```
