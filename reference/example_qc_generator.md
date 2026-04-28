# Create an example QC_generator object from simulated data

Builds a QC_generator R6 object with simulated peptide data and fake
FASTA-derived columns (protein_length, nr_tryptic_peptides) so that IBAQ
computation works. Useful for vignette defaults and testing.

## Usage

``` r
example_qc_generator(Nprot = 100)
```

## Arguments

- Nprot:

  number of simulated proteins (default 100)

## Value

a `QC_generator` R6 object

## Examples

``` r
pap <- example_qc_generator(Nprot = 10)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_children_experiment
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
pap$get_prot_data()
```
