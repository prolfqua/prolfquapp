# simulate peptdata and fitting protein annotation for testing

simulate peptdata and fitting protein annotation for testing

## Usage

``` r
sim_data_protAnnot(Nprot = 100, PROTEIN = FALSE)
```

## Arguments

- Nprot:

  number of proteins to simulate

- PROTEIN:

  if TRUE simulate protein-level data

## Examples

``` r
res <- sim_data_protAnnot()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_children_experiment
res <- sim_data_protAnnot(PROTEIN = TRUE)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_children_experiment
```
