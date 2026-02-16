# compute IBAQ values

compute IBAQ values

## Usage

``` r
compute_IBAQ_values(
  lfqdata,
  protein_annotation,
  protein_length = "protein_length",
  nr_tryptic_peptides = "nr_tryptic_peptides"
)
```

## Arguments

- lfqdata:

  LFQData object

- protein_annotation:

  ProteinAnnotation object

- protein_length:

  column name for protein length

- nr_tryptic_peptides:

  column name for number of tryptic peptides

## Examples

``` r
pAlf <- sim_data_protAnnot()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
xd <- compute_IBAQ_values(pAlf$lfqdata, pAlf$pannot)
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
xd$response()
#> [1] "IBAQValue"
pAlf <- sim_data_protAnnot(PROTEIN = TRUE)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
xd <- compute_IBAQ_values(pAlf$lfqdata, pAlf$pannot)
#> Warning: nothing to aggregate from, returning unchanged data.
xd$response()
#> [1] "IBAQValue"
```
