# ProteinDataPrep

ProteinDataPrep

ProteinDataPrep

## Details

Handles data preparation for differential expression analysis:
contaminant/decoy filtering, peptide-to-protein aggregation, and
normalization.

## Public fields

- `prolfq_app_config`:

  ProlfquAppConfig

- `lfq_data_peptide`:

  LFQData peptide level

- `lfq_data`:

  LFQData protein level (after aggregation)

- `lfq_data_transformed`:

  normalized LFQData

- `aggregator`:

  aggregator object

- `rowAnnot`:

  ProteinAnnotation

- `summary`:

  data.frame with contaminant/decoy summary

## Methods

### Public methods

- [`ProteinDataPrep$new()`](#method-ProteinDataPrep-new)

- [`ProteinDataPrep$cont_decoy_summary()`](#method-ProteinDataPrep-cont_decoy_summary)

- [`ProteinDataPrep$remove_cont_decoy()`](#method-ProteinDataPrep-remove_cont_decoy)

- [`ProteinDataPrep$aggregate()`](#method-ProteinDataPrep-aggregate)

- [`ProteinDataPrep$get_aggregation_plots()`](#method-ProteinDataPrep-get_aggregation_plots)

- [`ProteinDataPrep$write_aggregation_plots()`](#method-ProteinDataPrep-write_aggregation_plots)

- [`ProteinDataPrep$transform_data()`](#method-ProteinDataPrep-transform_data)

- [`ProteinDataPrep$clone()`](#method-ProteinDataPrep-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize ProteinDataPrep

#### Usage

    ProteinDataPrep$new(lfq_data_peptide, rowAnnot, prolfq_app_config)

#### Arguments

- `lfq_data_peptide`:

  LFQData object at peptide level

- `rowAnnot`:

  ProteinAnnotation object

- `prolfq_app_config`:

  ProlfquAppConfig object

------------------------------------------------------------------------

### Method `cont_decoy_summary()`

Count number of contaminants and decoys

#### Usage

    ProteinDataPrep$cont_decoy_summary()

------------------------------------------------------------------------

### Method `remove_cont_decoy()`

Remove contaminants and decoys from peptide data

#### Usage

    ProteinDataPrep$remove_cont_decoy()

------------------------------------------------------------------------

### Method [`aggregate()`](https://rdrr.io/r/stats/aggregate.html)

Aggregate peptide data to protein level

#### Usage

    ProteinDataPrep$aggregate()

------------------------------------------------------------------------

### Method `get_aggregation_plots()`

Get aggregation plots

#### Usage

    ProteinDataPrep$get_aggregation_plots(exp_nr_children = 2)

#### Arguments

- `exp_nr_children`:

  minimum number of peptides per protein; default 2

------------------------------------------------------------------------

### Method `write_aggregation_plots()`

Write aggregation plots to file

#### Usage

    ProteinDataPrep$write_aggregation_plots(exp_nr_children = 2)

#### Arguments

- `exp_nr_children`:

  minimum number of peptides per protein; default 2

------------------------------------------------------------------------

### Method `transform_data()`

Transform and normalize protein-level data

#### Usage

    ProteinDataPrep$transform_data()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ProteinDataPrep$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
pep <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
pep <- prolfqua::LFQData$new(pep$data, pep$config)
pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
pA <- prolfquapp::ProteinAnnotation$new(pep, row_annot = pA, description = "fasta.annot")
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
GRP2 <- prolfquapp::make_DEA_config_R6()
GRP2$processing_options$transform <- "robscale"

data_prep <- prolfquapp::ProteinDataPrep$new(pep, pA, GRP2)
data_prep$cont_decoy_summary()
#>   totalNrOfProteins percentOfContaminants percentOfFalsePositives
#> 1               100                     0                       0
#>   NrOfProteinsNoDecoys
#> 1                  100
data_prep$remove_cont_decoy()
#> Joining with `by = join_by(protein_Id)`
#> INFO [2026-03-11 08:17:12] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
data_prep$aggregate()
#> INFO [2026-03-11 08:17:12] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
#> INFO [2026-03-11 08:17:13] END OF PROTEIN AGGREGATION
data_prep$transform_data()
#> INFO [2026-03-11 08:17:13] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Warning: Expected 1 pieces. Additional pieces discarded in 1200 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(protein_Id, sampleName)`
#> INFO [2026-03-11 08:17:13] Transforming data : robscale.
```
