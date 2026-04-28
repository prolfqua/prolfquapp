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

- `lfq_data_peptide_transformed`:

  transformed peptide-level LFQData (for nested facades)

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

- [`ProteinDataPrep$transform_peptide_data()`](#method-ProteinDataPrep-transform_peptide_data)

- [`ProteinDataPrep$build_deanalyse()`](#method-ProteinDataPrep-build_deanalyse)

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

### Method `transform_peptide_data()`

Transform peptide-level data (for nested facades like lmer/ropeca)

#### Usage

    ProteinDataPrep$transform_peptide_data()

------------------------------------------------------------------------

### Method `build_deanalyse()`

Build a DEAnalyse object with the correct data for the chosen facade

#### Usage

    ProteinDataPrep$build_deanalyse(contrasts, default_model = NULL)

#### Arguments

- `contrasts`:

  named character vector of contrast definitions

- `default_model`:

  facade registry key, or NULL to read from config

#### Returns

DEAnalyse R6 object

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
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
pep <- prolfqua::LFQData$new(pep$data, pep$config)
pA <- data.frame(protein_Id = unique(pep$data_long()$protein_Id))
pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
pA <- prolfquapp::ProteinAnnotation$new(pep, row_annot = pA, description = "fasta.annot")
#> Warning: no exp_nr_children column specified, computing using nr_children_experiment
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
#> INFO [2026-04-28 19:51:25] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
data_prep$aggregate()
#> INFO [2026-04-28 19:51:25] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
#> INFO [2026-04-28 19:51:27] END OF PROTEIN AGGREGATION
data_prep$transform_data()
#> INFO [2026-04-28 19:51:27] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(protein_Id, sampleName, isotopeLabel)`
#> INFO [2026-04-28 19:51:27] Transforming data : robscale.
```
