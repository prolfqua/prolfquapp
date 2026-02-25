# will replace make_DEA_report

will replace make_DEA_report

will replace make_DEA_report

## Public fields

- `prolfq_app_config`:

  ProlfquAppConfig

- `lfq_data_peptide`:

  LFQData peptide level

- `lfq_data`:

  LFQData

- `lfq_data_transformed`:

  transformed LFQData

- `lfq_data_subset`:

  subset of LFQData

- `aggregator`:

  aggregator

- `rowAnnot`:

  ProteinAnnotation

- `contrasts`:

  vector with contrasts

- `FDR_threshold`:

  fdr threshold

- `diff_threshold`:

  diff_threshold

- `summary`:

  data.frame with contaminant/decoy summary

- `annotated_contrasts`:

  contrasts joined with row annotations

- `annotated_contrasts_signif`:

  significant annotated contrasts

- `reference_proteins`:

  reference proteins to use for internal normalization

- `formula`:

  model formula

- `formula_glm_peptide`:

  glm peptide formula

- `models`:

  list of fitted models

- `contrast_results`:

  list of contrast results

- `m1_linear`:

  linearModel

- `m2_missing`:

  imputedModel

- `m3_merged`:

  mergedModel

- `m4_glm_protein`:

  m4_glm_protein

- `m4_glm_peptide`:

  m4_glm_peptide

- `default_model`:

  default_model

## Methods

### Public methods

- [`DEAnalyse$new()`](#method-DEAnalyse-new)

- [`DEAnalyse$cont_decoy_summary()`](#method-DEAnalyse-cont_decoy_summary)

- [`DEAnalyse$remove_cont_decoy()`](#method-DEAnalyse-remove_cont_decoy)

- [`DEAnalyse$aggregate()`](#method-DEAnalyse-aggregate)

- [`DEAnalyse$get_aggregation_plots()`](#method-DEAnalyse-get_aggregation_plots)

- [`DEAnalyse$write_aggregation_plots()`](#method-DEAnalyse-write_aggregation_plots)

- [`DEAnalyse$transform_data()`](#method-DEAnalyse-transform_data)

- [`DEAnalyse$create_model_formula()`](#method-DEAnalyse-create_model_formula)

- [`DEAnalyse$build_model_linear_protein()`](#method-DEAnalyse-build_model_linear_protein)

- [`DEAnalyse$get_strategy_glm_prot()`](#method-DEAnalyse-get_strategy_glm_prot)

- [`DEAnalyse$build_model_glm_protein()`](#method-DEAnalyse-build_model_glm_protein)

- [`DEAnalyse$build_model_glm_peptide()`](#method-DEAnalyse-build_model_glm_peptide)

- [`DEAnalyse$get_contrasts_linear_protein()`](#method-DEAnalyse-get_contrasts_linear_protein)

- [`DEAnalyse$get_contrasts_glm_peptide()`](#method-DEAnalyse-get_contrasts_glm_peptide)

- [`DEAnalyse$get_contrasts_glm_protein()`](#method-DEAnalyse-get_contrasts_glm_protein)

- [`DEAnalyse$get_contrasts_missing_protein()`](#method-DEAnalyse-get_contrasts_missing_protein)

- [`DEAnalyse$get_contrasts_merged_protein()`](#method-DEAnalyse-get_contrasts_merged_protein)

- [`DEAnalyse$get_annotated_contrasts()`](#method-DEAnalyse-get_annotated_contrasts)

- [`DEAnalyse$filter_contrasts()`](#method-DEAnalyse-filter_contrasts)

- [`DEAnalyse$filter_data()`](#method-DEAnalyse-filter_data)

- [`DEAnalyse$get_boxplots()`](#method-DEAnalyse-get_boxplots)

- [`DEAnalyse$contrasts_to_Grob()`](#method-DEAnalyse-contrasts_to_Grob)

- [`DEAnalyse$get_boxplots_contrasts()`](#method-DEAnalyse-get_boxplots_contrasts)

- [`DEAnalyse$write_boxplots_contrasts()`](#method-DEAnalyse-write_boxplots_contrasts)

- [`DEAnalyse$clone()`](#method-DEAnalyse-clone)

------------------------------------------------------------------------

### Method `new()`

initialize DEAnalyse with data and configuration

#### Usage

    DEAnalyse$new(
      lfq_data_peptide,
      rowAnnot,
      prolfq_app_config,
      contrasts,
      default_model = "mergedModel"
    )

#### Arguments

- `lfq_data_peptide`:

  LFQData object at peptide level

- `rowAnnot`:

  ProteinAnnotation object

- `prolfq_app_config`:

  ProlfquAppConfig object

- `contrasts`:

  vector with contrasts

- `default_model`:

  default model to use

------------------------------------------------------------------------

### Method `cont_decoy_summary()`

count number of decoys

#### Usage

    DEAnalyse$cont_decoy_summary()

------------------------------------------------------------------------

### Method `remove_cont_decoy()`

remove contaminants and decoys

#### Usage

    DEAnalyse$remove_cont_decoy()

------------------------------------------------------------------------

### Method [`aggregate()`](https://rdrr.io/r/stats/aggregate.html)

aggregate peptide data

#### Usage

    DEAnalyse$aggregate()

------------------------------------------------------------------------

### Method `get_aggregation_plots()`

get aggregation plots

#### Usage

    DEAnalyse$get_aggregation_plots(exp_nr_children = 2)

#### Arguments

- `exp_nr_children`:

  nr children to filter; default \>=2

------------------------------------------------------------------------

### Method `write_aggregation_plots()`

write aggregation plots

#### Usage

    DEAnalyse$write_aggregation_plots(exp_nr_children = 2)

#### Arguments

- `exp_nr_children`:

  nr children to filter; default \>=2

------------------------------------------------------------------------

### Method `transform_data()`

transform data

#### Usage

    DEAnalyse$transform_data()

------------------------------------------------------------------------

### Method `create_model_formula()`

create model formula

#### Usage

    DEAnalyse$create_model_formula()

------------------------------------------------------------------------

### Method `build_model_linear_protein()`

fit linear model

#### Usage

    DEAnalyse$build_model_linear_protein()

------------------------------------------------------------------------

### Method `get_strategy_glm_prot()`

get strategy

#### Usage

    DEAnalyse$get_strategy_glm_prot()

------------------------------------------------------------------------

### Method `build_model_glm_protein()`

fit generalized linear model

#### Usage

    DEAnalyse$build_model_glm_protein()

------------------------------------------------------------------------

### Method `build_model_glm_peptide()`

fit generalized linear model

#### Usage

    DEAnalyse$build_model_glm_peptide()

------------------------------------------------------------------------

### Method `get_contrasts_linear_protein()`

compute contrasts linear

#### Usage

    DEAnalyse$get_contrasts_linear_protein()

------------------------------------------------------------------------

### Method `get_contrasts_glm_peptide()`

get contrasts from glm model

#### Usage

    DEAnalyse$get_contrasts_glm_peptide()

------------------------------------------------------------------------

### Method `get_contrasts_glm_protein()`

get contrasts from glm model for peptides

#### Usage

    DEAnalyse$get_contrasts_glm_protein()

------------------------------------------------------------------------

### Method `get_contrasts_missing_protein()`

compute missing contrasts

#### Usage

    DEAnalyse$get_contrasts_missing_protein()

------------------------------------------------------------------------

### Method `get_contrasts_merged_protein()`

merge contrasts

#### Usage

    DEAnalyse$get_contrasts_merged_protein()

------------------------------------------------------------------------

### Method `get_annotated_contrasts()`

compute annotated contrasts by joining row annotations with default
model contrasts

#### Usage

    DEAnalyse$get_annotated_contrasts()

------------------------------------------------------------------------

### Method `filter_contrasts()`

filter contrasts for threshold

#### Usage

    DEAnalyse$filter_contrasts()

------------------------------------------------------------------------

### Method `filter_data()`

filter transformed lfq data for significant proteins.

#### Usage

    DEAnalyse$filter_data()

------------------------------------------------------------------------

### Method `get_boxplots()`

create boxplots

#### Usage

    DEAnalyse$get_boxplots()

------------------------------------------------------------------------

### Method `contrasts_to_Grob()`

create boxplots

#### Usage

    DEAnalyse$contrasts_to_Grob()

------------------------------------------------------------------------

### Method `get_boxplots_contrasts()`

get box with contrast information

#### Usage

    DEAnalyse$get_boxplots_contrasts()

------------------------------------------------------------------------

### Method `write_boxplots_contrasts()`

write boxplots contrasts to file

#### Usage

    DEAnalyse$write_boxplots_contrasts(filename = "boxplots")

#### Arguments

- `filename`:

  filename to write to

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DEAnalyse$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# example code

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
GRP2$processing_options$diff_threshold <- 0.2

GRP2$processing_options$transform <- "robscale"
pep$factors()
#> # A tibble: 12 × 3
#>    sample  sampleName group_
#>    <chr>   <chr>      <chr> 
#>  1 A_V1    A_V1       A     
#>  2 A_V2    A_V2       A     
#>  3 A_V3    A_V3       A     
#>  4 A_V4    A_V4       A     
#>  5 B_V1    B_V1       B     
#>  6 B_V2    B_V2       B     
#>  7 B_V3    B_V3       B     
#>  8 B_V4    B_V4       B     
#>  9 Ctrl_V1 Ctrl_V1    Ctrl  
#> 10 Ctrl_V2 Ctrl_V2    Ctrl  
#> 11 Ctrl_V3 Ctrl_V3    Ctrl  
#> 12 Ctrl_V4 Ctrl_V4    Ctrl  
contrasts <- c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")
# DEAnalyse$debug("get_contrasts_glm_peptide")
# DEAnalyse$debug("build_model_glm_protein")
deanalyse <- prolfquapp::DEAnalyse$new(pep, pA, GRP2, contrasts)
deanalyse$lfq_data_peptide$hierarchy_counts()
#> # A tibble: 1 × 3
#>   isotopeLabel protein_Id peptide_Id
#>   <chr>             <int>      <int>
#> 1 light               100        350
deanalyse$cont_decoy_summary()
#>   totalNrOfProteins percentOfContaminants percentOfFalsePositives
#> 1               100                     0                       0
#>   NrOfProteinsNoDecoys
#> 1                  100
deanalyse$prolfq_app_config$processing_options$remove_cont <- TRUE
deanalyse$remove_cont_decoy()
#> Joining with `by = join_by(protein_Id)`
#> INFO [2026-02-25 16:35:42] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
deanalyse$aggregate()
#> INFO [2026-02-25 16:35:42] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
#> INFO [2026-02-25 16:35:43] END OF PROTEIN AGGREGATION
pl <- deanalyse$get_aggregation_plots(exp_nr_children = 10)
#> Joining with `by = join_by(protein_Id)`
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the prolfqua package.
#>   Please report the issue at <https://github.com/wolski/prolfqua/issues>.
print(pl$plots[[3]])
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_line()`).

deanalyse$transform_data()
#> INFO [2026-02-25 16:35:43] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Warning: Expected 1 pieces. Additional pieces discarded in 1200 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(protein_Id, sampleName)`
#> INFO [2026-02-25 16:35:43] Transforming data : robscale.
mod <- deanalyse$build_model_linear_protein()
#> INFO [2026-02-25 16:35:43] fitted model with formula : normalized_abundance ~ group_
#> Joining with `by = join_by(protein_Id)`
contlm <- deanalyse$get_contrasts_linear_protein()


merged <- deanalyse$get_contrasts_merged_protein()
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> completing cases
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
stopifnot(nrow(merged$get_contrasts()) == 200)
stopifnot(nrow(merged$get_contrasts()) == 200)
# deanalyse$create_model_formula()
# deanalyse$build_model_glm_protein()
# deanalyse$build_model_glm_peptide()
xprot <- deanalyse$get_contrasts_glm_protein()
#> completing cases
#> INFO [2026-02-25 16:35:44] fitted model with formula : binresp ~ group_
#> Joining with `by = join_by(protein_Id)`
if(FALSE){
xprot$get_contrasts()
xprot$get_Plotter()$volcano()
xpep <- deanalyse$get_contrasts_glm_peptide()
xpep$get_Plotter()$volcano()
sr <- deanalyse$lfq_data_peptide$get_Summariser()



deanalyse$filter_contrasts()

xd <- deanalyse$filter_data()
xd <- deanalyse$contrasts_to_Grob()
bb <- deanalyse$get_boxplots()
bx <- deanalyse$get_boxplots_contrasts()
dev.off()
grid::grid.draw(bx$bxpl_grobs[[1]])
# deanalyse$write_boxplots_contrasts("test.pdf")
}
```
