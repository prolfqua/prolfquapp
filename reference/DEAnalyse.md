# Differential expression analysis engine

Differential expression analysis engine

Differential expression analysis engine

## Details

Takes a prepared `ProteinDataPrep` object and runs statistical
modelling: fits linear and GLM models, computes moderated contrasts,
merges results.

## Public fields

- `prolfq_app_config`:

  ProlfquAppConfig

- `lfq_data_peptide`:

  LFQData peptide level

- `lfq_data`:

  LFQData protein level

- `lfq_data_transformed`:

  normalized LFQData

- `rowAnnot`:

  ProteinAnnotation

- `contrasts`:

  vector with contrasts

- `FDR_threshold`:

  FDR threshold

- `diff_threshold`:

  difference threshold

- `summary`:

  data.frame with contaminant/decoy summary

- `annotated_contrasts`:

  contrasts joined with row annotations

- `annotated_contrasts_signif`:

  significant annotated contrasts

- `formula`:

  model formula

- `models`:

  list of fitted models

- `contrast_results`:

  list of contrast results

- `m1_linear`:

  Linear_Model

- `m2_missing`:

  Imputed_Mean

- `m3_merged`:

  mergedModel

- `m4_glm_protein`:

  glmModel

- `m4_glm_peptide`:

  glmModelPeptide

- `default_model`:

  default model key

## Methods

### Public methods

- [`DEAnalyse$new()`](#method-DEAnalyse-new)

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

- [`DEAnalyse$clone()`](#method-DEAnalyse-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize DEAnalyse from a prepared ProteinDataPrep object

#### Usage

    DEAnalyse$new(data_prep, contrasts, default_model = "mergedModel")

#### Arguments

- `data_prep`:

  ProteinDataPrep object with aggregated and normalized data

- `contrasts`:

  named vector of contrast definitions

- `default_model`:

  which model to use for final results

------------------------------------------------------------------------

### Method `create_model_formula()`

Create model formula from transformed data config

#### Usage

    DEAnalyse$create_model_formula()

------------------------------------------------------------------------

### Method `build_model_linear_protein()`

Fit linear model at protein level

#### Usage

    DEAnalyse$build_model_linear_protein()

------------------------------------------------------------------------

### Method `get_strategy_glm_prot()`

Get GLM strategy for protein-level missingness model

#### Usage

    DEAnalyse$get_strategy_glm_prot()

------------------------------------------------------------------------

### Method `build_model_glm_protein()`

Fit generalized linear model at protein level

#### Usage

    DEAnalyse$build_model_glm_protein()

------------------------------------------------------------------------

### Method `build_model_glm_peptide()`

Fit generalized linear model at peptide level (not yet implemented)

#### Usage

    DEAnalyse$build_model_glm_peptide()

------------------------------------------------------------------------

### Method `get_contrasts_linear_protein()`

Compute moderated contrasts from linear model

#### Usage

    DEAnalyse$get_contrasts_linear_protein()

------------------------------------------------------------------------

### Method `get_contrasts_glm_peptide()`

Compute moderated contrasts from GLM peptide model

#### Usage

    DEAnalyse$get_contrasts_glm_peptide()

------------------------------------------------------------------------

### Method `get_contrasts_glm_protein()`

Compute moderated contrasts from GLM protein model

#### Usage

    DEAnalyse$get_contrasts_glm_protein()

------------------------------------------------------------------------

### Method `get_contrasts_missing_protein()`

Compute moderated contrasts from missing-value imputation model

#### Usage

    DEAnalyse$get_contrasts_missing_protein()

------------------------------------------------------------------------

### Method `get_contrasts_merged_protein()`

Merge linear and missing-value contrasts (or use linear only if
model_missing = FALSE)

#### Usage

    DEAnalyse$get_contrasts_merged_protein()

------------------------------------------------------------------------

### Method `get_annotated_contrasts()`

Join default model contrasts with protein row annotations

#### Usage

    DEAnalyse$get_annotated_contrasts()

------------------------------------------------------------------------

### Method `filter_contrasts()`

Return contrast rows passing FDR and difference thresholds

#### Usage

    DEAnalyse$filter_contrasts()

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

contrasts <- c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl")

data_prep <- prolfquapp::ProteinDataPrep$new(pep, pA, GRP2)
data_prep$cont_decoy_summary()
#>   totalNrOfProteins percentOfContaminants percentOfFalsePositives
#> 1               100                     0                       0
#>   NrOfProteinsNoDecoys
#> 1                  100
data_prep$remove_cont_decoy()
#> Joining with `by = join_by(protein_Id)`
#> INFO [2026-03-19 20:25:04] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
data_prep$aggregate()
#> INFO [2026-03-19 20:25:04] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> 

#> [================================>--------------------------------------]  46%
#> 

#> [================================>--------------------------------------]  47%
#> 

#> [=================================>-------------------------------------]  48%
#> 

#> [==================================>------------------------------------]  49%
#> 

#> [===================================>-----------------------------------]  50%
#> 

#> [===================================>-----------------------------------]  51%
#> 

#> [====================================>----------------------------------]  52%
#> 

#> [=====================================>---------------------------------]  53%
#> 

#> [=====================================>---------------------------------]  54%
#> 

#> [======================================>--------------------------------]  55%
#> 

#> [=======================================>-------------------------------]  56%
#> 

#> [=======================================>-------------------------------]  57%
#> 

#> [========================================>------------------------------]  58%
#> 

#> [=========================================>-----------------------------]  59%
#> 

#> [==========================================>----------------------------]  60%
#> 

#> [==========================================>----------------------------]  61%
#> 

#> [===========================================>---------------------------]  62%
#> 

#> [============================================>--------------------------]  63%
#> 

#> [============================================>--------------------------]  64%
#> 

#> [=============================================>-------------------------]  65%
#> 

#> [==============================================>------------------------]  66%
#> 

#> [===============================================>-----------------------]  67%
#> 

#> [===============================================>-----------------------]  68%
#> 

#> [================================================>----------------------]  69%
#> 

#> [=================================================>---------------------]  70%
#> 

#> [=================================================>---------------------]  71%
#> 

#> [==================================================>--------------------]  72%
#> 

#> [===================================================>-------------------]  73%
#> 

#> [====================================================>------------------]  74%
#> 

#> [====================================================>------------------]  75%
#> 

#> [=====================================================>-----------------]  76%
#> 

#> [======================================================>----------------]  77%
#> 

#> [======================================================>----------------]  78%
#> 

#> [=======================================================>---------------]  79%
#> 

#> [========================================================>--------------]  80%
#> 

#> [=========================================================>-------------]  81%
#> 

#> [=========================================================>-------------]  82%
#> 

#> [==========================================================>------------]  83%
#> 

#> [===========================================================>-----------]  84%
#> 

#> [===========================================================>-----------]  85%
#> 

#> [============================================================>----------]  86%
#> 

#> [=============================================================>---------]  87%
#> 

#> [=============================================================>---------]  88%
#> 

#> [==============================================================>--------]  89%
#> 

#> [===============================================================>-------]  90%
#> 

#> [================================================================>------]  91%
#> 

#> [================================================================>------]  92%
#> 

#> [=================================================================>-----]  93%
#> 

#> [==================================================================>----]  94%
#> 

#> [==================================================================>----]  95%
#> 

#> [===================================================================>---]  96%
#> 

#> [====================================================================>--]  97%
#> 

#> [=====================================================================>-]  98%
#> 

#> [=====================================================================>-]  99%
#> 

#> [=======================================================================] 100%
#> 
                                                                              
#> 

#> Column added : exp_medpolish
#> INFO [2026-03-19 20:25:04] END OF PROTEIN AGGREGATION
data_prep$transform_data()
#> INFO [2026-03-19 20:25:04] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(protein_Id, sampleName,
#> isotopeLabel)`
#> INFO [2026-03-19 20:25:04] Transforming data : robscale.

deanalyse <- prolfquapp::DEAnalyse$new(data_prep, contrasts)
mod <- deanalyse$build_model_linear_protein()
#> INFO [2026-03-19 20:25:04] fitted model with formula : normalized_abundance ~ group_
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
```
