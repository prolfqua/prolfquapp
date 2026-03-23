# Differential expression analysis engine using prolfqua facade classes

Differential expression analysis engine using prolfqua facade classes

Differential expression analysis engine using prolfqua facade classes

## Details

Takes prepared LFQData (at the correct hierarchy level for the chosen
facade) and runs statistical modelling via prolfqua's ContrastsFacade
classes.

The caller (e.g. `ProteinDataPrep$build_deanalyse()`) is responsible for
providing data at the right level: aggregated protein-level for most
facades, or nested peptide-level for `lmer`/`ropeca`.

## Public fields

- `prolfq_app_config`:

  ProlfquAppConfig

- `lfq_data`:

  LFQData to model (transformed, at correct hierarchy level)

- `lfq_data_raw`:

  raw (untransformed) LFQData for reporting

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

- `contrast_results`:

  named list of facade objects

- `default_model`:

  facade registry key for the default model

## Methods

### Public methods

- [`DEAnalyse$new()`](#method-DEAnalyse-new)

- [`DEAnalyse$build_facade()`](#method-DEAnalyse-build_facade)

- [`DEAnalyse$build_default()`](#method-DEAnalyse-build_default)

- [`DEAnalyse$get_annotated_contrasts()`](#method-DEAnalyse-get_annotated_contrasts)

- [`DEAnalyse$filter_contrasts()`](#method-DEAnalyse-filter_contrasts)

- [`DEAnalyse$clone()`](#method-DEAnalyse-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize DEAnalyse

#### Usage

    DEAnalyse$new(
      lfq_data,
      rowAnnot,
      prolfq_app_config,
      contrasts,
      default_model = "lm_missing",
      lfq_data_raw = NULL,
      summary = NULL
    )

#### Arguments

- `lfq_data`:

  LFQData to model (transformed, at correct hierarchy level)

- `rowAnnot`:

  ProteinAnnotation object

- `prolfq_app_config`:

  ProlfquAppConfig object

- `contrasts`:

  named vector of contrast definitions

- `default_model`:

  facade registry key (default "lm_missing")

- `lfq_data_raw`:

  raw (untransformed) LFQData for reporting (optional)

- `summary`:

  data.frame with contaminant/decoy summary (optional)

------------------------------------------------------------------------

### Method `build_facade()`

Build a facade by registry key

#### Usage

    DEAnalyse$build_facade(name, modelstr = NULL)

#### Arguments

- `name`:

  facade registry key (e.g. "lm", "lm_missing", "limma")

- `modelstr`:

  model formula string; auto-generated if NULL

#### Returns

the facade object (invisibly)

------------------------------------------------------------------------

### Method `build_default()`

Build the default facade (as set in default_model)

#### Usage

    DEAnalyse$build_default()

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
#> INFO [2026-03-23 20:32:17] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
data_prep$aggregate()
#> INFO [2026-03-23 20:32:17] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
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
#> INFO [2026-03-23 20:32:17] END OF PROTEIN AGGREGATION
data_prep$transform_data()
#> INFO [2026-03-23 20:32:17] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(protein_Id, sampleName, isotopeLabel)`
#> INFO [2026-03-23 20:32:18] Transforming data : robscale.

deanalyse <- data_prep$build_deanalyse(contrasts)
deanalyse$build_default()
#> INFO [2026-03-23 20:32:18] model formula: normalized_abundance ~ group_
#> Joining with `by = join_by(protein_Id)`
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
stopifnot(nrow(deanalyse$contrast_results[[deanalyse$default_model]]$get_contrasts()) == 200)
```
