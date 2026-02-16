# Get all processing functions from all packages

Get all processing functions from all packages

## Usage

``` r
get_procfuncs(
  function_name = "prolfqua_preprocess_functions",
  prefix = "prolfqua"
)
```

## Arguments

- function_name:

  The name of the function to get

- prefix:

  The prefix of the package names

## Value

A list of processing functions

## Examples

``` r
get_procfuncs()
#> $prolfquapp.DIANN
#> $prolfquapp.DIANN$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 1)"
#> 
#> $prolfquapp.DIANN$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.DIANN_PEPTIDE
#> $prolfquapp.DIANN_PEPTIDE$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN_PEPTIDE$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 2)"
#> 
#> $prolfquapp.DIANN_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.FP_TMT
#> $prolfquapp.FP_TMT$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 1)"
#> 
#> $prolfquapp.FP_TMT$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.FP_TMT_PEPTIDE
#> $prolfquapp.FP_TMT_PEPTIDE$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 2)"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.MAXQUANT
#> $prolfquapp.MAXQUANT$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MAXQUANT$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE
#> $prolfquapp.MAXQUANT_PEPTIDE$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MSSTATS
#> $prolfquapp.MSSTATS$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_PEPTIDE
#> $prolfquapp.MSSTATS_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA
#> $prolfquapp.MSSTATS_FP_DIA$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.BGS
#> $prolfquapp.BGS$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.BGS$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.BGS_PEPTIDE
#> $prolfquapp.BGS_PEPTIDE$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.BGS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.DUMMY
#> $prolfquapp.DUMMY$get_files
#> [1] "prolfquapp::get_dummy_files"
#> 
#> $prolfquapp.DUMMY$preprocess
#> [1] "prolfquapp::preprocess_dummy"
#> 
#> $prolfquapp.DUMMY$extra_args
#> [1] "list()"
#> 
#> 
#> $prolfquapp.MZMINE
#> $prolfquapp.MZMINE$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINE$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINE$extra_args
#> [1] "list(annotated = FALSE)"
#> 
#> 
#> $prolfquapp.MZMINEannot
#> $prolfquapp.MZMINEannot$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINEannot$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINEannot$extra_args
#> [1] "list(annotated = TRUE)"
#> 
#> 
get_procfuncs("prolfqua_preprocess_functions", "prolfqua")
#> $prolfquapp.DIANN
#> $prolfquapp.DIANN$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 1)"
#> 
#> $prolfquapp.DIANN$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.DIANN_PEPTIDE
#> $prolfquapp.DIANN_PEPTIDE$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN_PEPTIDE$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 2)"
#> 
#> $prolfquapp.DIANN_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.FP_TMT
#> $prolfquapp.FP_TMT$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 1)"
#> 
#> $prolfquapp.FP_TMT$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.FP_TMT_PEPTIDE
#> $prolfquapp.FP_TMT_PEPTIDE$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 2)"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.MAXQUANT
#> $prolfquapp.MAXQUANT$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MAXQUANT$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE
#> $prolfquapp.MAXQUANT_PEPTIDE$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MSSTATS
#> $prolfquapp.MSSTATS$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_PEPTIDE
#> $prolfquapp.MSSTATS_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA
#> $prolfquapp.MSSTATS_FP_DIA$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.BGS
#> $prolfquapp.BGS$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.BGS$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.BGS_PEPTIDE
#> $prolfquapp.BGS_PEPTIDE$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.BGS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.DUMMY
#> $prolfquapp.DUMMY$get_files
#> [1] "prolfquapp::get_dummy_files"
#> 
#> $prolfquapp.DUMMY$preprocess
#> [1] "prolfquapp::preprocess_dummy"
#> 
#> $prolfquapp.DUMMY$extra_args
#> [1] "list()"
#> 
#> 
#> $prolfquapp.MZMINE
#> $prolfquapp.MZMINE$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINE$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINE$extra_args
#> [1] "list(annotated = FALSE)"
#> 
#> 
#> $prolfquapp.MZMINEannot
#> $prolfquapp.MZMINEannot$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINEannot$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINEannot$extra_args
#> [1] "list(annotated = TRUE)"
#> 
#> 
get_procfuncs("prolfqua_preprocess_functions", "prolfquapp")
#> $prolfquapp.DIANN
#> $prolfquapp.DIANN$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 1)"
#> 
#> $prolfquapp.DIANN$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.DIANN_PEPTIDE
#> $prolfquapp.DIANN_PEPTIDE$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN_PEPTIDE$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 2)"
#> 
#> $prolfquapp.DIANN_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.FP_TMT
#> $prolfquapp.FP_TMT$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 1)"
#> 
#> $prolfquapp.FP_TMT$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.FP_TMT_PEPTIDE
#> $prolfquapp.FP_TMT_PEPTIDE$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 2)"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.MAXQUANT
#> $prolfquapp.MAXQUANT$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MAXQUANT$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE
#> $prolfquapp.MAXQUANT_PEPTIDE$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MSSTATS
#> $prolfquapp.MSSTATS$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_PEPTIDE
#> $prolfquapp.MSSTATS_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA
#> $prolfquapp.MSSTATS_FP_DIA$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.BGS
#> $prolfquapp.BGS$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.BGS$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.BGS_PEPTIDE
#> $prolfquapp.BGS_PEPTIDE$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.BGS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.DUMMY
#> $prolfquapp.DUMMY$get_files
#> [1] "prolfquapp::get_dummy_files"
#> 
#> $prolfquapp.DUMMY$preprocess
#> [1] "prolfquapp::preprocess_dummy"
#> 
#> $prolfquapp.DUMMY$extra_args
#> [1] "list()"
#> 
#> 
#> $prolfquapp.MZMINE
#> $prolfquapp.MZMINE$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINE$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINE$extra_args
#> [1] "list(annotated = FALSE)"
#> 
#> 
#> $prolfquapp.MZMINEannot
#> $prolfquapp.MZMINEannot$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINEannot$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINEannot$extra_args
#> [1] "list(annotated = TRUE)"
#> 
#> 
get_procfuncs("prolfqua_preprocess_functions", "prolfquapp")
#> $prolfquapp.DIANN
#> $prolfquapp.DIANN$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 1)"
#> 
#> $prolfquapp.DIANN$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.DIANN_PEPTIDE
#> $prolfquapp.DIANN_PEPTIDE$get_files
#> [1] "prolfquapp::get_DIANN_files"
#> 
#> $prolfquapp.DIANN_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_DIANN"
#> 
#> $prolfquapp.DIANN_PEPTIDE$extra_args
#> [1] "list(q_value = 0.01, hierarchy_depth = 2)"
#> 
#> $prolfquapp.DIANN_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_diann"
#> 
#> 
#> $prolfquapp.FP_TMT
#> $prolfquapp.FP_TMT$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 1)"
#> 
#> $prolfquapp.FP_TMT$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.FP_TMT_PEPTIDE
#> $prolfquapp.FP_TMT_PEPTIDE$get_files
#> [1] "prolfquapp::get_FP_PSM_files"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_FP_PSM"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$extra_args
#> [1] "list(purity_threshold = 0.5, PeptideProphetProb = 0.9, hierarchy_depth = 2)"
#> 
#> $prolfquapp.FP_TMT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_FP_TMT"
#> 
#> 
#> $prolfquapp.MAXQUANT
#> $prolfquapp.MAXQUANT$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MAXQUANT$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE
#> $prolfquapp.MAXQUANT_PEPTIDE$get_files
#> [1] "prolfquapp::get_MQ_peptide_files"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MQ_peptide"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MAXQUANT_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MAXQUANT"
#> 
#> 
#> $prolfquapp.MSSTATS
#> $prolfquapp.MSSTATS$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_PEPTIDE
#> $prolfquapp.MSSTATS_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA
#> $prolfquapp.MSSTATS_FP_DIA$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$get_files
#> [1] "prolfquapp::get_MSstats_files"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_MSstats_FPDIA"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.MSSTATS_FP_DIA_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_MSSTATS"
#> 
#> 
#> $prolfquapp.BGS
#> $prolfquapp.BGS$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS$extra_args
#> [1] "list(hierarchy_depth = 1)"
#> 
#> $prolfquapp.BGS$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.BGS_PEPTIDE
#> $prolfquapp.BGS_PEPTIDE$get_files
#> [1] "prolfquapp::get_BGS_files"
#> 
#> $prolfquapp.BGS_PEPTIDE$preprocess
#> [1] "prolfquapp::preprocess_BGS"
#> 
#> $prolfquapp.BGS_PEPTIDE$extra_args
#> [1] "list(hierarchy_depth = 2)"
#> 
#> $prolfquapp.BGS_PEPTIDE$dataset
#> [1] "prolfquapp::dataset_template_BGS"
#> 
#> 
#> $prolfquapp.DUMMY
#> $prolfquapp.DUMMY$get_files
#> [1] "prolfquapp::get_dummy_files"
#> 
#> $prolfquapp.DUMMY$preprocess
#> [1] "prolfquapp::preprocess_dummy"
#> 
#> $prolfquapp.DUMMY$extra_args
#> [1] "list()"
#> 
#> 
#> $prolfquapp.MZMINE
#> $prolfquapp.MZMINE$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINE$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINE$extra_args
#> [1] "list(annotated = FALSE)"
#> 
#> 
#> $prolfquapp.MZMINEannot
#> $prolfquapp.MZMINEannot$get_files
#> [1] "prolfquapp::get_mzMine_files"
#> 
#> $prolfquapp.MZMINEannot$preprocess
#> [1] "prolfquapp::preprocess_mzMine"
#> 
#> $prolfquapp.MZMINEannot$extra_args
#> [1] "list(annotated = TRUE)"
#> 
#> 
get_procfuncs("prolfqua_preprocess_functions", "xdx")
#> NULL
```
