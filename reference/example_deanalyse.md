# Create an example DEAnalyse object from simulated data

Builds a complete DEAnalyse R6 object using simulated peptide data.
Useful for vignette defaults, examples, and testing.

## Usage

``` r
example_deanalyse(Nprot = 100)
```

## Arguments

- Nprot:

  number of simulated proteins (default 100)

## Value

a `DEAnalyse` R6 object with contrasts computed and annotated

## Examples

``` r
dea <- example_deanalyse(Nprot = 10)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_children_experiment
#> Joining with `by = join_by(protein_Id)`
#> INFO [2026-06-26 06:37:59] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
#> INFO [2026-06-26 06:37:59] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> completing cases
#> Column added : exp_medpolish
#> INFO [2026-06-26 06:37:59] END OF PROTEIN AGGREGATION
#> INFO [2026-06-26 06:37:59] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id)`
#> INFO [2026-06-26 06:37:59] Transforming data : robscale.
#> INFO [2026-06-26 06:37:59] model formula: normalized_abundance ~ group_
#> Warning: ContrastsLMMissingFacade (method = 'lm_missing') is deprecated: its second leg uses ContrastsMissing (group-mean substitution, no model fit). Prefer 'lm_impute' which refits failed/singular proteins with LOD imputation and borrowed variance, flagging rescued rows as estimate_type 'lod_imputed'. See ?ContrastsLMMissingFacade for migration.
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id)`
#> Joining with `by = join_by(protein_Id)`
dea$contrast_results[[dea$default_model]]$get_contrasts()
#> # A tibble: 20 × 14
#>    modelName  estimate_type protein_Id  contrast     diff std.error avgAbd
#>    <chr>      <chr>         <chr>       <chr>       <dbl>     <dbl>  <dbl>
#>  1 lm_missing observed      0EfVhX~0087 AVsC     -0.0601     0.0540   4.38
#>  2 lm_missing observed      7cbcrd~5725 AVsC      0.761      0.117    4.61
#>  3 lm_missing observed      9VUkAq~4703 AVsC     -0.642      0.100    4.51
#>  4 lm_missing observed      BEJI92~5282 AVsC      0.213      0.0918   4.30
#>  5 lm_missing observed      CGzoYe~2147 AVsC     -0.293      0.0804   4.82
#>  6 lm_missing observed      DoWup2~5896 AVsC      0.287      0.0427   4.51
#>  7 lm_missing observed      Fl4JiV~8625 AVsC      0.0509     0.0435   4.39
#>  8 lm_missing observed      HvIpHG~9079 AVsC     -0.394      0.0450   4.45
#>  9 lm_missing observed      JcKVfU~9653 AVsC     -0.0212     0.0990   5.07
#> 10 lm_missing observed      SGIVBl~5782 AVsC      0.0154     0.0573   4.77
#> 11 lm_missing observed      0EfVhX~0087 BVsC      0.278      0.0515   4.54
#> 12 lm_missing observed      7cbcrd~5725 BVsC      0.0530     0.125    4.26
#> 13 lm_missing observed      9VUkAq~4703 BVsC     -0.852      0.108    4.41
#> 14 lm_missing observed      BEJI92~5282 BVsC      0.165      0.0864   4.28
#> 15 lm_missing observed      CGzoYe~2147 BVsC     -0.157      0.0869   4.89
#> 16 lm_missing observed      DoWup2~5896 BVsC     -0.299      0.0427   4.21
#> 17 lm_missing observed      Fl4JiV~8625 BVsC      0.0825     0.0418   4.40
#> 18 lm_missing observed      HvIpHG~9079 BVsC      0.00208    0.0419   4.65
#> 19 lm_missing observed      JcKVfU~9653 BVsC      0.0236     0.0990   5.09
#> 20 lm_missing observed      SGIVBl~5782 BVsC     -0.178      0.0554   4.68
#> # ℹ 7 more variables: statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>
```
