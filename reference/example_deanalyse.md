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
#> INFO [2026-04-28 19:51:32] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
#> INFO [2026-04-28 19:51:32] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
#> INFO [2026-04-28 19:51:32] END OF PROTEIN AGGREGATION
#> INFO [2026-04-28 19:51:32] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(protein_Id, sampleName, isotopeLabel)`
#> INFO [2026-04-28 19:51:32] Transforming data : robscale.
#> INFO [2026-04-28 19:51:32] model formula: normalized_abundance ~ group_
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
dea$contrast_results[[dea$default_model]]$get_contrasts()
#> # A tibble: 20 × 14
#>    facade     modelName  protein_Id contrast     diff std.error avgAbd statistic
#>    <chr>      <fct>      <chr>      <chr>       <dbl>     <dbl>  <dbl>     <dbl>
#>  1 lm_missing WaldTest_… 0EfVhX~00… AVsC     -0.0601     0.0540   4.38   -1.09  
#>  2 lm_missing WaldTest_… 7cbcrd~57… AVsC      0.761      0.117    4.61    6.97  
#>  3 lm_missing WaldTest_… 9VUkAq~47… AVsC     -0.642      0.100    4.51   -6.71  
#>  4 lm_missing WaldTest_… BEJI92~52… AVsC      0.213      0.0918   4.30    2.45  
#>  5 lm_missing WaldTest_… CGzoYe~21… AVsC     -0.293      0.0804   4.82   -3.57  
#>  6 lm_missing WaldTest_… DoWup2~58… AVsC      0.287      0.0427   4.51    4.87  
#>  7 lm_missing WaldTest_… Fl4JiV~86… AVsC      0.0509     0.0435   4.39    1.14  
#>  8 lm_missing WaldTest_… HvIpHG~90… AVsC     -0.394      0.0450   4.45   -7.48  
#>  9 lm_missing WaldTest_… JcKVfU~96… AVsC     -0.0212     0.0990   5.07   -0.250 
#> 10 lm_missing WaldTest_… SGIVBl~57… AVsC      0.0154     0.0573   4.77    0.294 
#> 11 lm_missing WaldTest_… 0EfVhX~00… BVsC      0.278      0.0515   4.54    5.30  
#> 12 lm_missing WaldTest_… 7cbcrd~57… BVsC      0.0530     0.125    4.26    0.454 
#> 13 lm_missing WaldTest_… 9VUkAq~47… BVsC     -0.852      0.108    4.41   -8.25  
#> 14 lm_missing WaldTest_… BEJI92~52… BVsC      0.165      0.0864   4.28    2.02  
#> 15 lm_missing WaldTest_… CGzoYe~21… BVsC     -0.157      0.0869   4.89   -1.77  
#> 16 lm_missing WaldTest_… DoWup2~58… BVsC     -0.299      0.0427   4.21   -5.07  
#> 17 lm_missing WaldTest_… Fl4JiV~86… BVsC      0.0825     0.0418   4.40    1.92  
#> 18 lm_missing WaldTest_… HvIpHG~90… BVsC      0.00208    0.0419   4.65    0.0425
#> 19 lm_missing WaldTest_… JcKVfU~96… BVsC      0.0236     0.0990   5.09    0.278 
#> 20 lm_missing WaldTest_… SGIVBl~57… BVsC     -0.178      0.0554   4.68   -3.51  
#> # ℹ 6 more variables: df <dbl>, p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>
```
