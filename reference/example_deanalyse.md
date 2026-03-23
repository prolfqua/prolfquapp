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
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
#> Joining with `by = join_by(protein_Id)`
#> INFO [2026-03-23 20:32:25] removing contaminants and reverse sequences with patterns: ^zz|^CON|Cont_^REV_|^rev_
#> INFO [2026-03-23 20:32:25] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
#> INFO [2026-03-23 20:32:25] END OF PROTEIN AGGREGATION
#> INFO [2026-03-23 20:32:25] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(protein_Id, sampleName, isotopeLabel)`
#> INFO [2026-03-23 20:32:25] Transforming data : robscale.
#> INFO [2026-03-23 20:32:25] model formula: normalized_abundance ~ group_
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
#> Joining with `by = join_by(protein_Id)`
dea$contrast_results[[dea$default_model]]$get_contrasts()
#> # A tibble: 20 × 14
#>    facade     modelName  protein_Id contrast     diff std.error avgAbd statistic
#>    <chr>      <fct>      <chr>      <chr>       <dbl>     <dbl>  <dbl>     <dbl>
#>  1 lm_missing WaldTest_… 0EfVhX~00… AVsC     -0.0424     0.0502   4.38    -0.745
#>  2 lm_missing WaldTest_… 7cbcrd~57… AVsC      0.761      0.117    4.61     7.74 
#>  3 lm_missing WaldTest_… 9VUkAq~47… AVsC     -0.642      0.100    4.51    -7.41 
#>  4 lm_missing WaldTest_… BEJI92~52… AVsC      0.266      0.0957   4.29     3.14 
#>  5 lm_missing WaldTest_… CGzoYe~21… AVsC     -0.293      0.0804   4.82    -3.95 
#>  6 lm_missing WaldTest_… DoWup2~58… AVsC      0.287      0.0427   4.51     5.40 
#>  7 lm_missing WaldTest_… Fl4JiV~86… AVsC      0.0470     0.0420   4.38     0.891
#>  8 lm_missing WaldTest_… HvIpHG~90… AVsC     -0.405      0.0429   4.44    -7.62 
#>  9 lm_missing WaldTest_… JcKVfU~96… AVsC     -0.0212     0.0990   5.07    -0.244
#> 10 lm_missing WaldTest_… SGIVBl~57… AVsC      0.0102     0.0574   4.78     0.168
#> 11 lm_missing WaldTest_… 0EfVhX~00… BVsC      0.286      0.0502   4.54     5.03 
#> 12 lm_missing WaldTest_… 7cbcrd~57… BVsC      0.0530     0.125    4.26     0.504
#> 13 lm_missing WaldTest_… 9VUkAq~47… BVsC     -0.852      0.108    4.41    -9.10 
#> 14 lm_missing WaldTest_… BEJI92~52… BVsC      0.199      0.0957   4.26     2.35 
#> 15 lm_missing WaldTest_… CGzoYe~21… BVsC     -0.157      0.0869   4.89    -1.96 
#> 16 lm_missing WaldTest_… DoWup2~58… BVsC     -0.299      0.0427   4.21    -5.63 
#> 17 lm_missing WaldTest_… Fl4JiV~86… BVsC      0.0951     0.0420   4.41     1.80 
#> 18 lm_missing WaldTest_… HvIpHG~90… BVsC      0.00593    0.0429   4.64     0.112
#> 19 lm_missing WaldTest_… JcKVfU~96… BVsC      0.0236     0.0990   5.09     0.271
#> 20 lm_missing WaldTest_… SGIVBl~57… BVsC     -0.185      0.0574   4.68    -3.04 
#> # ℹ 6 more variables: df <dbl>, p.value <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   sigma <dbl>, FDR <dbl>
```
