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
#> INFO [2026-07-05 20:58:40] AGGREGATING PEPTIDE DATA: medpolish.
#> Column added : log_abundance
#> starting aggregation
#> completing cases
#> Column added : exp_medpolish
#> INFO [2026-07-05 20:58:40] END OF PROTEIN AGGREGATION
#> INFO [2026-07-05 20:58:40] Transforming using robscale.
#> Column added : log2_exp_medpolish
#> data is : TRUE
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id)`
#> INFO [2026-07-05 20:58:40] Transforming data : robscale.
#> INFO [2026-07-05 20:58:40] model formula: normalized_abundance ~ group_
#> determine linear functions:
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
dea$contrast_results[[dea$default_model]]$get_contrasts()
#> # A tibble: 20 × 14
#>    modelName estimate_type protein_Id  contrast     diff std.error avgAbd
#>    <chr>     <chr>         <chr>       <chr>       <dbl>     <dbl>  <dbl>
#>  1 lm_impute observed      0EfVhX~0087 AVsC     -0.0601     0.0540   4.38
#>  2 lm_impute observed      7cbcrd~5725 AVsC      0.761      0.117    4.61
#>  3 lm_impute observed      9VUkAq~4703 AVsC     -0.642      0.100    4.51
#>  4 lm_impute observed      BEJI92~5282 AVsC      0.213      0.0918   4.30
#>  5 lm_impute observed      CGzoYe~2147 AVsC     -0.293      0.0804   4.82
#>  6 lm_impute observed      DoWup2~5896 AVsC      0.287      0.0427   4.51
#>  7 lm_impute observed      Fl4JiV~8625 AVsC      0.0509     0.0435   4.39
#>  8 lm_impute observed      HvIpHG~9079 AVsC     -0.394      0.0450   4.45
#>  9 lm_impute observed      JcKVfU~9653 AVsC     -0.0212     0.0990   5.07
#> 10 lm_impute observed      SGIVBl~5782 AVsC      0.0154     0.0573   4.77
#> 11 lm_impute observed      0EfVhX~0087 BVsC      0.278      0.0515   4.54
#> 12 lm_impute observed      7cbcrd~5725 BVsC      0.0530     0.125    4.26
#> 13 lm_impute observed      9VUkAq~4703 BVsC     -0.852      0.108    4.41
#> 14 lm_impute observed      BEJI92~5282 BVsC      0.165      0.0864   4.28
#> 15 lm_impute observed      CGzoYe~2147 BVsC     -0.157      0.0869   4.89
#> 16 lm_impute observed      DoWup2~5896 BVsC     -0.299      0.0427   4.21
#> 17 lm_impute observed      Fl4JiV~8625 BVsC      0.0825     0.0418   4.40
#> 18 lm_impute observed      HvIpHG~9079 BVsC      0.00208    0.0419   4.65
#> 19 lm_impute observed      JcKVfU~9653 BVsC      0.0236     0.0990   5.09
#> 20 lm_impute observed      SGIVBl~5782 BVsC     -0.178      0.0554   4.68
#> # ℹ 7 more variables: statistic <dbl>, df <dbl>, p.value <dbl>, conf.low <dbl>,
#> #   conf.high <dbl>, sigma <dbl>, FDR <dbl>
```
