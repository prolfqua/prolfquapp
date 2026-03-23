# dataset transform data

dataset transform data

## Usage

``` r
aggregate_data(lfqdata, agg_method = c("medpolish", "lmrob", "topN"), N = 3)
```

## Arguments

- lfqdata:

  LFQData object

- agg_method:

  aggregation method

- N:

  number of top peptides for topN aggregation

## Examples

``` r
xx <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- prolfqua::LFQData$new(xx$data, xx$config)
aggregated <- aggregate_data(lfqdata, agg_method = "medpolish")
#> Column added : log_abundance
#> starting aggregation
#> Column added : exp_medpolish
aggregated$response()
#> [1] "exp_medpolish"
aggregated <- aggregate_data(lfqdata, agg_method = "lmrob")
#> Column added : log_abundance
#> starting aggregation
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Column added : exp_lmrob
aggregated$response()
#> [1] "exp_lmrob"
aggregated <- aggregate_data(lfqdata, agg_method = "topN")
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
aggregated$response()
#> [1] "srm_sum_N"
```
