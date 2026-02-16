# nice plot for paired analysis

nice plot for paired analysis

## Usage

``` r
writeLinesPaired(bb, outpath)
```

## Arguments

- bb:

  LFQData object with paired data

- outpath:

  output directory for plots

## Examples

``` r

xd <- prolfqua::sim_lfq_data_protein_config(with_missing = FALSE, paired = TRUE)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
xd <- prolfqua::LFQData$new(xd$data, xd$config)
xa <- prolfquapp::writeLinesPaired(xd)
xa[[1]]

xa[[2]]

xa[[3]]

xa[[4]]

```
