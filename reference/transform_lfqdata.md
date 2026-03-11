# transform lfq data using robscale, vsn or log2, Assumes that data is not transformed (still needs log2 transformation)

transform lfq data using robscale, vsn or log2, Assumes that data is not
transformed (still needs log2 transformation)

## Usage

``` r
transform_lfqdata(lfqdata, method = c("robscale", "vsn", "none", "log2"))
```

## Arguments

- lfqdata:

  [`LFQData`](https://rdrr.io/pkg/prolfqua/man/LFQData.html)

- method:

  normalization method to use

## Examples

``` r
istar <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
config <- prolfqua:::old2new(istar$config)
tmp <- prolfqua::LFQData$new(istar$data, config)
tmp2 <- transform_lfqdata(tmp)
#> INFO [2026-03-11 08:07:00] Transforming using robscale.
#> Column added : log2_peptide.intensity
#> data is : TRUE
#> Warning: Expected 2 pieces. Additional pieces discarded in 25780 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(protein_Id, sampleName, peptide_Id)`
#> INFO [2026-03-11 08:07:01] Transforming data : robscale.
```
