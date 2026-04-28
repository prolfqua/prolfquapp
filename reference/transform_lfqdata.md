# Transform lfq data using robscale, vsn or log2

Assumes that data is not transformed (still needs log2 transformation)

## Usage

``` r
transform_lfqdata(lfqdata, method = c("robscale", "vsn", "none", "log2"))
```

## Arguments

- lfqdata:

  [`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.html)

- method:

  normalization method to use

## Examples

``` r
istar <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
tmp <- prolfqua::LFQData$new(istar$data, istar$config)
tmp2 <- transform_lfqdata(tmp)
#> INFO [2026-04-28 19:51:56] Transforming using robscale.
#> Column added : log2_peptide.intensity
#> data is : TRUE
#> Joining with `by = join_by(protein_Id, sampleName, isotope, peptide_Id)`
#> INFO [2026-04-28 19:51:57] Transforming data : robscale.
```
