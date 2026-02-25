# transform lfq data using robscale, vsn or log2, Assumes that data is not transformed (still needs log2 transformation)

It will also run internal but then robscale must be used.

## Usage

``` r
transform_lfqdata(
  lfqdata,
  method = c("robscale", "vsn", "none", "log2"),
  internal = NULL
)
```

## Arguments

- lfqdata:

  [`LFQData`](https://rdrr.io/pkg/prolfqua/man/LFQData.html)

- method:

  normalization method to use

- internal:

  a data.frame with protein ids to be used for internal calibration,
  column name must be the same as

## Examples

``` r
istar <- prolfqua::prolfqua_data('data_ionstar')$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
config <- prolfqua:::old2new(istar$config)
tmp <- prolfqua::LFQData$new(istar$data, config)
d <- istar$d
internal <- dplyr::filter(d, protein_Id %in% sample(unique(d$protein_Id), 3 )) |>
  dplyr::select(all_of(tmp$config$hierarchy_keys()[1])) |> dplyr::distinct()
tmp2 <- transform_lfqdata(tmp, internal = internal)
#> INFO [2026-02-25 16:36:17] Transforming using robscale.
#> Column added : log2_peptide.intensity
#> data is : TRUE
#> Warning: Expected 2 pieces. Additional pieces discarded in 25780 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(protein_Id, sampleName, peptide_Id)`
#> INFO [2026-02-25 16:36:18] Transforming using robscale,
#> INFO [2026-02-25 16:36:18] Transforming Attempt of internal calibration.
#> Joining with `by = join_by(protein_Id)`
#> Column added : log2_peptide.intensity
#> Column added : log2_peptide.intensity
#> data is : TRUE
#> Warning: Expected 2 pieces. Additional pieces discarded in 25780 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(protein_Id, sampleName, peptide_Id)`
#> INFO [2026-02-25 16:36:18] Transforming data : robscale.
tmp2 <- transform_lfqdata(tmp)
#> INFO [2026-02-25 16:36:18] Transforming using robscale.
#> Column added : log2_peptide.intensity
#> data is : TRUE
#> Warning: Expected 2 pieces. Additional pieces discarded in 25780 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(protein_Id, sampleName, peptide_Id)`
#> INFO [2026-02-25 16:36:19] Transforming data : robscale.
```
