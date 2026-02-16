# read combined_protein.tsv file for FragPipe Version 16 or newer

read combined_protein.tsv file for FragPipe Version 16 or newer

## Usage

``` r
tidy_FragPipe_combined_protein(
  combprot,
  as_list = FALSE,
  spcnames = c("Total Spectral Count", "Unique Spectral Count", "Razor Spectral Count"),
  intnames = c("Total Intensity", "Unique Intensity", "Razor Intensity"),
  maxlfqnames = c("MaxLFQ Total Intensity", "MaxLFQ Unique Intensity",
    "MaxLFQ Razor Intensity")
)
```

## Arguments

- combprot:

  path to combined_protein.tsv file

- as_list:

  return as list

## Value

tidy dataframe or list with df (e.g. total.spectral.count or
total.intensity etc).

## See also

Other FragPipe:
[`FragPipe`](https://prolfqua.github.io/prolfquapp/reference/FragPipe.md),
[`tidy_FragPipe_MSstats_csv()`](https://prolfqua.github.io/prolfquapp/reference/tidy_FragPipe_MSstats_csv.md),
[`tidy_FragPipe_combined_protein_deprec()`](https://prolfqua.github.io/prolfquapp/reference/tidy_FragPipe_combined_protein_deprec.md)
