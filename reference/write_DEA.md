# Write differential expression analysis results

Write differential expression analysis results

## Usage

``` r
write_DEA(
  GRP2,
  outpath,
  ORA = TRUE,
  GSEA = TRUE,
  xlsxname = "AnalysisResults",
  write = TRUE
)
```

## Arguments

- GRP2:

  return value of
  [`make_DEA_report`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_report.md)

- outpath:

  path to place output

- ORA:

  if TRUE write ORA gene lists

- GSEA:

  if TRUE write GSEA rank files

- xlsxname:

  file name for xlsx

- write:

  if TRUE write files to disk

## See also

Other workflow:
[`make_DEA_report()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_report.md),
[`make_SummarizedExperiment()`](https://prolfqua.github.io/prolfquapp/reference/make_SummarizedExperiment.md)
