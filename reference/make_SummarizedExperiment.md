# Convert prolfqua differential expression analysis results to SummarizedExperiment

Convert prolfqua differential expression analysis results to
SummarizedExperiment

## Usage

``` r
make_SummarizedExperiment(
  GRP2,
  colname = NULL,
  rowname = NULL,
  strip = "~lfq~light",
  .url_builder = bfabric_url_builder
)
```

## Arguments

- GRP2:

  return value of
  [`make_DEA_report`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_report.md)

- colname:

  column name for sample identifier

- rowname:

  column name for row identifier

- strip:

  regex pattern to strip from rownames

- .url_builder:

  function to build bfabric URLs

## Value

SummarizedExperiment

## See also

Other workflow:
[`make_DEA_report()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_report.md),
[`write_DEA()`](https://prolfqua.github.io/prolfquapp/reference/write_DEA.md)
