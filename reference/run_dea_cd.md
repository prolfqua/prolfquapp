# Run differential expression analysis for CompoundDiscoverer ZIP exports

Reads the embedded prolfqua sample annotation and long feature table
from a CompoundDiscoverer ZIP export, then runs the same DEA preparation
and model pipeline as
[`run_dea`](https://prolfqua.github.io/prolfquapp/reference/run_dea.md).

## Usage

``` r
run_dea_cd(input = NULL, config, files = NULL, subset_column = NULL)
```

## Arguments

- input:

  ZIP file path or directory containing one ZIP export

- config:

  a
  [`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md)
  object

- files:

  optional file list from
  [`get_CD_export_files`](https://prolfqua.github.io/prolfquapp/reference/get_CD_export_files.md)

- subset_column:

  optional long-table subset column to filter features

## Value

list with `deanalyse`, `xd`, `annotation`, and `files`
