# Generate differential expression analysis reports

Writes results of DEA see
[`generate_DEA_reports`](https://prolfqua.github.io/prolfquapp/reference/generate_DEA_reports.md)

## Usage

``` r
write_DEA_all(
  grp2,
  name = "",
  boxplot = TRUE,
  render = TRUE,
  ORA = TRUE,
  GSEA = TRUE,
  markdown = "_Grp2Analysis.Rmd",
  toc = TRUE
)
```

## Arguments

- grp2:

  ProlfquAppConfig configuration with results

- name:

  optional name prefix for output files

- boxplot:

  if TRUE generate boxplots

- render:

  if TRUE render HTML reports

- ORA:

  if TRUE write ORA gene lists

- GSEA:

  if TRUE write GSEA rank files

- markdown:

  Rmd template to render

- toc:

  if TRUE include table of contents
