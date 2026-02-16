# Create DEA report in html and write data to xlsx table

For use examples see run_scripts directory

## Usage

``` r
make_DEA_report(lfqdata, protAnnot, GRP2)

render_DEA(
  GRP2,
  outpath,
  htmlname = "Result2Grp",
  word = FALSE,
  toc = TRUE,
  markdown = "_Grp2Analysis.Rmd"
)
```

## Arguments

- lfqdata:

  LFQData object

- protAnnot:

  ProteinAnnotation object

- GRP2:

  return value of `make_DEA_report`

- outpath:

  path to place output

- htmlname:

  name for html file

- word:

  default FALSE, if true create word document.s

- toc:

  if TRUE include table of contents

- markdown:

  which file to render

## See also

Other workflow:
[`make_SummarizedExperiment()`](https://prolfqua.github.io/prolfquapp/reference/make_SummarizedExperiment.md),
[`write_DEA()`](https://prolfqua.github.io/prolfquapp/reference/write_DEA.md)

Other workflow:
[`make_SummarizedExperiment()`](https://prolfqua.github.io/prolfquapp/reference/make_SummarizedExperiment.md),
[`write_DEA()`](https://prolfqua.github.io/prolfquapp/reference/write_DEA.md)
