# Synchronize opt and config

Synchronize opt and config

## Usage

``` r
sync_opt_config(opt, config)
```

## Arguments

- opt:

  list of command-line options

- config:

  ProlfquAppConfig configuration object

## Examples

``` r
opt <- list()
config <- make_DEA_config_R6()
xx <- sync_opt_config(opt, config)
stopifnot(names(xx$opt) %in% c("software","outdir"))
stopifnot(xx$opt$software == "DIANN")
stopifnot(xx$opt$outdir == ".")
opt$software <- "FP_TMT"
opt$outdir <- "testdir"
xx <- sync_opt_config(opt, config)
stopifnot(xx$config$path == "testdir")
stopifnot(xx$config$software == "FP_TMT")
```
