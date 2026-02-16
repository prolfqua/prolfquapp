# read minimal yaml and convert to R6 object

read minimal yaml and convert to R6 object

## Usage

``` r
list_to_R6_app_config(dd)
```

## Arguments

- dd:

  list containing configuration data

## Value

ProlfquAppConfig R6 object

## Examples

``` r
DEAconfig <- make_DEA_config_R6(WORKUNITID = "3333")
configList <- prolfqua::R6_extract_values(DEAconfig)
stopifnot(class(configList) == "list")
old <- configList$zipdir_name
config <- list_to_R6_app_config(configList)
stopifnot(config$zipdir_name == old)
stopifnot("ProlfquAppConfig" %in% class(config))
stopifnot(config$zipdir_name == configList$zipdir_name)
```
