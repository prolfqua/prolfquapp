# read yaml file and convert to R6 configuration object

read yaml file and convert to R6 configuration object

## Usage

``` r
read_BF_yamlR6(ymlfile, application = "DIANN")
```

## Arguments

- ymlfile:

  path to yaml configuration file

- application:

  software application name

## Value

ProlfquAppConfig R6 object

## Examples

``` r
if (FALSE) {
  yfile <- prolfqua::find_package_file("prolfquapp", "application/DIANN/config.yaml")
  file.exists(yfile)
  config <- read_BF_yamlR6(yfile)
}
```
