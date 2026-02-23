# get configuration from yaml file or create default configuration

get configuration from yaml file or create default configuration

## Usage

``` r
get_config(yamlfile, WORKUNITID = "HelloWorld", ORDERID = "123")
```

## Arguments

- yamlfile:

  path to yaml configuration file (optional)

- WORKUNITID:

  workunit identifier for default configuration

- ORDERID:

  order identifier for default configuration

## Value

ProlfquAppConfig R6 object

## Examples

``` r

get_config()
#> <ProlfquAppConfig>
#>   Public:
#>     RES: list
#>     as_list: function () 
#>     clone: function (deep = FALSE) 
#>     ext_reader: ExternalReader, R6
#>     get_input_dir: function () 
#>     get_result_dir: function () 
#>     get_zipdir: function () 
#>     group: G_
#>     initialize: function (processing_options, project_spec, ext_reader, zipdir_name = ".", 
#>     path: .
#>     pop: list
#>     prefix: DEA
#>     processing_options: ProcessingOptions, R6
#>     project_spec: ProjectSpec, R6
#>     set_zipdir_name: function () 
#>     software: DIANN
#>     zipdir_name: DEA_20260223_PI123_O123_WUHelloWorld_none
```
