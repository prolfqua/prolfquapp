# build bfabric urls

build bfabric urls

## Usage

``` r
bfabric_url_builder(project_spec)
```

## Arguments

- project_spec:

  ProjectSpec R6 object with project, order, workunit IDs

## Examples

``` r

ps <- ProjectSpec$new()
ps$project_Id <- 32258
ps$order_Id <- 34628
ps$workunit_Id <- 302212
bfabric_url_builder(ps)
#> $orderURL
#> [1] "https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=34628&tab=details"
#> 
#> $projectURL
#> [1] "https://fgcz-bfabric.uzh.ch/bfabric/project/show.html?id=34628&tab=details"
#> 
#> $workunitURL
#> [1] "https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=34628&tab=details"
#> 

ps <- ProjectSpec$new()
ps$order_Id <- 34628
ps$workunit_Id <- 302212
bfabric_url_builder(ps)
#> $orderURL
#> [1] "https://fgcz-bfabric.uzh.ch/bfabric/order/show.html?id=34628&tab=details"
#> 
#> $projectURL
#> NULL
#> 
#> $workunitURL
#> [1] "https://fgcz-bfabric.uzh.ch/bfabric/workunit/show.html?id=34628&tab=details"
#> 
```
