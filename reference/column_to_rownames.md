# convert tibble to data.frame with rownames

convert tibble to data.frame with rownames

## Usage

``` r
column_to_rownames(.data, var = "rowname", sep = "~lfq~")
```

## Arguments

- .data:

  a tibble or data.frame

- var:

  name of the column with new row.names

- sep:

  separator for uniting columns

## Value

a data.frame with rownames

## Examples

``` r
ind <- tibble::tibble(a = 1:3, rowname = letters[1:3])
column_to_rownames(ind)
#>   a rowname
#> a 1       a
#> b 2       b
#> c 3       c
```
