# add vector of contrasts to annotation data frame

add vector of contrasts to annotation data frame

## Usage

``` r
add_contrasts_vec(xx, Contrasts)
```

## Arguments

- xx:

  annotation data frame

- Contrasts:

  character vector of contrasts

## Examples

``` r
annot <- data.frame(Group = rep(c("A","B","C"), each = 3))
annot$Name
#> NULL
```
