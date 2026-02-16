# extract contrast from annotation file

extract contrast from annotation file

## Usage

``` r
extract_contrasts(annot, prefix = "G_", group = "group")
```

## Arguments

- annot:

  annotation data frame

- prefix:

  prefix for group levels

- group:

  name of the group column

## Examples

``` r
annot <- data.frame(names = c("a1","b1"), group= c("a","b"), ddd = c("T","C"))
testthat::expect_error(extract_contrasts(annot))
#> INFO [2026-02-16 14:19:22] levels: c("a", "b")
annot$control <- annot$ddd
contrast <- extract_contrasts(annot)
#> INFO [2026-02-16 14:19:22] levels: c("a", "b") c("T", "C")
#> a b 
stopifnot(contrast == "G_a - G_b")

annot$Contrast <- c("G_a - G_b","G_b - G_a")
annot$ContrastName <- c("a_vs_b","b_vs_a")
annot$control <- NULL
ct <- extract_contrasts(annot)
#> INFO [2026-02-16 14:19:22] levels: c("a", "b")
stopifnot(length(ct) == 2)
```
