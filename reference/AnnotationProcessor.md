# AnnotationProcessor

AnnotationProcessor

AnnotationProcessor

## Public fields

- `QC`:

  is it a QC run

- `prefix`:

  name for one factor designs

- `repeated`:

  is it a repeated measurement

- `SAINT`:

  is it a AP MS experiment, then use Bait\_ as prefix

- `file_pattern`:

  colnames for file

- `grouping_pattern`:

  colnames grouping variable

- `subject_pattern`:

  colnames for pairing variable

- `control_pattern`:

  contrast specification columns

- `control_col_pattern`:

  columns which contains C or T.

- `sample_name_pattern`:

  sample name column

- `norm_value_pattern`:

  normalization value column (e.g., Creatinine)

- `strict`:

  should name check be strict

## Methods

### Public methods

- [`AnnotationProcessor$new()`](#method-AnnotationProcessor-new)

- [`AnnotationProcessor$check_annotation()`](#method-AnnotationProcessor-check_annotation)

- [`AnnotationProcessor$read_annotation()`](#method-AnnotationProcessor-read_annotation)

- [`AnnotationProcessor$extract_contrasts()`](#method-AnnotationProcessor-extract_contrasts)

- [`AnnotationProcessor$add_contrasts_vec()`](#method-AnnotationProcessor-add_contrasts_vec)

- [`AnnotationProcessor$clone()`](#method-AnnotationProcessor-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    AnnotationProcessor$new(
      QC = FALSE,
      prefix = "G_",
      repeated = TRUE,
      SAINT = FALSE
    )

#### Arguments

- `QC`:

  default FALSE

- `prefix`:

  default "G\_"

- `repeated`:

  default TRUE

- `SAINT`:

  default FALSE

------------------------------------------------------------------------

### Method `check_annotation()`

check annotation

#### Usage

    AnnotationProcessor$check_annotation(annot)

#### Arguments

- `annot`:

  annotation

------------------------------------------------------------------------

### Method [`read_annotation()`](https://prolfqua.github.io/prolfquapp/reference/read_annotation.md)

read annotation

#### Usage

    AnnotationProcessor$read_annotation(dsf)

#### Arguments

- `dsf`:

  either dataframe or file path.

------------------------------------------------------------------------

### Method [`extract_contrasts()`](https://prolfqua.github.io/prolfquapp/reference/extract_contrasts.md)

check annotation

#### Usage

    AnnotationProcessor$extract_contrasts(annot, group)

#### Arguments

- `annot`:

  annotation

- `group`:

  group column e.g. group

------------------------------------------------------------------------

### Method [`add_contrasts_vec()`](https://prolfqua.github.io/prolfquapp/reference/add_contrasts_vec.md)

add vector of contrasts to annot table

#### Usage

    AnnotationProcessor$add_contrasts_vec(annot, Contrasts)

#### Arguments

- `annot`:

  annotation

- `Contrasts`:

  vector with contrasts

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AnnotationProcessor$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# AnnotationProcessor$debug("read_annotation")
ap <- AnnotationProcessor$new(prefix = "G_")

annot <- data.frame(
file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
group = c("a","a","b","b"),
CONTROL = c("C","C","T","T"),
Subject = c("X","Y","X","Y"))
ap$check_annotation(annot)
#> Warning: column starting with :^name is missing.
af <- annot
af$file <- NULL
testthat::expect_error(ap$check_annotation(af), "column starting with :")
af <- annot
af$group <- NULL
testthat::expect_error(ap$check_annotation(af),"column starting with :")
#> Warning: column starting with :^name is missing.
aa <- ap$read_annotation(annot)
#> Warning: column starting with :^name is missing.
#> Registered S3 method overwritten by 'prolfqua':
#>   method         from    
#>   print.pheatmap pheatmap
#> INFO [2026-02-25 16:35:40] levels: c("a", "b") c("C", "T")
#> b a 
stopifnot(length(aa$atable$factor_keys_depth()) == 2)
stopifnot(all(c("atable", "annot", "contrasts") %in% names(aa)))
stopifnot(aa$contrasts == "G_b - G_a")
af <- annot
af$CONTROL <- NULL
testthat::expect_error(ap$check_annotation(af),"you must specify a CONTROL column")
#> Warning: column starting with :^name is missing.
af <- annot
af$Subject <- NULL
testthat::expect_warning(ap$check_annotation(af),"column starting with")


# should not throw exception since QC does not require group or subject
ap <- AnnotationProcessor$new(QC = TRUE)
af <- annot
# af$group <- NULL
af$CONTROL <- NULL
af$Subject <- NULL
ap$check_annotation(af)
#> Warning: column starting with :^name is missing.
aa <- ap$read_annotation(af)
#> Warning: column starting with :^name is missing.

stopifnot(aa$atable$factor_keys() == "G_")
stopifnot(aa$atable$factors == "group")
aa <- ap$read_annotation(annot)
#> Warning: column starting with :^name is missing.
aa$atable$fileName
#> [1] "file"
aa$atable$sampleName
#> [1] "sampleName"
as <- annot
as$sample <- c("s1","s2","s3","s4")
aa <- ap$read_annotation(annot)
#> Warning: column starting with :^name is missing.
aa$atable$sampleName
#> [1] "sampleName"
stopifnot(is.null(aa$annotation))

annot <- data.frame(
file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
Name = c("a1.raw","a2.raw","a3.raw","a4.raw"),
"Grouping Var" = c("a","a","b","b"),
CONTROL = c("C","C","T","T"),
Subject = c("X","Y","X","Y"))
ax <- ap$read_annotation(annot)
```
