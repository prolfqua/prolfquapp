# Differential expression analysis from peptide input to protein output

Differential expression analysis from peptide input to protein output

Differential expression analysis from peptide input to protein output

## Details

Runs facades that consume peptide-level measurements but emit
protein-level contrasts. The input `LFQData` keeps peptide hierarchy
columns, while its active hierarchy depth points to the protein level.

## Super class

[`prolfquapp::DEAnalyse`](https://prolfqua.github.io/prolfquapp/reference/DEAnalyse.md)
-\> `DEAnalysePeptideToProtein`

## Methods

### Public methods

- [`DEAnalysePeptideToProtein$build_default()`](#method-DEAnalysePeptideToProtein-build_default)

- [`DEAnalysePeptideToProtein$clone()`](#method-DEAnalysePeptideToProtein-clone)

Inherited methods

- [`prolfquapp::DEAnalyse$build_facade()`](https://prolfqua.github.io/prolfquapp/html/DEAnalyse.html#method-DEAnalyse-build_facade)
- [`prolfquapp::DEAnalyse$filter_contrasts()`](https://prolfqua.github.io/prolfquapp/html/DEAnalyse.html#method-DEAnalyse-filter_contrasts)
- [`prolfquapp::DEAnalyse$get_annotated_contrasts()`](https://prolfqua.github.io/prolfquapp/html/DEAnalyse.html#method-DEAnalyse-get_annotated_contrasts)
- [`prolfquapp::DEAnalyse$initialize()`](https://prolfqua.github.io/prolfquapp/html/DEAnalyse.html#method-DEAnalyse-initialize)

------------------------------------------------------------------------

### Method `build_default()`

Build the default peptide-to-protein facade.

#### Usage

    DEAnalysePeptideToProtein$build_default()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DEAnalysePeptideToProtein$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
