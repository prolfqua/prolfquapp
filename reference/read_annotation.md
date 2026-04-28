# read annotation files

read annotation files

## Usage

``` r
read_annotation(dsf, repeated = TRUE, SAINT = FALSE, prefix = "G_", QC = FALSE)
```

## Arguments

- dsf:

  annotation table

- repeated:

  is this a repeated measurement

- SAINT:

  is this a SAINTexpress analysis

- prefix:

  prefix for group levels

- QC:

  if TRUE, read as QC annotation

## Value

list with annot (annotation table), atable (analtysis table
configuration), contrasts list with contrasts.

## Examples

``` r
annot <- data.frame(
file = c("a1.raw","a2.raw","a3.raw","a4.raw"),
name = c("aa","ba","aa","ba"),
group = c("a","a","b","b"))
read_annotation(annot, QC = TRUE)
#> $atable
#> <AnalysisConfiguration>
#>   Public:
#>     annotation_vars: function () 
#>     bin_resp: 
#>     clone: function (deep = FALSE) 
#>     factor_depth: 1
#>     factor_keys: function () 
#>     factor_keys_depth: function () 
#>     factors: list
#>     file_name: file
#>     get_response: function () 
#>     hierarchy: list
#>     hierarchy_depth: 1
#>     hierarchy_keys: function (rev = FALSE) 
#>     hierarchy_keys_depth: function (names = TRUE) 
#>     id_required: function () 
#>     id_vars: function () 
#>     ident_q_value: qValue
#>     ident_score: 
#>     initialize: function () 
#>     is_response_transformed: FALSE
#>     isotope_label: isotopeLabel
#>     min_peptides_protein: 2
#>     norm_value: NULL
#>     nr_children: nr_children
#>     opt_mz: 
#>     opt_rt: 
#>     opt_se: 
#>     pop_response: function () 
#>     sample_name: name
#>     sep: ~
#>     set_response: function (col_name) 
#>     value_vars: function () 
#>     work_intensity: NULL
#> 
#> $annot
#>     file name group
#> 1 a1.raw aa_1     a
#> 2 a2.raw ba_1     a
#> 3 a3.raw aa_2     b
#> 4 a4.raw ba_2     b
#> 
```
