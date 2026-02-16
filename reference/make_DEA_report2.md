# make DEA

make DEA

## Usage

``` r
make_DEA_report2(lfqdata, protAnnot, GRP2)
```

## Arguments

- lfqdata:

  LFQData object

- protAnnot:

  ProteinAnnotation object

- GRP2:

  ProlfquAppConfig configuration object

## Examples

``` r
# example code

pep <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
pep <- prolfqua::LFQData$new(pep$data, pep$config)
pA <- data.frame(protein_Id = unique(pep$data$protein_Id))
pA <- pA |> dplyr::mutate(fasta.annot = paste0(pA$protein_Id, "_description"))
pA <- prolfquapp::ProteinAnnotation$new(pep, row_annot = pA, description = "fasta.annot")
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
GRP2 <- prolfquapp::make_DEA_config_R6()

pep$factors()
#> # A tibble: 12 × 3
#>    sample  sampleName group_
#>    <chr>   <chr>      <chr> 
#>  1 A_V1    A_V1       A     
#>  2 A_V2    A_V2       A     
#>  3 A_V3    A_V3       A     
#>  4 A_V4    A_V4       A     
#>  5 B_V1    B_V1       B     
#>  6 B_V2    B_V2       B     
#>  7 B_V3    B_V3       B     
#>  8 B_V4    B_V4       B     
#>  9 Ctrl_V1 Ctrl_V1    Ctrl  
#> 10 Ctrl_V2 Ctrl_V2    Ctrl  
#> 11 Ctrl_V3 Ctrl_V3    Ctrl  
#> 12 Ctrl_V4 Ctrl_V4    Ctrl  
GRP2$pop <- list(Contrasts = c("AVsC" = "group_A - group_Ctrl", BVsC = "group_B - group_Ctrl"))
grp <- make_DEA_report2(pep, pA, GRP2)
#> INFO [2026-02-16 15:06:17] Transforming using log2
#> Column added : log2_abundance
#> INFO [2026-02-16 15:06:17] Transforming data : none.
#> abundance  already in data : sample sampleName group_ isotopeLabel protein_Id abundance qValue nr_peptides .
#> FORMULA :normalized_abundance ~ group_
#> Joining with `by = join_by(protein_Id)`
#> INFO [2026-02-16 15:06:17] fitted model with formula : normalized_abundance ~ group_
#> determine linear functions:
#> Warning: Warn 'linfct_matrix_contrasts':In argument: `AVsC = group_A - group_Ctrl`.
#> Warning: Warn 'linfct_matrix_contrasts':In argument: `BVsC = group_B - group_Ctrl`.
#> Warning: Warn 'linfct_matrix_contrasts':In argument: `avg_AVsC = (group_A + group_Ctrl)/2`.
#> Warning: Warn 'linfct_matrix_contrasts':In argument: `avg_BVsC = (group_B + group_Ctrl)/2`.
#> get_contrasts -> contrasts_linfct
#> contrasts_linfct
#> Joining with `by = join_by(protein_Id, contrast)`
#> completing cases
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> AVsC=group_A - group_Ctrl
#> BVsC=group_B - group_Ctrl
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id, contrast)`
#> Joining with `by = join_by(protein_Id)`
```
