# FragPipe read FragPipe combined protein files up to Version 15

FragPipe read FragPipe combined protein files up to Version 15

## Usage

``` r
tidy_FragPipe_combined_protein_deprec(
  combprot,
  intnames = c("total.intensity", "unique.intensity", "razor.intensity",
    "total.ion.count", "unique.ion.count", "razor.ion.count", "total.spectral.count",
    "unique.spectral.count", "razor.spectral.count"),
  protIDcol = "protein.group",
  subgroup = "subgroup",
  as_list = FALSE
)
```

## Arguments

- combprot:

  path to combined_protein.tsv file

- intnames:

  intensity column prefix

- protIDcol:

  default protein.group

- subgroup:

  default subgroup

## See also

Other FragPipe:
[`FragPipe`](https://prolfqua.github.io/prolfquapp/reference/FragPipe.md),
[`tidy_FragPipe_MSstats_csv()`](https://prolfqua.github.io/prolfquapp/reference/tidy_FragPipe_MSstats_csv.md),
[`tidy_FragPipe_combined_protein()`](https://prolfqua.github.io/prolfquapp/reference/tidy_FragPipe_combined_protein.md)

## Examples

``` r
prottsv <- prolfqua::find_package_file("prolfquapp", "samples/FragPipe/combined_protein_small.tsv")

prot <- tidy_FragPipe_combined_protein_deprec(prottsv)
#> annotation columns : protein.group
#> subgroup
#> protein
#> protein.id
#> entry.name
#> gene.names
#> protein.length
#> coverage
#> organism
#> protein.existence
#> description
#> protein.probability
#> top.peptide.probability
#> unique.stripped.peptides
#> summarized.total.spectral.count
#> summarized.unique.spectral.count
#> summarized.razor.spectral.count
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup, raw.file)`
#> Joining with `by = join_by(protein.group, subgroup)`
stopifnot( dim(prot) ==c(19980,27))
```
