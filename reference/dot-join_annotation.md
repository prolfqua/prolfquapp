# Enrich a quant/result table with the protein annotation

Right-joins so every row of `x` (the quant or result table) is preserved
and never multiplied.

## Usage

``` r
.join_annotation(annotation, x, hierarchy_keys)
```

## Arguments

- annotation:

  protein annotation data frame

- x:

  quant/result table to annotate

- hierarchy_keys:

  config hierarchy keys (e.g. `lfqdata$hierarchy_keys()`); the join uses
  the subset present in both frames

## Value

`x` enriched with annotation columns; one row per row of `x`

## Details

Joins on the \*\*hierarchy keys\*\* shared by both frames — the config's
feature-identity columns (e.g. `protein_Id`; plus deeper keys like
`site` for a PTM `protein_Id` + `site` hierarchy) — never on a
coincidentally-shared value column. `hierarchy_keys` is intersected with
the columns actually present in both frames, so a protein-level
annotation joined to protein-level contrasts uses `protein_Id`, while a
site-level analysis uses `protein_Id` + `site`. Joining on only
`protein_Id` when both carry `site` would suffix it to `site.x`/`site.y`
and drop the bare `site` the report needs
(\`unite(all_of(hierarchy_keys))\`). Keeps the row-preserving
`right_join` and a uniqueness guard on the resolved key(s).
