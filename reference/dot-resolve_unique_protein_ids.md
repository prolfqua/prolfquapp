# Resolve duplicate protein IDs to one row each

Within each duplicated-id group: drop decoy rows when a forward exists
(keep the forward), then prefer reviewed `sp|` over `tr|`, else keep the
first. Guarantees one row per id and logs the counts. Standalone decoys
(no forward twin) are left untouched.

## Usage

``` r
.resolve_unique_protein_ids(row_annot, pID, full_id, pattern_decoys = NULL)
```

## Arguments

- row_annot:

  annotation data frame

- pID:

  name of the protein-id column

- full_id:

  name of the column carrying the raw, prefixed id

- pattern_decoys:

  optional configured decoy regex

## Value

`row_annot` with one row per `pID`
