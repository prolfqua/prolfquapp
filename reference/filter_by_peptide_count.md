# Drop proteins supported by fewer than \`nr_peptides\` distinct peptides

Reader-local minimum-peptides-per-protein filter. Counts the number of
distinct \`peptide_col\` values per \`protein_col\` in \`data\` and
keeps only the rows whose parent protein reaches \`nr_peptides\`. The
reader supplies the columns because only the reader knows which column
holds the stripped (unmodified) peptide sequence.

## Usage

``` r
filter_by_peptide_count(data, protein_col, peptide_col, nr_peptides = 1)
```

## Arguments

- data:

  long-format quant table (before \`setup_analysis\`)

- protein_col:

  name of the parent protein column (single column)

- peptide_col:

  name of the (stripped) peptide column to count (single column)

- nr_peptides:

  minimum number of distinct peptides per protein (\>= 1)

## Value

\`data\` with the rows of under-supported proteins removed

## Details

Filtering the quant/peptide table is sufficient: \`ProteinAnnotation\`
is right-joined onto the \`LFQData\` protein set (it seeds \`row_annot\`
from the \`LFQData\` proteins and left-joins annotation onto it), so the
annotation, the logged protein summary, and IBAQ all follow the filtered
quant automatically – see the "ProteinAnnotation is annotation, not a
filter" invariant in \`AGENTS.md\`. Apply this to the long table
\*before\* \`setup_analysis()\` / \`LFQData\$new()\`.

A no-op when \`nr_peptides \<= 1\` (or \`NULL\`), so default runs are
unaffected.

## Examples

``` r
d <- data.frame(
  prot = c("A", "A", "A", "B", "C", "C"),
  pep = c("p1", "p2", "p2", "p1", "p1", "p2")
)
# A has 2 distinct peptides, B has 1, C has 2
nrow(filter_by_peptide_count(d, "prot", "pep", 2)) # drops B
#> INFO [2026-07-03 08:59:30] nr_peptides filter (>= 2 distinct pep per prot): kept 2 / 3 proteins (dropped 1)
#> [1] 5
nrow(filter_by_peptide_count(d, "prot", "pep", 1)) # keeps all
#> [1] 6
```
