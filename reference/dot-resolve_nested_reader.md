# Resolve the reader for a (possibly nested) facade

Nested facades (e.g. `firth_nested`, `lmer`, `ropeca`) fit models on
peptide-level data and therefore require a peptide-level reader. Peptide
readers are registered as `"<reader>_PEPTIDE"` and differ from their
protein-level counterpart only by `hierarchy_depth`. When a nested
facade is paired with a protein-level reader, this transparently
switches the software key to the matching peptide-level reader (e.g.
`"prolfquapp.DIANN"` -\> `"prolfquapp.DIANN_PEPTIDE"`) instead of
failing. Non-nested facades, and readers that are already peptide-level,
are returned unchanged.

## Usage

``` r
.resolve_nested_reader(software, is_nested, available, facade = "")
```

## Arguments

- software:

  software key (e.g. "prolfquapp.DIANN")

- is_nested:

  logical; whether the facade needs peptide-level data

- available:

  character vector of registered software keys
  (`names(get_procfuncs())`)

- facade:

  facade name, used only for messages

## Value

the (possibly remapped) software key
