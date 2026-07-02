# Add \`nr_peptides\` to reader args only for readers that declare it

Each reader owns the minimum-peptides-per-protein count against its own
stripped-peptide column, so the threshold is forwarded only to readers
whose formals declare an \`nr_peptides\` argument. A reader without it
is left unfiltered; a warning is raised when the caller asked for
filtering (\`nr_peptides \> 1\`) so the skip is explicit rather than a
silent no-op.

## Usage

``` r
.forward_nr_peptides(base_args, preprocess_fn, nr_peptides)
```

## Arguments

- base_args:

  named list of arguments passed to the reader

- preprocess_fn:

  the resolved reader function

- nr_peptides:

  minimum distinct peptides per protein

## Value

\`base_args\`, with \`nr_peptides\` added iff the reader supports it
