# Validate and normalize an \`nr_peptides\` threshold

Central guard for the minimum-peptides-per-protein threshold, applied
where the value enters a \`ProlfquAppConfig\` (native YAML, programmatic
config, and CLI override). Accepts clean YAML numerics (\`2\`, \`2.0\`)
and numeric-looking strings (\`"2"\`); rejects values that are not a
single whole number \`\>= 1\`.

## Usage

``` r
.validate_nr_peptides(x)
```

## Arguments

- x:

  candidate value (numeric, integer, or character scalar); \`NULL\`
  defaults to \`1L\`

## Value

the value coerced to a positive whole-number integer
