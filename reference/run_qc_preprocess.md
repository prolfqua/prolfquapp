# Preprocess quantification data for QC reporting

Resolves configuration (from YAML or defaults), reads the annotation,
and runs software-specific preprocessing. Returns everything needed to
construct a
[`QC_generator`](https://prolfqua.github.io/prolfquapp/reference/QC_generator.md).

## Usage

``` r
run_qc_preprocess(
  indir,
  dataset,
  software,
  yaml_file = NULL,
  outdir = "qc_dir",
  project = "",
  order = "",
  workunit = ""
)
```

## Arguments

- indir:

  directory containing quantification output files

- dataset:

  path to annotation CSV/TSV/XLSX

- software:

  software key (e.g. "DIANN", "SIM"); looked up in
  [`prolfqua_preprocess_functions`](https://prolfqua.github.io/prolfquapp/reference/prolfqua_preprocess_functions.md)

- yaml_file:

  optional path to config YAML; if it exists, used instead of generating
  a default config

- outdir:

  output directory (used only when no yaml_file)

- project:

  project ID (used only when no yaml_file)

- order:

  order ID (used only when no yaml_file)

- workunit:

  workunit ID (used only when no yaml_file)

## Value

list with `xd` (preprocessed data), `files` (discovered file paths), and
`config` (ProlfquAppConfig)
