---
name: prolfquapp-dea
description: Run, set up, or troubleshoot the prolfquapp differential expression analysis CLI, especially prolfqua_dea.sh, prolfqua_dataset.sh, prolfqua_yaml.sh, prolfqua_contrasts.sh, prolfqua_qc.sh, dataset annotation design, CONTROL/ContrastName/Contrast columns, paired designs with subject/bioreplicate blocking, two-factor/factorial contrasts, and the prolfquapp R functions that copy or drive those scripts. Use when users ask how to get prolfqua_dea.sh into a working directory, prepare the annotation/YAML inputs, define contrasts, choose the software key, run DEA for DIA-NN, MaxQuant, FragPipe, Spectronaut/BGS, MSstats, MZMine, or debug missing scripts/config/report outputs.
---
# Prolfquapp DEA CLI

Use this skill for `prolfquapp` command-line differential expression analysis workflows.

The source of truth is the package `README.md`. In the `prolfqua_fml` ecosystem workspace, read
`prolfquapp/README.md` before giving detailed workflow advice. For code-level debugging, also inspect
`prolfquapp/inst/application/bin/prolfqua_dea.sh`, `prolfquapp/inst/application/CMD_DEA_V2.R`, and
`prolfquapp/R/copy_helpers.R`.

## Core Workflow

After quantification, put the search-engine output files and the FASTA used for the search into one input directory.
Then copy the CLI wrappers into the working directory:

```bash
R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
```

On Linux/macOS, make them executable:

```bash
chmod a+x prolfqua_*
```

This should create:

```text
prolfqua_dea.sh
prolfqua_dea_cd.sh
prolfqua_yaml.sh
prolfqua_qc.sh
prolfqua_dataset.sh
prolfqua_contrasts.sh
```

If Docker is being used, the README pattern is to prefix commands with `./prolfquapp_docker.sh`.

## Prepare Inputs

Generate an annotation template from the quantification output:

```bash
./prolfqua_dataset.sh -i data_dir/ -s prolfquapp.DIANN -d annotation.xlsx
```

Fill in the annotation before DEA. Without a design and at least one contrast definition, DEA cannot run. Typical
columns include:

- file identifier such as `Relative.Path`, `Path`, `raw.file`, or `channel`
- `name`
- main factor such as `group` or `experiment`
- optional blocking factor such as `subject` or `bioreplicate`
- optional `control`, where `C` is control and `T` is treatment

Create the YAML config:

```bash
./prolfqua_yaml.sh -y config.yaml
```

Edit the generated YAML for analysis parameters not exposed on the command line.

Optionally add contrast definitions:

```bash
./prolfqua_contrasts.sh annotation.xlsx --control WT -o annotation_with_control.xlsx
./prolfqua_contrasts.sh annotation.xlsx --f1 treatment --f2 time -o annotation_with_contrasts.xlsx
```

## Dataset Design And Contrasts

`prolfquapp::read_annotation()` detects columns by name pattern, not by exact spelling:

- file: starts with `channel`, `Relative`, `raw`, `file`, or `run`
- sample label: starts with `name`
- main group: starts with `group`, `bait`, or `Experiment`
- pairing/blocking: starts with `subject` or `BioReplicate`
- contrast definition: `ContrastName`, `Contrast`, or a column starting with `control`

The main group values are sanitized by removing whitespace and replacing `-+/*()` with `_`. Keep levels simple:
`WT`, `KO`, `T0`, `T150`, `MI_T0`, etc. The default model prefix is `G_`; contrast expressions must reference model
levels as `G_<level>`, for example `G_KO - G_WT`.

### Single-Factor Designs

For an unpaired one-factor design, use one grouping column and a `CONTROL` column. The CLI helper adds this for all
non-control levels:

```bash
./prolfqua_contrasts.sh annotation.xlsx --control WT -o annotation_with_control.xlsx
```

Example:

```text
file,name,group,CONTROL
s1.raw,WT_1,WT,C
s2.raw,WT_2,WT,C
s3.raw,KO_1,KO,T
s4.raw,OE_1,OE,T
```

This yields contrasts like `KO_vs_WT = G_KO - G_WT` and `OE_vs_WT = G_OE - G_WT`. If the group column is not named
`group`/`experiment`, pass it explicitly:

```bash
./prolfqua_contrasts.sh annotation.xlsx --control WT --group treatment -o annotation_with_control.xlsx
```

### Paired Or Blocked Designs

For paired analysis, keep the same group/CONTROL setup and add exactly one subject column named like `subject` or
`bioreplicate`. Each subject should have observations in the relevant groups. `read_annotation()` then adds the subject
factor and the model blocks on subject.

```text
file,name,group,subject,CONTROL
s1.raw,S01_WT,WT,S01,C
s2.raw,S01_KO,KO,S01,T
s3.raw,S02_WT,WT,S02,C
s4.raw,S02_KO,KO,S02,T
```

For unpaired designs, delete the `subject`/`bioreplicate` column instead of filling it with dummy values.

### Two-Factor / Factorial Designs

For a two-factor design, start with two factor columns and run:

```bash
./prolfqua_contrasts.sh annotation.xlsx --f1 treatment --f2 time -o annotation_with_contrasts.xlsx
```

The helper creates a united `Group` column plus `ContrastName` and `Contrast` columns. For `treatment = MI/MINOCA` and
`time = T0/T150`, the generated model levels are `G_MI_T0`, `G_MI_T150`, `G_MINOCA_T0`, and `G_MINOCA_T150`.

Typical generated rows look like:

```text
ContrastName,Contrast
MINOCA_vs_MI,( (G_MINOCA_T0 + G_MINOCA_T150)/2 - (G_MI_T0 + G_MI_T150)/2 )
MINOCA_vs_MI_at_T0,G_MINOCA_T0 - G_MI_T0
MINOCA_vs_MI_at_T150,G_MINOCA_T150 - G_MI_T150
interaction_MINOCA_vs_MI_at_T150_vs_T0,(G_MINOCA_T150 - G_MI_T150) - (G_MINOCA_T0 - G_MI_T0)
```

Use `--interactions FALSE` when only main-effect and level-specific contrasts are wanted:

```bash
./prolfqua_contrasts.sh annotation.xlsx --f1 treatment --f2 time --interactions FALSE -o annotation_with_contrasts.xlsx
```

### Manual Contrast Columns

If `ContrastName` and `Contrast` are present, prolfquapp uses those rows and does not derive contrasts from `CONTROL`.
Only rows with non-empty `Contrast` matter; remaining sample rows can be blank. This is the right route for custom
factorial comparisons.

Rules:

- `ContrastName` is the output/report name.
- `Contrast` is an R expression over `G_`-prefixed model levels.
- Every referenced level must exist in the annotation's main group column after sanitization.
- If a level name contains syntax-sensitive characters, fix the upstream group level; do not patch generated wrappers.

## Run DEA

The minimal DEA command is:

```bash
./prolfqua_dea.sh -i data_dir/ -d annotation.xlsx -y config.yaml -w NameOfAnalysis -s prolfquapp.DIANN
```

Inputs:

- `-i`, `--indir`: quantification output directory
- `-d`, `--dataset`: annotation file, usually from `prolfqua_dataset.sh` or `prolfqua_contrasts.sh`
- `-y`, `--yaml`: prolfquapp YAML configuration
- `-w`, `--workunit`: workunit or analysis name
- `-s`, `--software`: software key (see below)
- `-o`, `--outdir`: optional output directory
- `-m`, `--model`: optional model/facade override

Every wrapper supports `--help`. When unsure about a flag, run `./prolfqua_dea.sh --help` for the authoritative,
installed-version list rather than relying on the summary above, which can drift from the code.

Expected output is a folder starting with `DEA_` containing HTML reports, Excel tables, rank files, normalized data, and
`SummarizedExperiment.rds`.

## Software Keys

The `-s` value must be a key from `prolfquapp::get_procfuncs()`. These keys are **package-namespaced** — the CLI checks
`software %in% names(get_procfuncs())`, and `get_procfuncs()` prefixes each preprocessor with its providing package. So
the accepted form for DIA-NN is `prolfquapp.DIANN`, not bare `DIANN`; the bare name never matches and the run aborts
with `Software '...' not found. Available: ...`. The same namespaced lookup is used by both `prolfqua_dataset.sh` and
`prolfqua_dea.sh`.

Do not guess keys. List the exact accepted values from the installed package:

```bash
R --vanilla -e "print(names(prolfquapp::get_procfuncs()))"
```

The README-level workflow mentions DIA-NN, MaxQuant, FragPipe-TMT, FragPipe-DIA, FragPipe-LFQ, and Spectronaut/BGS, but
the authoritative, install-specific keys (and their package prefixes) come from `get_procfuncs()`. Downstream packages
that register their own preprocessors appear under their own prefix (for example `prolfquasaint.*`).

## If `prolfqua_dea.sh` Is Missing

First check whether `prolfquapp` is installed and where R sees it:

```bash
R --vanilla -e "cat(system.file(package = 'prolfquapp'), '\n')"
```

Then copy wrappers again:

```bash
R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
```

For a direct path check:

```bash
R --vanilla -e "cat(system.file('application/bin/prolfqua_dea.sh', package = 'prolfquapp'), '\n')"
```

`prolfqua_dea.sh` is only a wrapper. It resolves the installed package directory and runs
`application/CMD_DEA_V2.R`. If the wrapper exists but fails to find the R script, fix the package installation or the
installed package contents; do not patch the copied wrapper as a workaround.

## Debugging Rules

Fix the root cause in the upstream package or input file:

- Missing wrapper: use `prolfquapp::copy_shell_script()` or fix package installation.
- Missing `CMD_DEA_V2.R`: reinstall/build `prolfquapp` correctly.
- Bad software key: inspect `names(prolfquapp::get_procfuncs())` and use the full namespaced key (e.g.
  `prolfquapp.DIANN`); do not pass the bare name.
- Empty contrast group ("group '...' has 0 samples after matching the annotation to the quantification data", or, on
  older builds, a raw `subscript out of bounds` from `linfct_from_model`): a contrast references a model level whose
  samples are not present in the quantification report. The annotation and the report disagree — the report is missing
  the raw files for that group. Fix the input data (correct the annotation so every contrast level maps to runs that
  exist in the report, or re-run quantification to include the missing runs). Do not silently drop the contrast or the
  group.
- Bad annotation: fix the annotation file generated from the dataset step; do not add skip logic.
- Bad YAML: regenerate with `prolfqua_yaml.sh` or repair the config field that drives the failing R6 object.
- Missing report templates: rebuild/install `prolfquapp` with vignettes/docs or run from the source root if that is what
  the package fallback expects.

For local ecosystem development, prefer the repository Makefile install flow over ad hoc installs unless doing a
one-off diagnostic:

```bash
make install
```

from the `prolfquapp` package directory, or root-level:

```bash
make installs
```
