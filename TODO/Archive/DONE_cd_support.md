# TODO: CompoundDiscoverer DEA Support

## Summary

Add `inst/application/CMD_DEA_CD.R` for CompoundDiscoverer ZIP exports. The script should produce the same result family as `CMD_DEA_V2.R`: DEA/QC HTML reports, XLSX tables, Parquet normalized data, YAML metadata, `SummarizedExperiment.rds`, `index.html`, logs, and archived inputs.

## Input

- Support v1 input as a ZIP export containing exactly one `*_long.csv` and one `*_prolfqua_samples.csv`.
- Read long data with columns `Feature_ID`, `Sample`, `Intensity`, and `Group`.
- Read embedded sample annotation and normalize `Sample` to the package's existing annotation convention by renaming it to `file`.
- Use the embedded annotation as the authoritative sample design and contrast source.

## Implementation

- Add a CD-specific preprocessing helper that converts the long export into `LFQData` plus `ProteinAnnotation`.
- Use a one-level hierarchy with `metabolite_feature_Id` and response `Intensity`.
- Handle duplicated `Feature_ID x Sample` rows by appending a stable suffix to make affected feature IDs unique:
  - `Feature [duplicate 1]`
  - `Feature [duplicate 2]`
- Apply duplicate suffixes to all occurrences of features that need splitting, so the resulting matrix has stable rows across samples.
- Add a CD-specific DEA runner that reuses the existing `ProteinDataPrep`, `DEAnalyse`, and `DEAReportGenerator` pipeline.
- Treat columns after `Group` in the long table as optional subset definitions and run one separate report per subset column.
- Generate the full unfiltered report by default in addition to subset reports.
- Skip IBAQ for CD because it is protein/peptide specific.
- Make YAML optional:
  - use the supplied YAML if present
  - otherwise create a default DEA config with `vsn` and `lm_missing`

## CLI

- Add `inst/application/CMD_DEA_CD.R`.
- Support:
  - `-i/--input`: CD ZIP file or a directory containing one CD ZIP
  - `-y/--yaml`: optional config YAML
  - `-o/--outdir`: output directory
  - `-w/--workunit`: optional workunit override
  - `-m/--model`: optional model override
  - `--subset-columns`: `auto`, `none`, or a comma-separated list of long-table subset columns
  - `--include-full`: include the full unfiltered CD table alongside subset runs
  - `--libPath`: optional R library path

## Tests

- Unit-test duplicate ID handling.
- Integration-test the CD script with a small generated ZIP fixture.
- Verify expected outputs exist: result directory, XLSX, `index.html`, normalized parquet, normalized YAML, and `SummarizedExperiment.rds`.

## Assumptions

- Raw CompoundDiscoverer Excel files and the existing wide PCnorm CSV are out of scope for v1.
- ORA/GSEA files may still be emitted by the existing reporter, using metabolite feature IDs instead of protein IDs.
- Existing report wording may still say "protein" in places; v1 focuses on output generation and data compatibility.
