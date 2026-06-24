Subject: New employee-only DEA application A414

Hi all,

We now have a new DEA (A414) application for prolfqua differential expression analysis. For now it is intended for employees only.

It is very similar to the existing DEA_Spectronaut, DEA_FragPipe-DiaNN, DEA FragPipe-TMT, DEA - MaxQuant, and DEA FragPipe-DDA-RESOURCE applications. However, instead of separate DEA applications for different upstream tools, DEA now infers the input format from the selected upstream quantification result.

Currently this covers DIA-NN / FragPipe DIA, Spectronaut, and MaxQuant outputs.

The usual DEA settings remain available:

- annotation dataset
- normalization
- fold-change and FDR thresholds
- contaminant/decoy handling, now with a dropdown menu
- peptide-level mode
- model selection

Model selection is now more explicit; special models are reached by choosing `Model = Extra` and then selecting the desired `Model (extra)`.

Model options:

- `lm`: Linear model with variance shrinkage/moderation for standard differential expression analysis.
- `lm_impute`: DEFAULT Linear model with LOD-based missing-value imputation and borrowed variance for sparse data.
- `saint`: SAINT-style model for affinity-purification / protein-interaction experiments.
- `Extra`
    - `rfit`: Rank-based robust regression, useful when outliers may affect linear-model fits.
    - `firth`: Firth logistic model for presence/absence or separation-prone data.
    - `firth_nested`: Peptide-level nested Firth model summarized to protein-level results.
    - `ropeca_nested`: Peptide-level ROPECA ranking model, mainly for creating rank files for GSEA, not for differential-expression testing.
