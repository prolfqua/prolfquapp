# TODO: `nr_children` Consistency

## Original Request

The preprocess functions return an `LFQData` object, which has a field `nr_children` (for proteins these are the number of peptides per protein per sample quantified). Make sure that all the preprocess functions correctly fill this field.

Also make sure that this field is correctly represented in the AnnData object as a layer.

Furthermore, make sure that we correctly pass `nr_children` to the `ContrastFacade`. You might need to have a look at the `prolfqua_fml/prolfqua` package.

---

## Analysis & Plan: `nr_children` Consistency

### Current State

There are two distinct `nr_children` concepts:

1. **`config$nr_children` / `atable$nr_children`** — Per-sample column in `LFQData` tracking how many child observations were aggregated into each row (e.g., precursors -> peptide). Used by `value_vars()`, AnnData layers, `ContrastsDEqMSFacade`, and plotting.

2. **`exp_nr_children` in `ProteinAnnotation`** — Per-protein metadata: expected number of peptides across the experiment. Used for filtering (e.g., require >= 2 peptides).

### Issue 1: Preprocess functions not setting `atable$nr_children`

| Preprocessor   | Computes `nr_children` in data? | Sets `atable$nr_children`? | Status |
|----------------|--------------------------------|---------------------------|--------|
| DIANN          | Yes (`nr_children = n()`)       | Yes (`"nr_children"`)      | OK |
| MSstats        | Yes (`nr_children = dplyr::n()`)| Yes (`"nr_peptides"`)      | OK |
| MSstats_FPDIA  | Yes (`nr_children = dplyr::n()`)| Yes (`"nr_peptides"`)      | OK |
| CompoundDisc   | Yes (gap status encoding)       | No — relies on default `"nr_children"` matching column name | OK by accident |
| FP_PSM         | Computes `nr_psm` but doesn't map it | No                   | **MISSING** |
| MaxQuant       | No (input is already at peptide level) | No                  | OK — `setup_analysis()` adds `nr_children = 1`, which is correct for peptide-level input |
| BGS            | No                              | No                        | OK — data at finest granularity, `setup_analysis()` adds `nr_children = 1`, `.add_nr_children()` recomputes during aggregation |
| MzMine         | No                              | No                        | OK — metabolomics, single-level features, `nr_children = 1` is appropriate |

### Issue 2: AnnData layer representation

`lfq_to_ann_data.R` already correctly iterates over `config$value_vars()` (which includes `nr_children`) and creates AnnData layers. This works if `nr_children` is properly set in the config. The round-trip via `preprocess_DIANN_anndata.R` and `anndata_to_lfqdata.R` also handles it.

**Mechanism:** `nr_children` is included in `config$value_vars()` (defined in `prolfqua/R/AnalysisConfiguration.R:211-222`), and the layer creation loop in `lfq_to_ann_data.R:35-44` iterates over `value_vars()` — so `nr_children` automatically becomes an AnnData layer. The reverse path in `anndata_to_lfqdata.R:138-165` reads it back from the stored `layer_names` in `adata$uns$prolfquapp`.

**Requirement:** `nr_children` must always be present as a layer in the AnnData object, so that DEqMS weighting and downstream tools have access to peptide/precursor counts per protein per sample.

**Current status:** This works IF the preprocess functions properly populate the `nr_children` column in the data. If the column is missing from `lfqdata$data`, `to_wide()` will fail during layer creation.

**Note on `value_vars()` and empty strings:** Optional fields like `ident_Score`, `opt_mz`, `opt_rt` default to `character()` (zero-length vector). This is NOT a bug — R's `c()` automatically strips `character(0)` elements: `c("a", character(0), "b")` returns `c("a", "b")`, not `c("a", "", "b")`. So no empty strings enter the layer creation loop.

### Issue 3: `ContrastFacade` usage

`ContrastsDEqMSFacade` at `ContrastsFacades.R:591` uses `lfqdata$config$nr_children` to get the count column for DEqMS weighting. This works correctly after aggregation because `.add_nr_children()` in `tidyMS_aggregation.R:602` recomputes and updates the column name during medpolish/topN aggregation.

**No changes needed in prolfqua** — the aggregation step already handles this. The fix is in the preprocess functions.

---

## Proposed Changes

### Fix 1: `preprocess_FP_PSM.R` — Map `nr_psm` to `nr_children`

The FP_PSM preprocessor already computes `nr_psm` (PSMs per peptide per sample) during aggregation, but doesn't register it. Two options:

- **Option A:** Rename `nr_psm` to `nr_children` in the summarize call and set `atable$nr_children <- "nr_children"`
- **Option B:** Keep `nr_psm` and set `atable$nr_children <- "nr_psm"`

Recommendation: **Option B** — preserves the semantic meaning; `atable$nr_children` just points to the right column name.

Note: When `aggregate = FALSE`, there's no `nr_psm` column computed, so `setup_analysis()` will default to adding `nr_children = 1`, which is correct (each PSM row is a single observation).

### Fix 2: `preprocess_BGS_default.R` — Compute `nr_children` for precursors per peptide

BGS has a 3-level hierarchy (`protein_Id` / `peptide_Id` / `elution_group`). When `hierarchy_depth = 2` (peptide level), `nr_children` should count elution groups per peptide per sample. When `hierarchy_depth = 1`, `setup_analysis` defaults to 1, and aggregation will recompute.

Looking more closely: BGS data comes at the elution group level (`EG.ModifiedSequence` + `FG.Charge`), so at `hierarchy_depth = 2` (peptide), `setup_analysis()` adds `nr_children = 1` which is fine since each row is one elution group. The aggregation to protein level will then correctly count peptides. Similarly at `hierarchy_depth = 1`, each row is still one elution group, and the aggregation step handles it.

Note: `hierarchyDepth` controls which hierarchy keys are active, not which level the data starts at. The data is still at the finest granularity. So `nr_children = 1` is correct at the finest level, and `.add_nr_children()` handles the rest during aggregation.

---

## Revised Assessment

After deeper analysis:

- **FP_PSM with aggregation:** The `nr_psm` column IS computed but not mapped to `nr_children`. This means per-PSM count information is lost. **Fix:** set `atable$nr_children <- "nr_psm"` when aggregation occurs.
- **MaxQuant, BGS, MzMine:** Input data is at the finest hierarchy level, so `nr_children = 1` (default from `setup_analysis()`) is semantically correct. The aggregation functions (`.add_nr_children()`) properly recompute when rolling up.
- **CompoundDisc:** Works by accident because the column is named `nr_children` matching the default.

**The only real fix needed is in FP_PSM.**

---

## Action Items

1. **[DONE] Fix `preprocess_FP_PSM.R`:** Set `atable$nr_children <- "nr_psm"` conditionally when the `nr_psm` column exists (i.e., when aggregation occurred).
2. ~~Fix AnnData layer creation — filter empty strings from `value_vars()`:~~ **Not needed.** R's `c()` strips `character(0)` automatically — no empty strings enter the loop.
3. **AnnData `nr_children` layer:** Already works correctly — `nr_children` is in `value_vars()` and becomes a layer automatically. No code change needed.
4. **Verify** that MaxQuant/BGS work correctly with the default `nr_children = 1` since `.add_nr_children()` recomputes during protein aggregation.
