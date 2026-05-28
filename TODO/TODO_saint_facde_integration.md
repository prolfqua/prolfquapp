can we review, the saint_model.R file what is it for? Why do we need it? why is it in prolfquapp and not prolfquasaint?

Does the code depend on prolfquapp? 

To me this looks a bit like glue code... nothing wrong .. but It also might point to a design problem of the DEAnalyse class... 

But we did not had such a problem with the other facede members, 

So it might also point to a problem with the SaintExpress implementation.

I even thnink we a might have a prolfqua skill for integrating new models into prolfqua.

Can you investigate: 1. compatibility of saint express model with prolfqua modelling interface.
2. Properly integrating saintexpress into the prolfquafacade. 

3. integrating of saintexpress into the DEAnalyser through the facade

---

# Investigation & Plan

## Answers to the framing questions

**What is `prolfquapp/R/saint_model.R` for?** It is glue code that pulls
SAINT-specific bait/control/length/gene columns out of an `LFQData` object
+ row annotation, builds the SAINT input tables via
`prolfquasaint::protein_2localSaint()`, runs SAINTexpress via
`prolfquasaint::runSaint()`, and wraps the result in
`prolfquasaint::ContrastsSAINTexpress`. It is consumed from
`R6_DEAnalyse$build_facade()` when `name == "saint"`.

**Why is it in `prolfquapp` and not `prolfquasaint`?** Historical: the
glue resolves columns chosen by `prolfquapp`'s annotation reader
(`read_annotation(SAINT = TRUE)`), and accesses fields on the prolfquapp
`rowAnnot` object. `prolfquasaint` deliberately does not depend on
`prolfquapp` (the dependency direction is the opposite, per
[TODO/Archive/TODO_saint_model_path.md](Archive/TODO_saint_model_path.md)).
Pulling the glue into `prolfquasaint` is fine — the LFQData inputs it
needs come from `prolfqua` and the column resolution is generic.

**Does the code depend on `prolfquapp`?** Only weakly: it reads
`row_annot$row_annot` (a plain data.frame) and uses base/dplyr ops on
`LFQData`. No `prolfquapp` classes are referenced inside
`.build_saint_contrast_result()` itself, so it can move to
`prolfquasaint` once we standardise on passing a plain row-annotation
data.frame.

**Is this a design problem in `DEAnalyse`?** Yes. Every other facade
member is consumed uniformly:
`facade_class$new(lfqdata, modelstr, contrasts)` → `get_contrasts()` /
`get_Plotter()` / `to_wide()` / `get_missing()` with a stable column
schema. SAINT bypasses that pattern because (a) the registry in
`prolfqua::FACADE_REGISTRY` does not know about SAINT, and (b)
`ContrastsSAINTexpress` returns SAINT-native columns
(`Bait`/`log2_EFCs`/`SaintScore`/`BFDR`) instead of the standard
prolfqua schema (`contrast`/`diff`/`p.value`/`FDR`/`statistic`/...).
That forces every downstream consumer to special-case SAINT.

**Is it a problem with the SAINTexpress implementation?** Partially.
`ContrastsSAINTexpress` already inherits from
`prolfqua::ContrastsInterface`, but it only honours the interface in
spirit — the column names diverge from what the
`prolfqua-adding-models` skill specifies. Bringing SAINT into compliance
is a small adapter change inside `prolfquasaint`.

## Current state (what is in the tree today)

- `prolfquapp/R/facade_model.R`: `.is_saint_model()`,
  `.resolve_facade_model()`, `.valid_facade_models()` extend
  `prolfqua::FACADE_REGISTRY` with `"saint"` outside the registry.
- `prolfquapp/R/saint_model.R::.build_saint_contrast_result()`: SAINT
  glue, returns a list `{contrast, input, result, protein_id}`.
- `R6_DEAnalyse$build_facade()`: explicit `if (.is_saint_model(name))`
  branch that calls the glue and stores `saint_input`/`saint_result`
  on the `DEAnalyse` object.
- `R6_DEAnalyse$get_annotated_contrasts()`, `filter_contrasts()`,
  `filter_significant_contrasts()`: branch on column names
  (`BFDR`/`log2_EFCs` vs `FDR`/`diff`).
- `R6_DEAReportGenerator`: branches in `prep_result_list()`,
  `write_DEA()`, `contrasts_to_Grob()`, `write_DEA_all()` (skip QC
  HTML), `make_SummarizedExperiment()` (`contrast_column` choice).
- `R6_ProteinDataPrep` and `cmd_helpers`: branch on
  `.is_saint_model()` to control annotation reading (`SAINT = TRUE`).
- `prolfquasaint::ContrastsSAINTexpress`: implements `get_contrasts()`,
  `get_Plotter()`, `to_wide()`, plus SAINT-specific `get_ora()`,
  `get_rank()`. Columns: `Prey`, `modelName`, `Bait`, `avgAbd`,
  `log2_EFCs`, `SaintScore`, `BFDR`.

## Goal

Make SAINT a first-class facade so that:

```r
fa <- build_contrast_analysis(lfqdata, modelstr = NULL, contrasts = NULL,
                              method = "saint", row_annot = ...)
fa$get_contrasts()   # standard prolfqua schema + extra SAINT columns
fa$to_wide()
fa$get_Plotter()
fa$get_missing()
```

works through the same code path as `lm`, `limma`, `firth`, etc.
Downstream code in `prolfquapp` should drop every `.is_saint_model()`
branch except for two places: (a) `read_annotation(SAINT = TRUE)` (an
input concern, not a modelling concern), and (b) skipping the legacy
QC HTML report (a report capability, decided from the facade's class or
its registry entry, not a name match).

## Design

### 1. Add `ContrastsSAINTFacade` in `prolfquasaint`

New file `prolfquasaint/R/ContrastsSAINTFacade.R`, mirroring
`prolfqua::ContrastsLimmaFacade`:

- fields: `model`, `contrast`, `.lfqdata`, `.contrast_names`, plus
  SAINT artifacts `saint_input`, `saint_result` for the report
  generator.
- `initialize(lfqdata, modelstr = NULL, contrasts = NULL, row_annot,
  spc = FALSE, engine = "r")`:
  - validate aggregated input via the same assertion used by other
    aggregated facades (`subject_id == hierarchy_keys`).
  - move `.saint_prepare_data()`, `.saint_bait_col()`,
    `.saint_control_col()`, and `.saint_first_existing()` from
    `prolfquapp/R/saint_model.R` into `prolfquasaint` (private
    helpers next to the facade). They are SAINT-specific.
  - run `protein_2localSaint()` → `runSaint()` → wrap result in
    `ContrastsSAINTexpress`; remap `Prey` to the project hierarchy
    key (as `.build_saint_contrast_result()` already does).
  - accept `modelstr`/`contrasts` for facade uniformity and warn if
    a caller passes a non-`NULL` value.
- `get_contrasts()`: return a data.frame translating
  `ContrastsSAINTexpress$get_contrasts()` to the prolfqua standard:
  - `<subject_id>` (protein id)
  - `modelName` ← `"ContrastSaint"`
  - `contrast` ← `Bait`
  - `avgAbd` ← `avgAbd`
  - `diff` ← `log2_EFCs`
  - `FDR` ← `BFDR`
  - `statistic` ← `SaintScore` (documented meaning)
  - `p.value` ← `NA_real_` (SAINT is not a p-value test; documented)
  - `std.error`, `df`, `conf.low`, `conf.high`, `sigma` ← `NA_real_`
  - SAINT-native columns (`Bait`, `SaintScore`, `BFDR`, `log2_EFCs`,
    `FoldChange`, `AvgSpec`/`AvgIntensity`) retained as extras so
    the SAINT XLSX sheets keep rendering with familiar names.
  - prepend `facade = "saint"`.
- `to_wide()`, `get_Plotter()`: delegate to `ContrastsSAINTexpress`
  (already produces the right shape and named scores).
- `get_missing()`: reuse the `.compute_missing()` pattern from
  `prolfqua/R/ContrastsFacades.R`; subject ids and contrast names
  come from `self$.lfqdata$data_long()` and `get_contrast_sides()`.
- `get_ora(up, FDR_threshold, diff_threshold)` and
  `get_rank(score)`: delegate to the wrapped
  `ContrastsSAINTexpress`. These extend the facade API (other
  facades do not declare them); see § 4 for the consumption pattern
  in `prolfquapp`.

### 2. Make the facade registry extensible

`prolfqua::FACADE_REGISTRY` is a hard-coded named list. To register
SAINT without `prolfqua` learning about it, add in `prolfqua`:

- a private `.facade_registry` environment seeded with
  `FACADE_REGISTRY`.
- exported `register_facade(name, class, needs, package)` and
  `lookup_facade(name)`.
- rewrite `build_contrast_analysis()` to use `lookup_facade()`
  instead of indexing `FACADE_REGISTRY` directly.

Then in `prolfquasaint/R/zzz.R::.onLoad()`:

```r
prolfqua::register_facade("saint",
                          class = "ContrastsSAINTFacade",
                          needs = "aggregated",
                          package = "prolfquasaint")
```

If touching `prolfqua` is undesirable for this slice, fall back to a
local `prolfquapp` helper `.facade_lookup(name)` that first checks
`prolfqua::FACADE_REGISTRY` and then a small `prolfquapp`-owned
registry containing only `saint`. The downstream cleanup in § 4 is
the same either way.

### 3. Drop the SAINT short-circuit in `R6_DEAnalyse`

In `prolfquapp/R/R6_DEAnalyse.R::build_facade()`:

- remove the `.is_saint_model(name)` branch.
- always go through the registry: resolve `entry$class` and
  `entry$package`, fetch via
  `utils::getFromNamespace(entry$class, entry$package %||% "prolfqua")`.
- always call
  `facade_class$new(self$lfq_data, modelstr, self$contrasts, ...)`.
  For SAINT, pass `row_annot = self$rowAnnot` and `spc`/`engine`
  through named arguments — only the SAINT facade declares them.
- after `facade_class$new()` returns, if `facade$saint_input` is
  non-`NULL`, mirror it onto the `DEAnalyse` object so existing
  reporting code that reads `dea$saint_input`/`dea$saint_result`
  keeps working. Mark these fields deprecated in the roxygen.

Drop the `if (!.is_saint_model(default_model)) { stopifnot(length(contrasts) >= 1) }`
guard: make `contrasts` optional in the facade signature and let
the SAINT facade ignore it.

### 4. Replace `.is_saint_model()` branches with schema-driven code

Once `get_contrasts()` returns the standard schema, the branches in
`R6_DEAnalyse`, `R6_DEAReportGenerator`, and `R6_ProteinDataPrep`
collapse:

- `get_annotated_contrasts()`: drop the SAINT branch; filter on
  `FDR < FDR_threshold & |diff| > diff_threshold`. Use
  `facade$get_ora()` when the facade defines it (capability check
  via `is.function(facade$get_ora)`), else fall back to the generic
  threshold filter.
- `filter_contrasts()` / `filter_significant_contrasts()`: same.
- `contrasts_to_Grob()`: pick `contrast`, `diff`, `statistic`, `FDR`
  from the standardized schema. Optionally also include
  `SaintScore` / `BFDR` when present (gate on column presence, not
  on a model-name match).
- `make_SummarizedExperiment()`: drop the `contrast_column`
  branch — always `"contrast"`.
- `write_DEA()`'s `saint` flag becomes
  `inherits(contrast_obj, "ContrastsSAINTFacade")` or is replaced
  by polymorphism: `.write_ORA()` and `.write_GSEA()` ask the
  facade for `get_ora()` / `get_rank()` when available and use the
  defaults otherwise.
- `write_DEA_all()`'s QC-HTML skip becomes a facade capability flag
  (e.g. `facade$supports_dea_qc` defaulting to `TRUE`, `FALSE` on
  the SAINT facade) rather than a name match.
- `prep_result_list()`'s `saint_inter`/`saint_prey`/`saint_bait`/
  `saint_list` block: read them from the facade fields
  (`facade$saint_input`, `facade$saint_result`) and add them only
  when those fields are non-`NULL`.
- `R6_ProteinDataPrep` and `cmd_helpers`: keep the
  `read_annotation(SAINT = TRUE)` branch (an input concern). Drive
  it from a facade-registry attribute (e.g.
  `entry$needs_saint_annotation`) rather than a model-name match.

### 5. Remove `prolfquapp/R/facade_model.R` and `prolfquapp/R/saint_model.R`

After § 1–4 land:

- `.build_saint_contrast_result()` and the `.saint_*` helpers move
  to `prolfquasaint`.
- `.is_saint_model()`, `.resolve_facade_model()`,
  `.valid_facade_models()` are no longer needed. The
  `"prolfqua" → "lm"/"lm_impute"` legacy alias in
  `.resolve_facade_model()` moves to either an alias table in the
  registry or a one-line resolver in `R6_DEAnalyse$initialize`.
- delete both files.

### 6. Tests

In `prolfquasaint`:

- `tests/testthat/test-ContrastsSAINTFacade.R`:
  - construction from a small simulated `LFQData` with bait/control
    annotation columns.
  - `get_contrasts()` returns the documented standard columns plus
    the SAINT-native extras.
  - `to_wide()`, `get_Plotter()`, `get_missing()` work.
  - `get_ora()` and `get_rank()` return non-empty results on the
    fixture.
  - registration: `prolfqua::lookup_facade("saint")` resolves the
    class when `prolfquasaint` is loaded.

In `prolfquapp`:

- update `tests/testthat/test-saint-model.R` to drive everything
  through `build_contrast_analysis(method = "saint", ...)` and
  `R6_DEAnalyse$build_facade("saint")`.
- focused test that `R6_DEAnalyse$build_facade("saint")` populates
  `dea$saint_input` / `dea$saint_result` (transitional
  compatibility) and exposes the standard contrast schema.
- regression test that confirms ORA / GSEA file presence and QC
  HTML absence for the SAINT path through the generic capability
  checks (not `.is_saint_model()`).

### 7. Documentation

- Roxygen for `ContrastsSAINTFacade` must state that `p.value` is
  `NA` and that `FDR` is BFDR / `statistic` is SaintScore.
- Update `prolfqua-adding-models` skill examples if they reference
  `FACADE_REGISTRY` directly to use `register_facade()`.
- Note in `prolfquapp/CLAUDE.md` that SAINT is reached via the
  facade registry just like `lm`/`limma`.

## Execution order

1. Land the facade + registration in `prolfquasaint` (§ 1–2), with
   `prolfqua::register_facade()` if § 2's preferred path is chosen.
   Ship the new `prolfquasaint` tests.
2. Switch `R6_DEAnalyse$build_facade()` to the registry path (§ 3)
   but keep the `.is_saint_model()` helper as a transitional shim
   that reads off the facade. Verify all current `prolfquapp` tests
   still pass.
3. Strip `.is_saint_model()` from every reporting / data-prep site
   (§ 4) one file at a time, each behind a focused test.
4. Delete `facade_model.R` and `saint_model.R` (§ 5) and move the
   YAML `"prolfqua"` alias to a single resolver.
5. Update docs (§ 7) and run `make installs && make tests &&
   make check-fast` from the workspace root.

---

## Course correction (supersedes § 1 schema mapping and § 4 column gymnastics)

User feedback on the first draft:

> Why does `DEAnalyse` need to know about the columns in the Contrast
> classes? Contrast classes should implement methods like `get_ora()` /
> `get_rank()` so that `DEAnalyse` does not see the internal
> implementation. The ~15 `.is_saint_model()` branches point to code
> that does not use an interface.

This reframes the problem. The right fix is not "make SAINT return the
same column names as LM"; it is "make every consumer call interface
methods so the column names stop mattering." The branches survive only
because `DEAnalyse` / `DEAReportGenerator` read columns
(`BFDR`/`FDR`/`Bait`/`contrast`/`log2_EFCs`/`diff`) directly. Replace
the reads with method calls and the branches disappear.

### Revised contract for `ContrastsInterface` / facades

Every contrast/facade implementation must expose:

- `get_contrasts(all = FALSE)` — backend-specific tidy table (kept for
  XLSX export and debugging). Consumers should not select columns from
  it by name; they call the methods below.
- `filter_significant(FDR_threshold, diff_threshold, ...)` — return
  rows that pass the backend's notion of "significant." LM facades
  filter on `FDR` and `|diff|`; SAINT filters on `BFDR` and
  `|log2_EFCs|`. The threshold semantics stay the same; the column
  resolution is the backend's job.
- `get_ora(up, FDR_threshold, diff_threshold)` — backend-specific
  selection of features for ORA, already in `ContrastsSAINTexpress`.
  Add a default implementation on the prolfqua side that falls back to
  `filter_significant()` + `diff > 0` / `< 0`.
- `get_rank(score = NULL)` — backend-specific rank table for GSEA. LM
  facades default to `signed = sign(diff) * -log10(p.value)`; SAINT
  uses `log2_EFCs`. Already on `ContrastsSAINTexpress`; lift to the
  interface with a default.
- `get_Plotter()`, `to_wide()`, `get_missing()` — already polymorphic.
- `contrast_summary_table()` — small per-protein data.frame consumed
  by `contrasts_to_Grob()` for the per-protein boxplot side tables.
  Returns generically named columns (e.g. `contrast`, `effect`,
  `score`, `fdr`) already rounded; LM facade fills it from
  `diff`/`statistic`/`FDR`, SAINT facade fills it from
  `log2_EFCs`/`SaintScore`/`BFDR`. `DEAReportGenerator` stops knowing
  what those columns are.
- `supports_dea_qc()` — capability flag. Default `TRUE`; SAINT facade
  overrides to `FALSE` so `write_DEA_all()` skips the legacy
  differential-expression QC HTML without naming SAINT.
- `contrast_column_name()` — the column to split contrast tables by in
  `make_SummarizedExperiment()`. Default `"contrast"`; SAINT facade
  overrides to `"Bait"`. This is the one place where a column name
  legitimately crosses the interface, because SE rowData layout is
  itself a contract.

Where the method does not yet exist on `prolfqua` contrast classes,
add it with a default implementation that uses the standard schema.
Once each method is in place, the corresponding `.is_saint_model()`
branch in `prolfquapp` becomes a single method call against the
contrast object.

### Mapping `.is_saint_model()` branches to the new methods

| Site (prolfquapp) | Today | After |
| --- | --- | --- |
| `R6_DEAnalyse$get_annotated_contrasts()` | SAINT branch picks `get_ora()`; LM branch picks `filter_significant_contrasts()` | `contrast_obj$filter_significant(FDR_threshold, diff_threshold)` |
| `R6_DEAnalyse$filter_contrasts()` | Same pair | Same single call |
| `R6_DEAnalyse$private$filter_significant_contrasts()` | Column-aware filter | Deleted; the contrast object owns this |
| `R6_DEAReportGenerator$contrasts_to_Grob()` | Column branch on `Bait`/`log2_EFCs`/`SaintScore`/`BFDR` vs `contrast`/`diff`/`statistic`/`FDR` | `contrast_obj$contrast_summary_table(rounded = TRUE)` per protein |
| `R6_DEAReportGenerator$make_SummarizedExperiment()` | `contrast_column <- if (saint) "Bait" else "contrast"` | `contrast_column <- contrast_obj$contrast_column_name()` |
| `R6_DEAReportGenerator$write_DEA()` | `saint` flag → branch in `.write_ORA` / `.write_GSEA` | `.write_ORA(contrast_obj, ...)` calls `contrast_obj$get_ora(...)`; `.write_GSEA(contrast_obj, ...)` calls `contrast_obj$get_rank(...)`. The `saint` parameter is removed from both helpers |
| `R6_DEAReportGenerator$write_DEA_all()` | Skip QC HTML for SAINT | `if (!contrast_obj$supports_dea_qc()) qc_file <- NULL` |
| `R6_DEAReportGenerator$prep_result_list()` | Reach into `dea$saint_input` / `dea$saint_result` | `contrast_obj$extra_artifacts()` returns a (possibly empty) named list that is merged into the result list; SAINT facade fills it with `saint_inter`/`saint_prey`/`saint_bait`/`saint_list` |
| `R6_ProteinDataPrep` SAINT branch | Decides annotation-reader shape | Driven from a facade-registry attribute (`entry$needs_saint_annotation`), not from a model-name match |
| `cmd_helpers` SAINT branch | `read_annotation(..., SAINT = saint_model)` | Same registry attribute |

`extra_artifacts()` is the escape valve for backend-specific data that
genuinely needs to ride along (SAINT input tables in the XLSX, the
prey/bait/inter sheets). It is opt-in: default returns `list()`.

### Implications for the schema in § 1

The schema-translation work in § 1 (`diff` ← `log2_EFCs`,
`FDR` ← `BFDR`, etc.) becomes unnecessary. `ContrastsSAINTFacade`'s
`get_contrasts()` can return SAINT-native columns unchanged. The
prolfqua standard schema only matters when downstream code wants to
treat all backends uniformly — and after this refactor, downstream
code asks the facade instead of reading columns.

The two columns that do need to be unified across backends are the
ones the `ContrastsPlotter` already standardises (`diff`, `contrast`,
`modelName` plus a score whose name is passed in). The SAINT facade's
existing `get_Plotter()` already constructs a plotter with
`diff = "log2_EFCs"`, `contrast = "Bait"`, `score = "SaintScore"`,
so that path is fine.

### Revised execution order

1. Add the new interface methods (`filter_significant`,
   `get_ora` default, `get_rank` default,
   `contrast_summary_table`, `supports_dea_qc`,
   `contrast_column_name`, `extra_artifacts`) to
   `prolfqua::ContrastsInterface` with default implementations that
   use the existing standard schema. Tests in `prolfqua`.
2. Override the relevant methods on `ContrastsSAINTexpress` /
   `ContrastsSAINTFacade` so SAINT speaks the interface in its own
   column names. Tests in `prolfquasaint`.
3. Land `register_facade()` in `prolfqua` and register the SAINT
   facade from `prolfquasaint/R/zzz.R`.
4. Replace `.is_saint_model()` branches one site at a time, each
   behind a focused test. Goal: zero `.is_saint_model()` calls in
   `prolfquapp/R/R6_*.R` and `cmd_helpers.R` after the refactor.
5. Delete `prolfquapp/R/facade_model.R` and
   `prolfquapp/R/saint_model.R`.
6. Move the YAML `"prolfqua"` → `"lm"`/`"lm_impute"` alias to a
   single resolver in `R6_DEAnalyse$initialize` (or to an aliases
   slot on the registry).
7. Docs (`prolfquapp/CLAUDE.md`, `prolfqua-adding-models` skill).
8. `make installs && make tests && make check-fast` from the
   workspace root.

---

## Revision 3: `ContrastConfiguration` (supersedes the interface-method expansion above)

User feedback on revision 2:

> Maybe it would be good also to have a similar configuration object
> for the contrast classes as we have for LFQData?

This is the right shape. `AnalysisConfiguration` works because it
holds the column-role mapping while `LFQData` exposes operations on
those roles; consumers call `lfqdata$response()`, never
`"transformedIntensity"`. Mirror that split for contrasts: a
`ContrastConfiguration` owns the column-role mapping, the contrast
object owns the operations that genuinely differ between backends.

### Pure-config vs hybrid

A pure-config approach (every operation derived from the config plus
generic code) does not quite work:

- `filter_significant` and `get_ora` are generic if the config carries
  one extra flag (`effect_symmetric`). LM uses `|diff| > t`; SAINT
  uses `log2_EFCs > t` or `< -t` depending on direction. Same code,
  branching on `effect_symmetric`.
- `get_rank` is *not* generic. LM's GSEA score is
  `sign(diff) * -log10(p.value)`; SAINT's is `log2_EFCs`. Different
  computation, not different column names.
- `extra_artifacts` is *not* generic. SAINT carries
  `saint_inter`/`saint_prey`/`saint_bait`/`saint_list`; LM carries
  nothing.

So: **config for column roles + behaviour flags, methods only for
operations whose logic actually differs.**

### `ContrastConfiguration` schema

R6 class `prolfqua::ContrastConfiguration` (flat, no nesting). Fields:

| Field | Type | LM default | SAINT |
| --- | --- | --- | --- |
| `subject_id` | `character()` | hierarchy keys | hierarchy keys |
| `model_name_col` | `character(1)` | `"modelName"` | `"modelName"` |
| `contrast_col` | `character(1)` | `"contrast"` | `"Bait"` |
| `effect_col` | `character(1)` | `"diff"` | `"log2_EFCs"` |
| `effect_symmetric` | `logical(1)` | `TRUE` | `FALSE` |
| `score_col` | `character(1)` | `"statistic"` | `"SaintScore"` |
| `score_direction` | `"two_sided" / "high" / "low"` | `"two_sided"` | `"high"` |
| `pvalue_col` | `character(1)` or `NA` | `"p.value"` | `NA_character_` |
| `fdr_col` | `character(1)` | `"FDR"` | `"BFDR"` |
| `avg_abundance_col` | `character(1)` | `"avgAbd"` | `"avgAbd"` |
| `supports_dea_qc` | `logical(1)` | `TRUE` | `FALSE` |
| `needs_saint_annotation` | `logical(1)` | `FALSE` | `TRUE` |

`score_direction` answers "is bigger better, smaller better, or both
extremes interesting" for the rank/plotter helpers. `effect_symmetric`
answers "is `up` decided by sign of effect alone (asymmetric backend
like SAINT) or by sign and magnitude both (symmetric backend like LM
where ORA up/down filter is `|diff| > t & diff > 0`)."

Accessors mirror `AnalysisConfiguration`: `cfg$effect_col()`,
`cfg$fdr_col()`, `cfg$contrast_col()`, ... `cfg$has_pvalue()` returns
`!is.na(cfg$pvalue_col())`.

### Updated `ContrastsInterface`

Required fields:
- `contrast_result` — backend-native tidy data.frame
- `config` — `ContrastConfiguration`

Required methods (defaults provided by `ContrastsInterface`):

| Method | Default implementation | When to override |
| --- | --- | --- |
| `get_contrasts(all = FALSE)` | return `contrast_result` | rarely (LM, SAINT keep theirs) |
| `get_config()` | return `self$config` | never |
| `filter_significant(FDR_threshold, diff_threshold)` | generic, uses `cfg$fdr_col()`, `cfg$effect_col()`, `cfg$effect_symmetric` | never |
| `get_ora(up = TRUE, FDR_threshold, diff_threshold)` | generic, uses the same fields | never |
| `get_rank(score = NULL)` | LM default: `sign(diff) * -log10(p.value)` if `cfg$has_pvalue()`; else `cfg$effect_col()` | SAINT, ROPECA |
| `get_Plotter()` | construct `ContrastsPlotter` from config | rare |
| `to_wide(columns = NULL)` | call `pivot_model_contrasts_to_wide` with config-derived defaults | rare |
| `get_missing()` | generic `.compute_missing()` | rare |
| `contrast_summary_table(rounded = TRUE)` | generic select of `contrast`/`effect`/`score`/`fdr` columns, renamed to canonical names, optionally rounded | rare |
| `extra_artifacts()` | `list()` | SAINT |

The interface no longer leaks column names. Consumers either ask the
config for a column role (when they need to *touch* the table, e.g.
inside `to_wide()`) or call a method (when they need an operation).
`DEAnalyse` and `DEAReportGenerator` only call methods.

### Mapping every `.is_saint_model()` site to the new API

Same 15 sites as revision 2, but each maps to either a config
accessor or a method, with zero column-name knowledge in
`prolfquapp`:

| Site (prolfquapp) | Today | After |
| --- | --- | --- |
| `R6_DEAnalyse$get_annotated_contrasts()` SAINT branch | `get_ora(up = TRUE, ...)` | `contrast_obj$filter_significant(FDR_threshold, diff_threshold)` for the "significant" set; `contrast_obj$get_ora(up = TRUE, ...)` only when ORA is needed |
| `R6_DEAnalyse$get_annotated_contrasts()` LM branch | `private$filter_significant_contrasts()` | same `filter_significant()` call |
| `R6_DEAnalyse$filter_contrasts()` | same pair | `contrast_obj$filter_significant(...)` |
| `R6_DEAnalyse$private$filter_significant_contrasts()` | column-aware filter | deleted |
| `R6_DEAReportGenerator$contrasts_to_Grob()` | column branch | `contrast_obj$contrast_summary_table(rounded = TRUE)` per protein |
| `R6_DEAReportGenerator$make_SummarizedExperiment()` | `if (saint) "Bait" else "contrast"` | `contrast_obj$get_config()$contrast_col()` |
| `R6_DEAReportGenerator$write_DEA()` | `saint` flag in `.write_ORA` / `.write_GSEA` | helpers take a contrast object and call `get_ora()` / `get_rank()`; `saint` parameter removed from both |
| `R6_DEAReportGenerator$write_DEA_all()` | skip QC HTML for SAINT | `if (!contrast_obj$get_config()$supports_dea_qc()) qc_file <- NULL` |
| `R6_DEAReportGenerator$prep_result_list()` | reach into `dea$saint_input` / `dea$saint_result` | `resultList <- c(resultList, contrast_obj$extra_artifacts())` |
| `R6_ProteinDataPrep` SAINT branch | annotation-reader shape | `entry$needs_saint_annotation` on the facade registry |
| `cmd_helpers` SAINT branch | `read_annotation(..., SAINT = saint_model)` | same registry attribute |

After this refactor: zero `.is_saint_model()` calls in
`prolfquapp/R/`. The only place SAINT is named is the registry entry.

### Backwards compatibility

- `dea$saint_input` / `dea$saint_result` keep working for one
  release: `R6_DEAnalyse$build_facade()` continues to populate them
  from `facade$extra_artifacts()`, deprecated in roxygen.
- `ContrastsSAINTexpress` keeps its existing `get_contrasts()` /
  `get_ora()` / `get_rank()` shape (already there). It gains a
  `config` field initialised to a SAINT-flavoured
  `ContrastConfiguration` and methods that delegate to the config
  (`filter_significant`, `contrast_summary_table`, `to_wide`,
  `get_missing`, `extra_artifacts`).
- LM/limma/firth/etc. facades gain a `config` field initialised to
  an LM-flavoured `ContrastConfiguration`. Their existing
  `get_contrasts()` results already match those column names, so the
  config is descriptive, not transformative.

---

## Detailed implementation plan

Each step lists files touched, the change, and the test that pins it.
Steps within a phase can land in separate PRs; phases must land in
order.

### Phase A — Foundation in `prolfqua`

#### A1. `ContrastConfiguration` class

- File: `prolfqua/R/ContrastConfiguration.R` (new).
- Implement R6 class with the fields above, plus accessor methods for
  each field (e.g. `effect_col()`, `fdr_col()`, ...). Provide:
  - `ContrastConfiguration$new(subject_id, ...)`: arguments default
    to the LM-flavoured values listed in the table.
  - `has_pvalue()`: convenience wrapper.
  - `clone(deep = TRUE)` semantics that match
    `AnalysisConfiguration`.
- Roxygen: document every field, link to `ContrastsInterface`.
- Test: `tests/testthat/test-ContrastConfiguration.R`. Cover
  construction, defaults, deep clone, accessor return values.
- `make document && make check-fast`.

#### A2. Generic helpers on `ContrastsInterface`

- File: `prolfqua/R/ContrastsInterface.R`.
- Add a `config` field (default `NULL`) and a `get_config()` method
  returning `self$config` (or a no-arg constructed default if
  `NULL`, for backward compatibility with adapters that have not
  been migrated yet).
- Add **default implementations** (so subclasses get them for free):
  - `filter_significant(FDR_threshold, diff_threshold)`: dispatch on
    `cfg$effect_symmetric`.
  - `get_ora(up, FDR_threshold, diff_threshold)`: same predicate
    plus the `up` direction.
  - `get_rank(score = NULL)`: if `cfg$has_pvalue()` use
    `sign(effect) * -log10(p.value)`; else use `effect_col`.
  - `contrast_summary_table(rounded = TRUE)`: select the four
    config-driven columns and rename to canonical
    `contrast`/`effect`/`score`/`fdr`.
  - `extra_artifacts()`: returns `list()`.
- Add unit tests using a tiny fake contrast object (data.frame plus
  config) to cover both symmetric and asymmetric branches of the
  defaults.

#### A3. `register_facade()` / `lookup_facade()`

- File: `prolfqua/R/ContrastsFacades.R` (extend).
- Convert `FACADE_REGISTRY` to an internal environment plus exported
  `register_facade(name, class, needs, package = "prolfqua",
  needs_saint_annotation = FALSE)` and `lookup_facade(name)`.
- Keep `FACADE_REGISTRY` exported for inspection (read-only view) so
  downstream code that lists methods keeps working.
- Rewrite `build_contrast_analysis()` to call `lookup_facade()`.
- Test: `tests/testthat/test-register-facade.R`. Register a fake
  facade, look it up, build via `build_contrast_analysis()`,
  unregister.

#### A4. Wire existing facades to `ContrastConfiguration`

- Files: `prolfqua/R/Contrasts.R`, `ContrastsLimma.R`,
  `ContrastsModerated.R`, `ContrastsLimpa.R`, `ContrastsROPECA.R`,
  `ContrastFirth.R`, `ContrastsSimpleImpute.R`,
  `ContrastsModeratedDEqMS.R`, `ContrastsFacades.R`.
- In each facade `initialize()`, build an LM-flavoured config
  (override `effect_symmetric = TRUE`, `score_direction =
  "two_sided"`, `pvalue_col = "p.value"`, `fdr_col = "FDR"`,
  `effect_col = "diff"`, `contrast_col = "contrast"`,
  `score_col = "statistic"`, `supports_dea_qc = TRUE`,
  `needs_saint_annotation = FALSE`).
- `ContrastsROPECAFacade` overrides `score_col = "statistic"` and
  keeps its own `get_contrasts()` since it already normalises
  columns; just attach the config.
- Verify existing tests in
  `prolfqua/tests/testthat/test-ContrastsFacades.R` still pass.

### Phase B — SAINT in `prolfquasaint`

#### B1. SAINT-flavoured config and method overrides

- File: `prolfquasaint/R/ContrastSaintExpress.R`.
- Give `ContrastsSAINTexpress` a `config` field initialised in
  `initialize()` with: `effect_col = "log2_EFCs"`,
  `effect_symmetric = FALSE`, `score_col = "SaintScore"`,
  `score_direction = "high"`, `pvalue_col = NA_character_`,
  `fdr_col = "BFDR"`, `avg_abundance_col = "avgAbd"`,
  `contrast_col = "Bait"`, `supports_dea_qc = FALSE`,
  `needs_saint_annotation = TRUE`, `subject_id = self$subject_id`.
- Keep existing `get_contrasts()`, `get_ora()`, `get_rank()` —
  signatures already match. Update them to read thresholds from
  arguments (already so) and column names from `self$config`
  (refactor; behaviour unchanged).
- Override `get_rank()` because the default would compute
  `-log10(p.value)` and SAINT has no p-value.
- Override `extra_artifacts()` to return the SAINT input/result
  tables once they are stored on the facade (see B2).
- Override `to_wide()` to keep current SAINT defaults.

#### B2. `ContrastsSAINTFacade` and registration

- File: `prolfquasaint/R/ContrastsSAINTFacade.R` (new).
- Mirror `ContrastsLimmaFacade`. Fields: `model`, `contrast`,
  `.lfqdata`, `.contrast_names`, `saint_input`, `saint_result`.
- `initialize(lfqdata, modelstr = NULL, contrasts = NULL,
  row_annot = NULL, spc = FALSE, engine = "r")`:
  - assert aggregated input.
  - move `.saint_prepare_data` / `.saint_bait_col` /
    `.saint_control_col` / `.saint_first_existing` from
    `prolfquapp/R/saint_model.R` into a new
    `prolfquasaint/R/saint_prepare.R` (private).
  - run `protein_2localSaint()` → `runSaint()` → wrap in
    `ContrastsSAINTexpress`; remap `Prey` to the project hierarchy
    key.
  - store `saint_input` and `saint_result` on the facade for
    `extra_artifacts()`.
- Delegate every interface method to `self$contrast` (the wrapped
  `ContrastsSAINTexpress`). Override `extra_artifacts()` to surface
  the SAINT-specific tables.
- File: `prolfquasaint/R/zzz.R` (new or extend).
  - `.onLoad` calls
    `prolfqua::register_facade("saint", "ContrastsSAINTFacade",
    needs = "aggregated", package = "prolfquasaint",
    needs_saint_annotation = TRUE)`.
- Tests: `prolfquasaint/tests/testthat/test-ContrastsSAINTFacade.R`.
  - construction from a small simulated `LFQData` with
    bait/control annotation columns.
  - `get_contrasts()` returns the SAINT-native columns.
  - `filter_significant()`, `get_ora()`, `get_rank()`,
    `contrast_summary_table()`, `extra_artifacts()`, `to_wide()`,
    `get_Plotter()`, `get_missing()` all work without
    `prolfquapp`.
  - registration: `prolfqua::lookup_facade("saint")` resolves the
    class when `prolfquasaint` is loaded.

#### B3. DESCRIPTION updates

- `prolfquasaint/DESCRIPTION`: confirm `Imports: prolfqua`. (Likely
  already there.)
- `prolfquapp/DESCRIPTION`: already imports `prolfquasaint` from the
  earlier commit. No change.

### Phase C — Strip `.is_saint_model()` from `prolfquapp`

Each step lands behind a focused test. Goal: zero
`.is_saint_model()` calls in `prolfquapp/R/`.

#### C1. Route SAINT through the registry

- File: `prolfquapp/R/R6_DEAnalyse.R::build_facade()`.
- Replace the `if (.is_saint_model(name))` branch with a single
  registry lookup: `entry <- prolfqua::lookup_facade(name)`. Pull
  `entry$class` from `entry$package` (default `"prolfqua"`). Call
  `facade_class$new(self$lfq_data, modelstr, self$contrasts,
  row_annot = self$rowAnnot)` and ignore extra args on facades
  that do not declare `row_annot`.
- After construction, if `facade$saint_input` is non-`NULL`, mirror
  onto `self$saint_input` / `self$saint_result` for one-release
  back-compat; roxygen-deprecate those fields.
- Drop the `if (!.is_saint_model(default_model)) {
  stopifnot(length(contrasts) >= 1) }` guard. The SAINT facade
  ignores `contrasts`.
- Test: facade dispatch works for `"lm"`, `"lm_missing"`, and
  `"saint"` against a tiny fixture, with no `.is_saint_model()`
  reference in the test body.

#### C2. Significance filtering

- File: `prolfquapp/R/R6_DEAnalyse.R`.
- Delete `private$filter_significant_contrasts()`.
- `get_annotated_contrasts()` and `filter_contrasts()` call
  `contrast_obj$filter_significant(self$FDR_threshold,
  self$diff_threshold)`. The SAINT-only `get_ora()` branch in
  `get_annotated_contrasts()` becomes
  `contrast_obj$get_ora(up = TRUE, ...)` *only when ORA is what we
  need*; otherwise use `filter_significant()`.
- Test: significance filtering produces the same rows as today for
  both an LM fixture and a SAINT fixture.

#### C3. Per-protein summary grobs

- File: `prolfquapp/R/R6_DEAReportGenerator.R::contrasts_to_Grob()`.
- Replace the column branch with
  `contrast_obj$contrast_summary_table(rounded = TRUE)` per
  protein, then `gridExtra::tableGrob()`. The summary table is
  already rounded and uses canonical column names.
- Test: grob has the expected canonical column names for both
  backends.

#### C4. `make_SummarizedExperiment()`

- File: `prolfquapp/R/R6_DEAReportGenerator.R`.
- `contrast_column <- contrast_obj$get_config()$contrast_col()`.
- Test: SE rowData splits by `"contrast"` for LM and `"Bait"` for
  SAINT.

#### C5. ORA / GSEA write helpers

- File: `prolfquapp/R/R6_DEAReportGenerator.R`,
  `prolfquapp/R/report_helpers.R` (or wherever `.write_ORA` /
  `.write_GSEA` live).
- Helpers take a contrast object and call `contrast_obj$get_ora(...)`
  / `contrast_obj$get_rank(...)`. Remove `saint` parameter.
- Test: ORA files and rank files are written for both backends with
  the same file naming.

#### C6. QC HTML capability

- File: `prolfquapp/R/R6_DEAReportGenerator.R::write_DEA_all()`.
- `if (contrast_obj$get_config()$supports_dea_qc()) qc_file <-
  ... else qc_file <- NULL`.
- Test: QC HTML is rendered for LM and skipped for SAINT.

#### C7. Extra artifacts in result list

- File: `prolfquapp/R/R6_DEAReportGenerator.R::prep_result_list()`.
- After building `resultList`, merge in
  `contrast_obj$extra_artifacts()`. Remove the
  `saint_inter`/`saint_prey`/`saint_bait`/`saint_list` block.
- Test: SAINT extra tables appear in the result list; LM result
  list has no spurious SAINT keys.

#### C8. Annotation reading

- Files: `prolfquapp/R/R6_ProteinDataPrep.R`,
  `prolfquapp/R/cmd_helpers.R`.
- Replace `.is_saint_model(model)` with a registry-attribute lookup:
  `entry <- prolfqua::lookup_facade(model);
  saint_annot <- isTRUE(entry$needs_saint_annotation)`.
- Test: `read_annotation()` is called with `SAINT = TRUE` for the
  SAINT model and `FALSE` otherwise.

### Phase D — Cleanup

#### D1. Delete `prolfquapp/R/facade_model.R`

- Move the legacy YAML alias `"prolfqua" → "lm" / "lm_impute"` into
  a one-line resolver in `R6_DEAnalyse$initialize` (or to an
  `aliases` table on the registry; prefer the inline resolver to
  avoid scope creep).
- Confirm `.valid_facade_models()` and `.is_saint_model()` are
  unreferenced (`git grep` across `prolfquapp/R/` and
  `prolfquapp/inst/`).
- Delete the file.

#### D2. Delete `prolfquapp/R/saint_model.R`

- Confirm `.build_saint_contrast_result()` is unreferenced.
- Delete.

#### D3. Drop transitional `dea$saint_input` / `dea$saint_result`

(Only after one release with the deprecation warning.)

- Roxygen-remove the deprecated fields.
- Confirm reports read from `contrast_obj$extra_artifacts()`.

### Phase E — Docs

- `prolfqua/man/*` regenerated by `make document` after A1–A4.
- `prolfqua/vignettes/` (if any walks through facades): add a note
  on `ContrastConfiguration` and `register_facade()`.
- `prolfqua-adding-models` skill at
  `~/.claude/skills/prolfqua-adding-models/SKILL.md`: replace
  references to `FACADE_REGISTRY` with `register_facade()` and add
  a section on building a `ContrastConfiguration`.
- `prolfquapp/CLAUDE.md`: note that SAINT is reached through the
  facade registry like `lm` and `limma`.

### Phase F — Validation

- Per-package: `make document && make check-fast` in each of
  `prolfqua`, `prolfquasaint`, `prolfquapp`.
- Ecosystem: from `prolfqua_fml/` root, `make installs && make
  tests && make check-fast`.
- Integration: re-render the SAINT WU345302 example and verify the
  report has correct contrasts, ORA, GSEA, and no empty `NA`
  facets.

### Acceptance criteria

- `grep -r ".is_saint_model\|\\bsaint_model\\b" prolfquapp/R/` is
  empty.
- `grep -r "saint" prolfquapp/R/` matches only YAML keyword
  documentation and the `read_annotation` integration.
- `prolfqua::lookup_facade("saint")$class` resolves to
  `"ContrastsSAINTFacade"` when `prolfquasaint` is loaded.
- All existing `prolfquapp` SAINT tests pass without modification
  to their expected outputs.
- Adding a new backend now requires (a) one
  `ContrastConfiguration`, (b) a `Contrasts*` class with at most
  `get_rank` / `extra_artifacts` overrides, (c) a facade, and (d) a
  `register_facade()` call — no edits to `prolfquapp`.

## Open questions

- Should `prolfqua` learn about a `register_facade()` API, or
  should the registry stay closed and `prolfquapp` keep a local
  override? The former is cleaner and matches how
  `build_contrast_analysis()` is meant to be the single entry
  point; the latter avoids touching `prolfqua` at all in this
  slice.
- Should `p.value` for SAINT be `1 - SaintScore` to keep downstream
  plotting / histogram code happy, or stay `NA`? Bias is toward
  `NA` plus making the histogram code tolerate missing `p.value`.
  Confirm before implementation.
- Does the SAINT facade need to support `nested` LFQData (peptide
  level) as well, or only aggregated input as today? Current code
  only handles aggregated; keep that constraint unless a use case
  appears.
