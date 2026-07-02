# TODO: expose `nr_peptides` (minimum peptides per protein) end-to-end

> **Status: IMPLEMENTED 2026-07-02** in prolfquapp (reader-local design below).
> Full test suite green (388 pass); an adversarial review pass found one real
> bug (non-finite threshold slipped validation → helper crash), now fixed +
> regression-tested. The prolfqua cleanup and PTM stripped-peptide element
> remain as their own deferred TODOs.
>
> Revised 2026-07-02 after `REVIEW_expose_nr_peptides.md` (second review) and a
> design interview.
>
> **Design decision (supersedes the earlier "generic `LFQData`-side wrapper"
> sketch): the filter is READER-LOCAL.** Each reader applies the count + drop
> against *its own* known stripped-peptide column, to *its own* quant table and
> protein annotation — DIA-NN already does this correctly and is the template.
> The earlier idea of a single generic step reading `xd$lfqdata` was **wrong**:
> only the reader knows which column is the stripped (unmodified) peptide, and
> `peptide_Id` is not that column in every reader (see the mapping table below).

## Requirements

- A DEA run drops any parent protein — and its child peptides — whose parent
  protein is supported by fewer than **N distinct stripped (unmodified)
  peptides**. Default `N = 1` (no filtering; existing runs unaffected).
- The cut is applied **by the reader**, to both its quant table (→ `LFQData`)
  and the protein annotation it builds, so the returned `xd` is consistent by
  construction (`xd$lfqdata`, `xd$protein_annotation`, IBAQ, and the logged
  summary all agree). `ProteinAnnotation` follows the filtered quant; it never
  drives the drop.
- The count basis is **each reader's own stripped-peptide column** — the reader
  supplies it to a shared helper, because only the reader knows it.
- MZMine (metabolomics) is **exempt** — one feature per protein, so a `≥2` cut
  would wipe every feature; the reader accepts the option and ignores it.
- The value is carried on the **A414/A413 native prolfquapp YAML path**
  (`processing_options: {nr_peptides: N}`) and — full scope — also via a CLI flag
  and the YAML generator.

## Hard constraint: `ProteinAnnotation` is never the filter

`ProteinAnnotation` is always **right-joined** into the `LFQData`, so annotation
never drops quant rows (`AGENTS.md` invariant → "ProteinAnnotation is annotation,
not a filter"). The reader-local design honors this naturally: the reader
computes the kept-protein set from the quant/peptide table, filters the quant,
and **builds the annotation from the already-filtered set** (DIA-NN builds
`prot_annot` from the filtered `nrPEP`, `R/preprocess_DIANN.R:294`). Do **not**
use `ProteinAnnotation$filter_by_nr_children()` or any `ProteinAnnotation` method
to drive dropping.

## Current implementation, verified in code (2026-07-02)

The `nr_peptides` field exists, and DIA-NN fully implements the filter, but the
value never reaches any reader on the DEA path, so it is inert everywhere.

- **DIA-NN already implements it correctly and consistently.** `get_nr_pep()`
  counts distinct `Stripped.Sequence` per `Protein.Group` from the raw report
  (`R/preprocess_DIANN.R:1-9`); when `nr_peptides > 1` it filters *both* the
  `peptide` quant table and `nrPEP` (`R/preprocess_DIANN.R:268-271`); and
  `prot_annot` is then built from the filtered `nrPEP`
  (`R/preprocess_DIANN.R:294`). So `lfqdata` and `protein_annotation` come out in
  sync — the review's finding #1 (state drift) does not arise for this shape.
- **The value never reaches the reader.** `preprocess_software()` forwards only
  `pattern_contaminants` / `pattern_decoys` into the reader `do.call`
  (`R/preprocess_software.R:168-180`); `nr_peptides` is not forwarded, so DIA-NN's
  filter runs at its own default of 1 (dead code on the DEA path).
- **8 of 9 readers don't accept the parameter.** Only `preprocess_DIANN` has
  `nr_peptides` in its signature (`R/preprocess_DIANN.R:217`). MaxQuant, FP_PSM,
  MSstats, MSstats_FPDIA, BGS, SIM, dummy, and mzMine do not — so *uniformly*
  adding `nr_peptides` to the `do.call` arg list would raise "unused argument" in
  every one of them (the do.call hazard the first review flagged is real).
- **Only the reader knows the stripped-peptide column** — `peptide_Id` is not
  universally the stripped sequence:

  | Reader | `protein_Id` | `peptide_Id` | stripped seq? |
  |---|---|---|---|
  | DIANN | `Protein.Group` | `Stripped.Sequence` | yes |
  | MaxQuant | `leading.razor.protein` | `sequence` | yes (MQ `sequence` is unmodified) |
  | FP_PSM | `Protein` | `Peptide` | likely yes (`Peptide` vs `Modified Peptide`) — verify |
  | MSstats / FP_DIA | `ProteinName` | `PeptideSequence` | **unknown — may carry modifications** |
  | BGS | `PG.ProteinGroups` | `PEP.GroupingKey` | **unknown — may carry modifications** |
  | SIM | `protein_Id` | `peptide_Id` | yes (simulated) |
  | MZMine / dummy | — | (none) | n/a |

- **CLI / YAML generator omit it:** `CMD_DEA_V2.R` (`inst/application/CMD_DEA_V2.R`)
  has no flag and `sync_opt_config()` (`R/utils.R`) doesn't merge one;
  `CMD_MAKE_YAML.R` always writes `nr_peptides: 1.0`.
- **Config plumbing is fine:** `ProcessingOptions$nr_peptides = 1`
  (`R/R6_AppConfiguration.R:24-25`) is populated by `list_to_R6_app_config()`
  (`R/R6_AppConfiguration.R:321-348`) — the value simply is not consumed.

## Design (reader-local filter + one shared helper)

1. **Shared helper (prolfquapp), parameterized by the reader's columns.**
   e.g. `filter_by_peptide_count(data, protein_col, peptide_col, nr_peptides)`:
   count distinct `peptide_col` per `protein_col`, keep proteins with
   `>= nr_peptides`, return the filtered table (and/or the kept-protein set). It
   is a no-op for `nr_peptides <= 1`. The **reader passes its own** `protein_col`
   and stripped `peptide_col`, because only the reader knows them.
2. **Each supported reader calls the helper** with its columns and applies the
   result to both the quant table and the annotation it builds — exactly as
   DIA-NN does today. Refactor DIA-NN's inline block onto the shared helper so
   there is one implementation. `xd` stays consistent by construction.
3. **Forwarding via formals inspection (avoids the do.call hazard).**
   `preprocess_software()` forwards `nr_peptides` **only to readers whose
   `formals()` declare it**, so migrated readers pick it up and non-migrated /
   genuinely-N/A readers are untouched. If `nr_peptides > 1` and the reader does
   **not** declare it, **log a warning** ("reader X does not support nr_peptides
   filtering; ignoring") so the skip is explicit, never a silent no-op (review
   findings #3, #4). `run_dea()` passes `config$processing_options$nr_peptides`
   into `preprocess_software()`.
4. **MZMine is an explicit exemption:** `preprocess_mzMine` gains an
   `nr_peptides` param that it *accepts and ignores* (documented), so it is
   exempt without a warning.
5. **Validate the threshold centrally** (review #5): require a scalar whole
   number `>= 1`; coerce clean YAML numerics (`2.0` → `2`); fail clearly on `0`,
   negatives, non-integers, `NA`, or non-numeric. Absent → default 1.

Why reader-local (not the generic wrapper): only the reader knows the
stripped-peptide column; the wrapper reading `peptide_Id` would silently count
the wrong thing for MSstats/BGS; and building annotation from the filtered set in
the reader keeps `xd` consistent with no resync step.

### Adapter behavior

| Adapter family | stripped-peptide column | behavior |
|---|---|---|
| DIANN | `Stripped.Sequence` | filter (already implemented; refactor onto shared helper) |
| MAXQUANT | `sequence` | add param + call helper |
| FP_TMT (FP_PSM) | `Peptide` (verify vs modified) | add param + call helper |
| MSSTATS / MSSTATS_FP_DIA | `PeptideSequence` (verify stripped) | add param + call helper |
| BGS | `PEP.GroupingKey` (verify stripped) | add param + call helper |
| SIM | `peptide_Id` | add param + call helper — regression fixture |
| MZMINE / MZMINEannot | none | accept + **ignore** (explicit exemption) |
| DUMMY | none | accept + ignore (interface stub) |
| PTM / site (`prolfquappPTMreaders`) | not exposed yet (child = site) | **deferred** — accept the param now (so the forwarded do.call doesn't break), implement once a stripped-peptide element exists (handoff) |

## Implementation plan

- [ ] Add the shared `filter_by_peptide_count(data, protein_col, peptide_col, nr_peptides)`
      helper (location — see open questions) with a no-op at `nr_peptides <= 1`.
- [ ] Refactor `preprocess_DIANN` to call the shared helper (no behavior change).
- [ ] Add an `nr_peptides = 1` param to each remaining supported reader and call
      the helper with that reader's `(protein_col, peptide_col)`, applied to the
      quant table and the annotation build.
- [ ] Add `nr_peptides = 1` to `preprocess_mzMine` and `preprocess_dummy`,
      accepted and ignored.
- [ ] Forward `nr_peptides` in `preprocess_software()` via `formals()` inspection
      + a warning when unsupported and `nr_peptides > 1`; pass it from `run_dea()`
      (`config$processing_options$nr_peptides`).
- [ ] Central validation of `nr_peptides` (whole number `>= 1`).
- [ ] Full CLI scope: `--nr_peptides` (or `--min_nr_peptides`) in `CMD_DEA_V2.R`
      + `sync_opt_config()`; an `nr_peptides` option in `CMD_MAKE_YAML.R` /
      `run_make_yaml()`.
- Files: `R/preprocess_software.R`, `R/preprocess_DIANN.R`,
  `R/preprocess_MaxQuant.R`, `R/preprocess_FP_PSM.R`, `R/preprocess_MSstats.R`,
  `R/preprocess_BGS_default.R`, `R/preprocess_SIM.R`, `R/preprocess_MzMine.R`,
  `R/preprocess_dummy.R`, `R/cmd_helpers.R`, `R/utils.R`,
  `inst/application/CMD_DEA_V2.R`, `inst/application/CMD_MAKE_YAML.R`, plus the new
  helper file.

## Tests

- **Primary regression (proves the fix):** with `nr_peptides = 2`, a
  one-stripped-peptide protein is **absent** from both `xd$lfqdata` and
  `xd$protein_annotation$row_annot`; with `nr_peptides = 1` it **remains**.
  `--software SIM` (peptide-level) and a DIA-NN fixture are the first targets.
- **State consistency:** after a run at `nr_peptides = 2`, `result$xd$lfqdata`,
  `result$xd$protein_annotation$row_annot`, and `deanalyse$lfq_data_raw` agree on
  the kept proteins; IBAQ / log-summary use the filtered `xd`.
- **Shared-helper unit test:** distinct-count + drop correct at N = 2 / 3 on a
  small table, no-op at N = 1.
- **MZMine:** option accepted and ignored — no features dropped at N = 2.
- **Unsupported reader:** a reader without the param + `nr_peptides > 1` warns
  (not a silent no-op).
- **Validation:** `1`, `2`, YAML `2.0` pass; `0`, negative, non-integer,
  `NA`/missing (unless defaulting to 1), non-numeric fail clearly.
- **Full-scope:** `CMD_DEA_V2.R` `--nr_peptides` override via `sync_opt_config()`;
  `CMD_MAKE_YAML.R` writes the chosen value.
- A414 (**done**): `build_config()` writes `processing_options["nr_peptides"]`
  from the `min_nr_peptides` executable parameter
  (`slurmworker/config/A414_DEA/tests/test_process_prolfqua.py`).

## Cross-package

- **A414** (done): exposes `min_nr_peptides` (dropdown 1/2/3/4, default 1) and
  emits `processing_options: {nr_peptides: N}` in its native prolfquapp YAML;
  verified to reach `pop$nr_peptides`. Inert until this prolfquapp fix ships.
- **PTM readers** (deferred, bigger ask): site-level `LFQData` has no
  stripped-peptide child, so peptides-per-protein can't be counted there yet.
  **New sub-requirement from this design:** the PTM readers must *accept* the
  `nr_peptides` param (even if initially ignored) so that once
  `preprocess_software()` forwards it via formals inspection, nothing breaks and
  they can be migrated later. Handed off:
  `prolfquappPTMreaders/TODO/TODO_stripped_peptide_hierarchy_for_nr_peptides.md`.
- **prolfqua** cleanup (separate): `prolfqua::filter_proteins_by_peptide_count()`
  / `nr_B_in_A()` is depth-coupled and reads `config$min_peptides_protein`, so it
  is *not* reused here. Handed off:
  `prolfqua/TODO/TODO_filtering_by_childcount.md`.
- **Release:** A414 runs the published `prolfqua/prolfquapp:latest`, so the
  behavior reaches production only after a new prolfquapp release (version bump +
  tag → CI `publish_docker.yml` republishes `:latest`).

## How the review findings are resolved

- **#1 (state drift):** dissolved — the reader builds annotation from the
  filtered set (DIA-NN template), so `xd` is consistent by construction; no
  resync step in `run_dea()`.
- **#2 (existing prolfqua filter not reusable):** not reused; reader-local shared
  helper instead; prolfqua-side cleanup handed off to `TODO_filtering_by_childcount.md`.
- **#3 (unsupported downstream readers):** explicit — forward by `formals()`, warn
  when `nr_peptides > 1` and unsupported; no accidental site-counting.
- **#4 (explicit MZMine skip):** MZMine accepts + ignores the param (explicit
  exemption, no warning).
- **#5 (validation):** central whole-number `>= 1` validation added to the plan.

## Open questions

- **MSstats (`PeptideSequence`) / BGS (`PEP.GroupingKey`):** is the column the
  stripped/unmodified sequence, or does it carry modifications? If modified, the
  reader needs a stripping step (or a different source column) before the count
  is correct.
- **FP_PSM (`Peptide`):** confirm it is the unmodified peptide, not the modified one.
- **Shared helper name / location:** a `prolfquapp` reader utility (leaning here)
  vs later promotion to `prolfqua` (tracked in `TODO_filtering_by_childcount.md`).
- **CLI flag name:** `--nr_peptides` vs `--min_nr_peptides` (align with the A414
  parameter label `min_nr_peptides`).
- **Validation location:** in the shared helper's guard vs at config
  construction (`list_to_R6_app_config` / `sync_opt_config`).
