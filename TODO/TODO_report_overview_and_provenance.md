# Report overview and provenance

> Give every Quarto report a short visual abstract at the start and a consistent, compact provenance record at the end.

## Requirements

- Each report opens with a one-screen **Overview**. It replaces a pre-existing Introduction tab; reports without an
  Introduction gain it as their first tab or first section.
- The Overview uses a static visual abstract authored for that report type, a short explanation, and a concise input
  summary. The input summary begins with three compact cards for samples, groups, and quantified proteins (or the
  report's analysed feature type). It contains no run-specific result figures and does not embed data tables.
- The final **Session Info** area contains two subtabs: **Report provenance** and **R session info**.
- Report provenance records the workunit, order and project context, creator, creation time, input-data reference,
  software/model where available, and the package version. It remains a short field/value table, not embedded input
  data.

## Design

- Keep report content explicit in each QMD. The FGCZ Quarto template supplies styling and authoring guidance only; it
  does not infer a report type, choose figures, or inject report content.
- Store one AI-generated PNG visual abstract per report type in `vignettes/visual_abstracts/`. Package and stage that
  directory beside each runtime-rendered QMD so the report remains self-contained.
- Use internal provenance helpers to keep field names, creator/timestamp capture, and missing-value behaviour
  consistent across all five reports.

## Implementation plan

- [x] Add visual abstracts for primary DEA, tabbed DEA, DEA QC, protein abundances, and QC/sample-size estimation.
- [x] Package and stage the visual-abstract directory for runtime Quarto rendering.
- [x] Add the Overview and final two-subtab Session Info layout to all five QMD reports.
- [x] Add the compact samples/groups/proteins overview cards to all five QMD reports.
- [x] Update the FGCZ Quarto report-authoring skill with the agreed layout rules.
- [x] Render all five QMD reports and verify the Overview image and both Session Info subtabs.

## Open questions

- Review the first generated visual abstracts in the rendered reports and replace individual assets if their scientific
  emphasis or visual style needs refinement.
