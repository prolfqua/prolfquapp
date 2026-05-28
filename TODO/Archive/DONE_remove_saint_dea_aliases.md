# Remove SAINT DEA Aliases

## Summary

Remove the `.saint_add_dea_aliases()` compatibility layer so SAINT results keep
SAINT-native columns (`Bait`, `log2_EFCs`, `BFDR`, `SaintScore`) through
analysis and reporting.

## Implementation Plan

- In `prolfqua`, move `ContrastsPlotter` contrast data storage from public
  `contrast_df` to private state and keep all existing plot methods working.
- In `prolfquasaint`, make `ContrastsSAINTexpress$get_Plotter()` use
  `self$get_contrasts()` rather than direct internal data access.
- In `prolfquapp`, delete `.saint_add_dea_aliases()` and remove all call sites.
- For SAINT significance filtering, use
  `contrast_obj$get_ora(up = TRUE, FDR_threshold, diff_threshold)` and join row
  annotations by the subject identifier.
- Keep regular DEA filtering on `FDR` and `diff`.
- Update report helpers/templates so SAINT logic does not require synthetic
  `contrast`, `diff`, `FDR`, or `statistic` columns.

## Validation

- Run targeted tests for `ContrastsPlotter`, `ContrastSaintExpress`, and the
  `prolfquapp` SAINT/report paths.
- Reinstall local packages with package Makefiles after tests pass.
- Rerun the WU345302 SAINT example and verify HTML, XLSX, ORA, and RNK output.
