# Expose Normalization in CMD_DEA_CD

## Summary

Add `-n` / `--normalization` to `inst/application/CMD_DEA_CD.R` so CompoundDiscoverer DEA runs can select the
normalization method without requiring a YAML file.

## Scope

- Expose the CLI values requested for CD:
  - `vsn`
  - `robscale`
  - `none`
- Keep the default as `vsn`.
- If a YAML file is supplied, allow `--normalization` to override the YAML config, matching the existing `--model`
  override behavior.
- Validate invalid normalization values early with a clear error.
- Add or update subprocess coverage for `CMD_DEA_CD.R`.

## Downstream

`pd_metaboanalyst` will pass this option through `prolfqua_dea_cd.sh` from the app by adding a normalization dropdown.
