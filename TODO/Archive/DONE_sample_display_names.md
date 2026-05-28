# TODO: Derived Short Sample Display Names

## Goal

Avoid plot/report breakage from very long annotation `Name` values without overwriting the original annotation column.

## Plan

1. Add annotation-processor support for a derived display sample-name column.
   - Keep the original `Name` column unchanged.
   - When sample names are longer than the configured suffix length, derive a `sampleName` column from the right-most suffix.
   - Check whether derived labels are unique and make them unique if needed.
   - Also avoid destructive mutation for duplicate short names by creating a derived unique display column.

2. Keep downstream behavior stable for normal annotations.
   - If names are already short and unique, keep using the original `Name` column as `atable$sample_name`.
   - Use the derived `sampleName` only when shortening or de-duplication is needed.

3. Add focused regression tests for:
   - short unique names preserving `Name` as the sample-name column;
   - long names producing a short derived `sampleName` while preserving `Name`;
   - duplicate suffixes being made unique;
   - duplicate short names using a derived unique column instead of mutating `Name`.

4. Regenerate roxygen documentation and run targeted tests.
