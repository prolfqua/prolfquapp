# Fix index.html relative links

## Problem

`write_index_html()` can write Windows absolute paths into `href` attributes.
When served by the `pd_metaboanalyst` local HTTP server, a link like
`E:/CD_output/.../Results.../DE_....html` is resolved relative to the current
HTTP directory and becomes a broken URL containing both the DEA folder and the
Windows drive path.

## Root Cause

The link builder uses `gsub(result_dir, "", p)` to derive relative paths. That
is regex-based, separator-sensitive, and unsafe for Windows paths with drive
letters and backslashes.

## Plan

- [x] Add an internal path helper that converts Windows and Unix paths to URL
  separators.
- [x] Derive links relative to the index output directory with prefix matching,
  not regex replacement.
- [x] URL-encode path segments while preserving `/`.
- [x] Add tests that Windows-style absolute paths produce relative browser
  links such as `./Results_.../DE_....html` and never include `E:/...`.
- [x] Run targeted tests and `git diff --check`.

## Verification

- `Rscript -e "devtools::test(filter = '^report_helpers$')"`
