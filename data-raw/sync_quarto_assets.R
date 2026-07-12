## Keep the vendored Quarto assets used by package vignettes in step with the
## installed fgczquartotemplate package.
##
## `devtools::build_vignettes()` renders the QMD sources directly, so it does
## not call fgczquartotemplate::fgcz_render() to stage assets at runtime. The
## reports refer to the assets by bare filename; copy and verify them here
## before package or vignette builds.
##
##   Rscript data-raw/sync_quarto_assets.R

assets <- c(
  "_metadata.yml",
  "fgcz.scss",
  "fgcz_header_quarto.html",
  "fgcz-plot-finder.html"
)
target_dir <- normalizePath("vignettes", mustWork = TRUE)
source_files <- fgczquartotemplate::fgcz_quarto_dir(assets)
target_files <- file.path(target_dir, assets)

fgczquartotemplate::fgcz_copy_assets(target_dir)

if (
  !identical(
    unname(tools::md5sum(source_files)),
    unname(tools::md5sum(target_files))
  )
) {
  stop(
    "FGCZ Quarto assets differ after synchronization: ",
    paste(assets, collapse = ", ")
  )
}

message(
  "OK: synchronized FGCZ Quarto assets into vignettes/: ",
  paste(assets, collapse = ", ")
)
