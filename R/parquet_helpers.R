# Internal helpers for Parquet exports.

write_parquet_isolated <- function(
  data,
  sink,
  rscript = file.path(R.home("bin"), "Rscript")
) {
  dir.create(dirname(sink), showWarnings = FALSE, recursive = TRUE)

  input_file <- tempfile("prolfquapp_parquet_data_", fileext = ".rds")
  script_file <- tempfile("prolfquapp_parquet_writer_", fileext = ".R")
  on.exit(unlink(c(input_file, script_file)), add = TRUE)

  saveRDS(data, input_file)
  writeLines(
    c(
      "args <- commandArgs(trailingOnly = TRUE)",
      "data <- readRDS(args[[1]])",
      "arrow::write_parquet(data, sink = args[[2]])"
    ),
    con = script_file
  )

  output <- system2(
    rscript,
    c("--vanilla", script_file, input_file, sink),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- if (is.null(attr(output, "status"))) 0L else attr(output, "status")

  list(
    sink = sink,
    exit_status = status,
    output = output,
    success = identical(status, 0L) && file.exists(sink)
  )
}
