#' capture output of function to send it to log
#' @export
capture_output <- function(expr) {
  con <- textConnection("output", "w", local = TRUE)
  sink(con)
  on.exit({
    sink()
    close(con)
  })
  eval(expr)
  paste(output, collapse = "\n")
}
