#' Route base messages and warnings into the prolfquapp logger
#'
#' Installs global calling handlers so that any `message()` or `warning()`
#' signalled by prolfqua core, prolfquapp, or their dependencies is re-emitted
#' through [logger::log_info()] / [logger::log_warn()] with the standard
#' `INFO [timestamp]` layout and captured by the active log appender, instead of
#' printing as an untagged line on `stderr`. Also sets
#' `options(readr.show_col_types = FALSE)` so readr's column-specification output
#' (`Rows: .. Columns: ..`, delimiter, `Use spec()` hints) is silenced -- pure
#' noise in batch runs.
#'
#' Call once near the top of a command-line entry script, after the initial
#' [logger::log_appender()] setup. This is intended for the prolfquapp CLI entry
#' points; library and interactive callers keep the default `message()` /
#' `warning()` behaviour.
#'
#' @return Invisibly `NULL`; called for its side effects.
#' @export
route_messages_to_logger <- function() {
  options(readr.show_col_types = FALSE)
  globalCallingHandlers(
    message = function(m) {
      msg <- sub("\n$", "", conditionMessage(m))
      if (nzchar(msg)) {
        logger::log_info(.escape_braces(msg))
      }
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      logger::log_warn(.escape_braces(conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  invisible(NULL)
}

# logger formats messages through glue, which would try to evaluate any {..}
# braces present in a forwarded message/warning. Double them so glue emits them
# verbatim.
.escape_braces <- function(x) {
  gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
}
