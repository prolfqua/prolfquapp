.resolve_facade_model <- function(model, model_missing = FALSE) {
  if (identical(model, "prolfqua")) {
    if (isTRUE(model_missing)) {
      return("lm_impute")
    }
    return("lm")
  }
  model
}
