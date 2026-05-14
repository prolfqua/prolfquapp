.resolve_facade_model <- function(model, model_missing = FALSE) {
  if (identical(model, "prolfqua")) {
    if (isTRUE(model_missing)) {
      return("lm_impute")
    }
    return("lm")
  }
  model
}

.is_saint_model <- function(model) {
  identical(model, "saint")
}

.valid_facade_models <- function() {
  c(names(prolfqua::FACADE_REGISTRY), "saint")
}
