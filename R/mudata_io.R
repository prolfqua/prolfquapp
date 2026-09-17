#' Write a MuData container atomically
#'
#' Preserves complete AnnData modalities and shared sample metadata. Both axes
#' are shared at the container level; individual modalities retain their own
#' ordered axes. Temporary serialization files are private to this operation.
#'
#' @param modalities Named list of anndataR AnnData objects.
#' @param path Destination H5MU file.
#' @param obs Shared sample annotation, indexed by sample identifiers.
#' @param uns Container metadata as a portable named list.
#' @return Destination path, invisibly.
#' @export
write_h5mu <- function(modalities, path, obs, uns = list()) {
  .validate_mudata_inputs(modalities, obs)
  if (!dir.exists(dirname(path))) {
    stop("MuData output directory does not exist: ", dirname(path))
  }
  temporary <- tempfile(".h5mu-write-", tmpdir = dirname(path), fileext = ".h5mu")
  on.exit(unlink(temporary), add = TRUE)
  .write_mudata_container(modalities, obs, uns, temporary)
  restored <- read_h5mu(temporary)
  if (!identical(names(restored$modalities), names(modalities))) {
    stop("MuData round-trip changed modality order.")
  }
  for (name in names(modalities)) {
    original <- modalities[[name]]
    recovered <- restored$modalities[[name]]
    for (slot in c("X", "obs", "var", "layers", "varm", "obsm", "varp", "obsp", "uns")) {
      comparison <- .compare_mudata_slot(original, recovered, slot)
      if (!isTRUE(comparison)) {
        stop("MuData round-trip changed ", name, "/", slot, ": ", paste(comparison, collapse = "; "))
      }
    }
  }
  if (
    !isTRUE(all.equal(obs, restored$obs, check.attributes = FALSE)) ||
      !isTRUE(all.equal(.mudata_value(uns), .mudata_value(restored$uns), check.attributes = FALSE))
  ) {
    stop("MuData round-trip changed container metadata.")
  }
  if (!file.rename(temporary, path)) {
    stop("Could not publish MuData: ", path)
  }
  invisible(normalizePath(path, mustWork = TRUE))
}

.compare_mudata_slot <- function(original, recovered, slot) {
  left <- original[[slot]]
  right <- recovered[[slot]]
  if (slot %in% c("layers", "varm", "obsm", "varp", "obsp", "uns")) {
    left <- as.list(left)
    right <- as.list(right)
  }
  all.equal(.mudata_value(left), .mudata_value(right), check.attributes = FALSE)
}

.validate_mudata_inputs <- function(modalities, obs) {
  .validate_mudata_names(names(modalities))
  keys <- names(modalities)
  if (!is.data.frame(obs) || anyDuplicated(rownames(obs))) {
    stop("MuData obs must be a data frame with unique sample names.")
  }
  for (name in keys) {
    value <- modalities[[name]]
    if (!inherits(value, "AbstractAnnData")) {
      stop("MuData modality is not AnnData: ", name)
    }
    if (!setequal(value$obs_names, rownames(obs))) {
      stop("MuData modality sample set differs from shared obs: ", name)
    }
  }
}

.mudata_attribute <- function(handle, name, value, scalar = TRUE) {
  rhdf5::h5writeAttribute(
    value,
    handle,
    name,
    asScalar = scalar,
    encoding = "UTF-8",
    variableLengthString = TRUE
  )
}

.write_mudata_container <- function(modalities, obs, uns, path) {
  carrier_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(carrier_path), add = TRUE)
  feature_names <- Reduce(union, lapply(modalities, function(value) value$var_names))
  obsmap <- lapply(modalities, function(value) matrix(as.integer(match(rownames(obs), value$obs_names)), ncol = 1L))
  varmap <- lapply(modalities, function(value) {
    positions <- match(feature_names, value$var_names)
    positions[is.na(positions)] <- 0L
    matrix(as.integer(positions), ncol = 1L)
  })
  carrier <- anndataR::AnnData(
    obs = obs,
    var = data.frame(row.names = feature_names),
    uns = list(container = uns, obsmap = obsmap, varmap = varmap),
    obsm = lapply(obsmap, function(value) matrix(value > 0L, ncol = 1L)),
    varm = lapply(varmap, function(value) matrix(value > 0L, ncol = 1L))
  )
  .write_mudata_h5ad(carrier, carrier_path)
  rhdf5::h5createFile(path)
  destination <- rhdf5::H5Fopen(path)
  on.exit(rhdf5::H5Fclose(destination), add = TRUE)
  source <- rhdf5::H5Fopen(carrier_path, flags = "H5F_ACC_RDONLY")
  on.exit(rhdf5::H5Fclose(source), add = TRUE)
  .mudata_attribute(destination, "encoding-type", "MuData")
  .mudata_attribute(destination, "encoding-version", "0.1.0")
  .mudata_attribute(destination, "axis", -1L)
  .mudata_attribute(destination, "encoder", "prolfquapp")
  .mudata_attribute(destination, "encoder-version", as.character(utils::packageVersion("prolfquapp")))
  for (slot in c("obs", "var", "obsm", "varm", "obsp", "varp")) {
    rhdf5::H5Ocopy(source, slot, destination, slot)
  }
  for (slot in c("obsmap", "varmap")) {
    rhdf5::H5Ocopy(source, paste0("uns/", slot), destination, slot)
  }
  rhdf5::H5Ocopy(source, "uns/container", destination, "uns")
  rhdf5::h5createGroup(destination, "mod")
  group <- rhdf5::H5Gopen(destination, "mod")
  .mudata_attribute(group, "mod-order", names(modalities), scalar = FALSE)
  rhdf5::H5Gclose(group)
  for (name in names(modalities)) {
    .write_mudata_modality(modalities[[name]], destination, name)
  }
}

.write_mudata_modality <- function(adata, destination, name) {
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  .write_mudata_h5ad(adata, path)
  source <- rhdf5::H5Fopen(path, flags = "H5F_ACC_RDONLY")
  on.exit(rhdf5::H5Fclose(source), add = TRUE)
  rhdf5::H5Ocopy(source, "/", destination, paste0("mod/", name))
}

#' Read a complete MuData container
#'
#' @param path Path to an H5MU file.
#' @return List with named AnnData `modalities`, shared `obs`, and container `uns`.
#' @export
read_h5mu <- function(path) {
  attributes <- rhdf5::h5readAttributes(path, "/")
  if (!identical(attributes[["encoding-type"]], "MuData")) {
    stop("Expected a MuData container: ", path)
  }
  keys <- as.character(rhdf5::h5readAttributes(path, "mod")[["mod-order"]])
  if (length(keys) == 0L || anyDuplicated(keys)) {
    stop("MuData has no valid modality order: ", path)
  }
  modalities <- lapply(keys, function(name) .read_mudata_modality(path, name))
  names(modalities) <- keys
  shared <- .read_mudata_shared(path, modalities[[1L]]$obs)
  .validate_mudata_inputs(modalities, shared$obs)
  list(modalities = modalities, obs = shared$obs, uns = shared$uns)
}

.read_mudata_modality <- function(path, name) {
  temporary <- tempfile(fileext = ".h5ad")
  on.exit(unlink(temporary), add = TRUE)
  rhdf5::h5createFile(temporary)
  source <- rhdf5::H5Fopen(path, flags = "H5F_ACC_RDONLY")
  on.exit(rhdf5::H5Fclose(source), add = TRUE)
  destination <- rhdf5::H5Fopen(temporary)
  group <- rhdf5::H5Gopen(source, paste0("mod/", name))
  for (slot in rhdf5::h5ls(group, recursive = FALSE)$name) {
    rhdf5::H5Ocopy(group, slot, destination, slot)
  }
  .mudata_attribute(destination, "encoding-type", "anndata")
  .mudata_attribute(destination, "encoding-version", "0.1.0")
  rhdf5::H5Gclose(group)
  rhdf5::H5Fclose(destination)
  .read_mudata_h5ad(temporary)
}

.read_mudata_shared <- function(path, obs) {
  temporary <- tempfile(fileext = ".h5ad")
  on.exit(unlink(temporary), add = TRUE)
  carrier <- anndataR::AnnData(obs = obs, var = data.frame(row.names = ".metadata"))
  carrier$write_h5ad(temporary)
  source <- rhdf5::H5Fopen(path, flags = "H5F_ACC_RDONLY")
  on.exit(rhdf5::H5Fclose(source), add = TRUE)
  destination <- rhdf5::H5Fopen(temporary)
  for (slot in c("obs", "uns")) {
    rhdf5::h5delete(destination, slot)
    rhdf5::H5Ocopy(source, slot, destination, slot)
  }
  rhdf5::H5Fclose(destination)
  restored <- .read_mudata_h5ad(temporary)
  list(obs = as.data.frame(restored$obs), uns = restored$uns)
}

# anndataR currently encodes one-row dataframe columns (including the index)
# as scalars. AnnData requires arrays in these positions. Repair only that
# serialization defect; no data or installed dependency is modified.
.write_mudata_h5ad <- function(adata, path) {
  adata$write_h5ad(path, compression = "gzip", mode = "w")
  .repair_mudata_booleans(adata, path)
  frames <- c("obs", "var", .mudata_frame_paths(adata$uns, "uns"))
  for (group in frames) {
    attributes <- rhdf5::h5readAttributes(path, group)
    columns <- c(attributes[["_index"]], attributes[["column-order"]])
    for (column in columns) {
      .mudata_column_array(path, paste0(group, "/", column))
    }
  }
}

.mudata_column_array <- function(path, name) {
  encoding <- rhdf5::h5readAttributes(path, name)[["encoding-type"]]
  if (!encoding %in% c("string", "numeric-scalar")) {
    return(invisible(NULL))
  }
  value <- rhdf5::h5read(path, name)
  rhdf5::h5delete(path, name)
  rhdf5::h5write(array(value, dim = length(value)), path, name)
  handle <- rhdf5::H5Fopen(path)
  on.exit(rhdf5::H5Fclose(handle), add = TRUE)
  dataset <- rhdf5::H5Dopen(handle, name)
  on.exit(rhdf5::H5Dclose(dataset), add = TRUE)
  .mudata_attribute(dataset, "encoding-type", if (is.character(value)) "string-array" else "array")
  .mudata_attribute(dataset, "encoding-version", "0.2.0")
}

.validate_mudata_names <- function(keys) {
  if (length(keys) == 0L) {
    stop("MuData requires named modalities.")
  }
  invalid <- c(anyNA(keys), any(!nzchar(keys)), anyDuplicated(keys) > 0L, any(grepl("/", keys, fixed = TRUE)))
  if (any(invalid)) stop("MuData requires uniquely named modalities without slashes.")
}

# Python's newer nullable string encoding is not yet read by anndataR.
# Decode it in the private copy as categorical, then restore character columns.
.read_mudata_h5ad <- function(path) {
  nodes <- rhdf5::h5ls(path)
  candidates <- unique(nodes$group[nodes$name == "mask"])
  converted <- character()
  for (name in candidates) {
    if (identical(rhdf5::h5readAttributes(path, name)[["encoding-type"]], "nullable-string-array")) {
      .decode_mudata_strings(path, name)
      converted <- c(converted, name)
    }
  }
  adata <- anndataR::read_h5ad(path)
  for (name in converted) {
    keys <- strsplit(sub("^/", "", name), "/", fixed = TRUE)[[1L]]
    slot <- keys[[1L]]
    adata[[slot]] <- .restore_mudata_strings(adata[[slot]], keys[-1L])
  }
  adata
}

.decode_mudata_strings <- function(path, name) {
  values <- as.character(rhdf5::h5read(path, paste0(name, "/values")))
  mask <- as.logical(rhdf5::h5read(path, paste0(name, "/mask")))
  values[mask] <- NA_character_
  categories <- unique(values[!is.na(values)])
  codes <- match(values, categories) - 1L
  codes[mask] <- -1L
  rhdf5::h5delete(path, name)
  rhdf5::h5createGroup(path, name)
  rhdf5::h5write(array(categories, length(categories)), path, paste0(name, "/categories"))
  rhdf5::h5write(array(codes, length(codes)), path, paste0(name, "/codes"))
  handle <- rhdf5::H5Fopen(path)
  on.exit(rhdf5::H5Fclose(handle), add = TRUE)
  group <- rhdf5::H5Gopen(handle, name)
  .mudata_attribute(group, "encoding-type", "categorical")
  .mudata_attribute(group, "encoding-version", "0.2.0")
  .mudata_attribute(group, "ordered", FALSE)
  rhdf5::H5Gclose(group)
  for (child in c("categories", "codes")) {
    dataset <- rhdf5::H5Dopen(handle, paste0(name, "/", child))
    .mudata_attribute(dataset, "encoding-type", if (child == "categories") "string-array" else "array")
    .mudata_attribute(dataset, "encoding-version", "0.2.0")
    rhdf5::H5Dclose(dataset)
  }
}

.restore_mudata_strings <- function(value, keys) {
  if (!length(keys)) {
    return(as.character(value))
  }
  key <- keys[[1L]]
  if (is.data.frame(value) && !key %in% names(value)) {
    return(value)
  }
  value[[key]] <- .restore_mudata_strings(value[[key]], keys[-1L])
  value
}

.mudata_frame_paths <- function(value, path) {
  if (is.data.frame(value)) {
    return(path)
  }
  if (!is.list(value)) {
    return(character())
  }
  unlist(
    Map(function(element, key) .mudata_frame_paths(element, paste0(path, "/", key)), value, names(value)),
    use.names = FALSE
  )
}

# anndataR's boolean writer drops matrix dimensions when converting to integer.
.write_mudata_boolean <- function(path, name, value) {
  rhdf5::h5delete(path, name)
  handle <- rhdf5::H5Fopen(path)
  on.exit(rhdf5::H5Fclose(handle), add = TRUE)
  space <- rhdf5::H5Screate_simple(dim(value), maxdims = NULL, native = TRUE)
  on.exit(rhdf5::H5Sclose(space), add = TRUE)
  type <- rhdf5::H5Tenum_create("H5T_NATIVE_SCHAR")
  rhdf5::H5Tenum_insert(type, "FALSE", 0L)
  rhdf5::H5Tenum_insert(type, "TRUE", 1L)
  dataset <- rhdf5::H5Dcreate(handle, name, dtype_id = type, h5space = space)
  on.exit(rhdf5::H5Dclose(dataset), add = TRUE)
  rhdf5::H5Dwrite(dataset, as.raw(as.vector(t(value))), h5type = type)
  .mudata_attribute(dataset, "encoding-type", "array")
  .mudata_attribute(dataset, "encoding-version", "0.2.0")
}

.repair_mudata_booleans <- function(adata, path) {
  for (slot in c("obsm", "varm", "layers", "obsp", "varp")) {
    for (key in names(adata[[slot]])) {
      .repair_mudata_boolean(adata[[slot]][[key]], path, paste0(slot, "/", key))
    }
  }
}

.repair_mudata_boolean <- function(value, path, key) {
  if (is.matrix(value) && is.logical(value) && !anyNA(value)) {
    .write_mudata_boolean(path, key, value)
  }
}

.mudata_value <- function(value) {
  if (is.factor(value)) {
    return(as.character(value))
  }
  if (is.data.frame(value)) {
    return(structure(lapply(value, .mudata_value), class = "data.frame", row.names = rownames(value)))
  }
  if (is.list(value)) {
    if (!is.null(names(value))) {
      value <- value[sort(names(value))]
    }
    return(lapply(value, .mudata_value))
  }
  if (length(dim(value)) <= 1L) {
    return(unname(as.vector(value)))
  }
  value
}
