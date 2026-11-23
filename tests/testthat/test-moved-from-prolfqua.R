## script_copy_helper_vec() and separate_hierarchy() moved here from prolfqua 1.9.5

test_that("copying an absent package script warns", {
  expect_warning(
    script_copy_helper_vec("not-a-package-file", workdir = withr::local_tempdir()),
    "Could not copy script file"
  )
})

test_that("copying a script over an existing file without overwrite warns", {
  path <- withr::local_tempdir()
  file.create(file.path(path, "DESCRIPTION"))
  expect_warning(
    script_copy_helper_vec("DESCRIPTION", workdir = path, overwrite = FALSE),
    "Could not copy script file"
  )
})

test_that("copying a package file returns its new path", {
  path <- withr::local_tempdir()
  copied <- script_copy_helper_vec("DESCRIPTION", workdir = path)
  expect_identical(copied, file.path(path, "DESCRIPTION"))
  expect_true(file.exists(copied))
})

test_that("separate_hierarchy splits hierarchy keys at the separator into their configured columns", {
  bb <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 3)
  lfq <- prolfqua::LFQData$new(bb$data, bb$config)
  config <- lfq$get_config()$clone(deep = TRUE)
  data <- lfq$data_long()
  data$protein_Id[1] <- "noseparator"

  res <- separate_hierarchy(data, config)
  expect_identical(colnames(res)[4:6], c("protein_Id", "proteinID", "idtype2"))
  parts <- strsplit(data$protein_Id[-1], "~", fixed = TRUE)
  expect_identical(res$proteinID[-1], vapply(parts, `[`, character(1), 1))
  expect_identical(res$idtype2[-1], vapply(parts, `[`, character(1), 2))
  # too few parts fill from the start
  expect_identical(res$proteinID[1], "noseparator")
  expect_true(is.na(res$idtype2[1]))
  # keys below the hierarchy depth are left alone
  expect_false("peptideID" %in% colnames(res))

  config$hierarchy_depth <- 2
  deeper <- separate_hierarchy(data, config)
  expect_identical(deeper$peptideID, data$peptide_Id)
})
