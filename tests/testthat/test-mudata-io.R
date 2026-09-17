test_that("MuData preserves overlapping modalities and root metadata", {
  obs <- data.frame(Group = c("A", "B"), row.names = c("s1", "s2"))
  enriched <- anndataR::AnnData(
    X = matrix(c(1, NA, 3, 4), 2L),
    obs = obs,
    var = data.frame(protein_Id = c("p1", "p1"), row.names = c("site1", "site2")),
    varm = list(dea = matrix(c(0.1, 0.2), ncol = 1)),
    uns = list(columns = "diff")
  )
  total <- anndataR::AnnData(X = matrix(c(1, 2), ncol = 1), obs = obs, var = data.frame(row.names = "p1"))
  cf <- anndataR::AnnData(X = matrix(c(0, NA), ncol = 1), obs = obs, var = data.frame(row.names = "site1"))
  path <- tempfile(fileext = ".h5mu")
  on.exit(unlink(path))
  mods <- list(enriched = enriched, total = total, cf = cf)
  write_h5mu(mods, path, obs, list(stage = "PTM_statistics"))
  restored <- read_h5mu(path)
  expect_identical(names(restored$modalities), names(mods))
  expect_equal(restored$obs, obs)
  expect_equal(restored$uns$stage, "PTM_statistics")
  expect_equal(restored$modalities$enriched$X, enriched$X)
  expect_equal(restored$modalities$enriched$varm$dea, enriched$varm$dea)
  expect_equal(restored$modalities$cf$var_names, "site1")
  expect_equal(as.integer(rhdf5::h5read(path, "varmap/cf")), c(1L, 0L, 0L))
  expect_equal(rhdf5::h5readAttributes(path, "/")$axis, -1L)
})

test_that("invalid MuData cannot replace an existing artifact", {
  obs <- data.frame(row.names = c("s1", "s2"))
  adata <- anndataR::AnnData(X = matrix(1, 2, 1), obs = obs)
  path <- tempfile(fileext = ".h5mu")
  on.exit(unlink(path))
  write_h5mu(list(enriched = adata), path, obs)
  before <- tools::md5sum(path)
  expect_error(write_h5mu(list(enriched = adata), path, obs[1, , drop = FALSE]), "sample set")
  expect_identical(tools::md5sum(path), before)
})
