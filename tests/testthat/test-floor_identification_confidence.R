test_that("floor_identification_confidence raises confidences below 0.5 and keeps NA", {
  expect_equal(floor_identification_confidence(c(0, 0.01, 0.5, 0.9, 1, NA)), c(0.5, 0.5, 0.5, 0.9, 1, NA))
})

test_that("the APB reader floors 1 - pep at 0.5", {
  long <- data.frame(pep = c(0.9, 0.2, 0.001))
  expect_equal(.apb_identification_confidence(long, character()), c(0.5, 0.8, 0.999))
})
