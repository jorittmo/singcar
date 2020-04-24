context("ch.ttest produce accurate output")

test_that("output is correct class", {
  outp <- ch.ttest(-2, rnorm(15))
  expect_is(outp, "htest")
})


test_that("errors are occuring as they should", {

  controls <- rnorm(15)
  expect_error(ch.ttest(NA, controls), "Case is NA")

  controls[1] <- NA
  expect_error(ch.ttest(-2, controls), "controls contains NA, set na.rm = TRUE to proceed")

  expect_error(ch.ttest(-2, 1), "Controls do not have enough observations")

  expect_error(ch.ttest(rnorm(2), controls, na.rm = TRUE), "Case should only have 1 observation")

})
