context("ch.ttest produce accurate output")

test_that("output is correct class", {
  outp <- ch.ttest(-2, rnorm(15))
  expect_is(outp, "htest")
})


test_that("errors are occuring as they should", {

  controls <- rnorm(15)
  expect_error(ch.ttest(NA, controls), "Case is NA")

  controls[1] <- NA
  expect_error(ch.ttest(-2, controls), "Controls contains NA, set na.rm = TRUE to proceed")

  expect_error(ch.ttest(-2, 1), "Not enough obs. Set sd and n for input of controls to be treated as mean")

  expect_error(ch.ttest(rnorm(2), controls, na.rm = TRUE), "Case should only have 1 observation")

})

test_that("summary input works as expected", {

  t1 <- ch.ttest(-2, controls = 0, controls.sd = 1, controls.n = 15)[["statistic"]][["t"]]

  expect_equal(round(t1, 3), -1.936)

  t1 <- ch.ttest(2, controls = 0, controls.sd = 1, controls.n = 15)[["statistic"]][["t"]]

  expect_equal(round(t1, 3), 1.936)

})
