
test_that("output is correct class", {
  outp <- TD(-2, rnorm(15))
  expect_is(outp, "htest")
})


test_that("errors are occuring as they should", {

  controls <- rnorm(15)
  expect_error(TD(NA, controls), "Case is NA")

  controls[1] <- NA
  expect_error(TD(-2, controls), "Controls contains NA, set na.rm = TRUE to proceed")

  expect_error(TD(-2, 1), "Not enough obs. Set sd and sample size for input of controls to be treated as mean")

  expect_error(TD(-2, controls = 0, sd = 1), "Input sample size")

  expect_error(TD(rnorm(2), controls, na.rm = TRUE), "Case should only have 1 observation")

  expect_error(TD(2, controls, na.rm = TRUE, conf_level = 2), "Confident level must be between 0 and 0.9999999")

  expect_error(TD(2, controls, na.rm = TRUE, conf_level = -2), "Confident level must be between 0 and 0.9999999")

})

test_that("summary input works as expected", {

  t1 <- TD(-2, controls = 0, sd = 1, sample_size = 15)[["statistic"]][["t"]]

  expect_equal(round(t1, 3), -1.936)

  t2 <- TD(2, controls = 0, sd = 1, sample_size = 15)[["statistic"]][["t"]]

  expect_equal(round(t2, 3), 1.936)

})
