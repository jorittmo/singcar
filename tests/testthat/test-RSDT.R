test_that("output is correct class", {
  outp <- RSDT(-2, -3, rnorm(15), rnorm(15))
  expect_is(outp, "htest")
})


test_that("errors and warnings are occuring as they should", {

  conx <- rnorm(15)
  cony <- rnorm(15)

  cony[1] <- NA

  expect_error(RSDT(-2, -4, conx, cony),
               "Controls contains NA, set na.rm = TRUE to proceed")

  expect_warning(RSDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on controls.y resulted in removal of non-NAs on controls.x")
  expect_warning(RSDT(-2, -4, cony, rnorm(15), na.rm = TRUE),
                 "Removal of NAs on controls.x resulted in removal of non-NAs on controls.y")
  conx[2] <- NA

  expect_warning(RSDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on one control sample resulted in removal of non-NAs on the other")

  expect_error(RSDT(-2, -4, rnorm(15), rnorm(20)), "Sample sizes must be equal")

  expect_error(RSDT(c(-2,-4), 1, rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(RSDT(1, c(-2,-4), rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(RSDT(NA, -2, rnorm(15), rnorm(15)), "One or both case scores is NA")
  expect_error(RSDT(-2, NA, rnorm(15), rnorm(15)), "One or both case scores is NA")

  expect_error(RSDT(-2, -4, 0, 0),
               "Please give sd and n on task x if controls.x is to be treated as mean")
  expect_error(RSDT(-2, -4, rnorm(20), 0),
               "Please give sd and n on task y if controls.y is to be treated as mean")

  expect_error(RSDT(-2, -4, rnorm(20), 0, controls.y.sd = 1), "Please set sample size")

  expect_error(RSDT(-2, -4, rnorm(20), 0, controls.n = 20),
               "Please give sd and n on task y if controls.y is to be treated as mean")

  expect_error(RSDT(-2, -4, rnorm(20), 0, controls.n = 20),
               "Please give sd and n on task y if controls.y is to be treated as mean")

  expect_message(RSDT(-2, -4, rnorm(20), rnorm(20), controls.x.sd = 1),
               "Value on controls.x.sd will be ignored")

  expect_message(RSDT(-2, -4, rnorm(20), rnorm(20), controls.y.sd = 1),
                 "Value on controls.y.sd will be ignored")

  expect_message(RSDT(-2, -4, rnorm(20), rnorm(20), controls.n = 15),
                 "Value on controls.n will be ignored")


})
