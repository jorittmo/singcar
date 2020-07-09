test_that("output is correct class", {
  outp <- UDT(-2, -3, rnorm(15), rnorm(15))
  expect_is(outp, "htest")
})


test_that("errors and warnings are occuring as they should", {

  conx <- rnorm(15)
  cony <- rnorm(15)

  cony[1] <- NA

  expect_error(UDT(-2, -4, conx, cony),
               "Controls contains NA, set na.rm = TRUE to proceed")

  expect_warning(UDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on controls_b resulted in removal of non-NAs on controls_a")
  expect_warning(UDT(-2, -4, cony, rnorm(15), na.rm = TRUE),
                 "Removal of NAs on controls_a resulted in removal of non-NAs on controls_b")
  conx[2] <- NA

  expect_warning(UDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on one control sample resulted in removal of non-NAs on the other")

  expect_error(UDT(-2, -4, rnorm(15), rnorm(20)), "Sample sizes must be equal")

  expect_error(UDT(c(-2,-4), 1, rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(UDT(1, c(-2,-4), rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(UDT(NA, -2, rnorm(15), rnorm(15)), "One or both case scores is NA")
  expect_error(UDT(-2, NA, rnorm(15), rnorm(15)), "One or both case scores is NA")

  expect_error(UDT(-2, -4, 0, 0),
               "Please give sd and n on task A if controls_a is to be treated as mean")
  expect_error(UDT(-2, -4, rnorm(20), 0),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_error(UDT(-2, -4, rnorm(20), 0, sd_b = 1), "Please set sample size")

  expect_error(UDT(-2, -4, rnorm(20), 0, sample_size = 20),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_error(UDT(-2, -4, rnorm(20), 0, sample_size = 20),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_message(UDT(-2, -4, rnorm(20), rnorm(20), sd_a = 1),
                 "Value on sd_a will be ignored")

  expect_message(UDT(-2, -4, rnorm(20), rnorm(20), sd_b = 1),
                 "Value on sd_b will be ignored")

  expect_message(UDT(-2, -4, rnorm(20), rnorm(20), sample_size = 15),
                 "Value on sample_size will be ignored")


})
