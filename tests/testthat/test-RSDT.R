test_that("output is correct class", {
  outp <- RSDT(-2, -3, rnorm(15), rnorm(15))
  expect_is(outp, "htest")
})

test_that("all RSDT functionality works", {
  conx <- rnorm(15)
  cony <- rnorm(15)
  expect_error(
    RSDT(-2, 0, conx, cony, alternative = "g"),
    NA
  )
  expect_error(
    RSDT(-2, 0, conx, cony, alternative = "l"),
    NA
  )
  cony[1] <- NA
  expect_warning(
    RSDT(-2, 0, conx, cony, alternative = "t", na.rm = TRUE),
    "Removal of NAs on controls_b resulted in removal of non-NAs on controls_a"
  )

})

test_that("summary input gives same result as raw", {
  x <- MASS::mvrnorm(20, mu = c(100, 80),
                     Sigma = matrix(c(15^2, 108,
                                      108, 10^2),
                                    nrow = 2, byrow = T),
                     empirical = TRUE)

  expect_equal(
    RSDT(70, 80, 100, 80, sd_a = 15, sd_b = 10, sample_size = 20,
        r_ab = 0.72, alternative = "t"),
    RSDT(70, 80, x[ , 1], x[ , 2],  alternative = "t")
  )

})




test_that("errors and warnings are occuring as they should for RSDT", {

  conx <- rnorm(15)
  cony <- rnorm(15)

  cony[1] <- NA

  expect_error(RSDT(-2, -4, conx, cony),
               "Controls contains NA, set na.rm = TRUE to proceed")

  expect_warning(RSDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on controls_b resulted in removal of non-NAs on controls_a")
  expect_warning(RSDT(-2, -4, cony, rnorm(15), na.rm = TRUE),
                 "Removal of NAs on controls_a resulted in removal of non-NAs on controls_b")
  conx[2] <- NA

  expect_warning(RSDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on one control sample resulted in removal of non-NAs on the other")

  expect_error(RSDT(-2, -4, rnorm(15), rnorm(20)), "Sample sizes must be equal")

  expect_error(RSDT(c(-2,-4), 1, rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(RSDT(1, c(-2,-4), rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(RSDT(NA, -2, rnorm(15), rnorm(15)), "One or both case scores is NA")
  expect_error(RSDT(-2, NA, rnorm(15), rnorm(15)), "One or both case scores is NA")

  expect_error(RSDT(-2, -4, 0, 0),
               "Please give sd and n on task A if controls_a is to be treated as mean")
  expect_error(RSDT(-2, -4, rnorm(20), 0),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_error(RSDT(-2, -4, rnorm(20), 0, sd_b = 1), "Please set sample size")

  expect_error(RSDT(-2, -4, rnorm(20), 0, sample_size = 20),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_error(RSDT(-2, -4, rnorm(20), 0, sample_size = 20),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_message(RSDT(-2, -4, rnorm(20), rnorm(20), sd_a = 1),
               "Value on sd_a will be ignored")

  expect_message(RSDT(-2, -4, rnorm(20), rnorm(20), sd_b = 1),
                 "Value on sd_b will be ignored")

  expect_message(RSDT(-2, -4, rnorm(20), rnorm(20), sample_size = 15),
                 "Value on sample_size will be ignored")


})
