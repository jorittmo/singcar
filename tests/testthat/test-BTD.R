
test_that("output is correct class", {
  outp <- BTD(-2, rnorm(15))
  expect_is(outp, "htest")
})


test_that("errors are occuring as they should", {

  controls <- rnorm(15)
  expect_error(BTD(NA, controls), "Case is NA")

  controls[1] <- NA
  expect_error(BTD(-2, controls), "Controls contains NA, set na.rm = TRUE to proceed")

  expect_error(BTD(-2, 1), "Not enough obs. Set sd and n for input of controls to be treated as mean")

  expect_error(BTD(-2, controls = 0, sd = 1), "Input sample size")

  expect_error(BTD(rnorm(2), controls, na.rm = TRUE), "Case should only have 1 observation")

  expect_error(BTD(2, controls, na.rm = TRUE, int_level = 2), "Interval level must be between 0 and 1")

  expect_error(BTD(2, controls, na.rm = TRUE, int_level = -2), "Interval level must be between 0 and 1")

})

test_that("TD and BTD gives equal output", {

  p_bdt <- BTD(-2, controls = 0, sd = 1, sample_size = 15)[["p.value"]]
  p_td <- TD(-2, controls = 0, sd = 1, sample_size = 15)[["p.value"]]

  expect_equal(round(p_bdt, 2), round(p_td, 2))

  p_bdt <- BTD(-2, controls = 0, sd = 1, sample_size = 15, alternative = "g")[["p.value"]]
  p_td <- TD(-2, controls = 0, sd = 1, sample_size = 15, alternative = "g")[["p.value"]]

  expect_equal(round(p_bdt, 2), round(p_td, 2))

  p_bdt <- BTD(2, controls = 0, sd = 1, sample_size = 15, alternative = "t", iter = 1000)[["p.value"]]
  p_td <- TD(2, controls = 0, sd = 1, sample_size = 15, alternative = "t")[["p.value"]]

  expect_equal(round(p_bdt, 1), round(p_td, 1))

  p_bdt <- BTD(-2, controls = 0, sd = 1, sample_size = 15, alternative = "t", iter = 1000)[["p.value"]]
  p_td <- TD(-2, controls = 0, sd = 1, sample_size = 15, alternative = "t")[["p.value"]]

  expect_equal(round(p_bdt, 1), round(p_td, 1))



})
