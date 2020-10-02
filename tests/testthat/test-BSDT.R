test_that("check that summary stats give same output as non-summary for BSDT", {
  conx <- c(-0.49094434, -1.25947544, -0.40339230, -0.12673867,
            0.43919457, -0.11994170,  0.68001403, -1.12678679,
            -0.91987671,  0.49401765, 1.07773774,  0.09122404,
            -0.12151505, -0.76224593, -1.44645118,  1.05976883,
            2.59665103,  1.27920207, -0.75500952, -0.18543234)

  cony <- c(-0.23357273, -2.09724146, -2.13526191,  0.89081394,
            0.76778599,  0.61269488,  0.54535923, -1.30459451,
            -0.25811271,  0.64790674, 0.52118776, -0.33601990,
            -0.37181171, -0.31255830, -0.03610260,  0.85230736,
            2.06553374, -0.10455893,  0.33921656, -0.05297144)


set.seed(123)
p_nsum <- BSDT(case_a = -2, case_b = -3, controls_a = conx, controls_b = cony, iter = 1000)[["p.value"]]


set.seed(123)
p_sum <- BSDT(case_a = -2, case_b = -3, controls_a = mean(conx), controls_b = mean(cony),
              sd_a = stats::sd(conx), sd_b = stats::sd(cony),
              sample_size = 20, r_ab = cor(conx, cony), iter = 1000)[["p.value"]]

expect_equal(p_nsum, p_sum)

})

test_that("check that BSDT produces same as C and G", {
  conx <- c(-0.49094434, -1.25947544, -0.40339230, -0.12673867,
            0.43919457, -0.11994170,  0.68001403, -1.12678679,
            -0.91987671,  0.49401765, 1.07773774,  0.09122404,
            -0.12151505, -0.76224593, -1.44645118,  1.05976883,
            2.59665103,  1.27920207, -0.75500952, -0.18543234)

  cony <- c(-0.23357273, -2.09724146, -2.13526191,  0.89081394,
            0.76778599,  0.61269488,  0.54535923, -1.30459451,
            -0.25811271,  0.64790674, 0.52118776, -0.33601990,
            -0.37181171, -0.31255830, -0.03610260,  0.85230736,
            2.06553374, -0.10455893,  0.33921656, -0.05297144)

  set.seed(123)
  p_nsum <- BSDT(case_a = -2, case_b = -3, controls_a = conx, controls_b = cony, iter = 10000)[["p.value"]]

  CandGPval <- 0.27849

  expect_equal(p_nsum, CandGPval, tol = 0.05)

})

test_that("check that BSDT produce roughly same output for chol_sim", {

  conx <- c(-0.49094434, -1.25947544, -0.40339230, -0.12673867,
            0.43919457, -0.11994170,  0.68001403, -1.12678679,
            -0.91987671,  0.49401765, 1.07773774,  0.09122404,
            -0.12151505, -0.76224593, -1.44645118,  1.05976883,
            2.59665103,  1.27920207, -0.75500952, -0.18543234)

  cony <- c(-0.23357273, -2.09724146, -2.13526191,  0.89081394,
            0.76778599,  0.61269488,  0.54535923, -1.30459451,
            -0.25811271,  0.64790674, 0.52118776, -0.33601990,
            -0.37181171, -0.31255830, -0.03610260,  0.85230736,
            2.06553374, -0.10455893,  0.33921656, -0.05297144)

  set.seed(123)
  p_cholsim <- BSDT(case_a = -2, case_b = -3, controls_a = conx, controls_b = cony, iter = 10000, chol_sim = TRUE)[["p.value"]]
  set.seed(123)
  p_loopchol <- BSDT(case_a = -2, case_b = -3, controls_a = conx, controls_b = cony, iter = 10000, chol_sim = FALSE)[["p.value"]]


  expect_equal(p_cholsim, p_loopchol, tol = 0.01)

})


test_that("errors and warnings are occuring as they should for BSDT", {

  conx <- rnorm(15)
  cony <- rnorm(15)

  cony[1] <- NA

  expect_error(RSDT(-2, -4, conx, cony),
               "Controls contains NA, set na.rm = TRUE to proceed")

  expect_warning(BSDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on controls_b resulted in removal of non-NAs on controls_a")
  expect_warning(BSDT(-2, -4, cony, rnorm(15), na.rm = TRUE),
                 "Removal of NAs on controls_a resulted in removal of non-NAs on controls_b")
  conx[2] <- NA

  expect_warning(BSDT(-2, -4, conx, cony, na.rm = TRUE),
                 "Removal of NAs on one control sample resulted in removal of non-NAs on the other")

  expect_error(BSDT(-2, -4, rnorm(15), rnorm(20)), "Sample sizes must be equal")

  expect_error(BSDT(c(-2,-4), 1, rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(BSDT(1, c(-2,-4), rnorm(15), rnorm(15)), "Case scores should be single value")

  expect_error(BSDT(NA, -2, rnorm(15), rnorm(15)), "One or both case scores is NA")
  expect_error(BSDT(-2, NA, rnorm(15), rnorm(15)), "One or both case scores is NA")

  expect_error(BSDT(-2, -4, 0, 0),
               "Please give sd and n on task A if controls_a is to be treated as mean")
  expect_error(BSDT(-2, -4, rnorm(20), 0),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_error(BSDT(-2, -4, rnorm(20), 0, sd_b = 1), "Please set sample size")

  expect_error(BSDT(-2, -4, rnorm(20), 0, sample_size = 20),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_error(BSDT(-2, -4, rnorm(20), 0, sample_size = 20),
               "Please give sd and n on task B if controls_b is to be treated as mean")

  expect_message(BSDT(-2, -4, rnorm(20), rnorm(20), sd_a = 1),
                 "Value on sd_a will be ignored")

  expect_message(BSDT(-2, -4, rnorm(20), rnorm(20), sd_b = 1),
                 "Value on sd_b will be ignored")

  expect_message(BSDT(-2, -4, rnorm(20), rnorm(20), sample_size = 15),
                 "Value on sample_size will be ignored")


})

