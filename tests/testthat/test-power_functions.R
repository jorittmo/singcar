
###########################################

# TD_power

###########################################

test_that("all TD_power functionality works", {
  expect_error(
    TD_power(-2, power = 0.25),
    NA)
  expect_error(
    TD_power(-2, sample_size =  10),
    NA)
  expect_error(
    TD_power(-2, sample_size =  10, alternative = "l"),
    NA)
  expect_error(
    TD_power(-2, sample_size =  10, alternative = "g"),
    NA)
  expect_error(
    TD_power(-2, sample_size =  10, alternative = "t"),
    NA)

})


test_that("errors are working for TD_power", {
  expect_error(TD_power(-2),
               "Must supply either sample size or desired power")
  expect_error(TD_power(-2, sample_size = 20, power = 0.8),
               "Must supply only one of sample size or desired power")
  expect_error(TD_power(-2, power = 1.2),
               "Desired power must be between 0 and 1")
  expect_error(TD_power(-2, power = -1),
               "Desired power must be between 0 and 1")
  expect_error(TD_power(-2, sample_size = 20, alpha = 1.2),
               "Type I error rate must be between 0 and 1")
  expect_error(TD_power(-2, sample_size = 20, alpha = -1),
               "Type I error rate must be between 0 and 1")

})


test_that("output is approx the same as in Crawford & Garthwaite (2006), TD_power", {

  CG <- c(44.83, 54.63, 59.46, 62.06, 63.00,
          72.80, 83.87, 87.88, 90.11, 90.73)/100

  x <- c(TD_power(-2, sample_size = 5),
         TD_power(-2, sample_size = 10),
         TD_power(-2, sample_size = 20),
         TD_power(-2, sample_size = 50),
         TD_power(-2, sample_size = 100),
         TD_power(-3, sample_size = 5),
         TD_power(-3, sample_size = 10),
         TD_power(-3, sample_size = 20),
         TD_power(-3, sample_size = 50),
         TD_power(-3, sample_size = 100))

  expect_equal(x, CG, tolerance = 1e-3)



})

test_that("TD_power and UDT_power produce equal output when r = 0.5", {


  expect_equal(TD_power(-2, sample_size = 15, alternative = "less"),
               UDT_power(-2, 0, r_ab = 0.5, sample_size = 15, alternative = "less"))

  expect_equal(TD_power(-2, sample_size = 15, alternative = "two.sided"),
               UDT_power(-2, 0, r_ab = 0.5, sample_size = 15, alternative = "two.sided"))

  expect_equal(TD_power(-2, sample_size = 15, alternative = "greater"),
               UDT_power(-2, 0, r_ab = 0.5, sample_size = 15, alternative = "greater"))

  expect_equal(TD_power(-2, power = 0.5, alternative = "less"),
               UDT_power(-2, 0, r_ab = 0.5, power = 0.5, alternative = "less"))

  suppressMessages(
    expect_equal(TD_power(-2, power = 0.5, alternative = "two.sided"),
                 UDT_power(-2, 0, r_ab = 0.5, power = 0.5, alternative = "two.sided"))
  )

  suppressMessages(
  expect_equal(TD_power(-2, power = 0.5, alternative = "greater"),
               UDT_power(-2, 0, r_ab = 0.5, power = 0.5, alternative = "greater"))
  )


})


test_that("Output is dataframe when power is given and int when n is given, TD_power", {

  x <- TD_power(-4, power = 0.8)

  expect_equal(class(x), "data.frame")

  x <- TD_power(-4, sample_size = 20)

  expect_equal(class(x), "numeric")

})


test_that("when n is high, produce same power a z-test", {


  expect_equal(TD_power(-2, sample_size = 1000),
               pnorm(qnorm(0.05), mean = -2),
               tolerance = 0.001)


})

test_that("standardised and non-std scores give equal results TD", {

  expect_equal(TD_power(case = 70,
                        mean = 100,
                        sd = 15,
                        sample_size = 15,
                        alternative = "less",
                        alpha = 0.05),
               TD_power(case = -2,
                        mean = 0,
                        sd = 1,
                        sample_size = 15,
                        alternative = "less",
                        alpha = 0.05),
               tolerance = 0.0001)


})

###########################################

# BTD_power

###########################################

test_that("alt functionality works", {

  expect_error(
    BTD_power(-2, sample_size =  10, alternative = "l", iter = 20, nsim = 20),
    NA)
  expect_error(
    BTD_power(-2, sample_size =  10, alternative = "g", iter = 20, nsim = 20),
    NA)
  expect_error(
    BTD_power(-2, sample_size =  10, alternative = "t", iter = 20, nsim = 20),
    NA)

})

test_that("TD_power and BTD_power produce approx equal output", {

  set.seed(12345)
  expect_equal(TD_power(-2, sample_size = 15, alternative = "less"),
               BTD_power(-2, sample_size = 15, alternative = "less"), tolerance = 0.05)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "two.sided"),
               BTD_power(-2, sample_size = 15, alternative = "two.sided"), tolerance = 0.05)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "greater"),
               BTD_power(-2, sample_size = 15, alternative = "greater"), tolerance = 0.05)



})

test_that("errors are working for BTD_power", {

  expect_error(BTD_power(-2, sample_size = 1),
               "Sample size must be greater than 1")
  expect_error(BTD_power(-2, sample_size = 10, alpha = 2),
               "Type I error rate must be between 0 and 1")

})

test_that("standardised and non-std scores give equal results BTW_pow", {

  expect_equal(BTD_power(case = 70,
                        mean = 100,
                        sd = 15,
                        sample_size = 15,
                        alternative = "less",
                        alpha = 0.05,
                        iter = 500,
                        nsim = 500),
               BTD_power(case = -2,
                        mean = 0,
                        sd = 1,
                        sample_size = 15,
                        alternative = "less",
                        alpha = 0.05,
                        iter = 500,
                        nsim = 500),
               tolerance = 0.1)


})


###########################################

# BTD_cov_power

###########################################

test_that("alt functionality works", {

  expect_error(
    BTD_cov_power(-2, 0, sample_size =  10, alternative = "l", iter = 20, nsim = 20),
    NA)
  expect_error(
    BTD_cov_power(-2, 0, sample_size =  10, alternative = "g", iter = 20, nsim = 20),
    NA)
  expect_error(
    BTD_cov_power(-2, 0, sample_size =  10, alternative = "t", iter = 20, nsim = 20),
    NA)

})

test_that("TD_power and BTD_cov_power produce approx equal output", {

  set.seed(123456)
  expect_equal(TD_power(-2, sample_size = 15, alternative = "less"),
               BTD_cov_power(-2, 0, sample_size = 15, alternative = "less",
                             cor_mat = diag(2),
                             iter=50, nsim = 50), tolerance = 0.05)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "two.sided"),
               BTD_cov_power(-2, 0, sample_size = 15, alternative = "two.sided",
                             cor_mat = diag(2),
                             iter=50, nsim = 50), tolerance = 0.05)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "greater"),
               BTD_cov_power(-2, 0, sample_size = 15, alternative = "greater",
                             cor_mat = diag(2),
                             iter=50, nsim = 50), tolerance = 0.05)



})

test_that("errors are working for BTD_cov_power", {

  expect_error(BTD_cov_power(-2, sample_size = 1),
               "Sample size must be greater than 1")
  expect_error(BTD_cov_power(-2, sample_size = 10, alpha = 2),
               "Type I error rate must be between 0 and 1")
  expect_error(BTD_cov_power(-2, case_cov = c(0, 0, 0), control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
                             sample_size = 10, cor_mat = diag(4)- 0.7 + diag(0.7, nrow = 4), iter = 50, nsim = 50),
               "cor_mat is not positive definite")
  expect_error(BTD_cov_power(-2, case_cov = c(0, 0, 0), control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
                             sample_size = 10, cor_mat = diag(3), iter = 50, nsim = 50),
               "Number of variables and number of correlations do not match")


})

test_that("standardised and non-std scores give equal results BTD_pow_cov", {

  expect_equal(BTD_cov_power(case = 70,
                             case_cov = 70,
                             control_task = c(100, 15),
                             control_covar = c(100, 15),
                             cor_mat = diag(2),
                             sample_size = 15,
                             alternative = "less",
                             alpha = 0.05,
                             nsim = 50,
                             iter = 50),
               BTD_cov_power(case = -2,
                             case_cov = -2,
                             control_task = c(0, 1),
                             control_covar = c(0, 1),
                             cor_mat = diag(2),
                             sample_size = 15,
                             alternative = c("less", "greater", "two.sided"),
                             alpha = 0.05,
                             nsim = 50,
                             iter = 50),
               tolerance = 0.2)
})

###########################################

# UDT_power

###########################################


test_that("alt functionality works", {

  expect_error(
    UDT_power(-2, -1, sample_size =  20, alternative = "l"),
    NA)
  expect_error(
    UDT_power(-2, -1, sample_size =  20, alternative = "g"),
    NA)
  expect_error(
    UDT_power(-2, -1, sample_size =  20, alternative = "t"),
    NA)

})

test_that("Output is dataframe when power is given and int when n is given, UDT_power", {

  x <- UDT_power(-5, -1, power = 0.8)

  expect_equal(class(x), "data.frame")

  x <- UDT_power(-4, -1, sample_size = 20)

  expect_equal(class(x), "numeric")

})

test_that("output is approx the same as simulations, TD_power", {

  sim <- c(0.29013, 0.39963, 0.45752, 0.48926, 0.50659)

  x <- c(UDT_power(-2, 0, sample_size = 5),
         UDT_power(-2, 0, sample_size = 10),
         UDT_power(-2, 0, sample_size = 20),
         UDT_power(-2, 0, sample_size = 50),
         UDT_power(-2, 0, sample_size = 100))

  expect_equal(x, sim, tolerance = 1e-2)

})

test_that("errors are working for UDT_power", {
  expect_error(UDT_power(-2, 0),
               "Must supply either sample size or desired power")
  expect_error(UDT_power(-2, 0, sample_size = 20, power = 0.8),
               "Must supply only one of sample size or desired power")
  expect_error(UDT_power(-2, 0, power = 1.2),
               "Desired power must be between 0 and 1")
  expect_error(UDT_power(-2, 0, power = -1),
               "Desired power must be between 0 and 1")
  expect_error(UDT_power(-2, 0, sample_size = 1),
               "Sample size must be greater than 1")
  expect_error(UDT_power(-2, 0, sample_size = 20, alpha = 1.2),
               "Type I error rate must be between 0 and 1")
  expect_error(UDT_power(-2, 0,  sample_size = 20, alpha = -1),
               "Type I error rate must be between 0 and 1")
  expect_error(UDT_power(-2, 0,  sample_size = 20, r_ab = 1.2),
               "Correlation between task a and b must be between -1 and 1")


})


###########################################

# RSDT_power

###########################################



test_that("alt functionality works", {

  expect_error(
    RSDT_power(-2, -1, sample_size =  20, alternative = "l", nsim = 20),
    NA)
  expect_error(
    RSDT_power(-2, -1, sample_size =  20, alternative = "g", nsim = 20),
    NA)
  expect_error(
    RSDT_power(-2, -1, sample_size =  20, alternative = "t", nsim = 20),
    NA)

})

test_that("errors are working for RSDT_power", {

  expect_error(RSDT_power(-2, 0, sample_size = 1),
               "Sample size must be greater than 1")
  expect_error(RSDT_power(-2, 0, sample_size = 20, alpha = 1.2),
               "Type I error rate must be between 0 and 1")
  expect_error(RSDT_power(-2, 0,  sample_size = 20, alpha = -1),
               "Type I error rate must be between 0 and 1")
  expect_error(RSDT_power(-2, 0,  sample_size = 20, r_ab = 1.2),
               "Correlation between task a and b must be between -1 and 1")

})


test_that("TD_power and RSDR_power power estimations in the same range when r = 0.5", {

  set.seed(12344)
  expect_equal(TD_power(-2, sample_size = 15, alternative = "less"),
               RSDT_power(-2, 0, r_ab = 0.5, sample_size = 15,
                          alternative = "less", nsim = 100), tolerance = 0.06)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "two.sided"),
               RSDT_power(-2, 0, r_ab = 0.5, sample_size = 15,
                          alternative = "two.sided", nsim = 100), tolerance = 0.06)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "greater"),
               RSDT_power(-2, 0, r_ab = 0.5, sample_size = 15,
                          alternative = "greater", nsim = 100), tolerance = 0.06)



})

test_that("standardised and non-std scores give equal results RSDT_pow", {

  expect_equal(RSDT_power(case_a = 70,
                          case_b = 100,
                          mean_a = 100,
                          mean_b = 100,
                          sd_a = 15,
                          sd_b = 15,
                          r_ab = 0.5,
                          sample_size = 15,
                          alternative = "t",
                          alpha = 0.05,
                          nsim = 1000),
               RSDT_power(case_a = -2,
                          case_b = 0,
                          mean_a = 0,
                          mean_b = 0,
                          sd_a = 1,
                          sd_b = 1,
                          r_ab = 0.5,
                          sample_size = 15,
                          alternative = "t",
                          alpha = 0.05,
                          nsim = 1000),
               tolerance = 0.1)


})

###########################################

# BSDT_power

###########################################

test_that("alt functionality works", {

  expect_error(
    BSDT_power(-2, -1, sample_size =  20, alternative = "l", nsim = 20, iter = 20),
    NA)
  expect_error(
    BSDT_power(-2, -1, sample_size =  20, alternative = "g", nsim = 20, iter = 20),
    NA)
  expect_error(
    BSDT_power(-2, -1, sample_size =  20, alternative = "t", nsim = 20, iter = 20),
    NA)

})

test_that("choice of prior functionality works", {

  expect_error(
    BSDT_power(-2, -1, sample_size =  20,
               calibrated = TRUE, nsim = 20, iter = 20),
    NA)
  expect_error(
    BSDT_power(-2, -1, sample_size =  20,
               calibrated = FALSE, nsim = 20, iter = 20),
    NA)

})

test_that("errors are working for BSDT_power", {

  expect_error(BSDT_power(-2, 0, sample_size = 1),
               "Sample size must be greater than 1")
  expect_error(BSDT_power(-2, 0, sample_size = 20, alpha = 1.2),
               "Type I error rate must be between 0 and 1")
  expect_error(BSDT_power(-2, 0,  sample_size = 20, alpha = -1),
               "Type I error rate must be between 0 and 1")
  expect_error(BSDT_power(-2, 0,  sample_size = 20, r_ab = 1.2),
               "Correlation between task a and b must be between -1 and 1")

})


test_that("TD_power and BSDT_power power estimations in the same range when r = 0.5", {

  set.seed(123123)
  expect_equal(TD_power(-2, sample_size = 15, alternative = "less"),
               BSDT_power(-2, 0, r_ab = 0.5, sample_size = 15,
                          alternative = "less", nsim = 50, iter = 50), tolerance = 0.06)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "two.sided"),
               BSDT_power(-2, 0, r_ab = 0.5, sample_size = 15,
                          alternative = "two.sided", nsim = 50, iter = 50), tolerance = 0.06)

  expect_equal(TD_power(-2, sample_size = 15, alternative = "greater"),
               BSDT_power(-2, 0, r_ab = 0.5, sample_size = 15,
                          alternative = "greater", nsim = 50, iter = 50), tolerance = 0.06)



})

test_that("standardised and non-std scores give equal results BSDT_pow", {

  expect_equal(BSDT_power(case_a = 70,
                          case_b = 100,
                          mean_a = 100,
                          mean_b = 100,
                          sd_a = 15,
                          sd_b = 15,
                          r_ab = 0.5,
                          sample_size = 15,
                          alternative = "t",
                          alpha = 0.05,
                          nsim = 50),
               BSDT_power(case_a = -2,
                          case_b = 0,
                          mean_a = 0,
                          mean_b = 0,
                          sd_a = 1,
                          sd_b = 1,
                          r_ab = 0.5,
                          sample_size = 15,
                          alternative = "t",
                          alpha = 0.05,
                          nsim = 50),
               tolerance = 0.2)
})

###########################################

# BSDT_cov_power

###########################################

test_that("alt functionality works", {

  expect_error(
    BSDT_cov_power(c(-2, -2), -1, sample_size =  20, alternative = "l", nsim = 20, iter = 20),
    NA)
  expect_error(
    BSDT_cov_power(c(-2, -2), -1, sample_size =  20, alternative = "g", nsim = 20, iter = 20),
    NA)
  expect_error(
    BSDT_cov_power(c(-2, -2), -1, sample_size =  20, alternative = "t", nsim = 20, iter = 20),
    NA)
})

test_that("choice of prior functionality works", {

  expect_error(
    BSDT_cov_power(c(-2, -2), -1, sample_size =  20,
                    calibrated = TRUE, nsim = 20, iter = 20),
    NA)
  expect_error(
    BSDT_cov_power(c(-2, -2), -1, sample_size =  20,
                    calibrated = FALSE, nsim = 20, iter = 20),
    NA)

})


test_that("errors are working for BTD_cov_power", {

  expect_error(BSDT_cov_power(c(-2, 0), 0, sample_size = 1),
               "Sample size must be greater than 1")
  expect_error(BSDT_cov_power(c(-2, 0), 0, sample_size = 10, alpha = 2),
               "Type I error rate must be between 0 and 1")
  expect_error(BSDT_cov_power(c(-2, 0), case_cov = c(0, 0, 0), control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
                             sample_size = 10, cor_mat = diag(4)- 0.7 + diag(0.7, nrow = 4), iter = 50, nsim = 50),
               "cor_mat is not positive definite")
  expect_error(BSDT_cov_power(c(-2, 0), case_cov = c(0, 0, 0), control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
                             sample_size = 10, cor_mat = diag(4), iter = 50, nsim = 50),
               "Number of variables and number of correlations do not match")
  expect_error(BSDT_cov_power(c(-2, 0, 0), case_cov = c(0, 0, 0), control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
                              sample_size = 10, cor_mat = diag(4), iter = 50, nsim = 50),
               "case_tasks should be of length 2")

})

test_that("standardised and non-std scores give equal results BSDT_cov_pow", {

  expect_equal(
    BSDT_cov_power(case_tasks = c(70, 100),
                   case_cov = 100,
                   control_tasks = matrix(c(100, 100, 15, 15), ncol = 2),
                   control_covar = c(100, 15),
                   cor_mat = diag(3),
                   sample_size = 15,
                   alternative = c("two.sided", "greater", "less"),
                   alpha = 0.05,
                   nsim = 100,
                   iter = 100,
                   calibrated = TRUE),
    BSDT_cov_power(case_tasks = c(-2, 0),
                   case_cov = 0,
                   control_tasks = matrix(c(0, 0, 1, 1), ncol = 2),
                   control_covar = c(0, 1),
                   cor_mat = diag(3),
                   sample_size = 15,
                   alternative = c("two.sided", "greater", "less"),
                   alpha = 0.05,
                   nsim = 100,
                   iter = 100,
                   calibrated = TRUE),
    tolerance = 0.2)
})
