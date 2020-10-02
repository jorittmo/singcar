test_that("we get approx same results as C&G on BSDT_cov", {

  set.seed(123456)
  x <- MASS::mvrnorm(20, mu = c(100, 80, 50),
                     Sigma = matrix(c(15^2, 108, 78,
                                      108, 10^2, 12,
                                      78, 12, 10^2),
                                    nrow = 3, byrow = T),
                     empirical = TRUE)

  # p-values and intervals from C&G programs given these values

  calib_prior_tt <- c(0.16247, -2.8085, -0.465, 0.251, 32.114) # Values from C&G software
  jeff_prior_tt <- c(0.14009, -2.882, -0.559, 0.197, 28.816) # Values from C&G software


  sc_c <- BSDT_cov(c(70, 78), 30, x[ , 1:2], x[ , 3],
                          calibrated = T, iter = 10000)
  sc_calib_tt <- c(sc_c[["p.value"]],
                   sc_c[["interval"]][["Lower Z-DCCC CI"]],
                   sc_c[["interval"]][["Upper Z-DCCC CI"]],
                   sc_c[["interval"]][["Lower p CI"]],
                   sc_c[["interval"]][["Upper p CI"]])

  expect_equal(sc_calib_tt, calib_prior_tt, tolerance = 0.05)


  sc_j <- BSDT_cov(c(70, 78), 30, x[ , 1:2], x[ , 3],
                         calibrated = F, iter = 10000)

  sc_jeff_tt <- c(sc_j[["p.value"]],
                  sc_j[["interval"]][["Lower Z-DCCC CI"]],
                  sc_j[["interval"]][["Upper Z-DCCC CI"]],
                  sc_j[["interval"]][["Lower p CI"]],
                  sc_j[["interval"]][["Upper p CI"]])

  expect_equal(sc_jeff_tt, jeff_prior_tt, tolerance = 0.05)


})

test_that("alternative hypotheses direction", {

  set.seed(123456)
  x <- MASS::mvrnorm(20, mu = c(100, 80, 50),
                     Sigma = matrix(c(15^2, 108, 78,
                                      108, 10^2, 12,
                                      78, 12, 10^2),
                                    nrow = 3, byrow = T),
                     empirical = TRUE)



  pos_z <- BSDT_cov(c(78, 70), 30, x[ , 2:1], x[ , 3],
                   calibrated = T, iter = 1000, alternative = "less")[["p.value"]]
  expect_equal(pos_z > 0.5, TRUE)

  pos_z <- BSDT_cov(c(78, 70), 30, x[ , 2:1], x[ , 3],
                    calibrated = T, iter = 1000, alternative = "greater")[["p.value"]]
  expect_equal(pos_z < 0.5, TRUE)




  neg_z <- BSDT_cov(c(70, 78), 30, x[ , 1:2], x[ , 3],
                    calibrated = T, iter = 1000, alternative = "less")[["p.value"]]
  expect_equal(neg_z < 0.5, TRUE)

  neg_z <- BSDT_cov(c(70, 78), 30, x[ , 1:2], x[ , 3],
                    calibrated = T, iter = 1000, alternative = "greater")[["p.value"]]
  expect_equal(neg_z > 0.5, TRUE)



})

test_that("errors and warnings are occuring as they should for BSDT", {

  na_con <- rnorm(15)
  na_con[1] <- NA

  expect_error(BSDT_cov(c(-2, 0), 0, cbind(rnorm(15), na_con), rnorm(15)),
               "control_tasks contains NA")

  expect_error(BSDT_cov(c(-2, 0), 0, cbind(rnorm(15), rnorm(15)), na_con),
               "control_covar contains NA")

  expect_error(BSDT_cov(case_tasks = c(-2, 0), case_covar = 0,
                        control_tasks = matrix(c(0, 0, 1, 1), ncol = 2),
                        control_covar = matrix(c(0, 1), nrow = 1), sample_size = 20, cor_mat = diag(3)),
               "If input is summary data, set use_sumstats = TRUE")

  expect_error(BSDT_cov(case_tasks = c(-2, 0), case_covar = 0,
                        control_tasks = matrix(c(0, 0, 1, 1), ncol = 2),
                        control_covar = matrix(c(0, 1), nrow = 1), use_sumstats = TRUE),
               "Please supply both correlation matrix and sample size")

  expect_error(BSDT_cov(c(-2, 0, 0), 0, cbind(rnorm(15), rnorm(15)), rnorm(15)),
               "case_task should have length 2")

  expect_error(BSDT_cov(c(-2, 0), 0, cbind(rnorm(15), rnorm(15), rnorm(15)), rnorm(15)),
               "columns in control_tasks should be 2")

  expect_error(BSDT_cov(c(-2, 0), 0, cbind(rnorm(15), rnorm(15)), rnorm(16)),
               "Must supply equal number of observations for tasks and covariates")
  expect_error(BSDT_cov(c(-2, 0), 0, cbind(rnorm(15), rnorm(15)), cbind(rnorm(16), rnorm(16))),
               "Must supply equal number of observations for tasks and covariates")



})


test_that("summary data and raw gives equal results", {


  sigm <- matrix(c(1, 0.5, 0.1, 0.5, 1, 0.1, 0.1, 0.1, 1), ncol =3)

  con <- MASS::mvrnorm(15, c(0, 0, 0), Sigma = sigm, empirical = TRUE)

  set.seed(123)
  expect_equal(BSDT_cov(c(-2, 0), 0, con[ , 1:2], con[ , 3])[["p.value"]],
               BSDT_cov(case_tasks = c(-2, 0), case_covar = 0,
                        control_tasks = matrix(c(0, 0, 1, 1), ncol = 2),
                        control_covar = matrix(c(0, 1), nrow = 1),
                        sample_size = 15, cor_mat = sigm, use_sumstats = TRUE)[["p.value"]], tolerance = 0.01)



})

