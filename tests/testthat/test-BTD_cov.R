test_that("summary input yields same result as raw", {
  x <- MASS::mvrnorm(18, mu = c(100, 13),
                     Sigma = matrix(c(15^2, 0.65*15*3,
                                      0.65*15*3, 3^2),
                                    nrow = 2, byrow = T),
                     empirical = TRUE)
  set.seed(123456)
  sumstats <- BTD_cov(-2, 0, c(100, 15), c(13, 3), use_sumstats = TRUE,
                      cor_mat = matrix(c(1, 0.65, 0.65, 1), nrow=2),
                      sample_size = 18)[["p.value"]]
  set.seed(123456)
  raw <- BTD_cov(-2, 0, x[ , 1], x[ , 2])[["p.value"]]
  expect_equal(sumstats, raw, tol = 0.01)


})


test_that("we get approx same results as C&G on BTD_cov", {
  x <- MASS::mvrnorm(18, mu = c(100, 13),
                     Sigma = matrix(c(15^2, 0.65*15*3,
                                      0.65*15*3, 3^2),
                                    nrow = 2, byrow = T),
                     empirical = TRUE)


  # p-values and intervals from C&G programs given these values

  cg_ot <- c(0.04362 , -2.653, -1.071, 0.3987, 14.2189)


  set.seed(1234597)
  sc_ot <- BTD_cov(78, 13, x[ , 1], x[ , 2], iter = 10000)
  sc_ot <- c(sc_ot[["p.value"]],
             sc_ot[["interval"]][["Lower Z-CCC CI"]],
             sc_ot[["interval"]][["Upper Z-CCC CI"]],
             sc_ot[["interval"]][["Lower p CI"]],
             sc_ot[["interval"]][["Upper p CI"]])


  expect_equal(sc_ot, cg_ot, tolerance = 1e-2)


})

test_that("alternative hypotheses direction", {
  x <- MASS::mvrnorm(18, mu = c(100, 13),
                     Sigma = matrix(c(15^2, 0.65*15*3,
                                      0.65*15*3, 3^2),
                                    nrow = 2, byrow = T),
                     empirical = TRUE)


  set.seed(123456234)
  pos_z <- BTD_cov(105, 13, x[ , 1], x[ , 2],
                   iter = 1000, alternative = "less")[["p.value"]]
  expect_equal(pos_z > 0.5, TRUE)
  set.seed(123456234)
  pos_z <- BTD_cov(105, 13, x[ , 1], x[ , 2],
                   iter = 1000, alternative = "greater")[["p.value"]]
  expect_equal(pos_z < 0.5, TRUE)

  set.seed(123456234)
  neg_z <- BTD_cov(78, 13, x[ , 1], x[ , 2],
                   iter = 1000, alternative = "less")[["p.value"]]
  expect_equal(neg_z < 0.5, TRUE)
  set.seed(123456234)
  neg_z <- BTD_cov(78, 13, x[ , 1], x[ , 2],
                   iter = 1000, alternative = "greater")[["p.value"]]
  expect_equal(neg_z > 0.5, TRUE)


})

test_that("errors and warnings are occuring as they should for BTD", {

  expect_error(BTD_cov(1, 0, 0, 0, use_sumstats = TRUE, sample_size = NULL),
               "Please supply both correlation matrix and sample size")
  expect_error(BTD_cov(1, 0, 0, 0, use_sumstats = TRUE, sample_size = 20, cor_mat = NULL),
               "Please supply both correlation matrix and sample size")

  expect_error(BTD_cov(-2, 0, rnorm(15), rnorm(15), int_level = 1.1),
               "Interval level must be between 0 and 1")

  expect_error(BTD_cov(c(-2, 0), 0, rnorm(15), rnorm(15)),
               "case_task should be single value")

  expect_error(BTD_cov(-2, 0, c(0, 1), c(0, 1), use_sumstats = TRUE, sample_size = 20,
                       cor_mat = diag(c(-2, -2))),
               "cor_mat is not positive definite")

  expect_error(BTD_cov(-2, 0, c(0, 1), c(0, 1), use_sumstats = FALSE, sample_size = 20,
                       cor_mat = diag(2)),
               "If input is summary data, set use_sumstats = TRUE")

  expect_error(BTD_cov(-2, 0, c(0, 1), c(0, 1), use_sumstats = TRUE, sample_size = 20,
                       cor_mat = diag(3)),
               "Number of variables and number of correlations does not match")

})
