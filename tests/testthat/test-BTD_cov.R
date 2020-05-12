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
             sc_ot[["interval"]][["Lower zccc CI"]],
             sc_ot[["interval"]][["Upper zccc CI"]],
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
