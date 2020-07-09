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
                   sc_c[["interval"]][["Lower zdccc CI"]],
                   sc_c[["interval"]][["Upper zdccc CI"]],
                   sc_c[["interval"]][["Lower p CI"]],
                   sc_c[["interval"]][["Upper p CI"]])

  expect_equal(sc_calib_tt, calib_prior_tt, tolerance = 1e-2)


  sc_j <- BSDT_cov(c(70, 78), 30, x[ , 1:2], x[ , 3],
                         calibrated = F, iter = 10000)

  sc_jeff_tt <- c(sc_j[["p.value"]],
                  sc_j[["interval"]][["Lower zdccc CI"]],
                  sc_j[["interval"]][["Upper zdccc CI"]],
                  sc_j[["interval"]][["Lower p CI"]],
                  sc_j[["interval"]][["Upper p CI"]])

  expect_equal(sc_jeff_tt, jeff_prior_tt, tolerance = 1e-2)


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

