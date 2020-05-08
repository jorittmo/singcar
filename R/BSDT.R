#' Bayesian Standardised Difference Test
#'
#' A test on the difference between two tasks in a single case, by comparison to
#' the difference of the same two tasks in a control sample. Calculates a
#' standardised effects size of task difference as well as a point estimate of
#' the proportion of the control population that would be expected to show a
#' more extreme task difference.
#'
#' @param case.x Case's score on task X.
#' @param case.y Case's score on task Y.
#' @param controls.x Controls' scores on task X. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can
#'   supply a vector as input for task X while mean and SD for task Y.
#' @param controls.y Controls' scores on task Y. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can
#'   supply a vector as input for task Y while mean and SD for task X.
#' @param controls.x.sd If single value for task X is given as input you must
#'   supply the standard deviation of the sample.
#' @param controls.y.sd If single value for task Y is given as input you must
#'   supply the standard deviation of the sample.
#' @param controls.n If X or Y is given as mean and SD you must supply the
#'   sample size. If controls.x is given as vector and controls.y as mean and
#'   SD, controls.n must equal the number of observations in controls.x.
#' @param cor.x.y If X or Y is given as mean and SD you must supply the
#'   correlation between the tasks.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter. Since the direction
#'   of the expected effect depends on which task is set as X and which is set
#'   as Y, be very careful if changing this parameter.
#' @param int.level Level of confidence for credible intervals.
#' @param iter Number of iterations.
#' @param unstandardised Estimate z-value based on standardised or
#'   unstandardised task scores.
#' @param calibrated set to \code{TRUE} to use a calibrated prior distribution.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab if exact.method set to \code{TRUE},
#'   returns the value of an exact t-statistic, however, because of the
#'   underlying equation, it cannot be negative. Set exact.method to
#'   \code{FALSE} for an approximate t-value with correct sign. \cr\cr
#'   \code{parameter} \tab the degrees of freedom for the t-statistic.\cr\cr
#'   \code{p.value}    \tab the p-value for the test.\cr\cr \code{estimate} \tab
#'   case scores expressed as z-scores on task X and Y. Standardised effect size
#'   (Z-DCC) of task difference between case and controls and point estimate of
#'   the proportion of the control population estimated to show a more extreme
#'   task difference. \cr\cr \code{sample.size}   \tab the size of the control
#'   sample\cr\cr \code{null.value}   \tab the value of the difference under the
#'   null hypothesis.\cr\cr  \code{alternative}     \tab a character string
#'   describing the alternative hypothesis.\cr\cr \code{method} \tab a character
#'   string indicating what type of test was performed.\cr\cr \code{data.name}
#'   \tab a character string giving the name(s) of the data}
#'
#' @export
#'
#' @examples
#' BSDT(-3.857, -1.875, controls.x = 0, controls.y = 0, controls.x.sd = 1,
#' controls.y.sd = 1, controls.n = 20, cor.x.y = 0.68)
#' BSDT(-3.857, -1.875, controls.x = rnorm(20), controls.y = rnorm(20))
#'
#' @references
#' Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a
#' control or normative sample in neuropsychology: Development of a Bayesian
#' approach. \emph{Cognitive Neuropsychology, 24}(4), 343–372.
#' \url{https://doi.org/10.1080/02643290701290146}




BSDT <- function (case.x, case.y, controls.x, controls.y,
                  controls.x.sd = NULL, controls.y.sd = NULL,
                  controls.n = NULL, cor.x.y = NULL,
                  alternative = c("two.sided", "greater", "less"),
                  int.level = 0.95,
                  iter = 1000,
                  unstandardised = FALSE,
                  calibrated = FALSE,
                  na.rm = FALSE) {

  alternative <- match.arg(alternative)

  if (length(case.x) > 1 | length(case.y) > 1) stop("Case scores should be single value")
  if (length(controls.x) > 1 & length(controls.y) > 1) {
    if (length(controls.x) != length(controls.y)) stop("Sample sizes must be equal")
  }

  if (length(controls.x) > 1 & length(controls.y) > 1 & is.null(controls.n) == FALSE) message("Value on controls.n will be ignored")

  if (length(controls.x) > 1 & is.null(controls.x.sd) == FALSE) message("Value on controls.x.sd will be ignored")
  if (length(controls.y) > 1 & is.null(controls.y.sd) == FALSE) message("Value on controls.y.sd will be ignored")

  if (length(controls.x) == 1 & is.null(controls.x.sd) == TRUE) stop("Please give sd and n on task x if controls.x is to be treated as mean")
  if (length(controls.y) == 1 & is.null(controls.y.sd) == TRUE) stop("Please give sd and n on task y if controls.y is to be treated as mean")


  # Handling of NA use cases below
  if(is.na(case.x) == TRUE | is.na(case.y) == TRUE) stop("One or both case scores is NA")

  if (na.rm == TRUE) {
    if (sum(is.na(controls.x))  > 0 & sum(is.na(controls.y)) == 0 ) {
      controls.y <- controls.y[!is.na(controls.x)]
      controls.x <- controls.x[!is.na(controls.x)]
      warning("Removal of NAs on controls.x resulted in removal of non-NAs on controls.y")
    }

    if (sum(is.na(controls.y))  > 0 & sum(is.na(controls.x)) == 0 ) {
      controls.x <- controls.x[!is.na(controls.y)]
      controls.y <- controls.y[!is.na(controls.y)]
      warning("Removal of NAs on controls.y resulted in removal of non-NAs on controls.x")
    }

    if (sum(is.na(controls.y))  > 0 & sum(is.na(controls.x)) > 0 ) {

      if (identical(!is.na(controls.x), !is.na(controls.y)) == TRUE) {
        controls.x <- controls.x[!is.na(controls.x)]
        controls.y <- controls.y[!is.na(controls.y)]
      } else {
        conx <- controls.x[!is.na(controls.x) & !is.na(controls.y)]
        cony <- controls.y[!is.na(controls.x) & !is.na(controls.y)]

        controls.x <- conx
        controls.y <- cony

        warning("Removal of NAs on one control sample resulted in removal of non-NAs on the other")
      }

    }

  }
  if (sum(is.na(controls.x)) > 0 | sum(is.na(controls.y)) > 0) stop("Controls contains NA, set na.rm = TRUE to proceed")
  # End of NA use cases


  if (length(controls.x) > 1 & length(controls.y) > 1) {
    if (length(controls.x) != length(controls.y)) stop("Sample sizes must be equal")
  }

  con_m.x <- mean(controls.x) # Mean of the control sample on task x
  con_m.y <- mean(controls.y) # Mean of the control sample on task y

  con_sd.x <- stats::sd(controls.x) # Standard deviation of the control sample on task x
  if (length(controls.x) == 1 & is.null(controls.x.sd) == FALSE) con_sd.x <- controls.x.sd

  con_sd.y <- stats::sd(controls.y) # Standard deviation of the control sample on task y
  if (length(controls.y) == 1 & is.null(controls.y.sd) == FALSE) con_sd.y <- controls.y.sd


  # Since controls x and y need to be of equal length n is the length of any of them
  n <- length(controls.x)
  if (length(controls.x) == 1 | length(controls.y) == 1) {
    if (is.null(controls.n) == TRUE) stop("Please set sample size")
    n <- controls.n
    if (length(controls.x) > 1 & n != length(controls.x)) stop("Sample sizes must be equal")
    if (length(controls.y) > 1 & n != length(controls.y)) stop("Sample sizes must be equal")
  }

  if (is.null(cor.x.y) == TRUE & length(controls.x) == 1) stop("Please set correlation between tasks")
  if (is.null(cor.x.y) == TRUE & length(controls.y) == 1) stop("Please set correlation between tasks")

  if (is.null(cor.x.y) == FALSE){
    if (cor.x.y < -1 | cor.x.y > 1) stop("Correlation must be between -1 and 1")
  }

  r <- cor.x.y
  if (length(controls.x) > 1 & length(controls.y) > 1) r <- stats::cor(controls.x, controls.y)


  if (length(controls.x) > 1 & length(controls.y) > 1) {

    con_mat <- cbind(controls.x, controls.y)
    # Calculate SSCP matrix and call it A as in Crawford Garthwaite-notation - the scale matrix

    A <- (n - 1)* cov(con_mat)

  } else {

    sxx <- con_sd.x^2 * (n - 1)
    syy <- con_sd.y^2 * (n - 1)

    sxy <- con_sd.x*con_sd.y* r * (n - 1)

    A <- matrix(c(sxx, sxy, sxy, syy), nrow = 2)
  }

  if (calibrated == FALSE) {

    seed <- sample(1:10)

    set.seed(seed) # So that both the inverse wishart draws and the cholesky decomp on them are the same
    Sigma_hat <- CholWishart::rInvWishart(iter, n, A)

    set.seed(seed)
    Tchol <- CholWishart::rInvCholWishart(iter, n, A) # Simulates same as above but with cholesky decomp in C++
    Tchol <- aperm(Tchol, perm = c(2, 1, 3)) # Transposes each matrix to lower triangual instead of upper

    Mu_hat <- matrix(nrow = iter, ncol = 2)
    for (i in 1:iter) Mu_hat[i , ] <-  as.numeric(c(con_m.x, con_m.y) + (Tchol[ , , i]%*%stats::rnorm(2))/sqrt(n))


    ## For those that thinks apply() gives more readability, but it slows the code
    # Mu_hat <- t(apply(tc, 2, function(x) as.numeric((c(con_m.x, con_m.y) +
    #                                                   matrix(x, nrow = 2)%*%stats::rnorm(2))/sqrt(n))))

  } else { # if calibrated == TRUE

    A_ast <- ((n - 2)*A) / (n - 1)

    step_it <- iter
    Sigma_hat_acc_save <- array(dim = c(2, 2, 1))

    while(dim(Sigma_hat_acc_save)[3] < iter + 1) {

      # Sigma_hat <- CholWishart::rInvWishart(step_it, df = n - m - 2, solve(A_ast)) Är det verkligen A invers??
      # Invers enligt pappret, men får helt galna värden. Däremot med vanlig kryssprodukt blir det värden vi hade väntat oss
      Sigma_hat <- CholWishart::rInvWishart(step_it, df = n - 2, A_ast)

      rho_hat_pass <- Sigma_hat[1, 2, ] / sqrt(Sigma_hat[1, 1, ] * Sigma_hat[2, 2, ])

      u <- runif(step_it, min = 0, max = 1)

      Sigma_hat_acc <- Sigma_hat[ , , (u^2 <= (1 - rho_hat_pass^2))]

      Sigma_hat_acc_save <- abind::abind(Sigma_hat_acc_save, Sigma_hat_acc, along = 3)

      # Sigma_hat_acc_save <- array(c(Sigma_hat_acc_save, Sigma_hat_acc), # Bind the arrays together
      #                             dim = c(2, 2, (dim(Sigma_hat_acc_save)[3] + dim(Sigma_hat_acc)[3])))

      step_it <- iter - dim(Sigma_hat_acc_save)[3] + 1

    }

    Sigma_hat <- Sigma_hat_acc_save[ , , -1] # Remove the first matrix that is fild with NA
    rm(Sigma_hat_acc_save, step_it, u, Sigma_hat_acc, rho_hat_pass) # Remove all variables not needed

    Tchol <- array(dim = c(2, 2, iter))
    for (i in 1:iter) Tchol[ , , i] <- t(chol(Sigma_hat[ , , i]))

    Mu_hat <- matrix(nrow = iter, ncol = 2)
    for (i in 1:iter) Mu_hat[i , ] <-  as.numeric(c(con_m.x, con_m.y) + (Tchol[ , , i]%*%stats::rnorm(2))/sqrt(n))

  }



  if (unstandardised == FALSE) {

    zx <- (case.x - Mu_hat[ , 1]) / sqrt(Sigma_hat[1, 1, ])
    zy <- (case.y - Mu_hat[ , 2]) / sqrt(Sigma_hat[2, 2, ])

    rho_hat <- Sigma_hat[1, 2, ] / sqrt(Sigma_hat[1, 1, ] * Sigma_hat[2, 2, ])

    z_ast <- (zx - zy) / sqrt(2 - 2*rho_hat)

  } else {

    std.err <- sqrt(Sigma_hat[1, 1, ] + Sigma_hat[2, 2, ] - 2*Sigma_hat[1, 2, ])

    z_ast <- ((case.x - Mu_hat[ , 1]) - (case.y - Mu_hat[ , 2])) / std.err

  }


  if (alternative == "two.sided") {
    pval <- 2 * stats::pnorm(abs(z_ast), lower.tail = FALSE)
  } else if (alternative == "greater") {
    pval <- stats::pnorm(z_ast, lower.tail = FALSE)
  } else { # I.e. if alternative == "less"
    pval <- stats::pnorm(z_ast, lower.tail = TRUE)
  }


  alpha <- 1 - int.level

  z_ast_est <- mean(z_ast)
  names(z_ast_est) <- "est. z"

  zdcc_int <- stats::quantile(z_ast, c(alpha/2, (1 - alpha/2)))
  names(zdcc_int) <- c("Lower zdcc CI", "Upper zdcc CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100
  names(p_int) <- c("Lower p CI", "Upper p CI")

  std.x <- (case.x - con_m.x)/con_sd.x
  std.y <- (case.y - con_m.y)/con_sd.y

  zdcc <- (std.x - std.y) / sqrt(2 - 2*r) # Estimated effect size

  estimate <- c(std.x, std.y, zdcc, ifelse(alternative == "two.sided", (p_est/2*100), p_est*100))

  if (alternative == "two.sided") {
    alt.p.name <- "Proportion of control population with more extreme task difference, "
  } else if (alternative == "greater") {
    alt.p.name <- "Proportion of control population with more positive task difference, "
  } else {
    alt.p.name <- "Proportion of control population with more negative task difference, "
  }

  p.name <- paste0(alt.p.name,
                   100*int.level, "% credible interval [",
                   format(round(p_int[1], 2), nsmall = 2),", ",
                   format(round(p_int[2], 2), nsmall = 2),"]")

  zdcc.name <- paste0("Std. effect size (Z-DCC) for task diff. between case and controls, ",
                     100*int.level, "% credible interval [",
                     format(round(zdcc_int[1], 2), nsmall = 2),", ",
                     format(round(zdcc_int[2], 2), nsmall = 2),"]")

  names(estimate) <- c("Case score on task X as standard (z) score",
                       "Case score on task Y as standard (z) score",
                       zdcc.name,
                       p.name)


  typ.int <- 100*int.level
  names(typ.int) <- "Interval level (%)"
  interval <- c(typ.int, zdcc_int, p_int)


  # names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  names(con_m.x) <- "Mean (controls)"
  names(con_sd.x) <- "SD (controls)"
  names(con_m.y) <- "Mean (controls)"
  names(con_sd.y) <- "SD (controls)"
  names(n) <- "Sample size"
  dname <- paste0("Case score X: ", deparse(substitute(case.x)), ", ",
                  "Case score Y: ", deparse(substitute(case.y)), ", ",
                  "Controls score X: ", deparse(substitute(controls.x)), ", ",
                  "Controls score Y: ",deparse(substitute(controls.y)))

  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = z_ast_est,
                 #parameter = df,
                 p.value = p_est,
                 estimate = estimate,
                 null.value = null.value,
                 interval = interval,
                 desc = c(con_m.x, con_sd.x, con_m.y, con_sd.y, n),
                 alternative = alternative,
                 method = paste("Bayesian Standardised Difference Test"),
                 data.name = dname)

  class(output) <- "htest"
  output


}

