#' Bayesian Standardised Difference Test with Covariates
#'
#' @param case_tasks A vector of length 2. The case scores from the two tasks.
#' @param case_covar A vector containing the case scores on all covariates
#'   included.
#' @param control_tasks A matrix or dataframe with 2 columns and n rows
#'   containing the control scores for the two tasks.
#' @param control_covar A matrix or dataframe cointaining the control scores on
#'   the covariates included.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter. Since the direction
#'   of the expected effect depends on which task is set as X and which is set
#'   as Y, be very careful if changing this parameter.
#' @param int.level The probability level on the Bayesian credible intervals.
#' @param calibrated Whether or not to use the standard theory (Jeffreys) prior
#'   distribution (if set to \code{FALSE}) or a calibrated prior examined by
#'   Berger and Sun (2008) and sample size treated as n - 1. See Crawford et al.
#'   (2011) for further information. Calibrated prior is recommended.
#' @param iter Number of iterations to be performed. Greater number gives better
#'   estimation but takes longer to calculate.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab the average z-value over
#'   \code{iter} number of iterations. \cr\cr \code{p.value}    \tab the average
#'   p-value over \code{iter} number of iterations. \cr\cr \code{estimate} \tab
#'   case scores expressed as z-scores on task X and Y. Standardised effect size
#'   (Z-DCCC) of task difference between case and controls and point estimate of
#'   the proportion of the control population estimated to show a more extreme
#'   task difference. \cr\cr  \code{null.value}   \tab the value of the
#'   difference between tasks under the null hypothesis.\cr\cr \code{interval}
#'   \tab named numerical vector containing level of confidence and confidence
#'   intervals for both effect size and p-value.\cr\cr \code{desc}     \tab data
#'   frame containing means and standard deviations for controls as well as case
#'   scores. \cr\cr \code{cor.mat} \tab matrix giving the correlations between
#'   the tasks of interest and the covariates included. \cr\cr
#'   \code{sample.size}     \tab number of controls. Remember, if
#'   \code{calibrated} left as default (i.e. \code{TRUE}) sample size is treated
#'   as n - 1.\cr\cr \code{alternative}     \tab a character string describing
#'   the alternative hypothesis.\cr\cr \code{method} \tab a character string
#'   indicating what type of test was performed.\cr\cr \code{data.name} \tab a
#'   character string giving the name(s) of the data}
#' @export
#'
#' @examples
#'   controls <- MASS::mvrnorm(20, mu = c(100, 80, 50),
#'    Sigma = matrix(c(15^2, 108, 78, 108, 10^2, 12,78, 12, 10^2),
#'    nrow = 3, byrow = TRUE), empirical = TRUE)
#'
#'    BSDT_cov(c(70, 78), 30, controls[ , 1:2], controls[ , 3], iter = 1000)
#'
#'
BSDT_cov <- function (case_tasks, case_covar, control_tasks, control_covar,
                      alternative = c("two.sided", "greater", "less"),
                      int.level = 0.95,
                      calibrated = TRUE,
                      iter = 1000) {

  alternative <- match.arg(alternative)

  n <- nrow(control_tasks)

  m <- length(case_covar)

  k <- length(case_tasks)

  ## DATA ESTIMATION ##

  m_ct <- colMeans(control_tasks)
  sd_ct <- apply(control_tasks, 2, stats::sd)
  r <- stats::cor(control_tasks)[1, 2]

  X <- cbind(rep(1, n), control_covar)
  Y <- control_tasks

  B_ast <- solve(t(X) %*% X) %*% t(X) %*% Y

  Sigma_ast <- t(Y - X %*% B_ast) %*% (Y - X %*% B_ast)

  if (calibrated == TRUE) {

    ## PRIOR ##

    ## ACCEPT/REJECT STEP Berger and Sun (2008)

    A_ast <- ((n - m - 2)*Sigma_ast) / (n - m - 1)

    step_it <- iter
    Sigma_hat_acc_save <- array(dim = c(2, 2, 1))

    while(dim(Sigma_hat_acc_save)[3] < iter + 1) {

      Sigma_hat <- CholWishart::rInvWishart(step_it, df = (n - m - 2), A_ast)

      rho_hat_pass <- Sigma_hat[1, 2, ] / sqrt(Sigma_hat[1, 1, ] * Sigma_hat[2, 2, ])

      u <- stats::runif(step_it, min = 0, max = 1)

      Sigma_hat_acc <- array(Sigma_hat[ , , (u^2 <= (1 - rho_hat_pass^2))],
                             dim = c(2, 2, sum(u^2 <= (1 - rho_hat_pass^2))))

    #  Sigma_hat_acc_save <- abind::abind(Sigma_hat_acc_save, Sigma_hat_acc, along = 3)

      Sigma_hat_acc_save <- array(c(Sigma_hat_acc_save, Sigma_hat_acc), # Bind the arrays together
                                   dim = c(2, 2, (dim(Sigma_hat_acc_save)[3] + dim(Sigma_hat_acc)[3])))

      step_it <- iter - dim(Sigma_hat_acc_save)[3] + 1

    }

    Sigma_hat <- Sigma_hat_acc_save[ , , -1] # Remove the first matrix that is fild with NA
    rm(Sigma_hat_acc_save, step_it, u, Sigma_hat_acc, rho_hat_pass) # Remove all variables not needed


  } else {

    Sigma_hat <- CholWishart::rInvWishart(iter, df = (n - m + k - 2), Sigma_ast)

  }

  B_ast_vec <- c(B_ast)

  lazy <- Sigma_hat[ ,  , 1] %x% solve(t(X) %*% X) # Get the dimensions in a lazy way

  Lambda <- array(dim = c(nrow(lazy), ncol(lazy), iter))
  for(i in 1:iter) Lambda[ , , i] <- Sigma_hat[ , , i] %x% solve(t(X) %*% X) # %x% = Kronecker
  rm(lazy)

  B_vec <- matrix(ncol = (m+1)*k, nrow = iter)
  for(i in 1:iter) B_vec[i, ] <- MASS::mvrnorm(1, mu = B_ast_vec, Lambda[ , , i])

  mu_hat <- matrix(ncol = k, nrow = iter) # Ska det verkligen vara m +1 och inte k kolumner?!!!!!!!!!!!!!!!
  for (i in 1:iter) mu_hat[i, ] <- matrix(B_vec[i, ], ncol = (m + 1), byrow = TRUE) %*% c(1, case_covar) # THIS SEEMS CORRECT NOW

  # Each row indicates the conditional expected values
  # of the case on the tests. case_covar = values from the covariates

  rho_hat <- Sigma_hat[1, 2, ] / sqrt(Sigma_hat[1, 1, ] * Sigma_hat[2, 2, ])

  s1_hat <- sqrt(Sigma_hat[1, 1, ])

  s2_hat <- sqrt(Sigma_hat[2, 2, ])

  z1 <- (case_tasks[1] - mu_hat[ , 1]) / s1_hat ## IS THE MEAN THE INTERCEPT + COEF OR JUST COEF?
  z2 <- (case_tasks[2] - mu_hat[ , 2]) / s2_hat ## IS THE MEAN THE INTERCEPT + COEF OR JUST COEF?

  z_hat_dccc <- (z1-z2) / sqrt(2 - 2*rho_hat)

  if (alternative == "two.sided") {
    pval <- 2 * stats::pnorm(abs(z_hat_dccc), lower.tail = FALSE)
  } else if (alternative == "greater") {
    pval <- stats::pnorm(z_hat_dccc, lower.tail = FALSE)
  } else { # I.e. if alternative == "less"
    pval <- stats::pnorm(z_hat_dccc, lower.tail = TRUE)
  }

  alpha <- 1 - int.level

  z_ast_est <- mean(z_hat_dccc)
  names(z_ast_est) <- "est. z"

  zdccc_int <- stats::quantile(z_hat_dccc, c(alpha/2, (1 - alpha/2)))
  names(zdccc_int) <- c("Lower zdccc CI", "Upper zdccc CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100 # This seems wrong
  names(p_int) <- c("Lower p CI", "Upper p CI")

  std.y1 <- (case_tasks[1] - m_ct[1]) / sd_ct[1]
  std.y2 <- (case_tasks[2] - m_ct[2]) / sd_ct[2]

  zdccc <- std.y1 - std.y2 / sqrt(2 - 2*r) # THIS DOES NOT ALIGN WITH C&G. PERHAPS EMAIL GARTHWAITE.

  estimate <- c(std.y1, std.y2, zdccc, ifelse(alternative == "two.sided", (p_est/2*100), p_est*100))

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

  zdccc.name <- paste0("Std. effect size (Z-DCCC) for task diff. between case and controls, ",
                      100*int.level, "% credible interval [",
                      format(round(zdccc_int[1], 2), nsmall = 2),", ",
                      format(round(zdccc_int[2], 2), nsmall = 2),"]")

  names(estimate) <- c("Case score on task X as standard (z) score",
                       "Case score on task Y as standard (z) score",
                       zdccc.name,
                       p.name)



  typ.int <- 100*int.level
  names(typ.int) <- "Interval level (%)"
  interval <- c(typ.int, zdccc_int, p_int)


  colnames(control_tasks) <- c("Y1", "Y2")
  xname <- c()
  for (i in 1:length(case_covar)) xname[i] <- paste0("X", i)
  control_covar <- matrix(control_covar, ncol = length(case_covar),
                          dimnames = list(NULL, xname))
  cor.mat <- stats::cor(cbind(control_tasks, control_covar))

  desc <- data.frame(Means = colMeans(cbind(control_tasks, control_covar)),
                     SD = apply(cbind(control_tasks, control_covar), 2, stats::sd),
                     Case_score = c(case_tasks, case_covar))



  # names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0("Case score Y1: ", format(round(case_tasks[1], 2), nsmall = 2), ", ",
                  "Case score Y2: ", format(round(case_tasks[2], 2), nsmall = 2), ", ",
                  "Controls score Y1: ", format(round(m_ct[1], 2), nsmall = 2), ", ",
                  "Controls score Y2: ", format(round(m_ct[1], 2), nsmall = 2))

  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = z_ast_est,
                 #parameter = df,
                 p.value = p_est,
                 estimate = estimate,
                 null.value = null.value,
                 interval = interval,
                 desc = desc,
                 cor.mat = cor.mat,
                 sample.size = n,
                 alternative = alternative,
                 method = paste("Bayesian Standardised Difference Test with Covariates"),
                 data.name = dname)

  class(output) <- "htest"
  output

}


