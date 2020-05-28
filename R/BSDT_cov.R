#' Bayesian Standardised Difference Test with Covariates
#'
#' Tests whether the standardized difference between a case's scores
#' on two tasks (Y1 and Y2) is significantly different from the standardized
#' differences between these tasks in a control sample, while controlling for
#' the effects of covariates.  For example, it could be used to test whether a
#' case's standardized difference between two tasks is significantly greater
#' than the standardized differences in controls when controlling for a measure
#' of processing speed.  Under the null hypothesis the case's standardized
#' difference, conditioned on the covariate(s), is an observation from the
#' distribution of conditional standardized differences in the control
#' population. Returns (a) a signficance test, (b) point and
#' interval estimates of the effect size for the difference between the case and
#' controls, and (c) point and interval estimates of the abnormality of the
#' case's standardized difference (i.e., it estimates the percentage of controls
#' that would exhibit a more extreme standardized difference).
#'
#' Developed by Crawford, Garthwaite and Ryan (2011).
#'
#' @param case_tasks A vector of length 2. The case scores from the two tasks.
#' @param case_covar A vector containing the case scores on all covariates
#'   included.
#' @param control_tasks A matrix or dataframe with 2 columns and n rows
#'   containing the control scores for the two tasks. Or a matrix or dataframe
#'   containing summary statistics where the first column represents the means
#'   for each task and the second column represents the standard deviation.
#' @param control_covar A matrix or dataframe cointaining the control scores on
#'   the covariates included. Or a matrix or dataframe containing summary
#'   statistics where the first column represents the means for each covariate
#'   and the second column represents the standard deviation.
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
#' @param use_sumstats If set to \code{TRUE}, \code{control_tasks} and
#'   \code{control_covar} are treated as matrices with summary statistics. Where
#'   the first column represents the means for each variable and the second
#'   column represents the standard deviation.
#' @param cor_mat A correlation matrix of all variables included. NOTE: the two
#'   first variables should be the tasks of interest.
#' @param control_n An integer specifying the sample size of the controls.
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
#' @references
#'
#' Berger, J. O., & Sun, D. (2008). Objective Priors for the Bivariate Normal
#' Model. \emph{The Annals of Statistics, 36}(2), 963-982. JSTOR.
#'
#' Crawford, J. R., Garthwaite, P. H., & Ryan, K. (2011). Comparing a single
#' case to a control sample: Testing for neuropsychological deficits and
#' dissociations in the presence of covariates. \emph{Cortex, 47}(10),
#' 1166-1178. https://doi.org/10.1016/j.cortex.2011.02.017

BSDT_cov <- function (case_tasks, case_covar, control_tasks, control_covar,
                      alternative = c("two.sided", "greater", "less"),
                      int.level = 0.95, calibrated = TRUE, iter = 1000,
                      use_sumstats = FALSE, cor_mat = NULL, control_n = NULL) {

  alternative <- match.arg(alternative)

  # ERRORS BELOW

  if (use_sumstats & (is.null(cor_mat) | is.null(control_n))) stop("Please supply both correlation matrix and sample size")
  if (int.level < 0 | int.level > 1) stop("Interval level must be between 0 and 1")
  if (length(case_tasks) != 2) stop("case_task should have lenght 2")
  if (ncol(control_tasks) != 2) stop("columns in control_tasks should be 2")



  if (use_sumstats) {

    sum_stats <- rbind(control_tasks, control_covar)

    cov_mat <- diag(sum_stats[ , 2]) %*% cor_mat %*% diag(sum_stats[ , 2])

    lazy_gen <- MASS::mvrnorm(control_n, mu = sum_stats[ , 1], Sigma = cov_mat, empirical = TRUE)

    control_tasks <- lazy_gen[ , 1:2]

    control_covar <- lazy_gen[ , -c(1, 2)]

  }

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

  Sigma_ast <- t(Y - X %*% B_ast) %*% (Y - X %*% B_ast) # This should be multiplied by (1/n - m - 1) but the I get a discrepancy with C&G

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

  mu_hat <- matrix(ncol = k, nrow = iter)
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
  names(z_ast_est) <- "ave. z"

  zdccc_int <- stats::quantile(z_hat_dccc, c(alpha/2, (1 - alpha/2)))
  names(zdccc_int) <- c("Lower zdccc CI", "Upper zdccc CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100 # This seems wrong
  names(p_int) <- c("Lower p CI", "Upper p CI")



  z.y1 <- (case_tasks[1] - m_ct[1]) / sd_ct[1]
  z.y2 <- (case_tasks[2] - m_ct[2]) / sd_ct[2]

  mu_ast <- t(B_ast) %*% c(1, case_covar)

  cov_ast <- Sigma_ast/(n - 1) # This is aligns witht the program but in the paper it should be (n - m - 1)

  rho_ast <- cov_ast[1, 2]/sqrt(cov_ast[1, 1] * cov_ast[2, 2])

  std.y1 <- (case_tasks[1] - mu_ast[1]) / sqrt(cov_ast[1, 1])
  std.y2 <- (case_tasks[2] - mu_ast[2]) / sqrt(cov_ast[2, 2])

  zdccc <-  (std.y1 - std.y2) / sqrt(2 - 2*rho_ast)

  estimate <- c(z.y1, z.y2, zdccc, ifelse(alternative == "two.sided", (p_est/2*100), p_est*100))

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


