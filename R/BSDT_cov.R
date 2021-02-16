#' Bayesian Standardised Difference Test with Covariates
#'
#' Takes two single observations from a case on two variables (A and B) and
#' compares their standardised discrepancy to the discrepancies of the variables
#' in a control sample, while controlling for the effects of covariates, using
#' Bayesian methodology. This test is used when assessing a case conditioned on
#' some other variable, for example, assessing abnormality of discrepancy when
#' controlling for years of education or sex. Under the null hypothesis the case
#' is an observation from the distribution of discrepancies between the tasks of
#' interest coming from observations having the same score as the case on the
#' covariate(s). Returns a significance test, point and interval estimates of
#' difference between the case and the mean of the controls as well as point and
#' interval estimates of abnormality, i.e. an estimation of the proportion of
#' controls that would exhibit a more extreme conditioned score. This test is
#' based on random number generation which means that results may vary between
#' runs. This is by design and the reason for not using \code{set.seed()} to
#' reproduce results inside the function is to emphasise the randomness of the
#' test. To get more accurate and stable results please increase the number of
#' iterations by increasing \code{iter} whenever feasible. Developed by
#' Crawford, Garthwaite and Ryan (2011).
#'
#' Uses random generation of inverse wishart distributions from the
#' CholWishart package (Geoffrey Thompson, 2019).
#'
#' @param case_tasks A vector of length 2. The case scores from the two tasks.
#' @param case_covar A vector containing the case scores on all covariates
#'   included.
#' @param control_tasks A matrix or dataframe with 2 columns and n rows
#'   containing the control scores for the two tasks. Or if \code{use_sumstats}
#'   is set to \code{TRUE} a 2x2 matrix or dataframe
#'   containing summary statistics where the first column represents the means
#'   for each task and the second column represents the standard deviation.
#' @param control_covar A matrix or dataframe containing the control scores on
#'   the covariates included. Or if \code{use_sumstats}
#'   is set to \code{TRUE} a matrix or dataframe containing summary
#'   statistics where the first column represents the means for each covariate
#'   and the second column represents the standard deviation.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter. Since the direction
#'   of the expected effect depends on which task is set as A and which is set
#'   as B, be very careful if changing this parameter.
#' @param int_level The probability level on the Bayesian credible intervals, defaults to 95\%.
#' @param calibrated Whether or not to use the standard theory (Jeffreys) prior
#'   distribution (if set to \code{FALSE}) or a calibrated prior examined by
#'   Berger and Sun (2008). The sample estimation of the covariance matrix is based on
#'   the sample size being n - 1 when the calibrated prior is used. See Crawford
#'   et al. (2011) for further information. Calibrated prior is recommended.
#' @param iter Number of iterations to be performed. Greater number gives better
#'   estimation but takes longer to calculate. Defaults to 10000.
#' @param use_sumstats If set to \code{TRUE}, \code{control_tasks} and
#'   \code{control_covar} are treated as matrices with summary statistics. Where
#'   the first column represents the means for each variable and the second
#'   column represents the standard deviation.
#' @param cor_mat A correlation matrix of all variables included. NOTE: the two
#'   first variables should be the tasks of interest. Only needed if \code{use_sumstats}
#'   is set to \code{TRUE}.
#' @param sample_size An integer specifying the sample size of the controls.
#'   Only needed if \code{use_sumstats}
#'   is set to \code{TRUE}.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab the average z-value over
#'   \code{iter} number of iterations. \cr\cr \code{parameter} \tab the degrees
#'   of freedom used to specify the posterior distribution. \cr\cr
#'   \code{p.value}    \tab the average p-value over \code{iter} number of
#'   iterations. \cr\cr \code{estimate} \tab case scores expressed as z-scores
#'   on task A and B. Standardised effect size (Z-DCCC) of task difference
#'   between case and controls and point estimate of the proportion of the
#'   control population estimated to show a more extreme task difference. \cr\cr
#'   \code{null.value}   \tab the value of the difference between tasks under
#'   the null hypothesis.\cr\cr \code{interval} \tab named numerical vector
#'   containing level of confidence and confidence intervals for both effect
#'   size and p-value.\cr\cr \code{desc}     \tab data frame containing means
#'   and standard deviations for controls as well as case scores. \cr\cr
#'   \code{cor.mat} \tab matrix giving the correlations between the tasks of
#'   interest and the covariates included. \cr\cr \code{sample.size}     \tab
#'   number of controls. \cr\cr \code{alternative} \tab a character string
#'   describing the alternative hypothesis.\cr\cr \code{method} \tab a character
#'   string indicating what type of test was performed.\cr\cr \code{data.name}
#'   \tab a character string giving the name(s) of the data}
#' @export
#'
#' @examples
#' BSDT_cov(case_tasks = c(size_weight_illusion[1, "V_SWI"],
#'                         size_weight_illusion[1, "K_SWI"]),
#'          case_covar = size_weight_illusion[1, "YRS"],
#'          control_tasks = cbind(size_weight_illusion[-1, "V_SWI"],
#'                                size_weight_illusion[-1, "K_SWI"]),
#'          control_covar = size_weight_illusion[-1, "YRS"], iter = 100)
#'
#' @references
#'
#' Berger, J. O., & Sun, D. (2008). Objective Priors for the Bivariate Normal
#' Model. \emph{The Annals of Statistics, 36}(2), 963-982. JSTOR.
#'
#' Crawford, J. R., Garthwaite, P. H., & Ryan, K. (2011). Comparing a single
#' case to a control sample: Testing for neuropsychological deficits and
#' dissociations in the presence of covariates. \emph{Cortex, 47}(10),
#' 1166-1178. \doi{10.1016/j.cortex.2011.02.017}
#'
#' #' Geoffrey Thompson (2019). CholWishart: Cholesky Decomposition of the Wishart
#' Distribution. R package version 1.1.0.
#' \url{https://CRAN.R-project.org/package=CholWishart}

BSDT_cov <- function (case_tasks, case_covar, control_tasks, control_covar,
                      alternative = c("two.sided", "greater", "less"),
                      int_level = 0.95, calibrated = TRUE, iter = 10000,
                      use_sumstats = FALSE, cor_mat = NULL, sample_size = NULL) {

  ###
  # Set up of error and warning messages
  ###

  if (use_sumstats & (is.null(cor_mat) | is.null(sample_size))) stop("Please supply both correlation matrix and sample size")
  if (int_level < 0 | int_level > 1) stop("Interval level must be between 0 and 1")
  if (length(case_tasks) != 2) stop("case_task should have length 2")
  if (ncol(control_tasks) != 2) stop("columns in control_tasks should be 2")
  if (!is.null(cor_mat)) if (sum(eigen(cor_mat)$values > 0) < length(diag(cor_mat))) stop("cor_mat is not positive definite")
  if (!is.null(cor_mat) | !is.null(sample_size)) if (use_sumstats == FALSE) stop("If input is summary data, set use_sumstats = TRUE")
  if (sum(is.na(control_tasks)) > 0) stop("control_tasks contains NA")
  if (sum(is.na(control_covar)) > 0) stop("control_covar contains NA")

  if (use_sumstats == FALSE){
    cov_obs <- ifelse(is.vector(control_covar), length(control_covar), nrow(control_covar))
    if (nrow(control_tasks) != cov_obs) stop("Must supply equal number of observations for tasks and covariates")
    rm(cov_obs)
  }

  if (is.data.frame(control_tasks)) control_tasks <- as.matrix(control_tasks)
  if (is.data.frame(control_covar)) control_covar <- as.matrix(control_covar)

  ###
  # Extract relevant statistics
  ###

  alternative <- match.arg(alternative)

  if (use_sumstats) {

    # If summary statistics are used as input this is a lazy way of
    # generating corresponding data

    sum_stats <- rbind(control_tasks, control_covar)

    cov_mat <- diag(sum_stats[ , 2]) %*% cor_mat %*% diag(sum_stats[ , 2])

    lazy_gen <- MASS::mvrnorm(sample_size, mu = sum_stats[ , 1], Sigma = cov_mat, empirical = TRUE)

    control_tasks <- lazy_gen[ , 1:2]

    control_covar <- lazy_gen[ , -c(1, 2)]

  }

  n <- nrow(control_tasks)

  m <- length(case_covar)

  k <- length(case_tasks)

  ###
  # Data estimation
  ###

  m_ct <- colMeans(control_tasks)
  sd_ct <- apply(control_tasks, 2, stats::sd)
  r <- stats::cor(control_tasks)[1, 2]

  X <- cbind(rep(1, n), control_covar) # Design matrix
  Y <- control_tasks # Response matrix

  B_ast <- solve(t(X) %*% X) %*% t(X) %*% Y # Regression coefficients

  Sigma_ast <- t(Y - X %*% B_ast) %*% (Y - X %*% B_ast)

  ###
  # Get parameter distributions of Sigma depending on prior chosen
  ###

  if (calibrated == TRUE) {

    ###
    # Below follows the rejection sampling of Sigma
    # for the "calibrated" prior detailed in the vignette
    ###

    A_ast <- ((n - m - 2)*Sigma_ast) / (n - m - 1)
    df <- (n - m - 2)

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

    Sigma_hat <- Sigma_hat_acc_save[ , , -1] # Remove the first matrix that is filled with NA
    rm(Sigma_hat_acc_save, step_it, u, Sigma_hat_acc, rho_hat_pass) # Flush all variables not needed


  } else {

    # If the "standard theory" prior is requested

    df <- (n - m + k - 2)
    Sigma_hat <- CholWishart::rInvWishart(iter, df = (n - m + k - 2), Sigma_ast)

  }

  ###
  # Get parameter distributions as described in vignette
  ###

  B_ast_vec <- c(B_ast)

  lazy <- Sigma_hat[ ,  , 1] %x% solve(t(X) %*% X) # Get the dimensions in a lazy way

  Lambda <- array(dim = c(nrow(lazy), ncol(lazy), iter))
  for(i in 1:iter) Lambda[ , , i] <- Sigma_hat[ , , i] %x% solve(t(X) %*% X) # %x% = Kronecker
  rm(lazy)

  B_vec <- matrix(ncol = (m+1)*k, nrow = iter)
  for(i in 1:iter) B_vec[i, ] <- MASS::mvrnorm(1, mu = B_ast_vec, Lambda[ , , i])

  mu_hat <- matrix(ncol = k, nrow = iter)
  for (i in 1:iter) mu_hat[i, ] <- matrix(B_vec[i, ], ncol = (m + 1), byrow = TRUE) %*% c(1, case_covar)

  # Each row indicates the conditional case scores. case_covar = values from the covariates

  ###
  # Get conditional effect and p distributions
  ###

  rho_hat <- Sigma_hat[1, 2, ] / sqrt(Sigma_hat[1, 1, ] * Sigma_hat[2, 2, ])

  s1_hat <- sqrt(Sigma_hat[1, 1, ])

  s2_hat <- sqrt(Sigma_hat[2, 2, ])

  z1 <- (case_tasks[1] - mu_hat[ , 1]) / s1_hat
  z2 <- (case_tasks[2] - mu_hat[ , 2]) / s2_hat

  z_hat_dccc <- (z1-z2) / sqrt(2 - 2*rho_hat)

  if (alternative == "two.sided") {
    pval <- 2 * stats::pnorm(abs(z_hat_dccc), lower.tail = FALSE)
  } else if (alternative == "greater") {
    pval <- stats::pnorm(z_hat_dccc, lower.tail = FALSE)
  } else { # I.e. if alternative == "less"
    pval <- stats::pnorm(z_hat_dccc, lower.tail = TRUE)
  }

  ###
  # Get and name intervals
  ###

  alpha <- 1 - int_level

  zdccc_int <- stats::quantile(z_hat_dccc, c(alpha/2, (1 - alpha/2)))
  names(zdccc_int) <- c("Lower Z-DCCC CI", "Upper Z-DCCC CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100
  names(p_int) <- c("Lower p CI", "Upper p CI")

  ###
  # Get the point estimates for the conditional effects
  ###

  z.y1 <- (case_tasks[1] - m_ct[1]) / sd_ct[1]
  z.y2 <- (case_tasks[2] - m_ct[2]) / sd_ct[2]

  mu_ast <- t(B_ast) %*% c(1, case_covar)

  cov_ast <- Sigma_ast/(n - 1)

  rho_ast <- cov_ast[1, 2]/sqrt(cov_ast[1, 1] * cov_ast[2, 2])

  std.y1 <- (case_tasks[1] - mu_ast[1]) / sqrt(cov_ast[1, 1])
  std.y2 <- (case_tasks[2] - mu_ast[2]) / sqrt(cov_ast[2, 2])

  zdccc <-  (std.y1 - std.y2) / sqrt(2 - 2*rho_ast)

  ###
  # Name and store estimates
  ###

  prop <- ifelse(alternative == "two.sided", round((p_est/2*100), 2), p_est*100)
  estimate <- round(c(z.y1, z.y2, zdccc, prop), 6)

  if (alternative == "two.sided") {
    if (zdccc < 0) {
      alt.p.name <- "Proportion below case (%), "
    } else {
      alt.p.name <- "Proportion above case (%), "
    }
  } else if (alternative == "greater") {
    alt.p.name <- "Proportion above case (%), "
  } else {
    alt.p.name <- "Proportion below case (%), "
  }

  # The below sets the intervals to be shown in the print() output

  p.name <- paste0(alt.p.name,
                   100*int_level, "% CI [",
                   format(round(p_int[1], 2), nsmall = 2),", ",
                   format(round(p_int[2], 2), nsmall = 2),"]")

  zdccc.name <- paste0("Std. discrepancy (Z-DCCC), ",
                      100*int_level, "% CI [",
                      format(round(zdccc_int[1], 2), nsmall = 2),", ",
                      format(round(zdccc_int[2], 2), nsmall = 2),"]")

  names(estimate) <- c("Std. case score, task A (Z-CC)",
                       "Std. case score, task B (Z-CC)",
                       zdccc.name,
                       p.name)

  ###
  # Set names and type of intervals for output object
  ###

  typ.int <- 100*int_level
  names(typ.int) <- "Credible (%)"
  interval <- c(typ.int, zdccc_int, p_int)


  ###
  # Set names for each covariate (and the variates of interest)
  ###

  colnames(control_tasks) <- c("A", "B")
  xname <- c()
  for (i in 1:length(case_covar)) xname[i] <- paste0("COV", i)
  control_covar <- matrix(control_covar, ncol = length(case_covar),
                          dimnames = list(NULL, xname))
  cor.mat <- stats::cor(cbind(control_tasks, control_covar))

  ###
  # Create object with descriptives
  ###

  desc <- data.frame(Means = colMeans(cbind(control_tasks, control_covar)),
                     SD = apply(cbind(control_tasks, control_covar), 2, stats::sd),
                     Case_score = c(case_tasks, case_covar))



  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0("Case A = ", format(round(case_tasks[1], 2), nsmall = 2), ", ",
                  "B = ", format(round(case_tasks[2], 2), nsmall = 2), ", ",
                  "Ctrl. (m, sd) A: (", format(round(m_ct[1], 2), nsmall = 2),",", format(round(sd_ct[1], 2), nsmall = 2) , "), ",
                  "B: (", format(round(m_ct[2], 2), nsmall = 2),",", format(round(sd_ct[2], 2), nsmall = 2) , ")")

  # Build output to be able to set class as "htest" object for S3 methods.
  # See documentation for "htest" class for more info

  output <- list(parameter = df,
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

