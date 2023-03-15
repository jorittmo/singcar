#' Bayesian Test of Deficit with Covariates
#'
#' Takes a single observation and compares it to a distribution estimated by a
#' control sample, while controlling for the effect of covariates, using
#' Bayesian methodology. This test is used when assessing a case conditioned on
#' some other variable, for example, assessing abnormality when controlling for
#' years of education or sex. Under the null hypothesis the case is an
#' observation from the distribution of scores from the task of interest coming
#' from observations having the same score as the case on the covariate(s).
#' Returns a significance test, point and interval estimates of difference
#' between the case and the mean of the controls as well as point and interval
#' estimates of abnormality, i.e. an estimation of the proportion of controls
#' that would exhibit a more extreme conditioned score. This test is based on
#' random number generation which means that results may vary between runs. This
#' is by design and the reason for not using \code{set.seed()} to reproduce
#' results inside the function is to emphasise the randomness of the test. To
#' get more accurate and stable results please increase the number of iterations
#' by increasing \code{iter} whenever feasible. Developed by Crawford,
#' Garthwaite and Ryan (2011).
#'
#' Uses random generation of inverse wishart distributions from the
#' CholWishart package (Geoffrey Thompson, 2019).
#'
#' @param case_task The case score from the task of interest. Must be a single
#'   value.
#' @param case_covar A vector containing the case scores on all covariates
#'   included. Can be of any length except 0, in that case use
#'   \code{\link{BTD}}.
#' @param control_task A vector containing the scores from the controls on the
#'   task of interest. Or a vector of length 2 containing the mean and standard
#'   deviation of the task. In that order.
#' @param control_covar A vector, matrix or dataframe containing the control
#'   scores on the covariates included. If matrix or dataframe each column
#'   represents a covariate. Or a matrix or dataframe containing summary
#'   statistics where the first column represents the means for each covariate
#'   and the second column represents the standard deviation.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter.
#' @param int_level The probability level on the Bayesian credible intervals, defaults to 95\%.
#' @param iter Number of iterations to be performed. Greater number gives better
#'   estimation but takes longer to calculate. Defaults to 10000.
#' @param use_sumstats If set to \code{TRUE}, \code{control_tasks} and
#'   \code{control_covar} are treated as matrices with summary statistics. Where
#'   the first column represents the means for each variable and the second
#'   column represents the standard deviation.
#' @param cor_mat A correlation matrix of all variables included. NOTE: the
#'   first variable should be the task of interest.
#' @param sample_size An integer specifying the sample size of the controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab the average z-value over
#'   \code{iter} number of iterations. \cr\cr \code{parameter} \tab the degrees
#'   of freedom used to specify the posterior distribution. \cr\cr
#'   \code{p.value}    \tab the average p-value over \code{iter} number of
#'   iterations. \cr\cr \code{estimate} \tab case scores expressed as z-scores
#'   on task X and Y. Standardised effect size (Z-CCC) of task difference
#'   between case and controls and point estimate of the proportion of the
#'   control population estimated to show a more extreme task difference. \cr\cr
#'   \code{null.value}   \tab the value of the difference between tasks under
#'   the null hypothesis.\cr\cr \code{interval} \tab named numerical vector
#'   containing level of confidence and confidence intervals for both effect
#'   size and p-value.\cr\cr \code{desc}     \tab data frame containing means
#'   and standard deviations for controls as well as case scores. \cr\cr
#'   \code{cor.mat} \tab matrix giving the correlations between the task of
#'   interest and the covariates included. \cr\cr \code{sample.size} \tab number
#'   of controls..\cr\cr \code{alternative}     \tab a character string
#'   describing the alternative hypothesis.\cr\cr \code{method} \tab a character
#'   string indicating what type of test was performed.\cr\cr \code{data.name}
#'   \tab a character string giving the name(s) of the data}
#' @export
#'
#' @examples
#'
#' BTD_cov(case_task = size_weight_illusion[1, "V_SWI"],
#'          case_covar = size_weight_illusion[1, "YRS"],
#'          control_task = size_weight_illusion[-1, "V_SWI"],
#'          control_covar = size_weight_illusion[-1, "YRS"], iter = 100)
#'
#' @references
#'
#' Crawford, J. R., Garthwaite, P. H., & Ryan, K. (2011). Comparing
#' a single case to a control sample: Testing for neuropsychological deficits
#' and dissociations in the presence of covariates. \emph{Cortex, 47}(10),
#' 1166-1178. \doi{10.1016/j.cortex.2011.02.017}
#'
#' Geoffrey Thompson (2019). CholWishart: Cholesky Decomposition of the Wishart
#' Distribution. R package version 1.1.0.
#' \url{https://CRAN.R-project.org/package=CholWishart}

BTD_cov <- function (case_task, case_covar, control_task, control_covar,
                     alternative = c("less", "two.sided", "greater"),
                     int_level = 0.95, iter = 10000,
                     use_sumstats = FALSE, cor_mat = NULL, sample_size = NULL) {

  ###
  # Set up of error and warning messages
  ###

  case_task <- as.numeric(unlist(case_task))
  case_covar <- as.numeric(unlist(case_covar))

  control_task <- as.numeric(unlist(control_task))

  if (use_sumstats & (is.null(cor_mat) | is.null(sample_size))) stop("Please supply both correlation matrix and sample size")
  if (int_level < 0 | int_level > 1) stop("Interval level must be between 0 and 1")
  if (length(case_task) != 1) stop("case_task should be single value")
  if (!is.null(cor_mat)) if (sum(eigen(cor_mat)$values > 0) < length(diag(cor_mat))) stop("cor_mat is not positive definite")
  if (!is.null(cor_mat) | !is.null(sample_size)) if (use_sumstats == FALSE) stop("If input is summary data, set use_sumstats = TRUE")

  if (!is.matrix(control_covar) & !is.vector(control_covar)) control_covar <- as.matrix(control_covar)

  ###
  # Extract relevant statistics
  ###

  alternative <- match.arg(alternative)

  if (use_sumstats) {

    ###
    # If summary statistics are used as input this is a lazy way of
    # generating corresponding data
    ###

    sum_stats <- rbind(control_task, control_covar)

    if (length(sum_stats[ , 2]) != nrow(cor_mat)) stop("Number of variables and number of correlations does not match")

    cov_mat <- diag(sum_stats[ , 2]) %*% cor_mat %*% diag(sum_stats[ , 2])

    lazy_gen <- MASS::mvrnorm(sample_size, mu = sum_stats[ , 1], Sigma = cov_mat, empirical = TRUE)

    control_task <- lazy_gen[ , 1]

    control_covar <- lazy_gen[ , -1]

  }

  alpha <- 1 - int_level

  n <- length(control_task)

  m <- length(case_covar)

  k <- length(case_task)

  m_ct <- mean(control_task)

  sd_ct <- stats::sd(control_task)

  df <- (n - m + k - 2)

  ###
  # Data estimation
  ###

  X <- cbind(rep(1, n), control_covar) # Design matrix
  Y <- control_task # Response matrix

  B_ast <- solve(t(X) %*% X) %*% t(X) %*% Y # Regression coefficients

  Sigma_ast <- t(Y - X %*% B_ast) %*% (Y - X %*% B_ast) # Now sums of squares

  ###
  # Get parameter distributions as described in vignette
  ###

  Sigma_hat <- c(CholWishart::rInvWishart(iter, df = (n - m + k - 2), Sigma_ast))

  B_ast_vec <- c(B_ast)

  lazy <- Sigma_hat[1] %x% solve(t(X) %*% X) # Get the dimensions in a lazy way

  Lambda <- array(dim = c(nrow(lazy), ncol(lazy), iter))
  for(i in 1:iter) Lambda[ , , i] <- Sigma_hat[i] %x% solve(t(X) %*% X) # %x% = Kronecker product
  rm(lazy)

  B_vec <- matrix(ncol = (m+1)*k, nrow = iter)
  for(i in 1:iter) B_vec[i, ] <- MASS::mvrnorm(1, mu = B_ast_vec, Lambda[ , , i])


  mu_hat <- vector(length = iter)
  for (i in 1:iter) mu_hat[i] <- matrix(B_vec[i, ], ncol = (m+1), byrow = TRUE) %*% c(1, case_covar)

  ###
  # Get distributions for effect and p estimates
  ###

  z_hat_ccc <- (case_task - mu_hat) / sqrt(Sigma_hat)

  if (alternative == "less") {

    pval <- stats::pnorm(z_hat_ccc)

  } else if (alternative == "greater") {

    pval <- stats::pnorm(z_hat_ccc, lower.tail = FALSE)

  } else if (alternative == "two.sided") {

    pval <- 2 * stats::pnorm(abs(z_hat_ccc), lower.tail = FALSE)

  }

  mu_ast <- t(B_ast) %*% c(1, case_covar) # Conditional means

  var_ast <- Sigma_ast/(n - 1)

  zccc <- (case_task - mu_ast) / sqrt(var_ast) # Calculate point estimate based on conditional sample means and sd

  ###
  # Get intervals for effect and p estimates
  ###

  zccc_int <- stats::quantile(z_hat_ccc, c(alpha/2, (1 - alpha/2)))
  names(zccc_int) <- c("Lower Z-CCC CI", "Upper Z-CCC CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100
  names(p_int) <- c("Lower p CI", "Upper p CI")

  estimate <- round(c(zccc, p_est*100), 6)
  if (alternative == "two.sided") estimate <- c(zccc, (p_est/2)*100)

  ###
  # Set names for estimates
  ###

  zccc.name <- paste0("Std. case score (Z-CCC), ",
                     100*int_level, "% CI [",
                     format(round(zccc_int[1], 2), nsmall = 2),", ",
                     format(round(zccc_int[2], 2), nsmall = 2),"]")

  if (alternative == "less") {
    alt.p.name <- "Proportion below case (%), "
  } else if (alternative == "greater") {
    alt.p.name <- "Proportion above case (%), "
  } else {
    if (zccc < 0) {
      alt.p.name <- "Proportion below case (%), "
    } else {
      alt.p.name <- "Proportion above case (%), "
    }
  }

  p.name <- paste0(alt.p.name,
                   100*int_level, "% CI [",
                   format(round(p_int[1], 2), nsmall = 2),", ",
                   format(round(p_int[2], 2), nsmall = 2),"]")


  names(estimate) <- c(zccc.name, p.name)

  ###
  # Set names for each covariate
  ###

  control_task <- matrix(control_task, ncol = 1, dimnames = list(NULL, "Task"))
  xname <- c()
  for (i in 1:length(case_covar)) xname[i] <- paste0("COV", i)
  control_covar <- matrix(control_covar, ncol = length(case_covar),
                          dimnames = list(NULL, xname))
  cor.mat <- stats::cor(cbind(control_task, control_covar))

  ###
  # Create object with descriptives
  ###

  desc <- data.frame(Means = colMeans(cbind(control_task, control_covar)),
                     SD = apply(cbind(control_task, control_covar), 2, stats::sd),
                     Case_score = c(case_task, case_covar))

  ###
  # Set objects to use in the output object
  ###

  typ.int <- 100*int_level
  names(typ.int) <- "Interval level (%)"
  interval <- c(typ.int, zccc_int, p_int)

  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "diff. between case and controls"



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
                 method = paste("Bayesian Test of Deficit with Covariates"),
                 data.name = paste0("Case = ", format(round(case_task, 2), nsmall = 2),
                                    ", Controls (m = ", format(round(m_ct, 2), nsmall = 2),
                                    ", sd = ", format(round(sd_ct, 2), nsmall = 2),
                                    ", n = ", n, ")"))

  class(output) <- "htest"
  output

}

