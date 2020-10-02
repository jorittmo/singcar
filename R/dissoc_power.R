#' Power calculator for UDT
#'
#' Calculates exact power given sample size or sample size given power, using
#' analytical methods for the frequentist test of deficit for a specified case
#' scores, means and standard deviations for the control sample. The means and
#' standard deviations defaults to 0 and 1 respecitvely, so if no other values
#' are given, the case scores are interpreted as deviations from the mean in
#' standard deviations. The returned value will approximate the power for the
#' given parameters.
#'
#' @param case_a A single value from the expected case observation on task A.
#' @param case_b A single value from the expected case observation on task B.
#' @param mean_a The expected mean from the control sample on task A. Defaults
#'   to 0.
#' @param mean_b The expected mean from the control sample on task B. Defaults
#'   to 0.
#' @param sd_a The expected standard deviation from the control sample on task
#'   A. Defaults to 1.
#' @param sd_b The expected standard deviation from the control sample on task
#'   B. Defaults to 1.
#' @param r_ab The expected correlation between the tasks. Defaults to 0.5
#' @param sample_size The size of the control sample, vary this parameter to see
#'   how the sample size affects power. One of sample size or power must be
#'   specified, not both.
#' @param power A single value between 0 and 1 specifying desired power for
#'   calculating necessary sample size. One of sample size or power must be
#'   specified, not both.
#' @param alternative The alternative hypothesis. A string of either "two.sided"
#'   (default) or "one.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power. Defaults to 0.05.
#' @param spec A single value between 0 and 1. If desired power is given as
#'   input the function will utilise a search algorithm to find the sample size
#'   needed to reach the desired power. However, if the power specified is
#'   greater than what is actually possible to achieve the algorithm could
#'   search forever. Hence, when power does not increase substantially for any
#'   additional participant in the sample, the algorthm stops. By default the
#'   algorithm stops when power does not increase more than 0.5% for any added
#'   participant, but by varying \code{spec}, this specificity can be changed.
#'
#' @return Either a single value of the exact power, if sample size is given. Or
#'   a dataframe consisting of both the sample size and the exact power such
#'   size would yield.
#' @export
#'
#' @examples
#' UDT_power(case_a = -3, case_b = -1, mean_a = 0, mean_b = 0,
#'           sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20)
#' UDT_power(case_a = -3, case_b = -1, power = 0.8)
UDT_power <- function(case_a, case_b, mean_a = 0, mean_b = 0,
                      sd_a = 1, sd_b = 1, r_ab = 0.5,
                      sample_size = NULL, power = NULL,
                      alternative = c("two.sided", "greater", "less"),
                      alpha = 0.05, spec = 0.005) {

  if (!is.null(sample_size) & !is.null(power)) stop("Must supply only one of sample size or desired power")
  if (is.null(sample_size) & is.null(power)) stop("Must supply either sample size or desired power")
  if (!is.null(power)) if (power > 1 | power < 0) stop("Desired power must be between 0 and 1")
  if (!is.null(sample_size)) if (sample_size < 2) stop("Sample size must be greater than 1")
  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")
  if (r_ab < -1 | r_ab > 1) stop("Correlation between task a and b must be between -1 and 1")

  alternative <- match.arg(alternative)

  n = sample_size

  if (is.null(power)) {

    if (alternative == "two.sided") {

      tstat <- ((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

      power = stats::pt(stats::qt(alpha/2, df = n-1,
                                  lower.tail = T),
                        ncp = tstat,
                        df = n-1,
                        lower.tail = T) - stats::pt(-stats::qt(alpha/2, df = n-1,
                                                               lower.tail = T),
                                                    ncp = tstat,
                                                    df = n-1,
                                                    lower.tail = T) + 1

      return(power)

    }

    if (alternative == "less") {

      tstat <- ((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

      power = stats::pt(stats::qt(alpha, df = n-1,
                                       lower.tail = T),
                             ncp = tstat,
                             df = n-1,
                             lower.tail = T
      )
      return(power)
    }

    if (alternative == "greater") {

      tstat <- ((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

      power = stats::pt(stats::qt(alpha, df = n-1,
                                       lower.tail = F),
                             ncp = tstat,
                             df = n-1,
                             lower.tail = F
      )
      return(power)

    }

  }

  if (is.null(sample_size)) {

    n = 2

    search_pwr = 0

    prev_search_pwr = 1

    keep_going = TRUE



      while(keep_going == TRUE){

        if (alternative == "two.sided") {
          tstat <- ((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

          search_pwr = stats::pt(stats::qt(alpha/2, df = n-1,
                                           lower.tail = T),
                                 ncp = tstat,
                                 df = n-1,
                                 lower.tail = T) - stats::pt(-stats::qt(alpha/2, df = n-1,
                                                                        lower.tail = T),
                                                             ncp = tstat,
                                                             df = n-1,
                                                             lower.tail = T) + 1
        }




        if (alternative == "less") {

          tstat <- ((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

          search_pwr = stats::pt(stats::qt(alpha, df = n-1,
                                           lower.tail = T),
                                 ncp = tstat,
                                 df = n-1,
                                 lower.tail = T
          )

        }

        if (alternative == "greater") {

          tstat <- ((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

          search_pwr = stats::pt(stats::qt(alpha, df = n-1,
                                           lower.tail = F),
                                 ncp = tstat,
                                 df = n-1,
                                 lower.tail = F
          )


        }

        if (abs(prev_search_pwr - search_pwr) < spec) {
          keep_going = FALSE
          message(paste0("Power (", format(round(search_pwr, 5), nsmall = 5), ") will not increase more than ", spec*100, "%",
                       " for any additional participant over n = ", n))
        }

        prev_search_pwr = search_pwr

        n = n + 1

        if (search_pwr > power) keep_going = FALSE

      }

      return(data.frame(n = n - 1, power = search_pwr))

    }

}


#' Power calculator for RSDT
#'
#' Calculates approximate power, given sample size, using Monte Carlo
#' simulation, for specified case scores, means and standard deviations for the
#' control sample. The means and standard deviations defaults to 0 and 1
#' respectively, so if no other values are given the case scores are interpreted
#' as deviations from the mean in standard deviations. Hence, the effect size of
#' the dissociation (Z-DCC) would in that case be the difference between the two
#' case scores.
#'
#' @param case_a A single value from the expected case observation on task A.
#' @param case_b A single value from the expected case observation on task B.
#' @param mean_a The expected mean from the control sample on task A. Defaults
#'   to 0.
#' @param mean_b The expected mean from the control sample on task B. Defaults
#'   to 0.
#' @param sd_a The expected standard deviation from the control sample on task
#'   A. Defaults to 1.
#' @param sd_b The expected standard deviation from the control sample on task
#'   B. Defaults to 1.
#' @param r_ab The expected correlation between the tasks. Defaults to 0.5
#' @param sample_size The size of the control sample, vary this parameter to see
#'   how the sample size affects power.
#' @param alternative The alternative hypothesis. A string of either "two.sided"
#'   (default) or "one.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power. Defaults to 0.05.
#' @param nsim The number of simulations to run. Higher number gives better
#'   accuracy, but low numbers such as 10000 or even 1000 are usually sufficient
#'   for the purposes of this calculator.
#'
#' @return Returns a single value approximating the power of the test for the
#'   given parameters.
#' @export
#'
#' @examples
#' RSDT_power(case_a = -3, case_b = -1, mean_a = 0, mean_b = 0,
#'            sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20, nsim = 1000)


RSDT_power <- function(case_a, case_b, mean_a = 0, mean_b = 0,
                       sd_a = 1, sd_b = 1, r_ab = 0.5,
                       sample_size,
                       alternative = c("two.sided", "greater", "less"),
                       alpha = 0.05, nsim = 10000) {

  if (!is.null(sample_size)) if (sample_size < 2) stop("Sample size must be greater than 1")
  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")
  if (r_ab < -1 | r_ab > 1) stop("Correlation between task a and b must be between -1 and 1")

  alternative <- match.arg(alternative)

  n = sample_size

  rsdt_p_sim <- function(case_a, case_b,mean_a, mean_b,
                         sd_a, sd_b, r_ab) {

    cor_mat <- matrix(c(1, r_ab, r_ab, 1), nrow = 2)

    cov_mat <- diag(c(sd_a, sd_b)) %*% cor_mat %*% diag(c(sd_a, sd_b))

    controls <- MASS::mvrnorm(n + 1, mu = c(mean_a, mean_b), Sigma = cov_mat)

    case_a_gen <- controls[1, 1] + case_a
    case_b_gen <- controls[1, 2] + case_b

    controls <- controls[-1, ]

    pval <- singcar::RSDT(case_a_gen, case_b_gen,
                          controls[ , 1], controls[, 2],
                          alternative = alternative)[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- rsdt_p_sim(case_a, case_b, mean_a, mean_b,
                          sd_a, sd_b, r_ab)

  }

  power = sum(pval < alpha)/length(pval)


  return(power)
}

#' Power calculator for BSDT
#'
#' Calculates approximate power, given sample size, using Monte Carlo
#' simulation, for specified case scores, means and standard deviations for the
#' control sample. The means and standard deviations default to 0 and 1
#' respectively, so if no other values are given the case scores are interpreted
#' as deviations from the mean in standard deviations. Hence, the effect size of
#' the dissociation (Z-DCC) would in that case be the difference between the two
#' case scores. Is computationally heavy and might therefore take a few seconds.
#'
#' @param case_a A single value from the expected case observation on task A.
#' @param case_b A single value from the expected case observation on task B.
#' @param mean_a The expected mean from the control sample on task A. Defaults
#'   to 0.
#' @param mean_b The expected mean from the control sample on task B. Defaults
#'   to 0.
#' @param sd_a The expected standard deviation from the control sample on task
#'   A. Defaults to 1.
#' @param sd_b The expected standard deviation from the control sample on task
#'   B. Defaults to 1.
#' @param r_ab The expected correlation between the tasks. Defaults to 0.5
#' @param sample_size The size of the control sample, vary this parameter to see
#'   how the sample size affects power.
#' @param alternative The alternative hypothesis. A string of either "two.sided"
#'   (default) or "one.sided".
#' @param alpha The specified Type I error rate. This can be varied, with
#'   effects on power. Defaults to 0.05.
#' @param nsim The number of simulations to run. Higher number gives better
#'   accuracy, but low numbers such as 10000 or even 1000 are usually sufficient
#'   for the purposes of this calculator. Defaults to 1000 due to the
#'   computationally intense \code{BSTD}.
#' @param iter The number simulations used by \code{BSTD}. Defaults to 1000
#'
#' @return Returns a single value approximating the power of the test for the
#'   given parameters.
#' @export
#'
#' @examples
#' BSDT_power(case_a = -3, case_b = -1, mean_a = 0, mean_b = 0,
#'            sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20, nsim = 100, iter = 100)


BSDT_power <- function(case_a, case_b, mean_a = 0, mean_b = 0,
                       sd_a = 1, sd_b = 1, r_ab = 0.5,
                       sample_size,
                       alternative = c("two.sided", "greater", "less"),
                       alpha = 0.05, nsim = 1000, iter = 1000) {

  if (!is.null(sample_size)) if (sample_size < 2) stop("Sample size must be greater than 1")
  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")
  if (r_ab < -1 | r_ab > 1) stop("Correlation between task a and b must be between -1 and 1")

  alternative <- match.arg(alternative)

  n = sample_size

  bsdt_p_sim <- function() {

    cor_mat <- matrix(c(1, r_ab, r_ab, 1), nrow = 2)

    cov_mat <- diag(c(sd_a, sd_b)) %*% cor_mat %*% diag(c(sd_a, sd_b))

    controls <- MASS::mvrnorm(n + 1, mu = c(mean_a, mean_b), Sigma = cov_mat)

    case_a_gen <- controls[1, 1] + case_a
    case_b_gen <- controls[1, 2] + case_b

    controls <- controls[-1, ]

    pval <- singcar::BSDT(case_a_gen, case_b_gen,
                          controls[ , 1], controls[, 2],
                          alternative = alternative, iter = iter)[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- bsdt_p_sim()

  }

  power = sum(pval < alpha)/length(pval)


  return(power)
}




#' Power calculator for BSDT_cov
#'
#' Computationally intense. Lower \code{iter} and/or \code{nsim} for faster but
#' less precise calculations. Calculates approximate power, given sample size,
#' using Monte Carlo simulation for BSDT with covariates
#' for specified (expected) case score, means and standard deviations for the
#' control sample on the task of interest and included covariates. The number of
#' covariates defaults to 1, means and standard deviations for the tasks and
#' covariate default to 0 and 1, so if no other values are given the case scores
#' is interpreted as deviation from the mean in standard deviations for both tasks
#' and covariates.
#'
#' @param case_tasks A vector of length 2. The expected case scores from the
#'   tasks of interest.
#' @param case_cov A vector containing the expected case scores on all
#'   covariates included.
#' @param control_tasks A 2x2 matrix or dataframe containing the expected means
#'   (first column) and standard deviations (second column). Defaults to two
#'   variables with means 0 and sd = 1.
#' @param control_covar A px2 matrix or dataframe containing the expected means
#'   (first column) and standard deviations (second column), p being the number
#'   of covariates. Defaults to one covariate with mean 0 and sd = 1.
#' @param cor_mat A correlation matrix containing the correlations of the tasks
#'   of interest and the coviariate(s). The first two variables are treated as
#'   the tasks of interest. Defaults to no correlation between any.
#' @param sample_size Single value giving the size of the control sample for which you wish
#'   to calculate power.
#' @param alternative The alternative hypothesis. A string of either "less",
#'   "greater" or "two.sided" (default).
#' @param alpha The specified Type I error rate, default is 0.05. This can be
#'   varied, with effects on power.
#' @param nsim The number of simulations for the power calculation. Defaults to
#'   1000 due to BSDT already being computationally intense. Increase for better
#'   accuracy.
#' @param iter The number of simulations used by the BSDT_cov, defaults to 1000.
#'   Increase for better accuracy.
#' @param calibrated Whether or not to use the standard theory (Jeffreys) prior
#'   distribution (if set to \code{FALSE}) or a calibrated prior. See Crawford
#'   et al. (2011) for further information. Calibrated prior is recommended.
#'
#' @return Returns a single value approximating the power of the test for the
#'   given parameters.
#' @export
#'
#' @examples
#' BSDT_cov_power(c(-2, 0), case_cov = c(0, 0, 0),
#' control_covar = matrix(c(0, 0, 0, 1, 1, 1), ncol= 2),
#' sample_size = 10, cor_mat = diag(5), iter = 20, nsim = 20)

BSDT_cov_power <- function(case_tasks, case_cov, control_tasks = matrix(c(0, 0, 1, 1), ncol= 2), control_covar = c(0, 1),
                          cor_mat = diag(3),
                          sample_size,
                          alternative = c("two.sided", "greater", "less"),
                          alpha = 0.05,
                          nsim = 1000, iter = 1000,
                          calibrated = TRUE) {

  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")
  if (sum(eigen(cor_mat)$values > 0) < length(diag(cor_mat))) stop("cor_mat is not positive definite")
  if (sample_size < 2) stop("Sample size must be greater than 1")
  if (length(case_tasks) != 2) stop("case_tasks should be of length 2")

  alternative <- match.arg(alternative)
  n = sample_size

  sum_stats <- rbind(control_tasks, control_covar)

  if (length(sum_stats[ , 2]) != nrow(cor_mat)) stop("Number of variables and number of correlations do not match")

  Sigma <- diag(sum_stats[ , 2]) %*% cor_mat %*% diag(sum_stats[ , 2])
  mu <- sum_stats[ , 1]

  BSDT_cov_p_sim <- function() {


    con <- MASS::mvrnorm(n+1, mu = mu, Sigma = Sigma)

    case_scores <- c(case_tasks, case_cov)
    case_score_emp <- vector(length = length(case_scores))

    for (i in 1:length(case_scores)) {
      case_score_emp[i] <- con[1, i] + case_scores[i]
    }

    con <- con [-1, ]

    pval <- singcar::BSDT_cov(case_tasks = case_score_emp[c(1, 2)], case_covar = case_score_emp[-c(1, 2)],
                             control_task = con[ , c(1, 2)], control_covar = con[ , -c(1, 2)],
                             iter = iter, alternative = alternative,
                             calibrated = calibrated)[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- BSDT_cov_p_sim()

  }

  power = sum(pval < alpha)/length(pval)

  return(power)


}
