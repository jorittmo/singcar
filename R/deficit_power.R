#' Power calculator for TD
#'
#' Calculates exact power given sample size or sample size given power, using
#' analytical methods for the frequentist test of deficit for a specified case
#' score and mean and standard deviation for the control sample. The mean and
#' standard deviation defaults to 0 and 1, so if no other values are given the
#' case score is interpreted as deviation from the mean in standard deviations.
#'
#' @param case A single value from the expected case observation.
#' @param mean The expected mean of the control sample.
#' @param sd The expected standard deviation of the control sample.
#' @param sample_size The size of the control sample, vary this parameter to see
#'   how the sample size affects power. One of sample size or power must be
#'   specified, not both.
#' @param power A single value between 0 and 1 specifying desired power for
#'   calculating necessary sample size. One of sample size or power must be
#'   specified, not both.
#' @param alternative The alternative hypothesis. A string of either "less" (default),
#'   "greater" or "two.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power.
#' @param spec A single value between 0 and 1. If desired power is given as
#'   input the function will utilise a search algorithm to find the sample size
#'   needed to reach the desired power. However, if the power specified is
#'   greater than what is actually possible to achieve the algorithm could
#'   search forever. Hence, when power does not increase substantially for
#'   any additional participant in the sample, the algorithm stops.
#'   By default the algorithm stops when power does not increase more
#'   than 0.5\% for any added participant, but by varying \code{spec},
#'   this specificity can be changed.
#'
#' @return Either a single value of the exact power, if sample size is given. Or
#'   a dataframe consisting of both the sample size and the exact power such
#'   size would yield.
#' @export
#'
#' @examples
#' TD_power(case = -2, mean = 0, sd = 1, sample_size = 20)
#' TD_power(case = -2, mean = 0, sd = 1, power = 0.8)

TD_power <- function(case, mean = 0, sd = 1,
                     sample_size = NULL,
                     power = NULL,
                     alternative = c("less", "greater", "two.sided"),
                     alpha = 0.05, spec = 0.005) {

  if (!is.null(sample_size) & !is.null(power)) stop("Must supply only one of sample size or desired power")
  if (is.null(sample_size) & is.null(power)) stop("Must supply either sample size or desired power")
  if (!is.null(power)) if (power > 1 | power < 0) stop("Desired power must be between 0 and 1")
  if (!is.null(sample_size)) if (sample_size < 2) stop("Sample size must be greater than 1")
  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")


  alternative <- match.arg(alternative)
  n = sample_size

  if (is.null(power)) {

    if (alternative == "two.sided") {

      power = stats::pt(stats::qt(alpha/2, df = n-1,
                                  lower.tail = T),
                        ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                        df = n-1,
                        lower.tail = T) - stats::pt(-stats::qt(alpha/2, df = n-1,
                                                               lower.tail = T),
                                                    ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                                                    df = n-1,
                                                    lower.tail = T) + 1

      return(power)

    }

    if (alternative == "less") {

      power = stats::pt(stats::qt(alpha, df = n-1,
                                  lower.tail = T),
                        ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                        df = n-1,
                        lower.tail = T
      )



      return(power)

    }

    if (alternative == "greater") {

      power = stats::pt(stats::qt(alpha, df = n-1,
                                  lower.tail = F),
                        ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
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

        search_pwr = stats::pt(stats::qt(alpha/2, df = n-1,
                                         lower.tail = T),
                               ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                               df = n-1,
                               lower.tail = T) - stats::pt(-stats::qt(alpha/2, df = n-1,
                                                                      lower.tail = T),
                                                           ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                                                           df = n-1,
                                                           lower.tail = T) + 1

      }

      if (alternative == "less") {

        search_pwr = stats::pt(stats::qt(alpha, df = n-1,
                                    lower.tail = T),
                          ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                          df = n-1,
                          lower.tail = T
        )

      }

      if (alternative == "greater") {

        search_pwr = stats::pt(stats::qt(alpha, df = n-1,
                                    lower.tail = F),
                          ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                          df = n-1,
                          lower.tail = F
        )


      }


      if (abs(prev_search_pwr - search_pwr) < spec) {
        keep_going = FALSE
        message(paste0("Power (", format(round(search_pwr, 3), nsmall = 3), ") will not increase more than ", spec*100, "%",
                     " per participant for n > ", n))
      }

      prev_search_pwr = search_pwr

      n = n + 1

      if (search_pwr > power) keep_going = FALSE

    }

    return(data.frame(n = n - 1, power = search_pwr))
  }
}





#' Power calculator for BTD
#'
#' Calculates approximate power, given sample size, using Monte Carlo simulation for the
#' Bayesian test of deficit for a specified case score, mean and standard
#' deviation for the control sample. The mean and standard deviation defaults to
#' 0 and 1, so if no other values are given the case score is interpreted as
#' deviation from the mean in standard deviations.
#'
#' @param case A single value from the expected case observation.
#' @param mean The expected mean of the control sample.
#' @param sd The expected standard deviation of the control sample.
#' @param sample_size The size of the control sample, vary this parameter to see
#'   how the sample size affects power.
#' @param alternative The alternative hypothesis. A string of either "less" (default),
#'   "greater" or "two.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power.
#' @param nsim The number of simulations for the power calculation. Defaults to
#'   1000 due to BTD already being computationally intense.
#' @param iter The number of simulations used by the BTD. Defaults to 1000.
#'
#' @return Returns a single value approximating the power of the test for the
#'   given parameters.
#' @export
#'
#' @examples
#' BTD_power(case = -2, mean = 0, sd = 1, sample_size = 20)

BTD_power <- function(case, mean = 0, sd = 1,
                      sample_size,
                      alternative = c("less", "greater", "two.sided"),
                      alpha = 0.05,
                      nsim = 1000, iter = 1000) {

  if (sample_size < 2) stop("Sample size must be greater than 1")
  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")



  alternative <- match.arg(alternative)
  n = sample_size


  BTD_p_sim <- function() {


    con <- stats::rnorm(n, mean = mean, sd = sd)

    case <- stats::rnorm(1, mean = mean, sd = sd) + (case - mean)

    pval <- singcar::BTD(case, con, iter = iter, alternative = alternative)[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- BTD_p_sim()

  }

  power = sum(pval < alpha)/length(pval)

  return(power)

}



#' Power calculator for BTD_cov
#'
#' Computationally intense. Lower \code{iter} and/or \code{nsim} for less exact
#' but faster calculations. Calculates approximate power, given sample size,
#' using Monte Carlo simulation for the Bayesian test of deficit with covariates
#' for specified (expected) case score, means and standard deviations for the
#' control sample on the task of interest and included covariates. The number of
#' covariates defaults to 1, means and standard deviations for the task and
#' covariate defaults to 0 and 1, so if no other values are given the case score
#' is interpreted as deviation from the mean in standard deviations for both task
#' and covariate.
#'
#' @param case A single value from the expected case observation on the task of
#'   interest.
#' @param case_cov A vector of expected case observations from covariates of
#'   interest.
#' @param control_task A vector of length 2 containing the expected mean and standard
#'   deviation of the task of interest. In that order.
#' @param control_covar A matrix with 2 columns containing expected means (in the 1st
#'   column) and standard deviations (in the 2nd column) of the included
#'   covariates.
#' @param cor_mat A correlation matrix containing the correlations of the
#'   task of interest and the coviariate(s). The first variable is treated as
#'   the task of interest. Defaults to no correlation between any.
#' @param sample_size Single value of the size of the sample for which you wish
#'   to calculate power.
#' @param alternative The alternative hypothesis. A string of either "less" (default),
#'   "greater" or "two.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power.
#' @param nsim The number of simulations for the power calculation. Defaults to
#'   1000 due to BTD_cov already being computationally intense.
#' @param iter The number of simulations used by the BTD_cov. Defaults to 1000.
#'
#' @return Returns a single value approximating the power of the test for the
#'   given parameters.
#' @export
#'
#' @examples
#' cor_mat = matrix(c(1, 0.2, 0.3, 0.2, 1, 0.4, 0.3, 0.4, 1), ncol = 3)
#'
#' BTD_cov_power(case = -2, case_cov = c(105, 30), control_task = c(0, 1),
#' control_covar = matrix(c(100, 40, 15, 10), ncol = 2), sample_size = 15,
#' cor_mat = cor_mat, iter = 20, nsim = 20)

BTD_cov_power <- function(case, case_cov, control_task = c(0, 1), control_covar = c(0, 1),
                          cor_mat = diag(2),
                          sample_size,
                          alternative = c("less", "greater", "two.sided"),
                          alpha = 0.05,
                          nsim = 1000, iter = 1000) {

  if (alpha < 0 | alpha > 1) stop("Type I error rate must be between 0 and 1")

  if (sum(eigen(cor_mat)$values > 0) < length(diag(cor_mat))) stop("cor_mat is not positive definite")

  if (sample_size < 2) stop("Sample size must be greater than 1")

  alternative <- match.arg(alternative)
  n = sample_size

  sum_stats <- rbind(control_task, control_covar)

  if (length(sum_stats[ , 2]) != nrow(cor_mat)) stop("Number of variables and number of correlations do not match")

  Sigma <- diag(sum_stats[ , 2]) %*% cor_mat %*% diag(sum_stats[ , 2])
  mu <- sum_stats[ , 1]

  BTD_cov_p_sim <- function() {


    con <- MASS::mvrnorm(n+1, mu = mu, Sigma = Sigma)


    case_scores <- c(case, case_cov)
    case_score_emp <- vector(length = length(case_scores))

    for (i in 1:length(case_scores)) {
      case_score_emp[i] <- con[1, i] + (case_scores[i] - mu[i])
    }

    con <- con [-1, ]

    pval <- singcar::BTD_cov(case_task = case_score_emp[1], case_covar = case_score_emp[-1],
                             control_task = con[ , 1], control_covar = con[ , -1],
                             iter = iter, alternative = alternative)[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- BTD_cov_p_sim()

  }

  power = sum(pval < alpha)/length(pval)

  return(power)

}
