#' Power calculator for the unstandardised difference test (UDT)
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
                      alternative = c("two.sided", "one.sided"),
                      alpha = 0.05, spec = 0.005) {

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

    if (alternative == "one.sided") {

      tstat <- abs((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

      power = stats::pt(stats::qt(alpha, df = n-1,
                                  lower.tail = F),
                        ncp = tstat,
                        df = n-1,
                        lower.tail = F)

      return(power)

    }


  }

  if (is.null(sample_size)) {

    n = 2

    search_pwr = 0

    prev_search_pwr = 1

    keep_going = TRUE

    if (alternative == "two.sided") {

      while(keep_going == TRUE){

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

        if (abs(prev_search_pwr - search_pwr) < spec) {
          keep_going = FALSE
          print(paste0("Power (", format(round(search_pwr, 5), nsmall = 5), ") will not increase more than ", spec*100, "%",
                       " for any additional participant over n = ", n))
        }

        prev_search_pwr = search_pwr

        n = n + 1

        if (search_pwr > power) keep_going = FALSE

      }

      return(data.frame(n = n - 1, power = search_pwr))
    }

    if (alternative == "one.sided") {

      while(keep_going == TRUE){

        tstat <- abs((case_a - mean_a) - (case_b - mean_b)) / sqrt((sd_a^2 + sd_b^2 - 2*sd_a*sd_b*r_ab) * ((n + 1)/n))

        search_pwr = stats::pt(stats::qt(alpha, df = n-1,
                                         lower.tail = F),
                               ncp = tstat,
                               df = n-1,
                               lower.tail = F)

        if (abs(prev_search_pwr - search_pwr) < spec) {
          keep_going = FALSE
          print(paste0("Power (", format(round(search_pwr, 5), nsmall = 5), ") will not increase more than ", spec*100, "%",
                       " for any additional participant over n = ", n))
        }

        prev_search_pwr = search_pwr

        n = n + 1

        if (search_pwr > power) keep_going = FALSE

      }

      return(data.frame(n = n - 1, power = search_pwr))

    }
    }

}


#' Power calculator for the revised standardised difference test (RSDT)
#'
#' Calculates approximate power, given sample size, using Monte Carlo
#' simulation, for specified case scores, means and standard deviations for the
#' control sample. The means and standard deviations defaults to 0 and 1
#' respectively, so if no other values are given the case scores are interpreted
#' as deviations from the mean in standard deviations. Hence, the effect size of
#' the dissociation (zdcc) would in that case be the difference between the two
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
#'            sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20)


RSDT_power <- function(case_a, case_b, mean_a, mean_b, sd_a, sd_b, r_ab,
                       sample_size,
                       alternative = c("two.sided", "one.sided"),
                       alpha = 0.05, nsim = 10000) {

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
                          alternative = "two.sided")[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- rsdt_p_sim(case_a, case_b, mean_a, mean_b,
                          sd_a, sd_b, r_ab)

  }

  power = switch(alternative,
                 two.sided = sum(pval < alpha)/length(pval),
                 one.sided = sum(pval/2 < alpha)/length(pval)
  )

  return(power)
}

#' Power calculator for the Bayesian standardised difference test (BSDT)
#'
#' Calculates approximate power, given sample size, using Monte Carlo
#' simulation, for specified case scores, means and standard deviations for the
#' control sample. The means and standard deviations defaults to 0 and 1
#' respectively, so if no other values are given the case scores are interpreted
#' as deviations from the mean in standard deviations. Hence, the effect size of
#' the dissociation (zdcc) would in that case be the difference between the two
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
#' @param alpha The specified Type I error rate. This can also be varied, with
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
#'            sd_a = 1, sd_b = 1, r_ab = 0.5, sample_size = 20)


BSDT_power <- function(case_a, case_b, mean_a, mean_b, sd_a, sd_b, r_ab,
                       sample_size,
                       alternative = c("two.sided", "one.sided"),
                       alpha = 0.05, nsim = 1000, iter = 1000) {

  alternative <- match.arg(alternative)

  n = sample_size

  bsdt_p_sim <- function(case_a, case_b,mean_a, mean_b,
                         sd_a, sd_b, r_ab) {

    cor_mat <- matrix(c(1, r_ab, r_ab, 1), nrow = 2)

    cov_mat <- diag(c(sd_a, sd_b)) %*% cor_mat %*% diag(c(sd_a, sd_b))

    controls <- MASS::mvrnorm(n + 1, mu = c(mean_a, mean_b), Sigma = cov_mat)

    case_a_gen <- controls[1, 1] + case_a
    case_b_gen <- controls[1, 2] + case_b

    controls <- controls[-1, ]

    pval <- singcar::BSDT(case_a_gen, case_b_gen,
                          controls[ , 1], controls[, 2],
                          alternative = "two.sided", iter = iter)[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- bsdt_p_sim(case_a, case_b, mean_a, mean_b,
                          sd_a, sd_b, r_ab)

  }

  power = switch(alternative,
                 two.sided = sum(pval < alpha)/length(pval),
                 one.sided = sum(pval/2 < alpha)/length(pval)
  )

  return(power)
}

