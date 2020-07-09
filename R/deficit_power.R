#' Power calculator for the test of deficit
#'
#' Calculates exact power given sample size or sample size given power, using
#' analytical methods for the frequentist test of deficit for a specified case
#' score and, mean and standard deviation for the control sample. The mean and
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
#' @param alternative The alternative hypothesis. A string of either "two.sided"
#'   (default) or "one.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power.
#' @param spec A single value between 0 and 1. If desired power is given as
#'   input the function will utilise a search algorithm to find the sample size
#'   needed to reach the desired power. However, if the power specified is
#'   greater than what is actually possible to achieve the algorithm could
#'   search forever. Hence, when power does not increase substantially for
#'   any additional participant in the sample, the algorthm stops.
#'   By default the algorithm stops when power does not increase more
#'   than 0.5% for any added participant, but by varying \code{spec},
#'   this specificity can be changed.
#'
#' @return Either a single value of the exact power, if sample size is given. Or
#'   a dataframe consisting of both the sample size and the exact power such
#'   size would yield.
#' @export
#'   TD_power(case = )
#'
#' @examples

TD_power <- function(case, mean = 0, sd = 1,
                     sample_size = NULL,
                     power = NULL,
                     alternative = c("one.sided", "two.sided"),
                     alpha = 0.05, spec = 0.005) {

  if (!is.null(sample_size) & !is.null(power)) stop("Must supply only one of sample size or desired power")
  if (is.null(sample_size) & is.null(power)) stop("Must supply either sample size or desired power")

  alternative <- match.arg(alternative)
  n = sample_size

  if (is.null(power)) {

    if (alternative == "two.sided") {

      power = pt(qt(alpha/2, df = n-1,
                    lower.tail = T),
                 ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                 df = n-1,
                 lower.tail = T) - pt(-qt(alpha/2, df = n-1,
                                          lower.tail = T),
                                      ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                                      df = n-1,
                                      lower.tail = T) + 1

      return(power)

    }

    if (alternative == "one.sided") {

      power = pt(qt(0.05, df = n-1,
                    lower.tail = T),
                 ncp = ((case - mean)/(sd*sqrt((n+1)/n))),
                 df = n-1,
                 lower.tail = T
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

      search_pwr = pt(qt(0.05, df = n-1,
                         lower.tail = T),
                      ncp = ((zcc)/(sigma*sqrt((n+1)/n))),
                      df = n-1,
                      lower.tail = T
      )

      if (abs(prev_search_pwr - search_pwr) < spec) {
        keep_going = FALSE
        print(paste0("Power(", format(round(search_pwr, 5), nsmall = 5), ") will not increase more than ", spec*100, "%",
                     " for any additional participant over n = ", n))
      }

      prev_search_pwr = search_pwr

      n = n + 1

      if (search_pwr > power) keep_going = FALSE

    }

    return(data.frame(n = n - 1, power = search_pwr))
  }
}





#' Power calculator for Bayesian test of deficit
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
#' @param alternative The alternative hypothesis. A string of either "two.sided"
#'   (default) or "one.sided".
#' @param alpha The specified Type I error rate. This can also be varied, with
#'   effects on power.
#' @param nsim The number of simulations for the power calculation. Defaults to
#'   1000 due to BTD already being computationally intense. Defaults to 1000.
#' @param iter The number of simulations used by the BTD.
#'
#' @return Returns a single value approximating the power of the test for the
#'   given parameters.
#' @export
#'
#' @examples

BTD_power <- function(case, mean = 0, sd = 1,
                      sample_size,
                      alternative = c("one.sided", "two.sided"),
                      alpha = 0.05,
                      nsim = 1000, iter = 1000) {

  if (is.null(sample_size) & is.null(power)) stop("Must supply either sample size or desired power")

  alternative <- match.arg(alternative)
  n = sample_size


  BTD_p_sim <- function(case, mean, sd, n) {


    con <- rnorm(n, mean = mean, sd = sd)

    case <- rnorm(1, mean = mean, sd = sd) + case

    pval <- singcar::BTD(case, con, iter = iter, alternative = "two.sided")[["p.value"]]

    pval
  }


  pval <- vector(length = nsim)

  for(i in 1:nsim) {

    pval[i] <- BTD_p_sim(case = case, mean = mean, sd = sd,  n = n)

  }


  power = switch(alternative,
                 two.sided = sum(pval < alpha)/length(pval),
                 one.sided = sum(pval/2 < alpha)/length(pval)
  )

  return(power)

}

