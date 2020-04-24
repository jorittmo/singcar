#' Crawford and Howell's (1998) modified t-test
#'
#' Takes a single observation and compares it to a distribution estimated by a
#' control sample. Calculates standardised difference between case and controls
#' or proportions falling above or below the case score as well as associated
#' confidence intervals.
#'
#' Returns either the point estimate of the standardised difference
#' between the case score and the mean of the controls or the point estimate
#' of the p-value (i.e. the percentage of the population that would be
#' expected to obtain a lower or higher score, depending on the alternative
#' hypothesis).
#'
#' @section Note of caution:
#' Calculating the confidence intervals relies on finding non-centrality
#' parameters for non-central t-distributions. Depending on the degrees of
#' freedom, the confidence level and the effect size exact accuracy from the
#' \code{stats::qt()} function used can not be guaranteed. However, the
#' approximations should be good enough for most cases.
#' See \url{https://stat.ethz.ch/pipermail/r-help/2008-June/164843.html}.
#'
#' @param case Case observation, can only be a single value.
#' @param controls Numeric vector of observations from the control sample.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"less"} (default), \code{"greater"} or
#'   \code{"two.sided"}. You can specify just the initial letter.
#' @param estimate.p Set to \code{TRUE} to estimate the proportion of population
#'   falling above or below case score instead of standardised difference.
#' @param conf.int Calculate confidence intervals for desired estimate. Uses
#'   iterative method, set to \code{FALSE} for faster calculation (e.g. for
#'   simulations).
#' @param conf.level Level of confidence for intervals.
#' @param conf.int.spec The size of iterative steps for calculating confidence
#'   intervals. Smaller values gives more precise intervals but takes longer to
#'   calculate.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \tabular{llll}{
#' \code{statistic}   \tab the value of the t-statistic.\cr\cr  \code{parameter}
#' \tab the degrees of freedom for the t-statistic.\cr\cr \code{p.value}    \tab
#' the p-value for the test.\cr\cr \code{estimate}    \tab estimated standardised
#' difference or proporiton as well as mean and sd of controls.\cr\cr
#' \code{null.value}   \tab the value of the difference under the null
#' hypothesis.\cr\cr \code{conf.int}     \tab a confidence interval for the
#' specified estimate and alternative hypothesis.\cr\cr \code{interval}     \tab
#' string indicating the type of interval (proportion or difference). \cr\cr
#' \code{stderr}     \tab the standard error of the mean (difference), used as
#' denominator in the t-statistic formula.\cr\cr \code{alternative}     \tab a
#' character string describing the alternative hypothesis.\cr\cr \code{method}
#' \tab a character string indicating what type of t-test was performed.\cr\cr
#' \code{data.name}     \tab a character string giving the name(s) of the data
#' as well as sum summaries.
#' }
#'
#' @export
#'
#' @examples
#' ch.ttest(-2, rnorm(15), alternative = "l", estimate.p = TRUEuse)


ch.ttest <- function (case, controls, alternative = c("less", "greater", "two.sided"), estimate.p = FALSE,
                      conf.int = TRUE, conf.level = 0.95, conf.int.spec = 0.01, na.rm = FALSE) {

  if (length(case)>1) stop("Case should only have 1 observation")
  if (length(controls)<2) stop("Controls do not have enough observations")
  if (is.na(case)==TRUE) stop("Case is NA")

  if (na.rm == TRUE) controls <- controls[!is.na(controls)]
  if (sum(is.na(controls)) > 0) stop("controls contains NA, set na.rm = TRUE to proceed")

  alternative <- match.arg(alternative)

  con_m <- mean(controls) # Mean of the control sample
  con_sd <- stats::sd(controls) # Standard deviation of the control sample
  n <- length(controls)
  stderr <- con_sd * sqrt((n + 1)/n) # Standard error by C&H (1998) method

  std_fx <- (case-con_m)/con_sd
  tstat <- (case - con_m)/stderr # C&H's modified t
  df <- length(controls) - 1 # The degrees of freedom

  # Get p-value depending on alternative hypothesis

  if (alternative == "less") {

    pval <- stats::pt(tstat, df = df)

  } else if (alternative == "greater") {

    pval <- stats::pt(tstat, df = df, lower.tail = FALSE)

  } else if (alternative == "two.sided") {

    pval <- 2 * (1 - stats::pt(abs(tstat), df = df))

  }

  # Calculate the CIs with method described in Crawford and Garthwaite (2002) and Cumming and Finch (2001)

  # Below is a search algorithm to find the non-centrality parameter of two non-central t-distributions
  # which have their alpha/2 and 1-alpha/2 percentile at the std effect size * sqrt(n), respectively.
  # These non-centrality paramters / sqrt(n) are then taken as the limits of the CIs.

  if (conf.int == T) {

    alph <- 1 - conf.level

    stop_ci_lo <- FALSE
    ncp_lo <- std_fx*sqrt(n)
    perc_lo <- 1 - (alph/2)
    while (stop_ci_lo == FALSE) {

      # Here we search downwards with each step being as big as specified in conf.int.spec
      ncp_lo <- ncp_lo - conf.int.spec

      suppressWarnings( # Depending on ncp and percentile qt gives approximations, which produces warnings
        quant <- stats::qt(perc_lo, df = df, ncp = ncp_lo)
      )

      if (quant <= std_fx*sqrt(n)) { # Wen the specified quantile reaches std_fx*sqrt(n) the search stops
        stop_ci_lo <- TRUE
      }
    }

    stop_ci_up <- FALSE
    ncp_up <- std_fx*sqrt(n)
    perc_up <- (alph/2)
    while (stop_ci_up == FALSE) {

      # Here we search upwards with each step being as big as specified in conf.int.spec
      ncp_up <- ncp_up + conf.int.spec

      suppressWarnings( # Depending on ncp and percentile qt gives approximations, which produces warnings
        quant <- stats::qt(perc_up, df = df, ncp = ncp_up)
      )

      if (quant >= std_fx*sqrt(n)) { # Wen the specified quantile reaches std_fx*sqrt(n) the search stops
        stop_ci_up <- TRUE
      }
    }

    ci_lo <- ncp_lo/sqrt(n)
    ci_up <- ncp_up/sqrt(n)
    cint <- c(ci_lo, ci_up)

    if (estimate.p == TRUE) {
      if (alternative == "less") {

        ci_lo <- stats::pnorm(ci_lo)*100
        ci_up <- stats::pnorm(ci_up)*100
        cint <- c(ci_lo, ci_up)

      } else if (alternative == "greater") {

        ci_lo <- (1 - stats::pnorm(ci_lo))*100
        ci_up <- (1 - stats::pnorm(ci_up))*100

        # NOTE (!): Because of right side of dist, lower and upper CI must switch to
        # be consistent with lower CI to the left and upper to the right in output
        cint <- c(ci_up, ci_lo)

      } else {
        if (tstat < 0) {

          ci_lo <- stats::pnorm(ci_lo)*100
          ci_up <- stats::pnorm(ci_up)*100
          cint <- c(ci_lo, ci_up)

        } else {

          ci_lo <- (1 - stats::pnorm(ci_lo))*100
          ci_up <- (1 - stats::pnorm(ci_up))*100

          # NOTE (!): Because of right side of dist, lower and upper CI must switch to
          # be consistent with lower CI to the left and upper to the right in output
          cint <- c(ci_up, ci_lo)
        }
      }

    }

    attr(cint,"conf.level") <- conf.level # Give the CIs an attribute based in specified conf level

  } else {
    cint <- NULL # If conf.int set to FALSE no intervals will be calculated or given in output
  }

  if (estimate.p == FALSE) {
    # Set the estimates and names of estimates given in output
    estimate <- c(std_fx, con_m, con_sd)
    names(estimate) <- c("Standardised case difference", "mean (controls)", "SD (controls)")
  } else {

    estimate <- c(pval*100, con_m, con_sd)

    if (alternative == "two.sided") estimate <- c((pval/2)*100, con_m, con_sd)

    if (alternative == "less") {
      names(estimate) <- c("Proportion below case (%)", "mean (controls)", "SD (controls)")
    } else if (alternative == "greater") {
      names(estimate) <- c("Proportion above case (%)", "mean (controls)", "SD (controls)")
    } else {
      names(estimate) <- c(paste("Proportion", ifelse(tstat < 0, "below", "above"), "case (%)"),
                           "mean (controls)", "SD (controls)")
    }

  }

  # Set names for objects in output
  names(tstat) <- "t"
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between case and controls"


  # Set name for estimates in output
  if (estimate.p == F) {
    interval <- "Interval for standardised difference"
  } else {
    interval <- "Interval for proportion (%)"
  }


  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = tstat, parameter = df, p.value = pval,
                 estimate = estimate, null.value = null.value,
                 conf.int = cint,
                 interval = interval,
                 stderr = stderr,
                 alternative = alternative,
                 method = paste("Crawford-Howell (1998) t-test"),
                 data.name = paste0("case = ", format(round(case, 2), nsmall = 2),
                                    " and controls (M = ", format(round(con_m, 2), nsmall = 2),
                                    ", SD = ", format(round(con_sd, 2), nsmall = 2), ")"))

  class(output) <- "htest"
  output

}
