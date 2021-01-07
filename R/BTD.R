#' Bayesian Test of Deficit
#'
#' Takes a single observation and compares it to a distribution estimated by a
#' control sample using Bayesian methodology. Calculates standardised difference
#' between the case score and the mean of the controls and proportions falling
#' above or below the case score, as well as associated credible intervals. This
#' approach was developed by Crawford and Garthwaite (2007) but converge to the
#' results of \code{\link{TD}()}, which is faster. Returns the point estimate of
#' the standardised difference between the case score and the mean of the
#' controls and the point estimate of the p-value (i.e. the percentage of the
#' population that would be expected to obtain a lower or higher score,
#' depending on the alternative hypothesis). This test is based on random number
#' generation which means that results may vary between runs. This is by design
#' and the reason for not using \code{set.seed()} to reproduce results inside
#' the function is to emphasise the randomness of the test. To get more accurate
#' and stable results please increase the number of iterations by increasing
#' \code{iter} whenever feasible.
#'
#'
#' @param case Case observation, can only be a single value.
#' @param controls Numeric vector of observations from the control sample. If
#'   single value, treated as mean.
#' @param sd If input of controls is single value, the standard
#'   deviation of the sample must be given as well.
#' @param sample_size If input of controls is single value, the size of the
#'   sample must be given as well.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"less"} (default), \code{"greater"} or
#'   \code{"two.sided"}. You can specify just the initial letter.
#' @param int_level Level of confidence for credible intervals, defaults to 95\%.
#' @param iter Number of iterations. Set to higher for more accuracy, set to
#'   lower for faster calculations.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \tabular{llll}{
#' \code{statistic}   \tab the mean z-value over \code{iter} number of
#' iterations \cr\cr  \code{parameter} \tab the degrees of freedom used to
#' specify the posterior distribution. \cr\cr \code{p.value}    \tab the mean p-value
#' for all simulated Z-scores.\cr\cr \code{estimate}    \tab estimated standardised difference
#' (Z-CC) and point estimate of p-value. \cr\cr \code{null.value}   \tab the
#' value of the difference under the null hypothesis.\cr\cr \code{interval}
#' \tab named numerical vector containing credibility level and intervals for
#' both Z-CC and estimated proportion. \cr\cr \code{desc}     \tab named
#' numerical containing descriptive statistics: mean and standard deviations of
#' controls as well as sample size. \cr\cr \code{alternative}     \tab a
#' character string describing the alternative hypothesis.\cr\cr \code{method}
#' \tab a character string indicating what type of test was performed.\cr\cr
#' \code{data.name} \tab a character string giving the name(s) of the data as
#' well as summaries. }
#'
#' @export
#'
#' @examples
#' BTD(case = -2, controls = 0, sd = 1, sample_size = 20, iter = 1000)
#'
#' BTD(case = size_weight_illusion[1, "V_SWI"],
#'     controls = size_weight_illusion[-1, "V_SWI"], alternative = "l", iter = 1000)
#'
#' @references
#'
#' Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a
#' control or normative sample in neuropsychology: Development of a Bayesian
#' approach. \emph{Cognitive Neuropsychology, 24}(4), 343-372.
#' \doi{10.1080/02643290701290146}





BTD <- function (case, controls, sd = NULL, sample_size = NULL,
                alternative = c("less", "greater", "two.sided"),
                int_level = 0.95, iter = 10000, na.rm = FALSE) {

  if (length(case)>1) stop("Case should only have 1 observation")
  if (length(controls)<2 & is.null(sd) == TRUE) {
    stop("Not enough obs. Set sd and n for input of controls to be treated as mean")
  }

  if (length(controls)<2 & is.null(sd) == FALSE & is.null(sample_size) == TRUE) stop("Input sample size")
  if (is.na(case)==TRUE) stop("Case is NA")

  if (na.rm == TRUE) controls <- controls[!is.na(controls)]
  if (sum(is.na(controls)) > 0) stop("Controls contains NA, set na.rm = TRUE to proceed")

  if (int_level < 0 | int_level > 1) stop("Interval level must be between 0 and 1")

  alternative <- match.arg(alternative)

  con_m <- mean(controls) # Mean of the control sample

  con_sd <- stats::sd(controls) # Standard deviation of the control sample (returns NA if summary stats used)
  if (length(controls)<2 & is.null(sd) == FALSE) con_sd <- sd

  n <- length(controls)
  if (length(controls)<2 & is.null(sd) == FALSE & is.null(sample_size) == FALSE) n <- sample_size

  df <- n - 1 # The degrees of freedom
  alpha <- 1 - int_level

  ## BAYESIAN PROCESS AS EXPLAINED BY CRAWFORD AND GARTHWAITE (2007)##

  # By their notation I will use theta for the unknown variance instead of sigma

  theta_hat <- ((n - 1)*con_sd^2) /  stats::rchisq(iter, df = df)

  z <- stats::rnorm(iter)

  mu_hat <- con_m + (z * sqrt(theta_hat/n))

  z_ast <- (case - mu_hat)/sqrt(theta_hat)


  if (alternative == "less") {

    pval <- stats::pnorm(z_ast)

  } else if (alternative == "greater") {

    pval <- stats::pnorm(z_ast, lower.tail = FALSE)

  } else if (alternative == "two.sided") {

    pval <- 2 * stats::pnorm(abs(z_ast), lower.tail = FALSE)

  }

  zcc <- (case - con_m)/con_sd


  zcc_int <- stats::quantile(z_ast, c(alpha/2, (1 - alpha/2)))
  names(zcc_int) <- c("Lower Z-CC CI", "Upper Z-CC CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100
  names(p_int) <- c("Lower p CI", "Upper p CI")

  estimate <- c(zcc, p_est*100)
  if (alternative == "two.sided") estimate <- c(zcc, (p_est/2)*100)

  zcc.name <- paste0("Std. case score (Z-CC), ",
                     100*int_level, "% CI [",
                     format(round(zcc_int[1], 2), nsmall = 2),", ",
                     format(round(zcc_int[2], 2), nsmall = 2),"]")

  if (alternative == "less") {
    alt.p.name <- "Proportion below case (%), "
  } else if (alternative == "greater") {
    alt.p.name <- "Proportion above case (%), "
  } else {
    if (zcc < 0) {
      alt.p.name <- "Proportion below case (%), "
    } else {
      alt.p.name <- "Proportion above case (%), "
    }
  }

  p.name <- paste0(alt.p.name,
                   100*int_level, "% CI [",
                   format(round(p_int[1], 2), nsmall = 2),", ",
                   format(round(p_int[2], 2), nsmall = 2),"]")

  names(estimate) <- c(zcc.name, p.name)


  typ.int <- 100*int_level
  names(typ.int) <- "Interval level (%)"
  interval <- c(typ.int, zcc_int, p_int)

  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between case and controls"
  names(con_m) <- "Mean (controls)"
  names(con_sd) <- "SD (controls)"
  names(n) <- "Sample size"


  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(parameter = df,
                 p.value = p_est,
                 estimate = estimate,
                 null.value = null.value,
                 interval = interval,
                 desc = c(con_m, con_sd, n),
                 alternative = alternative,
                 method = paste("Bayesian Test of Deficit"),
                 data.name = paste0("Case = ", format(round(case, 2), nsmall = 2),
                                    ", Controls (m = ", format(round(con_m, 2), nsmall = 2),
                                    ", sd = ", format(round(con_sd, 2), nsmall = 2),
                                    ", n = ", n, ")"))

  class(output) <- "htest"
  output

}
