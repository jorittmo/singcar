#' Bayesian Test of Deficit
#'
#' Takes a single observation and
#' compares it to a distribution estimated by a control sample. Calculates
#' standardised difference between case and controls or proportions falling
#' above or below the case score as well as associated credible intervals.
#'
#' Returns the point estimate of the standardised difference
#' between the case score and the mean of the controls and the point estimate
#' of the p-value (i.e. the percentage of the population that would be
#' expected to obtain a lower or higher score, depending on the alternative
#' hypothesis).
#'
#'
#' @param case Case observation, can only be a single value.
#' @param controls Numeric vector of observations from the control sample. If
#'   single value, treated as mean.
#' @param controls.sd If input of controls is single value, the standard
#'   deviation of the sample must be given as well.
#' @param controls.n If input of controls is single value, the size of the
#'   sample must be given as well.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"less"} (default), \code{"greater"} or
#'   \code{"two.sided"}. You can specify just the initial letter.
#' @param int.level Level of confidence for credible intervals.
#' @param iter Number of iterations.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \tabular{llll}{
#' \code{statistic}   \tab the value of the t-statistic.\cr\cr  \code{parameter}
#' \tab the degrees of freedom for the t-statistic.\cr\cr \code{p.value}    \tab
#' the p-value for the test.\cr\cr \code{estimate}    \tab estimated
#' standardised difference (zcc) and point estimate of p-value. \cr\cr
#' \code{null.value}   \tab the value of the difference under the null
#' hypothesis.\cr\cr \code{interval}     \tab named numerical vector containing
#' credibility level and intervals for both zcc and estimated proportion. \cr\cr
#' \code{desc}     \tab named numerical containing descriptive statistics: mean
#' and standard deviations of controls as well as sample size.
#'  \cr\cr \code{alternative}     \tab a character string
#' describing the alternative hypothesis.\cr\cr \code{method} \tab a character
#' string indicating what type of t-test was performed.\cr\cr \code{data.name}
#' \tab a character string giving the name(s) of the data as well as sum
#' summaries. }
#'
#' @export
#'
#' @examples
#' BTD(case = -2, controls = 0, controls.sd = 1, controls.n = 20)
#'
#' @references
#' Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a
#' control or normative sample in neuropsychology: Development of a Bayesian
#' approach. \emph{Cognitive Neuropsychology, 24}(4), 343–372.
#' https://doi.org/10.1080/02643290701290146





BTD <- function (case, controls, controls.sd = NULL, controls.n = NULL,
                alternative = c("less", "greater", "two.sided"),
                int.level = 0.95, iter = 1000, na.rm = FALSE) {

  if (length(case)>1) stop("Case should only have 1 observation")
  if (length(controls)<2 & is.null(controls.sd) == TRUE) {
    stop("Not enough obs. Set sd and n for input of controls to be treated as mean")
  }

  if (length(controls)<2 & is.null(controls.sd) == FALSE & is.null(controls.n) == TRUE) stop("Input sample size")
  if (is.na(case)==TRUE) stop("Case is NA")

  if (na.rm == TRUE) controls <- controls[!is.na(controls)]
  if (sum(is.na(controls)) > 0) stop("Controls contains NA, set na.rm = TRUE to proceed")

  if (int.level < 0 | int.level > 1) stop("Interval level must be between 0 and 1")

  alternative <- match.arg(alternative)

  con_m <- mean(controls) # Mean of the control sample

  con_sd <- stats::sd(controls) # Standard deviation of the control sample (returns NA if summary stats used)
  if (length(controls)<2 & is.null(controls.sd) == FALSE) con_sd <- controls.sd

  n <- length(controls)
  if (length(controls)<2 & is.null(controls.sd) == FALSE & is.null(controls.n) == FALSE) n <- controls.n

  df <- n - 1 # The degrees of freedom
  alpha <- 1 - int.level

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
  z_ast_est <- mean(z_ast)

  zcc_int <- stats::quantile(z_ast, c(alpha/2, (1 - alpha/2)))
  names(zcc_int) <- c("Lower zcc CI", "Upper zcc CI")

  p_est <- mean(pval)

  p_int <- stats::quantile(pval, c(alpha/2, (1 - alpha/2)))*100
  if (alternative == "two.sided") p_int <- stats::quantile(pval/2, c(alpha/2, (1 - alpha/2)))*100
  names(p_int) <- c("Lower p CI", "Upper p CI")

  estimate <- c(zcc, p_est*100)
  if (alternative == "two.sided") estimate <- c(zcc, (p_est/2)*100)

  zcc.name <- paste0("Std. case difference (Z-CC), ",
                     100*int.level, "% credible interval [",
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
                   100*int.level, "% credible interval [",
                   format(round(p_int[1], 2), nsmall = 2),", ",
                   format(round(p_int[2], 2), nsmall = 2),"]")

  names(estimate) <- c(zcc.name, p.name)


  typ.int <- 100*int.level
  names(typ.int) <- "Interval level (%)"
  interval <- c(typ.int, zcc_int, p_int)

  names(z_ast_est) <- "est. z"
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between case and controls"
  names(con_m) <- "Mean (controls)"
  names(con_sd) <- "SD (controls)"
  names(n) <- "Sample size"


  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = z_ast_est,
                 parameter = df,
                 p.value = p_est,
                 estimate = estimate,
                 null.value = null.value,
                 interval = interval,
                 desc = c(con_m, con_sd, n),
                 alternative = alternative,
                 method = paste("Bayesian Test of deficit by Crawford and Garthwaite (2007)"),
                 data.name = paste0("case = ", format(round(case, 2), nsmall = 2),
                                    " and controls (M = ", format(round(con_m, 2), nsmall = 2),
                                    ", SD = ", format(round(con_sd, 2), nsmall = 2),
                                    ", N = ", n, ")"))

  class(output) <- "htest"
  output

}