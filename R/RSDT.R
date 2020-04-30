#' Revised Standardised Difference Test
#'
#' A test used to determine dissociation between two functions as assessed by
#' two tasks.Takes two single values from a case (from two tasks) and compares
#' them to two a distribution of scores for both tasks estimated by a control
#' sample. Calculates a standardised effects size of task difference as well as
#' a point estimate of the proportion of the control population that would be
#' expected to show a more extreme task difference.
#'
#' @param case.x Case's score on task X.
#' @param case.y Case's score on task Y.
#' @param controls.x Controls' scores on task X. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can supply
#'   a vector as input for task X while mean and SD for task Y.
#' @param controls.y Controls' scores on task Y. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can supply
#'   a vector as input for task Y while mean and SD for task X.
#' @param controls.x.sd If single value for task X is given as input you must
#'   supply the standard deviation of the sample.
#' @param controls.y.sd If single value for task Y is given as input you must
#'   supply the standard deviation of the sample.
#' @param controls.n If X or Y is given as mean and SD you must supply the
#'   sample size. If controls.x is given as vector and controls.y as mean and
#'   SD, controls.n must equal the number of observations in controls.x.
#' @param cor.x.y If X or Y is given as mean and SD you must supply the
#'   correlation between the tasks.
#' @param alpha Chosen risk of Type I errors.
#' @param exact.method Method for deriving the test statistic. Should rarely be
#'   changed. See Crawford and Garthwaite (2005) for further information.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab the value of the t-statistic.\cr\cr
#'   \code{parameter} \tab the degrees of freedom for the t-statistic.\cr\cr
#'   \code{p.value}    \tab the p-value for the test.\cr\cr \code{estimate}
#'   \tab case scores expressed as z-scores on task X and Y. Standardised effect
#'   size (Z-DCC) of task difference between case and controls and
#'   point estimate of the proportion of the control population estimated to show a
#'   more extreme task difference. \cr\cr \code{sample.size}   \tab the size of
#'   the control sample\cr\cr \code{null.value}   \tab the value of the
#'   difference under the null hypothesis.\cr\cr  \code{alternative}     \tab a
#'   character string describing the alternative hypothesis.\cr\cr \code{method}
#'   \tab a character string indicating what type of test was performed.\cr\cr
#'   \code{data.name}     \tab a character string giving the name(s) of the data}
#' @export
#'
#' @examples
#' RSDT(-3.857, -1.875, controls.x = 0, controls.y = 0, controls.x.sd = 1,
#' controls.y.sd = 1, controls.n = 20, cor.x.y = 0.68)
#' RSDT(-3.857, -1.875, controls.x = rnorm(20), controls.y = rnorm(20))


RSDT <- function (case.x, case.y, controls.x, controls.y,
                  controls.x.sd = NULL, controls.y.sd = NULL,
                  controls.n = NULL, cor.x.y = NULL,
                  alpha = 0.05, exact.method = T) {

  if (length(case.x) > 1 | length(case.y) > 1) stop("Case scores should be single value")

  if (length(controls.x) > 1 & length(controls.y) > 1) {
    if (length(controls.x) != length(controls.y)) stop("Sample sizes must be equal")
  }

  if (length(controls.x) > 1 & is.null(controls.x.sd) == FALSE) message("Value on controls.x.sd will be ignored")
  if (length(controls.y) > 1 & is.null(controls.y.sd) == FALSE) message("Value on controls.y.sd will be ignored")
  if (length(controls.x) == 1 & is.null(controls.x.sd) == TRUE) stop("Please give sd on task x")
  if (length(controls.y) == 1 & is.null(controls.y.sd) == TRUE) stop("Please give sd on task y")

  con_m.x <- mean(controls.x) # Mean of the control sample on task x
  con_m.y <- mean(controls.y) # Mean of the control sample on task y

  con_sd.x <- stats::sd(controls.x) # Standard deviation of the control sample on task x
  if (length(controls.x) == 1 & is.null(controls.x.sd) == FALSE) con_sd.x <- controls.x.sd

  con_sd.y <- stats::sd(controls.y) # Standard deviation of the control sample on task y
  if (length(controls.y) == 1 & is.null(controls.y.sd) == FALSE) con_sd.y <- controls.y.sd


  # Since controls x and y need to be of equal length n is the length of any of them
  n <- length(controls.x)
  if (length(controls.x) == 1 | length(controls.y) == 1) {
    if (is.null(controls.n) == TRUE) stop("Please set sample size")
    n <- controls.n
    if (length(controls.x) > 1 & n != length(controls.x)) stop("Sample sizes must be equal")
    if (length(controls.y) > 1 & n != length(controls.y)) stop("Sample sizes must be equal")
  }

  if (is.null(cor.x.y) == FALSE){
    if (cor.x.y < -1 | cor.x.y > 1) stop("Correlation must be between -1 and 1")
  }

  r <- cor.x.y

  if (length(controls.x) > 1 & length(controls.y) > 1) r <- stats::cor(controls.x, controls.y)


  std.x <- (case.x - con_m.x)/con_sd.x
  std.y <- (case.y - con_m.y)/con_sd.y


  df <- n - 1

  t.crit <- abs(stats::qt(alpha/2, df= df))

  zdcc <- (std.x - std.y) / sqrt(2 - 2*r)

  if (exact.method == T) {

    # Exact probability - point estimate - alternative method for RSDT

    a <- (1 + r)*(1 - r^2)

    b <- (1 - r)*(
      4*(n - 1)^2 + 4*(1 + r)*(n - 1) + (1 + r)*(5 + r)
    )

    c <- -2*(
      ((
        (case.x - con_m.x)/con_sd.x -
          (case.y - con_m.y)/con_sd.y
      )^2) * (
        (n*(n - 1)^2)/(n + 1)
      )
    )

    tstat <- sqrt(
      (-b + sqrt(b^2 - (4*a*c))) / (2*a)
    )
    names(tstat) <- "t"

    pval <- 2 * pt(abs(tstat), df = df, lower.tail = F)

  } else {

    # Method from which the above is derived

    denom <- sqrt(
      ((n + 1)/n) *
        (
          (2 - 2*r) + (
            2*(1 - r^2)/(n - 1)
          ) + (
            ((5 + t.crit^2) * (1 - r^2))/(2*(n - 1)^2)
          ) + (
            (r*(1 + t.crit^2)*(1 - r^2))/((2*(n - 1)^2))
          )
        )
    )

    psi <- (std.x - std.y)/denom

    tstat <- psi

    names(tstat) <- "approx. t"

    pval <- 2 * pt(abs(psi), df = df, lower.tail = F)

  }


  estimate <- c(std.x, std.y, zdcc, (pval/2*100))

  # Set names for objects in output
  names(estimate) <- c("Case score on task X as standard (z) score",
                       "Case score on task Y as standard (z) score",
                       "Std. effect size (Z-DCC) for task diff. between case and controls",
                       "Percentage of control population with more extreme task difference")
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0(deparse(substitute(case.x)), " - ",
                  deparse(substitute(case.y)), " and ",
                  deparse(substitute(controls.x)), " - ",
                  deparse(substitute(controls.y)))

  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = tstat,
                 parameter = df,
                 p.value = pval,
                 estimate = estimate,
                 sample.size = n,
                 null.value = null.value,
                 alternative = "two.sided",
                 method = paste("Revised Standardised Difference Test"),
                 data.name = dname)

  class(output) <- "htest"
  output


}
