#' Unstandardised Difference Test
#'
#' A test on the difference between two tasks in a single case, by comparison to
#' a control sample. Use only when the two tasks are measured on the \emph{same}
#' scale.
#'
#' @param case.x Case's score on task X.
#' @param case.y Case's score on task Y.
#' @param controls.x Controls' scores on task X. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can
#'   supply a vector as input for task X while mean and SD for task Y.
#' @param controls.y Controls' scores on task Y. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can
#'   supply a vector as input for task Y while mean and SD for task X.
#' @param controls.x.sd If single value for task X is given as input you must
#'   supply the standard deviation of the sample.
#' @param controls.y.sd If single value for task Y is given as input you must
#'   supply the standard deviation of the sample.
#' @param controls.n If X or Y is given as mean and SD you must supply the
#'   sample size. If controls.x is given as vector and controls.y as mean and
#'   SD, controls.n must equal the number of observations in controls.x.
#' @param cor.x.y If X or Y is given as mean and SD you must supply the
#'   correlation between the tasks.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter. Since the direction
#'   of the expected effect depends on which task is set as X and which is set
#'   as Y, be very careful if changing this parameter.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab the t-statistic. \cr\cr
#'   \code{parameter} \tab the degrees of freedom for the t-statistic.\cr\cr
#'   \code{p.value}    \tab the p-value of the test.\cr\cr \code{estimate} \tab
#'   unstandardised case scores, task difference and pont estimate of proportion
#'   control population expected to above or below the observed task difference.
#'   \cr\cr \code{control.desc}   \tab named numerical with descriptive
#'   statistics of the control samples. \cr\cr \code{null.value}   \tab the
#'   value of the difference under the null hypothesis.\cr\cr
#'   \code{alternative}     \tab a character string describing the alternative
#'   hypothesis.\cr\cr \code{method} \tab a character string indicating what
#'   type of test was performed.\cr\cr \code{data.name} \tab a character string
#'   giving the name(s) of the data}
#' @export
#'
#' @examples
#' UDT(-3.857, -1.875, controls.x = 0, controls.y = 0, controls.x.sd = 1,
#' controls.y.sd = 1, controls.n = 20, cor.x.y = 0.68)
#' UDT(-3.857, -1.875, controls.x = rnorm(20), controls.y = rnorm(20))
#'
#' @references {Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
#' Suspected Impairments and Dissociations in Single-Case Studies in
#' Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
#' Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
#' \url{https://doi.org/10.1037/0894-4105.19.3.318}}



UDT <- function (case.x, case.y, controls.x, controls.y,
                  controls.x.sd = NULL, controls.y.sd = NULL,
                  controls.n = NULL, cor.x.y = NULL,
                  alternative = c("two.sided", "greater", "less"),
                  na.rm = FALSE) {

  alternative <- match.arg(alternative)

  if (length(case.x) > 1 | length(case.y) > 1) stop("Case scores should be single value")
  if (length(controls.x) > 1 & length(controls.y) > 1) {
    if (length(controls.x) != length(controls.y)) stop("Sample sizes must be equal")
  }

  if (length(controls.x) > 1 & length(controls.y) > 1 & is.null(controls.n) == FALSE) message("Value on controls.n will be ignored")

  if (length(controls.x) > 1 & is.null(controls.x.sd) == FALSE) message("Value on controls.x.sd will be ignored")
  if (length(controls.y) > 1 & is.null(controls.y.sd) == FALSE) message("Value on controls.y.sd will be ignored")
  if (length(controls.x) == 1 & is.null(controls.x.sd) == TRUE) stop("Please give sd and n on task x if controls.x is to be treated as mean")
  if (length(controls.y) == 1 & is.null(controls.y.sd) == TRUE) stop("Please give sd and n on task y if controls.y is to be treated as mean")


  # Handling of NA use cases below
  if(is.na(case.x) == TRUE | is.na(case.y) == TRUE) stop("One or both case scores is NA")

  if (na.rm == TRUE) {
    if (sum(is.na(controls.x))  > 0 & sum(is.na(controls.y)) == 0 ) {
      controls.y <- controls.y[!is.na(controls.x)]
      controls.x <- controls.x[!is.na(controls.x)]
      warning("Removal of NAs on controls.x resulted in removal of non-NAs on controls.y")
    }

    if (sum(is.na(controls.y))  > 0 & sum(is.na(controls.x)) == 0 ) {
      controls.x <- controls.x[!is.na(controls.y)]
      controls.y <- controls.y[!is.na(controls.y)]
      warning("Removal of NAs on controls.y resulted in removal of non-NAs on controls.x")
    }

    if (sum(is.na(controls.y))  > 0 & sum(is.na(controls.x)) > 0 ) {

      if (identical(!is.na(controls.x), !is.na(controls.y)) == TRUE) {
        controls.x <- controls.x[!is.na(controls.x)]
        controls.y <- controls.y[!is.na(controls.y)]
      } else {
        conx <- controls.x[!is.na(controls.x) & !is.na(controls.y)]
        cony <- controls.y[!is.na(controls.x) & !is.na(controls.y)]

        controls.x <- conx
        controls.y <- cony

        warning("Removal of NAs on one control sample resulted in removal of non-NAs on the other")
      }

    }

  }
  if (sum(is.na(controls.x)) > 0 | sum(is.na(controls.y)) > 0) stop("Controls contains NA, set na.rm = TRUE to proceed")
  # End of NA use cases


  if (length(controls.x) > 1 & length(controls.y) > 1) {
    if (length(controls.x) != length(controls.y)) stop("Sample sizes must be equal")
  }

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

  if (is.null(cor.x.y) == TRUE & length(controls.x) == 1) stop("Please set correlation between tasks")
  if (is.null(cor.x.y) == TRUE & length(controls.y) == 1) stop("Please set correlation between tasks")

  if (is.null(cor.x.y) == FALSE){
    if (cor.x.y < -1 | cor.x.y > 1) stop("Correlation must be between -1 and 1")
  }

  r <- cor.x.y

  if (length(controls.x) > 1 & length(controls.y) > 1) r <- stats::cor(controls.x, controls.y)

  df <- n - 1

  def.x <- (case.x - con_m.x)
  def.y <- (case.y - con_m.y)

  dif <- (def.x - def.y)

  std.er <- sqrt(
    (con_sd.x^2 + con_sd.y^2 - 2*con_sd.x*con_sd.y*r) * ((n + 1) / n)
  )


  tstat <- dif/std.er
  names(tstat) <- "t"

  if (alternative == "two.sided") {
    pval <- 2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE)
  } else if (alternative == "greater") {
    pval <- stats::pt(tstat, df = df, lower.tail = FALSE)
  } else { # I.e. if alternative == "less"
    pval <- stats::pt(tstat, df = df, lower.tail = TRUE)
  }

  estimate <- c(def.x, def.y, dif, ifelse(alternative == "two.sided", (pval/2*100), pval*100))

  if (alternative == "two.sided") {
    p.name <- "Proportion of control population with more extreme task difference"
  } else if (alternative == "greater") {
    p.name <- "Proportion of control population with more positive task difference"
  } else {
    p.name <- "Proportion of control population with more negative task difference"
  }


  # Set names for objects in output
  names(estimate) <- c("Case score on task X",
                       "Case score on task Y",
                       "Task difference",
                       p.name)
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0("Case score X: ", deparse(substitute(case.x)), ", ",
                  "Case score Y: ", deparse(substitute(case.y)), ", ",
                  "Controls score X: ", deparse(substitute(controls.x)), ", ",
                  "Controls score Y: ", deparse(substitute(controls.y)))


  names(con_m.x) <- "Mean X"
  names(con_m.y) <- "Mean Y"
  names(con_sd.x) <- "SD X"
  names(con_sd.y) <- "SD Y"
  names(n) <- "Sample size"
  control.desc <- c(con_m.x, con_m.y, con_sd.x, con_sd.y, n)

  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = tstat,
                 parameter = df,
                 p.value = pval,
                 estimate = estimate,
                 control.desc = control.desc,
                 null.value = null.value,
                 alternative = alternative,
                 method = paste("Unstandardised Difference Test"),
                 data.name = dname)

  class(output) <- "htest"
  output


}



