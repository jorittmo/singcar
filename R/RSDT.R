#' Revised Standardised Difference Test
#'
#' A test on the difference between two tasks in a single case, by comparison to
#' a control sample. Calculates a standardised effects size of task difference
#' as well as a point estimate of the proportion of the control population that
#' would be expected to show a more extreme task difference.
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
#' @param exact.method If set to \code{FALSE} generates an approximate
#'   t-statistic used to derive the exact. The exact method can only generate an
#'   absolute t-statistic. The approximate method also requires a pre-set
#'   alpha-value, see Crawford and Garthwate (2005) for more information.
#' @param alpha Chosen risk of Type I errors. This is only relevant if setting
#'   \code{exact.method = FALSE} due to the test statistic depending on the
#'   critical value chosen.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab if exact.method set to \code{TRUE},
#'   returns the value of an exact t-statistic, however, because of the
#'   underlying equation, it cannot be negative. Set exact.method to
#'   \code{FALSE} for an approximate t-value with correct sign. \cr\cr
#'   \code{parameter} \tab the degrees of freedom for the t-statistic.\cr\cr
#'   \code{p.value}    \tab the p-value for the test.\cr\cr \code{estimate} \tab
#'   case scores expressed as z-scores on task X and Y. Standardised effect size
#'   (Z-DCC) of task difference between case and controls and point estimate of
#'   the proportion of the control population estimated to show a more extreme
#'   task difference. \cr\cr \code{sample.size}   \tab the size of the control
#'   sample\cr\cr \code{null.value}   \tab the value of the difference under the
#'   null hypothesis.\cr\cr  \code{alternative}     \tab a character string
#'   describing the alternative hypothesis.\cr\cr \code{method} \tab a character
#'   string indicating what type of test was performed.\cr\cr \code{data.name}
#'   \tab a character string giving the name(s) of the data}
#' @export
#'
#' @examples
#' RSDT(-3.857, -1.875, controls.x = 0, controls.y = 0, controls.x.sd = 1,
#' controls.y.sd = 1, controls.n = 20, cor.x.y = 0.68)
#' RSDT(-3.857, -1.875, controls.x = rnorm(20), controls.y = rnorm(20))
#'
#' @references {Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
#' Suspected Impairments and Dissociations in Single-Case Studies in
#' Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
#' Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
#' \url{https://doi.org/10.1037/0894-4105.19.3.318}}



RSDT <- function (case.x, case.y, controls.x, controls.y,
                  controls.x.sd = NULL, controls.y.sd = NULL,
                  controls.n = NULL, cor.x.y = NULL,
                  alternative = c("two.sided", "greater", "less"),
                  exact.method = T, alpha = 0.05, na.rm = FALSE) {

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

  std.x <- (case.x - con_m.x)/con_sd.x
  std.y <- (case.y - con_m.y)/con_sd.y

  zdcc <- (std.x - std.y) / sqrt(2 - 2*r) # Estimated effect size

  if (exact.method == T) {

    # Exact probability - point estimate

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
    names(tstat) <- "absolute t"

    if (alternative == "two.sided") {
      pval <- 2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE)
    } else if (alternative == "greater") {
      # Since equation (7) from Crawford and Garthwaite (the exact method)
      # cannot return a negative t-value we have to use zdcc to see
      # in which direction the effect it pointing and the impose the correct sign.
      if (zdcc > 0) {
        pval <- stats::pt(tstat, df = df, lower.tail = FALSE)
      } else {
        pval <- stats::pt(-tstat, df = df, lower.tail = FALSE)
      }

    }

     else { # I.e. if alternative == "less"

      if (zdcc < 0) {
        pval <- stats::pt(-tstat, df = df, lower.tail = TRUE)
      } else {
        pval <- stats::pt(tstat, df = df, lower.tail = TRUE)
      }

    }

  } else {

    # Method from which the above is derived

    if (alternative == "two.sided") {
      t.crit <- abs(stats::qt(alpha/2, df= df))
    } else if (alternative == "greater") {
      t.crit <- stats::qt(alpha, df= df, lower.tail = FALSE)
    } else {
      t.crit <- stats::qt(alpha, df= df)
    }

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

    if (alternative == "two.sided") {
      pval <- 2 * stats::pt(abs(psi), df = df, lower.tail = FALSE)
    } else if (alternative == "greater") {
      pval <- stats::pt(psi, df = df, lower.tail = FALSE)
    } else { # I.e. if alternative == "less"
      pval <- stats::pt(psi, df = df, lower.tail = TRUE)
    }

  }



  estimate <- c(std.x, std.y, zdcc, ifelse(alternative == "two.sided", (pval/2*100), pval*100))

  if (alternative == "two.sided") {
    p.name <- "Proportion of control population with more extreme task difference"
  } else if (alternative == "greater") {
    p.name <- "Proportion of control population with more positive task difference"
  } else {
    p.name <- "Proportion of control population with more negative task difference"
  }


  # Set names for objects in output
  names(estimate) <- c("Case score on task X as standard (z) score",
                       "Case score on task Y as standard (z) score",
                       "Std. effect size (Z-DCC) for task diff. between case and controls",
                       p.name)
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0("Case score X: ", deparse(substitute(case.x)), ", ",
                  "Case score Y: ", deparse(substitute(case.y)), ", ",
                  "Controls score X: ", deparse(substitute(controls.x)), ", ",
                  "Controls score Y: ",deparse(substitute(controls.y)))

  # Build output to be able to set class as "htest" object. See documentation for "htest" class for more info
  output <- list(statistic = tstat,
                 parameter = df,
                 p.value = pval,
                 estimate = estimate,
                 sample.size = n,
                 null.value = null.value,
                 alternative = alternative,
                 method = paste("Revised Standardised Difference Test"),
                 data.name = dname)

  class(output) <- "htest"
  output


}



