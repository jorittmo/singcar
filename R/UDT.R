#' Unstandardised Difference Test
#'
#' A test on the difference between two tasks in a single case, by comparison to
#' a control sample. Use only when the two tasks are measured on the \emph{same}
#' scale.
#'
#' @param case_a Case's score on task A.
#' @param case_b Case's score on task B.
#' @param controls_a Controls' scores on task A. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can
#'   supply a vector as input for task A while mean and SD for task B.
#' @param controls_b Controls' scores on task B. Takes either a vector of
#'   observations or a single value interpreted as mean. \emph{Note}: you can
#'   supply a vector as input for task B while mean and SD for task A.
#' @param sd_a If single value for task A is given as input you must
#'   supply the standard deviation of the sample.
#' @param sd_b If single value for task B is given as input you must
#'   supply the standard deviation of the sample.
#' @param sample_size If A or B is given as mean and SD you must supply the
#'   sample size. If controls_a is given as vector and controls_b as mean and
#'   SD, sample_size must equal the number of observations in controls_a.
#' @param r_ab If A and/or B is given as mean and SD you must supply the
#'   correlation between the tasks.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter. Since the direction
#'   of the expected effect depends on which task is set as A and which is set
#'   as B, be very careful if changing this parameter.
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
#' UDT(-3.857, -1.875, controls_a = 0, controls_b = 0, sd_a = 1,
#' sd_b = 1, sample_size = 20, r_ab = 0.68)
#'
#' UDT(case_a = size_weight_illusion[1, "V_SWI"], case_b = size_weight_illusion[1, "K_SWI"],
#'  controls_a = size_weight_illusion[-1, "V_SWI"], controls_b = size_weight_illusion[-1, "K_SWI"])
#'
#' @references {Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
#' Suspected Impairments and Dissociations in Single-Case Studies in
#' Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
#' Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
#' \url{https://doi.org/10.1037/0894-4105.19.3.318}}



UDT <- function (case_a, case_b, controls_a, controls_b,
                  sd_a = NULL, sd_b = NULL,
                  sample_size = NULL, r_ab = NULL,
                  alternative = c("two.sided", "greater", "less"),
                  na.rm = FALSE) {

  alternative <- match.arg(alternative)

  if (length(case_a) > 1 | length(case_b) > 1) stop("Case scores should be single value")
  if (length(controls_a) > 1 & length(controls_b) > 1) {
    if (length(controls_a) != length(controls_b)) stop("Sample sizes must be equal")
  }

  if (length(controls_a) > 1 & length(controls_b) > 1 & is.null(sample_size) == FALSE) message("Value on sample_size will be ignored")

  if (length(controls_a) > 1 & is.null(sd_a) == FALSE) message("Value on sd_a will be ignored")
  if (length(controls_b) > 1 & is.null(sd_b) == FALSE) message("Value on sd_b will be ignored")
  if (length(controls_a) == 1 & is.null(sd_a) == TRUE) stop("Please give sd and n on task A if controls_a is to be treated as mean")
  if (length(controls_b) == 1 & is.null(sd_b) == TRUE) stop("Please give sd and n on task B if controls_b is to be treated as mean")


  # Handling of NA use cases below
  if(is.na(case_a) == TRUE | is.na(case_b) == TRUE) stop("One or both case scores is NA")

  if (na.rm == TRUE) {
    if (sum(is.na(controls_a))  > 0 & sum(is.na(controls_b)) == 0 ) {
      controls_b <- controls_b[!is.na(controls_a)]
      controls_a <- controls_a[!is.na(controls_a)]
      warning("Removal of NAs on controls_a resulted in removal of non-NAs on controls_b")
    }

    if (sum(is.na(controls_b))  > 0 & sum(is.na(controls_a)) == 0 ) {
      controls_a <- controls_a[!is.na(controls_b)]
      controls_b <- controls_b[!is.na(controls_b)]
      warning("Removal of NAs on controls_b resulted in removal of non-NAs on controls_a")
    }

    if (sum(is.na(controls_b))  > 0 & sum(is.na(controls_a)) > 0 ) {

      if (identical(!is.na(controls_a), !is.na(controls_b)) == TRUE) {
        controls_a <- controls_a[!is.na(controls_a)]
        controls_b <- controls_b[!is.na(controls_b)]
      } else {
        con_a <- controls_a[!is.na(controls_a) & !is.na(controls_b)]
        con_b <- controls_b[!is.na(controls_a) & !is.na(controls_b)]

        controls_a <- con_a
        controls_b <- con_b

        warning("Removal of NAs on one control sample resulted in removal of non-NAs on the other")
      }

    }

  }
  if (sum(is.na(controls_a)) > 0 | sum(is.na(controls_b)) > 0) stop("Controls contains NA, set na.rm = TRUE to proceed")
  # End of NA use cases


  if (length(controls_a) > 1 & length(controls_b) > 1) {
    if (length(controls_a) != length(controls_b)) stop("Sample sizes must be equal")
  }

  con_m_a <- mean(controls_a) # Mean of the control sample on task x
  con_m_b <- mean(controls_b) # Mean of the control sample on task y

  con_sd_a <- stats::sd(controls_a) # Standard deviation of the control sample on task x
  if (length(controls_a) == 1 & is.null(sd_a) == FALSE) con_sd_a <- sd_a

  con_sd_b <- stats::sd(controls_b) # Standard deviation of the control sample on task y
  if (length(controls_b) == 1 & is.null(sd_b) == FALSE) con_sd_b <- sd_b


  # Since controls x and y need to be of equal length n is the length of any of them
  n <- length(controls_a)
  if (length(controls_a) == 1 | length(controls_b) == 1) {
    if (is.null(sample_size) == TRUE) stop("Please set sample size")
    n <- sample_size
    if (length(controls_a) > 1 & n != length(controls_a)) stop("Sample sizes must be equal")
    if (length(controls_b) > 1 & n != length(controls_b)) stop("Sample sizes must be equal")
  }

  if (is.null(r_ab) == TRUE & length(controls_a) == 1) stop("Please set correlation between tasks")
  if (is.null(r_ab) == TRUE & length(controls_b) == 1) stop("Please set correlation between tasks")

  if (is.null(r_ab) == FALSE){
    if (r_ab < -1 | r_ab > 1) stop("Correlation must be between -1 and 1")
  }

  r <- r_ab

  if (length(controls_a) > 1 & length(controls_b) > 1) r <- stats::cor(controls_a, controls_b)

  df <- n - 1

  def_a <- (case_a - con_m_a)
  def_b <- (case_b - con_m_b)

  dif <- (def_a - def_b)

  std.er <- sqrt(
    (con_sd_a^2 + con_sd_b^2 - 2*con_sd_a*con_sd_b*r) * ((n + 1) / n)
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

  estimate <- c(def_a, def_b, dif, ifelse(alternative == "two.sided", (pval/2*100), pval*100))

  if (alternative == "two.sided") {
    p.name <- "Proportion of control population with more extreme task difference"
  } else if (alternative == "greater") {
    p.name <- "Proportion of control population with more positive task difference"
  } else {
    p.name <- "Proportion of control population with more negative task difference"
  }


  # Set names for objects in output
  names(estimate) <- c("Case score on task A",
                       "Case score on task B",
                       "Task difference",
                       p.name)
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0("Case score A: ", deparse(substitute(case_a)), ", ",
                  "Case score B: ", deparse(substitute(case_b)), ", ",
                  "Controls score A: ", deparse(substitute(controls_a)), ", ",
                  "Controls score B: ", deparse(substitute(controls_b)))


  names(con_m_a) <- "Mean A"
  names(con_m_b) <- "Mean B"
  names(con_sd_a) <- "SD A"
  names(con_sd_b) <- "SD B"
  names(n) <- "Sample size"
  control.desc <- c(con_m_a, con_m_b, con_sd_a, con_sd_b, n)

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



