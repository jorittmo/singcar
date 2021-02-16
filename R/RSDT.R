#' Revised Standardised Difference Test
#'
#' A test on the discrepancy between two tasks in a single case, by comparison
#' to the discrepancy of means in the same two tasks in a control sample.
#' Standardises task scores as well as task discrepancy, so the tasks do not
#' need to be measured on the same scale. Calculates a standardised effect size
#' (Z-DCC) of task discrepancy as well as a point estimate of the proportion of
#' the control population that would be expected to show a more extreme
#' discrepancy. Developed by Crawford and Garthwaite (2005).
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
#' @param r_ab If A or B is given as mean and SD you must supply the
#'   correlation between the tasks.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or
#'   \code{"less"}. You can specify just the initial letter. Since the direction
#'   of the expected effect depends on which task is set as A and which is set
#'   as B, be very careful if changing this parameter.
#' @param na.rm Remove \code{NA}s from controls.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab Returns the value of a approximate
#'   t-statistic, however, because of the underlying equation, it cannot be
#'   negative. See effect direction from Z-DCC. \cr\cr \code{parameter} \tab the
#'   degrees of freedom for the t-statistic.\cr\cr \code{p.value}    \tab the
#'   p-value for the test.\cr\cr \code{estimate} \tab case scores expressed as
#'   z-scores on task A and Y. Standardised effect size (Z-DCC) of task
#'   difference between case and controls and point estimate of the proportion
#'   of the control population estimated to show a more extreme task
#'   discrepancy. \cr\cr \code{sample.size}   \tab the size of the control
#'   sample\cr\cr \code{null.value}   \tab the value of the discrepancy under
#'   the null hypothesis.\cr\cr  \code{alternative}     \tab a character string
#'   describing the alternative hypothesis.\cr\cr \code{method} \tab a character
#'   string indicating what type of test was performed.\cr\cr \code{data.name}
#'   \tab a character string giving the name(s) of the data}
#' @export
#'
#' @examples
#' RSDT(-3.857, -1.875, controls_a = 0, controls_b = 0, sd_a = 1,
#' sd_b = 1, sample_size = 20, r_ab = 0.68)
#'
#' RSDT(case_a = size_weight_illusion[1, "V_SWI"], case_b = size_weight_illusion[1, "K_SWI"],
#'  controls_a = size_weight_illusion[-1, "V_SWI"], controls_b = size_weight_illusion[-1, "K_SWI"])
#'
#' @references
#'
#' Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
#' Suspected Impairments and Dissociations in Single-Case Studies in
#' Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
#' Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
#' \doi{10.1037/0894-4105.19.3.318}



RSDT <- function (case_a, case_b, controls_a, controls_b,
                  sd_a = NULL, sd_b = NULL,
                  sample_size = NULL, r_ab = NULL,
                  alternative = c("two.sided", "greater", "less"),
                  na.rm = FALSE) {

  ###
  # Set up of error and warning messages
  ###

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

  ###
  # Extract relevant statistics and set up further errors
  ###

  alternative <- match.arg(alternative)

  con_m_a <- mean(controls_a) # Mean of the control sample on task A
  con_m_b <- mean(controls_b) # Mean of the control sample on task B

  con_sd_a <- stats::sd(controls_a) # Standard deviation of the control sample on task A
  if (length(controls_a) == 1 & is.null(sd_a) == FALSE) con_sd_a <- sd_a

  con_sd_b <- stats::sd(controls_b) # Standard deviation of the control sample on task B
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

  df <- n - 1 # Degrees of freedom

  ###
  # Calculate standardised effect sizes
  ###

  std_a <- (case_a - con_m_a)/con_sd_a
  std_b <- (case_b - con_m_b)/con_sd_b

  zdcc <- (std_a - std_b) / sqrt(2 - 2*r)

  ###
  # Get approximate t-value, see Garthwaite & Crawford (2004)
  ###

    a <- (1 + r)*(1 - r^2)

    b <- (1 - r)*(
      4*(n - 1)^2 + 4*(1 + r)*(n - 1) + (1 + r)*(5 + r)
    )

    c <- -2*(
      ((
        (case_a - con_m_a)/con_sd_a -
          (case_b - con_m_b)/con_sd_b
      )^2) * (
        (n*(n - 1)^2)/(n + 1)
      )
    )

    tstat <- sqrt(
      (-b + sqrt(b^2 - (4*a*c))) / (2*a)
    )
    names(tstat) <- "approx. abs. t"

    ###
    # Get p-value
    ###

    if (alternative == "two.sided") {
      pval <- 2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE)

      if (zdcc < 0) {
        p.name <- "Proportion below case (%)"
      } else {
        p.name <- "Proportion above case (%)"
      }



    } else if (alternative == "greater") {
      # Since equation (7) from Crawford and Garthwaite (the exact method)
      # cannot return a negative t-value we have to use zdcc to see
      # in which direction the effect it pointing and the impose the correct sign.
      if (zdcc > 0) {
        pval <- stats::pt(tstat, df = df, lower.tail = FALSE)
      } else {
        pval <- stats::pt(-tstat, df = df, lower.tail = FALSE)
      }
      p.name <- "Proportion above case (%)"
    }

     else { # I.e. if alternative == "less"

      if (zdcc < 0) {
        pval <- stats::pt(-tstat, df = df, lower.tail = TRUE)
      } else {
        pval <- stats::pt(tstat, df = df, lower.tail = TRUE)
      }
       p.name <- "Proportion below case (%)"
    }


  #   # Method from which the above is derived
  #   # Not included as functionality for now.
  #
  #   if (alternative == "two.sided") {
  #     t.crit <- abs(stats::qt(alpha/2, df= df))
  #   } else if (alternative == "greater") {
  #     t.crit <- stats::qt(alpha, df= df, lower.tail = FALSE)
  #   } else {
  #     t.crit <- stats::qt(alpha, df= df)
  #   }
  #
  #   denom <- sqrt(
  #     ((n + 1)/n) *
  #       (
  #         (2 - 2*r) + (
  #           2*(1 - r^2)/(n - 1)
  #         ) + (
  #           ((5 + t.crit^2) * (1 - r^2))/(2*(n - 1)^2)
  #         ) + (
  #           (r*(1 + t.crit^2)*(1 - r^2))/((2*(n - 1)^2))
  #         )
  #       )
  #   )
  #
  #   psi <- (std_a - std_b)/denom
  #
  #   tstat <- psi
  #
  #   names(tstat) <- "approx. t"
  #
  #   if (alternative == "two.sided") {
  #     pval <- 2 * stats::pt(abs(psi), df = df, lower.tail = FALSE)
  #   } else if (alternative == "greater") {
  #     pval <- stats::pt(psi, df = df, lower.tail = FALSE)
  #   } else { # I.e. if alternative == "less"
  #     pval <- stats::pt(psi, df = df, lower.tail = TRUE)
  #   }


  estimate <- c(std_a, std_b, zdcc, ifelse(alternative == "two.sided", (pval/2*100), pval*100))

  ###
  # Set names for objects in output
  ###

  names(estimate) <- c("Std. case score, task A (Z-CC)",
                       "Std. case score, task B (Z-CC)",
                       "Std. task discrepancy (Z-DCC)",
                       p.name)
  names(df) <- "df"
  null.value <- 0 # Null hypothesis: difference = 0
  names(null.value) <- "difference between tasks"
  dname <- paste0("Case A: ", format(round(case_a, 2), nsmall = 2), ", ",
                  "B: ", format(round(case_b, 2), nsmall = 2), ", ",
                  "Ctrl. A (m, sd): (", format(round(con_m_a, 2), nsmall = 2), ", ",format(round(con_sd_a, 2), nsmall = 2), "), ",
                  "B: (", format(round(con_m_b, 2), nsmall = 2), ", ",format(round(con_sd_b, 2), nsmall = 2), ")")
  names(pval) <- NULL

  # Build output to be able to set class as "htest" object for S3 methods.
  # See documentation for "htest" class for more info

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



