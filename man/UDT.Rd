% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UDT.R
\name{UDT}
\alias{UDT}
\title{Unstandardised Difference Test}
\usage{
UDT(
  case_a,
  case_b,
  controls_a,
  controls_b,
  sd_a = NULL,
  sd_b = NULL,
  sample_size = NULL,
  r_ab = NULL,
  alternative = c("two.sided", "greater", "less"),
  conf_int = TRUE,
  conf_level = 0.95,
  conf_int_spec = 0.01,
  na.rm = FALSE
)
}
\arguments{
\item{case_a}{Case's score on task A.}

\item{case_b}{Case's score on task B.}

\item{controls_a}{Controls' scores on task A. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task A while mean and SD for task B.}

\item{controls_b}{Controls' scores on task B. Takes either a vector of
observations or a single value interpreted as mean. \emph{Note}: you can
supply a vector as input for task B while mean and SD for task A.}

\item{sd_a}{If single value for task A is given as input you must
supply the standard deviation of the sample.}

\item{sd_b}{If single value for task B is given as input you must
supply the standard deviation of the sample.}

\item{sample_size}{If A or B is given as mean and SD you must supply the
sample size. If controls_a is given as vector and controls_b as mean and
SD, sample_size must equal the number of observations in controls_a.}

\item{r_ab}{If A and/or B is given as mean and SD you must supply the
correlation between the tasks.}

\item{alternative}{A character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter. Since the direction
of the expected effect depends on which task is set as A and which is set
as B, be very careful if changing this parameter.}

\item{conf_int}{Initiates a search algorithm for finding confidence
intervals. Defaults to \code{TRUE}, set to \code{FALSE} for faster
calculation (e.g. for simulations).}

\item{conf_level}{Level of confidence for intervals, defaults to 95\%.}

\item{conf_int_spec}{The size of iterative steps for calculating confidence
intervals. Smaller values gives more precise intervals but takes longer to
calculate. Defaults to a specificity of 0.01.}

\item{na.rm}{Remove \code{NA}s from controls.}
}
\value{
A list with class \code{"htest"} containing the following components:
  \tabular{llll}{ \code{statistic}   \tab the t-statistic. \cr\cr
  \code{parameter} \tab the degrees of freedom for the t-statistic.\cr\cr
  \code{p.value}    \tab the p-value of the test.\cr\cr \code{estimate} \tab
  unstandardised case scores, task difference and pont estimate of proportion
  control population expected to above or below the observed task difference.
  \cr\cr \code{control.desc}   \tab named numerical with descriptive
  statistics of the control samples. \cr\cr \code{null.value}   \tab the
  value of the difference under the null hypothesis.\cr\cr
  \code{alternative}     \tab a character string describing the alternative
  hypothesis.\cr\cr \code{method} \tab a character string indicating what
  type of test was performed.\cr\cr \code{data.name} \tab a character string
  giving the name(s) of the data}
}
\description{
A test on the discrepancy between two tasks in a single case, by comparison
to the mean of discrepancies of the same two tasks in a control sample. Use
\emph{only} when the two tasks are measured on the same scale with the same
underlying distribution because no standardisation is performed on task
scores. As a rule-of-thumb, the UDT may be applicable to pairs of tasks for
which it would be sensible to perform a paired t-test within the control
group. Calculates however a standardised effect size in the same manner as
\code{\link{RSDT}()}. This is original behaviour from Crawford and Garthwaite
(2005) but might not be appropriate. So use this standardised effect size
with caution. Calculates a standardised effect size of task discrepancy as
well as a point estimate of the proportion of the control population that
would be expected to show a more extreme discrepancy and respective
confidence intervals.
}
\details{
Running  \code{UDT} is equivalent to running \code{TD} on discrepancy scores
making it possible to run unstandardised tests with covariates by applying
\code{BTD_cov} to discrepancy scores.
}
\examples{
UDT(-3.857, -1.875, controls_a = 0, controls_b = 0, sd_a = 1,
sd_b = 1, sample_size = 20, r_ab = 0.68)

UDT(case_a = size_weight_illusion[1, "V_SWI"], case_b = size_weight_illusion[1, "K_SWI"],
 controls_a = size_weight_illusion[-1, "V_SWI"], controls_b = size_weight_illusion[-1, "K_SWI"])

}
\references{
Crawford, J. R., & Garthwaite, P. H. (2005). Testing for
Suspected Impairments and Dissociations in Single-Case Studies in
Neuropsychology: Evaluation of Alternatives Using Monte Carlo Simulations and
Revised Tests for Dissociations. \emph{Neuropsychology, 19}(3), 318 - 331.
\doi{10.1037/0894-4105.19.3.318}
}
