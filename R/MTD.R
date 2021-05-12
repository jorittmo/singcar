#' Multivariate Test of deficit
#'
#' @param case Vector of case scores
#' @param controls Matrix or data frame with scores from the control sample, each column representing a variable
#' @param conf_level Level of confidence for the confidence intervals
#' @param method One out of "pd", "pchi", "pf" and "pmd". Use "pmd" if the Mahalanobi's distance seems suspiciously small
#' @param mahalanobis_dist Mahalanobi's distance of the case if summary statistics are used
#' @param k The number of dimensions, if summary statistics are used
#' @param n The size of the control sample
#'
#' @return A list with class \code{"htest"} containing the following components:
#'   \tabular{llll}{ \code{statistic}   \tab Hotelling's T^2 statistic for the
#'   case's Mahalanobi's distance \cr\cr
#'   \code{p.value}    \tab The p value associated with the Hotelling statistic
#'    \cr\cr \code{estimate} \tab Estimates of the case Mahalanobis distance
#'    and index as well as abnormality \cr\cr \code{interval} \tab  List
#'    of interval measure for the estimates \cr\cr \code{sample.size} \tab number
#'   of controls.\cr\cr \code{method} \tab a character
#'   string indicating what type of test was performed and which abnormality measure
#'   used}
#' @export
#'
#' @examples
#'
MTD <- function(case, controls, conf_level = 0.95, method = c("pd", "pchi", "pf", "pmd"),
                mahalanobis_dist = NULL, k = NULL, n = NULL) {

  if (!is.null(mahalanobis_dist)){
    case <- NULL
    controls <- NULL
    delta_hat <- mahalanobis_dist
  }

  if (!is.null(case)){
    xbar <- colMeans(controls)
    S <- solve(cov(controls))
    k <- length(case)
    n <- nrow(controls)

    delta_hat = sqrt(
      t(case-xbar) %*% S %*% (case-xbar)
    )
  }

  if (delta_hat == 0) delta_hat <- 1e-10

  nu1 <- k
  nu2 <- n - k
  lamhat <- delta_hat^2
  method <- match.arg(method)



  if (method == "pmd") {

    likeli <- function(r, lam){
      m <- n/(n-1)
      ((((n^r)/2)*((lam/2)^((k/2)+r-1))*exp(-(n+1)*lam/2)) /
          (beta((n-k)/2, k/2+r)*gamma(k/2)*factorial(r)))*(m^((k/2)+r))*((1/(1+(m*lamhat)))^((n/2)+r))*(lamhat^((k/2)+r-1))
    }
    likeli <- Vectorize(likeli)

    post <- function(lam){
      x <- likeli(0:170, lam)
      x <- x[!is.nan(x) & !is.infinite(x)]
      sum(x)
    }
    post <- Vectorize(post)



    c <- integrate(post, 0, sqrt(lamhat)/2)$value +
      integrate(post, sqrt(lamhat)*0.5, sqrt(lamhat))$value +
      integrate(post, sqrt(lamhat), sqrt(lamhat)*1.5)$value +
      integrate(post, sqrt(lamhat)*1.5, sqrt(lamhat)*2)$value +
      integrate(post, sqrt(lamhat)*2, Inf)$value

    post_norm <- function(lam) {
      post(lam)/c
    }
    post_norm <- Vectorize(post_norm)


    postquant <- function(lim, quantile) {
      abs(integrate(post_norm, 0, lim)$value - quantile)
    }
    postquant <- Vectorize(postquant)

    lammed <- nlminb(delta_hat, postquant, quantile = 0.5, lower = 0)$par
    lamu <- nlminb(delta_hat, postquant, quantile = (1+conf_level)/2, lower = 0)$par
    laml <- nlminb(delta_hat, postquant, quantile = (1-conf_level)/2, lower = 0)$par


    pmd <- (1 - pchisq(lammed, df = nu1))*100
    pmdl <- (1 - pchisq(lamu, df = nu1))*100
    pmdu <- (1 - pchisq(laml, df = nu1))*100

    p_int <- c(pmdl, pmdu)
    dist_int <- c(sqrt(laml), sqrt(lamu))
    index_int <- c(laml, lamu)


    estimate <- c(sqrt(lammed), lammed, pmd)

    dist.name <- paste0("Median Mahalanobi's distance ",
                        100*conf_level, "% CI [",
                        format(round(sqrt(laml), 2), nsmall = 2),", ",
                        format(round(sqrt(lamu), 2), nsmall = 2),"]")
    ind.name <- paste0("Median Mahalanobi's index ",
                       100*conf_level, "% CI [",
                       format(round(laml, 2), nsmall = 2),", ",
                       format(round(lamu, 2), nsmall = 2),"]")
    p.name <- paste0("Proportion above case ",
                     100*conf_level, "% CI [",
                     format(round(pmdl, 2), nsmall = 2),", ",
                     format(round(pmdu, 2), nsmall = 2),"]")

    names(estimate) <- c(dist.name, ind.name, p.name)


  }

  if (method == "pd") {
    F0 <- (n*nu2)/((n-1)*nu1)*lamhat


    fquant <- function(nlam, qu) {
      abs(pf(F0, nu1, nu2, ncp = nlam) - qu)
    }
    fquant <- Vectorize(fquant)

    lammed <- nlminb(F0, fquant, qu = 0.5, lower = 0)$par/n
    laml <- nlminb(F0, fquant, qu = (1+conf_level)/2, lower = 0)$par/n
    lamu <- nlminb(F0, fquant, qu = (1-conf_level)/2, lower = 0)$par/n


    pd <- (1 - pchisq(lammed, nu1))*100
    pl <- (1 - pchisq(lamu, nu1))*100
    pu <- (1 - pchisq(laml, nu1))*100

    p_int <- c(pl, pu)
    dist_int <- c(sqrt(laml), sqrt(lamu))
    index_int <- c(laml, lamu)


    estimate <- c(sqrt(lammed), lammed, pd)

    dist.name <- paste0("Median Mahalanobi's distance ",
                        100*conf_level, "% CI [",
                        format(round(sqrt(laml), 2), nsmall = 2),", ",
                        format(round(sqrt(lamu), 2), nsmall = 2),"]")
    ind.name <- paste0("Median Mahalanobi's index ",
                       100*conf_level, "% CI [",
                       format(round(laml, 2), nsmall = 2),", ",
                       format(round(lamu, 2), nsmall = 2),"]")
    p.name <- paste0("Proportion above case ",
                     100*conf_level, "% CI [",
                     format(round(pl, 2), nsmall = 2),", ",
                     format(round(pu, 2), nsmall = 2),"]")

    names(estimate) <- c(dist.name, ind.name, p.name)
  }


  if (method == "pchi") {

    q <- (n/(n-1))*lamhat
    pchi <- pchisq(q, df = nu1, lower.tail = FALSE)*100

    F0 <- (n*nu2)/((n-1)*nu1)*lamhat
    fquant <- function(nlam, qu) {
      abs(pf(F0, nu1, nu2, ncp = nlam) - qu)
    }
    fquant <- Vectorize(fquant)

    laml <- nlminb(F0, fquant, qu = (1+conf_level)/2, lower = 0)$par/n
    lamu <- nlminb(F0, fquant, qu = (1-conf_level)/2, lower = 0)$par/n

    pl <- (1 - pchisq(lamu, nu1))*100
    pu <- (1 - pchisq(laml, nu1))*100

    p_int <- c(pl, pu)
    dist_int <- c(sqrt(laml), sqrt(lamu))
    index_int <- c(laml, lamu)


    estimate <- c(sqrt(lamhat), lamhat, pchi)

    dist.name <- paste0("Median Mahalanobi's distance ",
                        100*conf_level, "% CI [",
                        format(round(sqrt(lamu), 2), nsmall = 2),", ",
                        format(round(sqrt(laml), 2), nsmall = 2),"]")
    ind.name <- paste0("Median Mahalanobi's index ",
                       100*conf_level, "% CI [",
                       format(round(lamu, 2), nsmall = 2),", ",
                       format(round(laml, 2), nsmall = 2),"]")
    p.name <- paste0("Proportion above case ",
                     100*conf_level, "% CI [",
                     format(round(pl, 2), nsmall = 2),", ",
                     format(round(pu, 2), nsmall = 2),"]")

    names(estimate) <- c(dist.name, ind.name, p.name)
  }

  if (method == "pf") {

    T2 <- (n/(n+1))*lamhat
    q <- (nu2/((n-1)*nu1))*T2
    pF <- pf(q, nu1, nu2, lower.tail = FALSE)*100

    F0 <- (n*nu2)/((n-1)*nu1)*lamhat
    fquant <- function(nlam, qu) {
      abs(pf(F0, nu1, nu2, ncp = nlam) - qu)
    }
    fquant <- Vectorize(fquant)

    laml <- nlminb(F0, fquant, qu = (1+conf_level)/2, lower = 0)$par/n
    lamu <- nlminb(F0, fquant, qu = (1-conf_level)/2, lower = 0)$par/n

    pl <- (1 - pchisq(lamu, nu1))*100
    pu <- (1 - pchisq(laml, nu1))*100

    p_int <- c(pl, pu)
    dist_int <- c(sqrt(laml), sqrt(lamu))
    index_int <- c(laml, lamu)

    estimate <- c(sqrt(lamhat), lamhat, pF)

    dist.name <- paste0("Median Mahalanobi's distance ",
                        100*conf_level, "% CI [",
                        format(round(sqrt(laml), 2), nsmall = 2),", ",
                        format(round(sqrt(lamu), 2), nsmall = 2),"]")
    ind.name <- paste0("Median Mahalanobi's index ",
                       100*conf_level, "% CI [",
                       format(round(laml, 2), nsmall = 2),", ",
                       format(round(lamu, 2), nsmall = 2),"]")
    p.name <- paste0("Proportion above case ",
                     100*conf_level, "% CI [",
                     format(round(pl, 2), nsmall = 2),", ",
                     format(round(pu, 2), nsmall = 2),"]")

    names(estimate) <- c(dist.name, ind.name, p.name)
  }

  T2 <- (n/(n+1))*lamhat
  q <- (nu2/((n-1)*nu1))*T2
  p.val <- pf(q, nu1, nu2, lower.tail = FALSE)
  names(T2) <- "Hotelling's T^2"

  typ.int <- 100*conf_level
  names(typ.int) <- "Interval level (%)"
  interval <- list(Interval.level = typ.int, Distance = dist_int, Index = index_int, Abnormality = p_int)


  output <- list(statistic = T2, p.value = p.val,
                 estimate = estimate,
                 interval = interval,
                 sample.size = n,
                 method = paste0("Multivariate Test of deficit, abnormality estimate ", method))

  class(output) <- "htest"
  output

}
