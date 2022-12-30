#' Random non normal
#'
#' Random number generation for a non normal standard distribution using contamination
#' degree based on Fleishman coefficients (altering skewness and kurtosis)
#'
#' @param n number of observations.
#' @param f.coef Fleishman coefficients indicating normal contamination degree
#'
#' @return
#' a vector of length \code{n} containing a non normal sample
#'
#' @details
#' Being \eqn{X} a normal standard,  linear combination
#' \eqn{Z = a + bX + cX^2 + dX^3} follows an unknown distribution with skewness
#' \eqn{\gamma_1} and \eqn{\gamma_2} defining contamination degree from normal.
#' Coefficients \eqn{a, b, c, d} are known as Fleishman coefficients.
#' For example for  \eqn{a = 0, b = 1, c = 0, d = 0} does not exist contamination
#'
#' @export
#'
#' @examples
#' # No contaminated sample:
#' ncs <- rnonorm(n = 100, f.coef = c(0, 1, 0, 0))
#' hist(ncs, probability = TRUE)
#' curve(dnorm(x), add = TRUE, col = "red")
#'
#' # Contaminated sample:
#' cs <- rnonorm(n = 100, f.coef = c(-0.3137, 0.8632, 0.3137, 0.0227))
#' hist(cs, probability = TRUE)
#' curve(dnorm(x), add = TRUE, col = "red")
#'
#' @importFrom stats rnorm
rnonorm <- function(n, f.coef){
  z <- rnorm(n)
  f.coef[1] + (f.coef[2] * z) + (f.coef[3] * z^2) + (f.coef[4] * z^3)
}

