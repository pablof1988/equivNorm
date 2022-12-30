#' Non arbitrary Epsilon (Simulation Function)
#'
#' Simulation function to compute a non arbitrary epsilon value for an
#' equivalence normality test of lack of fit.
#'
#' @param nSim number of simulations.
#' @param n1 Size of sample 1.
#' @param n2 Size of sample 2.
#' @param alpha Significance level of test.
#' @param seed simulation seed
#'
#' @return
#' A non arbitrary value of equivalence/irrelevance for a lack of fit test to
#' prove normality
#'
#' @details
#' epsilon resulting from this function should be used only for the purpose
#' of pretesting normality assumption when testing means difference
#' @export
#'
#' @examples
#' eps <- epsilon(nSim = 10000, n1 = 7, n2 = 3)
#' eps
#'
#' # If the epsilon value is adequate, then for a normal sample, the frequency
#' # with which an equivalence of normality test fails
#' # (it is not declared as normal) should be low.
#'
#' nfails <- replicate(10000, {
#'   data <- rnorm(10)
#'   test <- normequiv(data, epsilon = eps)
#'   test["upper bound"] >= test["eps^2"]
#' })
#' sum(nfails) / 10000
#'
#' @importFrom BinNonNor fleishman.coef

epsilon <- function(nSim, n1, n2, alpha = 0.05, seed = 123){
  n <- n1 + n2
  cochran <- alpha + c(-0.2 * alpha, 0.2 * alpha)
  ic <- rep(alpha, 2)
  z <- qnorm(alpha / 2, lower.tail = FALSE)
  i <- 0
  while (ic[1] >= cochran[1] & ic[2] <= cochran[2] & i <= 3) {
    coef <- try(as.vector(fleishman.coef(1, skewness.vec = i,
                                     kurtosis.vec = 3 * (i + 0.25))),
                silent = TRUE)
    if(is.numeric(coef)){
      sample <- rnonorm(n, f.coef = coef)
      tiep <- tiepT(nSim, n1, n2, f.coef = coef, alpha, seed)
      se <- z * sqrt((tiep * (1 - tiep)) / nSim)
      ic <- tiep + c(-se, se)
    } else{
      i <- i + 0.25
    }
    i <- i + 0.25
  }
  test <- normequiv(sample, alpha = alpha)
  as.numeric(sqrt(test[2]))
}
