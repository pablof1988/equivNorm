#' H0 rejection proportion  for t-student (Simulation Function)
#'
#' Simulation function to compute TIEP estimation in a t-student test for
#' means comparison with normal or non-normal samples
#'
#' @param nSim number of simulations.
#' @param n1 Size of sample 1.
#' @param n2 Size of sample 2.
#' @param f.coef Fleishman coefficients indicating normal contamination degree.
#' @param alpha Significance level of test.
#' @param seed simulation seed
#'
#' @return
#' An estimation of Type I Error Probability TIEP for t-student test to
#' means comparison.
#'
#' @export
#'
#' @examples
#' # TIEP in a non contamination sample (Normal):
#' tiepT(nSim = 1000, n1 = 5, n2 = 3)
#'
#' # TIEP in a non normal sample (contaminated)
#' tiepT(nSim = 1000, n1 = 5, n2 = 3,
#'      f.coef = c(-0.3137, 0.8632, 0.3137, 0.0227))
#'
#' @importFrom stats t.test

tiepT <- function(nSim, n1, n2, f.coef = c(0, 1 ,0, 0), alpha = 0.05,
                  seed = 123){
  n <- n1 + n2
  facts <- factor(rep(1:2, c(n1, n2)))
  set.seed(seed)
  pvals <- replicate(nSim, {
    muestra <- rnonorm(n, f.coef)
    t.test(muestra ~ facts, var.equal = T, mu = 0, conf.level = 1 - alpha)$p.value
  })
  sum(pvals < alpha) / nSim
}
