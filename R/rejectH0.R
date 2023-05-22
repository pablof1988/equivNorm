#' H0 rejection proportion for all tests (Simulation Function).
#'
#' Simulation function to compute the proportion of rejections of H0 for four
#' comparisson mean tests: i) t-student direct (without pretest), ii) wilcoxon
#' direct (without pretest), iii) t-student or wilcoxon depending of chi-square
#' (goodnes of fit) pretest for normality, iv) t-student or wilcoxon depending
#' of Shapiro-Wilk (goodnes of fit) pretest for normality and v) t-student or
#' wilcoxon depending of equivalence (lack of fit) pretest for normality
#'
#' @param nSim number of simulations.
#' @param n1 Size of sample 1.
#' @param n2 Size of sample 2.
#' @param epsilon Irrelevance/Equivalence limit
#' @param diffMean Theoretical difference of means (= 0 by default)
#' @param f.coef Fleishman coefficients indicating normal contamination degree.
#' @param alpha Significance level of test.
#' @param seed simulation seed
#'
#' @return
#' An estimation of Type I Error Probability TIEP (if diffMean = 0) or Power
#' (if diffMean > 0) of tests.
#'
#' @details
#' when  \code{diffMean} is 0, function computes rejections proportions of
#' \eqn{H_0:\mu_1-\mu_2=0} when it is true i.e. TIEP,  while if \code{diffMean}
#' is different to zero function computes rejections proportions of
#' \eqn{H_0:\mu_1-\mu_2=0} when it is false i.e. Power.
#' @export
#'
#' @examples
#' # TIEP in a non contamination sample (Normal):
#' rejectH0(nSim = 1000, n1 = 5, n2 = 3, epsilon = 0.4, diffMean = 0)
#'
#' # TIEP in a non normal sample (contaminated)
#' rejectH0(nSim = 1000, n1 = 5, n2 = 3, epsilon = 0.4, diffMean = 0,
#'          f.coef = c(-0.3137, 0.8632, 0.3137, 0.0227))
#'
#' # POWER in a non normal sample (contaminated)
#' rejectH0(nSim = 1000, n1 = 5, n2 = 3, epsilon = 0.4, diffMean = 3,
#'          f.coef = c(-0.3137, 0.8632, 0.3137, 0.0227))
#'
#' @importFrom stats t.test wilcox.test

rejectH0 <- function(nSim, n1, n2, epsilon, diffMean = 0,
                     f.coef = c(0, 1, 0, 0), alpha = 0.05, seed = 123){
  n <- n1 + n2
  facts <- factor(rep(1:2, c(n1, n2)))
  mu <- rep(c(0, diffMean), c(n1, n2))
  set.seed(seed)
  tests <- replicate(nSim, {
    data <- rnonorm(n, f.coef = f.coef) + mu
    pretest <- bothTests(data, epsilon = epsilon)
    tstud <- t.test(data ~ facts, var.equal = T, conf.level = 1 - alpha)$p.value
    wilc <- wilcox.test(data ~ facts, conf.level = 1 - alpha)$p.value
    c('t_direct' = tstud, 'wilcx_direct' = wilc,
      ifelse(pretest == T, tstud, wilc))

  })
  rowSums(tests < alpha, na.rm = T) / nSim
}
