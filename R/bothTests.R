#' Both, Traditionl (Chi square & Shapiro Wilk) and equivalence tests to detect normality
#'
#' Chi square & Shapiro Wilk (goodnes of fit) and equivalence test (lack of fit) to detect
#' normality
#'
#' @param x A numeric vector with sample data.
#' @param mu Theoretical mean from normal distribution.
#' @param sigma Theoretical standard deviation from normal distribution.
#' @param alpha Significance level of test.
#' @param epsilon Irrelevance/Equivalence limit.
#'
#'
#' @return
#' logical (TRUE if normality is detected, FALSE if not) value to, chi-square,
#' Shapiro - Wilk and equivalence normality test.
#' In addition, p-value (for chi square) and Upper
#' limit with \eqn{\epsilon^2} (for equivalence) are shown as attribute.
#'
#' @details
#' Since we are fitting to a normal distribution with
#' parameters \eqn{\mu}  and \eqn{\sigma} known, chi-square statistic
#' uses k-1 degree of freedom. It is a fact only validate for
#' purposes of paper in reference, in real world parameters of a data set are
#' always unknown and chi-square statistics to prove normality uses k-3 degree
#' of freedom.
#'
#' Normality is detected in chi-square test when p-value is greater than the
#' significance level \code{alpha} and for the equivalence test we declare
#' normality when upper limit of euclidean square distance is less than
#' irrelevance limit \eqn{\epsilon^2} \code{epsilon}^2
#'
#' @export
#'
#' @examples
#' # Proving if data comes from a normal standard
#' set.seed(123)
#' data <- rnorm(30)
#' bothTests(x = data)
#'
#' # or from a normal N(mu = 5, sigma = 2) changing the default parameters
#' set.seed(123)
#' data <- rnorm(30, mean = 5, sd = 2)
#' bothTests(x = data, alpha = 0.01, epsilon = 0.15,
#'           mu = 5, sigma = 2)
#'
#' @importFrom stats qnorm pnorm pchisq
#' @importFrom graphics hist

bothTests <- function(x, mu = 0, sigma = 1, alpha = 0.05, epsilon = 0.4){
  n <- length(x)
  h <- hist(x, plot = F)
  breaks <- h$breaks
  nbreaks <- length(breaks)
  theoretical <- NULL
  for(i in 1:(nbreaks-1)){
    theoretical[i] <- pnorm(breaks[i + 1], mean = mu, sd = sigma) - pnorm(breaks[i], mean = mu, sd = sigma)
  }
  theoretical <- c(pnorm(breaks[1], mean = mu, sd = sigma), theoretical,
                   pnorm(breaks[nbreaks], lower.tail = F, mean = mu, sd = sigma))
  afr <- c(0, h$counts, 0)
  k <- length(afr)

  chi <- sum((afr - (theoretical*n))^2 / (theoretical*n))
  pvalchi <- pchisq(chi, k - 1, lower.tail = F)

  pvalsw <- shapiro.test(x)$p.val

  de <- afr / n
  deuc <- sum((de - theoretical)^2)

  u <- qnorm(alpha, lower.tail = F)
  fac1 <- sum(((de - theoretical)^2) * de)
  fac2 <- 0
  for(j1 in 1:k){
    for(j2 in 1:k){
      fac2 <- fac2 + (de[j1] - theoretical[j1]) * (de[j2] - theoretical[j2]) * de[j1] * de[j2]
    }
  }
  vn2 <- 4 * (fac1 - fac2)
  ee <- sqrt(vn2 / n)
  lupp <- deuc + (u * ee)

  res <- c('X2Test' = ifelse(pvalchi >= alpha, T, F),
           'swTest' = ifelse(pvalsw >= alpha, T, F),
           'equivTest' = ifelse(lupp < epsilon^2, T, F))
  attr(res, "Pval-Chi2, Pval-SW, Upper Limit, eps^2") <- c(pvalchi, pvalsw, lupp, epsilon^2)
  res
}
