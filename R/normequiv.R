#' Lack of fit test to prove normality
#'
#' Lack of fit test based on equivalence approach to prove normality
#'
#' @param x a numeric vector with sample data.
#' @param mu theoretical mean from normal distribution.
#' @param sigma theoretical standard deviation from normal distribution.
#' @param alpha significance level of test.
#' @param epsilon irrelevance/equivalence limit.
#'
#'
#' @return
#' a vector with Euclidean distance, upper bound of confidence interval
#' and irrelevance limit \eqn{\epsilon^2}. Decision (Normality or non normality) over the
#' is showed as an attribute
#'
#' @details
#' In \eqn{H_0:d^2(\pi, \pi^0) \geq \epsilon^2} (lack of fit) Vs
#' \eqn{H_1:d^2(\pi, \pi^0) < \epsilon^2}, null hypothesis \eqn{H_0} is
#' rejected (proving normality) when upper bound of confidence interval
#' is less than \eqn{\epsilon^2}
#'
#' @export
#'
#' @examples
#' # Proving if data comes from a normal standard
#' set.seed(123)
#' data <- rnorm(30)
#' normequiv(x = data)
#'
#' # or from a normal N(mu = 5, sigma = 2) changing the default parameters
#' set.seed(123)
#' data <- rnorm(30, mean = 5, sd = 2)
#' normequiv(x = data, alpha = 0.01, epsilon = 0.15,
#'          mu = 5, sigma = 1)
#'
#' @importFrom stats qnorm pnorm
#' @importFrom graphics hist
normequiv <- function(x, mu = 0, sigma = 1, alpha = 0.05, epsilon = 0.4){
  n <- length(x)
  h <- hist(x, plot = F)
  breaks <- h$breaks
  nbreaks <- length(breaks)
  theoretical <- NULL
  for(i in 1:(nbreaks-1)){
    theoretical[i] <- pnorm(breaks[i + 1], mean = mu, sd = sigma) - pnorm(breaks[i], mean = mu, sd = sigma)
  }
  theoretical <- c(pnorm(breaks[1],  mean = mu, sd = sigma), theoretical,
                   pnorm(breaks[nbreaks], lower.tail = F,  mean = mu, sd = sigma))
  afr <- c(0, h$counts, 0)
  k <- length(afr)
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
  res <- c("Euclidean Distance" = deuc, "upper bound" = lupp,  "eps^2" = epsilon^2)
  dec <- ifelse(lupp < epsilon^2, "Reject H0 - Normality", "Not Reject H0 - Non Normality")
  attr(res, "Decision") <- dec
  res
}
