#' Lack of fit test
#'
#' Lack of fit test based on equivalence approach to detect if a sample comes
#' from a target distribution
#'
#' @param probs Theoretical probabilities from a target distribution (for example: binomial, normal, etc).
#' @param afr Absolute frequencies obtained from a sample
#' @param alpha Significance level of test
#' @param epsilon Irrelevance/Equivalence limit
#'
#' @return
#' a vector with Euclidean distance and upper bound of confidence interval.
#' Decision (Reject or not the null hypothesis) is showed as an attribute
#'
#' @details
#' In \eqn{H_0:d^2(\pi, \pi^0) \geq \epsilon^2} (lack of fit) Vs \eqn{H_1:d^2(\pi, \pi^0) < \epsilon^2},
#' null hypothesis \eqn{H_0} is rejected (proving fit to target distribution)
#' when upper bound of confidence interval is less than \eqn{\epsilon^2}.
#'
#' @export
#'
#' @examples
#' # Proving if x comes from a binomial distribution (n = 3, p = 0.45)
#' set.seed(123)
#' x <- rbinom(100, size = 3, prob = 0.45)
#' frequencies <- table(x)
#' binomProbs <- dbinom(0:3, size = 3, prob = 0.45)
#' lackFitTest(probs = binomProbs, afr = frequencies)
#'
#' # or changing the default values
#' lackFitTest(probs = binomProbs, afr = frequencies,
#'             alpha = 0.01, epsilon = 0.15)
#'
#' @importFrom stats qnorm
lackFitTest <- function(probs, afr, alpha = 0.05, epsilon = 0.4){
  n <- sum(afr)
  k <- length(afr)
  de <- afr / n
  deuc <- sum((de - probs)^2)
  u <- qnorm(alpha, lower.tail = FALSE)
  fac1 <- sum(((de - probs)^2) * de)
  fac2 <- 0
  for(j1 in 1:k){
    for(j2 in 1:k){
      fac2 <- fac2 + (de[j1] - probs[j1]) * (de[j2] - probs[j2]) * de[j1] * de[j2]
    }
  }
  vn2 <- as.numeric(4 * (fac1 - fac2))
  ee <- sqrt(vn2 / n)
  lupp <- deuc + (u * ee)
  res <- c("Euclidean Distance" = deuc, "upper bound" = lupp)
  dec <- ifelse(lupp >= epsilon^2, "Not Reject H0", "Reject H0")
  attr(res, "Decision") <- dec
  res
}
