#' Find (one-sided) length of Kolmogorov-Smirnov bounds for a CDF
#' 
#' @param n The number of samples used to compute Fhat.
#' @param alpha The significance level.
#' @return The KS threshold
#' @export get.ks.threshold
get.ks.threshold = function(n, alpha = 0.05) {
  ks.thresh = uniroot(
    function(D) { kolmim::pkolmim(D, n) - 1 + alpha},
    c(0.3/sqrt(n), 2/sqrt(n))
  )$root
  ks.thresh
}

#' @export print.bounds
print.bounds = function(obj) {
  mu.int = obj[[1]]
  print(paste0("Identification interval for mu: (", signif(mu.int[1], 3), ", ", signif(mu.int[2], 3), ")"))
}