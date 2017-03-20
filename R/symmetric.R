#' Computes upper identification interval with symmetry constraints.
#' 
#' @param X The observed data.
#' @param sampling.ratio Bound on the sampling weights gamma.
#' @param xmin Used to construct histogram representation.
#' @param xmax Used to construct histogram representation.
#' @param buckets Used to construct histogram representation.
#' @param alpha Significance level used for KS bounds.
#' 
#' @return mu.bound The upper bound for mu(x).
#' 
#' @return Fhat Unweighted empirical CDF of the data.
#' @return xvals Points at which Fhat is evaluated.
#' @return Fhat.weighted Weighted version of Fhat that maximizes mu, subject to symmetry.
#' 
#' @export bounds.symmetric.internal
bounds.symmetric.internal = function(X, sampling.ratio = 5,
                                 xmin = NULL, xmax = NULL, buckets = 1000, alpha = 1/sqrt(length(X))) {
  
  n = length(X)
  
  if(is.null(xmin)) { xmin = min(X) }
  if(is.null(xmax)) { xmax = max(X) }
  
  if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
  
  xvals = seq(xmin, xmax, length.out = buckets + 1)
  
  Fhat = ecdf(X)(xvals)
  delta = qnorm(1 - alpha) * sqrt((1 + sampling.ratio) * (1 + 1/sampling.ratio) / 4 / n)
  delta = delta / 2 # TODO: is this justified?
  center.candidates = quantile(X, seq(1/sampling.ratio/2, 1 - 1/sampling.ratio/2, length.out = 10))
  
  Fhat.candidates = lapply(center.candidates, function(center) {
    hajek.constrained.symmetric(Fhat, xvals, center, delta, sampling.ratio)
  })
  
  mu.bound = sapply(Fhat.candidates, function(Fhat) {
    sum(xvals * (Fhat - c(0, Fhat[-length(Fhat)])))
  })
  
  opt.idx = which.max(mu.bound)
  
  ret = list(mu.bound=mu.bound[opt.idx],
             raw=data.frame(
               xvals=xvals,
               Fhat=Fhat,
               Fhat.weighted=Fhat.candidates[[opt.idx]]
             ))
  
  return(ret)
}