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
#' @return Xhat Unweighted empirical CDF of the data.
#' @return xvals Points at which Xhat is evaluated.
#' @return Xhat.weighted Weighted version of Xhat that maximizes mu, subject to symmetry.
#' 
#' @export bounds.symmetric.internal
bounds.symmetric.internal = function(X, sampling.ratio = 5,
                                 xmin = NULL, xmax = NULL, buckets = 1000, alpha = 1/sqrt(length(X))) {
  
  n = length(X)
  
  if(is.null(xmin)) { xmin = min(X) }
  if(is.null(xmax)) { xmax = max(X) }
  
  if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
  
  xvals = seq(xmin, xmax, length.out = buckets + 1)
  
  Xhat = ecdf(X)(xvals)
  delta = qnorm(1 - alpha) * sqrt((1 + sampling.ratio) * (1 + 1/sampling.ratio) / 4 / n)
  center.candidates = quantile(X, seq(1/sampling.ratio/2, 1 - 1/sampling.ratio/2, length.out = 10))
  
  Xhat.candidates = lapply(center.candidates, function(center) {
    hajek.constrained.symmetric(Xhat, xvals, sampling.ratio, center, delta)
  })
  
  mu.bound = sapply(Xhat.candidates, function(Xhat) {
    sum(xvals * (Xhat - c(0, Xhat[-length(Xhat)])))
  })
  
  opt.idx = which.max(mu.bound)
  
  ret = list(mu.bound=mu.bound[opt.idx],
             raw=data.frame(
               xvals=xvals,
               Xhat=Xhat,
               Xhat.weighted=Xhat.candidates[[opt.idx]]
             ))
  
  return(ret)
}