#' Find the smallest log-concave upper envelope for a function
#' 
#' @param H The function H(x), specified at equally spaced x-values.
#' @return Lhat.upper The log-concave upper bound.
#' @export get.logconc.upper
get.logconc.upper = function(H) {
  
  # append a "1" to H, to make CDF = 1 outside of support
  H.plus = c(H, 1)
  
  zero.idx = sum(H.plus == 0)
  lH.short = log(H.plus[(zero.idx + 1) : length(H.plus)])
  
  cvx.hull.all = chull(1:length(lH.short), lH.short)
  one.idx = which(cvx.hull.all == 1)
  if(one.idx == 1) {
    cvx.hull.reord = cvx.hull.all
  } else {
    cvx.hull.reord = c(cvx.hull.all[one.idx:length(cvx.hull.all)], cvx.hull.all[1:(one.idx - 1)])
  }
  
  cvx.hull.pos = cvx.hull.reord[1:which(cvx.hull.reord == length(lH.short))]
  
  logconc.short = approx(cvx.hull.pos, lH.short[cvx.hull.pos], 1:length(lH.short))$y
  
  Lhat.upper = c(rep(0, zero.idx), exp(logconc.short))
  
  # remove artificially appended "1"
  Lhat.upper[-length(H.plus)]
}

#' @export get.logconc.upper.threshold
get.logconc.upper.threshold = function(Xhat.upper, sampling.ratio, threshold.idx) {
  K = length(Xhat.upper)
  fhat.upper = c(Xhat.upper[-1], 1) - Xhat.upper
  fhat.upper = fhat.upper / sum(fhat.upper)
  w = c(rep(1, threshold.idx), rep(sampling.ratio, K - threshold.idx))
  hhat.upper = fhat.upper * w / sum(fhat.upper * w)
  Hhat.upper = cumsum(hhat.upper)
  get.logconc.upper(Hhat.upper)
}

#' Computes upper identification interval under the assumption that F is log-concave.
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
#' @return Xhat.upper KS bound for Xhat.
#' @return Lhat.upper Lhat function from paper.
#' @return Weighted version of Xhat that maximizes mu, subject to log-concavity
#' 
#' @export bounds.logconc.internal
bounds.logconc.internal = function(X, sampling.ratio = 5,
                                   xmin = NULL, xmax = NULL, buckets = 1000,
                                   alpha = 1/sqrt(length(X))) {
  
  n = length(X)
  
  if(is.null(xmin)) { xmin = min(X) }
  if(is.null(xmax)) { xmax = max(X) }
  
  if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
  
  xvals = seq(xmin, xmax, length.out = buckets + 1)
  
  Xhat = ecdf(X)(xvals)
  Xhat.upper = pmax(Xhat - get.ks.threshold(n, alpha), 0)

  thresholds = quantile(X, seq(1/sampling.ratio/2, 1 - 1/sampling.ratio/2, length.out = 20))
  Xhat.candidates = lapply(thresholds, function(threshold) {
    Lhat.upper = get.logconc.upper.threshold(Xhat.upper, sampling.ratio, sum(xvals <= threshold))
    hajek.constrained(Xhat, xvals, sampling.ratio, lower.bound = Lhat.upper)
  })
  
  mu.bound = sapply(Xhat.candidates, function(Xhat) {
    sum(xvals * (Xhat - c(0, Xhat[-length(Xhat)])))
  })
  
  opt.idx = which.max(mu.bound)
  
  ret = list(mu.bound=mu.bound[opt.idx],
             raw=data.frame(
               xvals=xvals,
               Xhat=Xhat,
               Xhat.upper=Xhat.upper,
               Xhat.weighted=Xhat.candidates[[opt.idx]]
             ))
  
  return(ret)
}