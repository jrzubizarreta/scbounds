library(kolmim)
library(lpSolve)

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

#' Find level-alpha Kolmogorov-Smirnov bounds for a CDF F
#' 
#' @param Fhat The empirical CDF.
#' @param n The number of samples used to compute Fhat.
#' @param alpha The significance level.
#' @return The ''upper'' KS bound (i.e., the one that lies below Fhat).
#' @export get.ks.upper
get.ks.upper = function(Fhat, n, alpha = 0.05) {
  ks.thresh = uniroot(
    function(D) { kolmim::pkolmim(D, n) - 1 + alpha},
    c(0.3/sqrt(n), 2/sqrt(n))
  )$root
  
  pmax(Fhat - ks.thresh, 0)
}

#' Main workhorse function. Maximizes the value of hat{mu} among all Hajek
#' ratio estimators that satisfy the following properties:
#' - The ratio of the sampling weights gamma is bounded
#' - The resulting weighted CDF lies "above" lower.bound (as a function)
#' 
#' @param Fhat The raw, unweighted, empirical CDF.
#' @param xvals The points at which Fhat is specified.
#' @param lower.bound The lower bound for the Hajek-weighted empirical CDF.
#' @param sampling.ratio The bound on the sampling ratio.
#' 
#' @return Fhat.weighted A weighted version of Fhat that maximizes hat{mu}
#' 
#' @export hajek.constrained
hajek.constrained = function(Fhat, xvals, lower.bound, sampling.ratio = 5) {
  
  K = length(Fhat)
  fhat = Fhat - c(0, Fhat[-K])
  
  objective.in = c(xvals * fhat, 0)
  const.mat = rbind(
    cbind(outer(1:K, 1:K, function(x, y) as.numeric(x >= y) * fhat[y]), 0),
    cbind(diag(1, K, K), -1),
    cbind(diag(-1, K, K), sampling.ratio),
    c(rep(0, K), 1),
    c(fhat, 0)
  )
  const.rhs = c(lower.bound, rep(0, 2*K + 1), 1)
  
  lp.out = lpSolve::lp(direction="max",
                       objective.in= objective.in,
                       const.mat=const.mat,
                       const.dir=c(rep(">=", length(const.rhs) - 1), "=="),
                       const.rhs=const.rhs)
  
  Fhat.weighted = cumsum(lp.out$solution[-(K+1)] * fhat)
  Fhat.weighted
}

#' Finds the "Uhat" function from the paper, i.e., the pointwise minimum of all log-concave
#' function that lie above some weighted version of Fhat.upper.
#' 
#' @param Fhat.upper The unweighted function.
#' @param sampling.ratio Bound on the sampling weights.
#' @return Lhat.upper The pointwise minimal bound.
#' 
#' @export get.logconc.upper.scan
get.logconc.upper.scan = function(Fhat.upper, sampling.ratio) {
  
  K = length(Fhat.upper)
  fhat.upper = Fhat.upper - c(0, Fhat.upper[-K])
  
  all.idx = 1:(K - 1)
  good.idx = which((Fhat.upper[-K] > 0) & (all.idx %% 5 == 0))
  
  Lhat.candidates = sapply(good.idx, function(thresh) {
    w = c(rep(1, thresh), rep(sampling.ratio, K - thresh))
    hhat.upper = fhat.upper * w / sum(fhat.upper * w)
    Hhat.upper = cumsum(hhat.upper)
    get.logconc.upper(Hhat.upper)
  })
  Lhat.upper = apply(Lhat.candidates, 1, min)
  Lhat.upper
}

#' Computes upper identification interval under the assumption that F is log-concave.
#' 
#' @param X The observed data.
#' @param alpha Significance level used for KS bounds.
#' @param sampling.ratio Bound on the sampling weights gamma.
#' @param xmin Used to construct histogram representation.
#' @param xmax Used to construct histogram representation.
#' @param buckets Used to construct histogram representation.
#' 
#' @return mu.bound The upper bound for mu(x).
#' 
#' @return Fhat Unweighted empirical CDF of the data.
#' @return xvals Points at which Fhat is evaluated.
#' @return Fhat.upper KS bound for Fhat.
#' @return Lhat.upper Lhat function from paper.
#' @return Weighted version of Fhat that maximizes mu, subject to log-concavity
#' 
#' @export bounds.logconc.internal
bounds.logconc.internal = function(X, alpha = 1/length(X), sampling.ratio = 5,
                                   xmin = NULL, xmax = NULL, buckets = 1000) {
  
  n = length(X)
  
  if(is.null(xmin)) { xmin = min(X) }
  if(is.null(xmax)) { xmax = max(X) }
  
  if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
  
  xvals = seq(xmin, xmax, length.out = buckets + 1)
  
  Fhat = ecdf(X)(xvals)
  Fhat.upper = get.ks.upper(Fhat, n, alpha)
  Lhat.upper = get.logconc.upper.scan(Fhat.upper, sampling.ratio)
  Fhat.weighted = hajek.constrained(Fhat, xvals, Lhat.upper, sampling.ratio)
  
  mu.bound = sum(xvals * (Fhat.weighted - c(0, Fhat.weighted[-length(Fhat.weighted)])))
  
  ret = list(mu.bound=mu.bound,
             raw=data.frame(
               xvals=xvals,
               Fhat=Fhat,
               Fhat.upper=Fhat.upper,
               Lhat.upper=Lhat.upper,
               Fhat.weighted=Fhat.weighted
             ))
  
  return(ret)
}

#' Computes upper identification interval without shape constraints, as
#' in Aronow & Lee (2013).
#' 
#' @param X The observed data.
#' @param sampling.ratio Bound on the sampling weights gamma.
#' @param xmin Used to construct histogram representation.
#' @param xmax Used to construct histogram representation.
#' @param buckets Used to construct histogram representation.
#' 
#' @return mu.bound The upper bound for mu(x).
#' 
#' @return Fhat Unweighted empirical CDF of the data.
#' @return xvals Points at which Fhat is evaluated.
#' @return Weighted version of Fhat that maximizes mu, subject to sampling ratio constraint.
#' 
#' @export bounds.plain.internal
bounds.plain.internal = function(X, sampling.ratio = 5,
                                   xmin = NULL, xmax = NULL, buckets = 1000) {
  
  n = length(X)
  
  if(is.null(xmin)) { xmin = min(X) }
  if(is.null(xmax)) { xmax = max(X) }
  
  if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
  
  xvals = seq(xmin, xmax, length.out = buckets + 1)
  
  Fhat = ecdf(X)(xvals)
  Fhat.AL = hajek.constrained(Fhat, xvals, 0 * Fhat, sampling.ratio)
  
  mu.bound = sum(xvals * (Fhat.AL - c(0, Fhat.AL[-length(Fhat.AL)])))
  
  ret = list(mu.bound=mu.bound,
             raw=data.frame(
               xvals=xvals,
               Fhat=Fhat,
               Fhat.weighted=Fhat.AL
             ))
  
  return(ret)
}
  
  
#' Main function. Computes identification interval.
#' 
#' @param X The observed data.
#' @param sampling.ratio Bound on the sampling weights gamma.
#' @param constraint Optional constraint on the shape of F.
#' @param alpha Significance level used for KS bounds.
#' @param xmin Used to construct histogram representation.
#' @param xmax Used to construct histogram representation.
#' @param buckets Used to construct histogram representation.
#' 
#' @return mu.interval The identification interval for mu(x).
#' @return lower.bound.internal Internal data used to create lower bound.
#' @return upper.bound.internal Internal data used to create upper bound.
#' 
#' @export bounds
bounds = function(X, sampling.ratio = 5, constraint = c("none", "logconcave", "symmetric", "gaussian"),
                  alpha = 1/length(X), xmin = min(X), xmax = max(X), buckets = 1000) {
  
  constraint = match.arg(constraint)
  
  if (constraint == "none") {

    low.obj = bounds.plain.internal(-X, sampling.ratio, -xmax, -xmin, buckets)
    high.obj = bounds.plain.internal(X, sampling.ratio, xmin, xmax, buckets)
    interval = c(-low.obj$mu.bound, high.obj$mu.bound)
    
  } else if (constraint == "logconcave") {
    
    low.obj = bounds.logconc.internal(-X, alpha, sampling.ratio, -xmax, -xmin, buckets)
    high.obj = bounds.logconc.internal(X, alpha, sampling.ratio, xmin, xmax, buckets)
    interval = c(-low.obj$mu.bound, high.obj$mu.bound)
    
  } else if (constraint == "symmetric") {
    stop()
  } else if (constraint == "gaussian") {
    stop()
  }
  
  ret = list(mu.interval=interval,
             constraint=constraint,
             lower.bound.internal=low.obj,
             upper.bound.internal=high.obj)
  class(ret) = "bounds"
  return(ret)
}

#' @export print.bounds
print.bounds = function(obj) {
  mu.int = obj[[1]]
  print(paste0("Identification interval for mu: (", signif(mu.int[1], 3), ", ", signif(mu.int[2], 3), ")"))
}
