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