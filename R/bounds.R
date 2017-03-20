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
    
    low.obj = bounds.logconc.internal(-X, sampling.ratio, -xmax, -xmin, buckets, alpha)
    high.obj = bounds.logconc.internal(X, sampling.ratio, xmin, xmax, buckets, alpha)
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
