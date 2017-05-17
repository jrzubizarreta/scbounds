#' Computational workhorse function. Maximizes the value of hat{mu} among all Hajek
#' ratio estimators that satisfy the following properties:
#' - The ratio of the sampling weights gamma is bounded
#' - The resulting weighted CDF lies "between" lower.bound and upper.bound (as a function)
#' 
#' @param Xcdf The raw, unweighted, empirical CDF.
#' @param xvals The points at which Xcdf is specified.
#' @param sampling.ratio The bound on the sampling ratio.
#' @param lower.bound A lower bound for the Hajek-weighted empirical CDF.
#' @param upper.bound An upper bound for the Hajek-weighted empirical CDF.
#' 
#' @return Xcdf.weighted A weighted version of Xcdf that maximizes hat{mu}
#' 
#' @export hajek.constrained
hajek.constrained = function(Xcdf, xvals, sampling.ratio, lower.bound = NULL, upper.bound = NULL) {
  
  K = length(Xcdf)
  X.density = Xcdf - c(0, Xcdf[-K])
  
  # Maximize the sum with weights wi subject to t <= wi <= gamma *t
  objective.in = c(xvals * X.density, 0)
  const.mat = rbind(
    cbind(diag(1, K, K), -1),
    cbind(diag(-1, K, K), sampling.ratio),
    c(X.density, 0)
  )
  const.rhs = c(rep(0, 2*K), 1)
  const.dir=c(rep(">=", 2*K), "==")
  
  if (!is.null(lower.bound)) {
    const.mat = rbind(const.mat,
                      cbind(outer(1:K, 1:K, function(x, y) as.numeric(x >= y) * X.density[y]), 0))
    const.rhs = c(const.rhs, lower.bound)
    const.dir = c(const.dir, rep(">=", length(lower.bound)))
  }
  
  if (!is.null(upper.bound)) {
    const.mat = rbind(const.mat,
                      cbind(outer(1:K, 1:K, function(x, y) as.numeric(x >= y) * X.density[y]), 0))
    const.rhs = c(const.rhs, upper.bound)
    const.dir = c(const.dir, rep("<=", length(upper.bound)))
  }
  
  lp.out = lpSolve::lp(direction="max",
                       objective.in= objective.in,
                       const.mat=const.mat,
                       const.dir=const.dir,
                       const.rhs=const.rhs)
  
  Xcdf.weighted = cumsum(lp.out$solution[-(K+1)] * X.density)
  Xcdf.weighted
}

#' Computational workhorse function. Maximizes the value of hat{mu} among all Hajek
#' ratio estimators that satisfy the following properties:
#' - The ratio of the sampling weights gamma is bounded
#' - The resulting CDF is symmetric to within tolerance delta.
#' 
#' @param Xcdf The raw, unweighted, empirical CDF.
#' @param xvals The points at which Xcdf is specified.
#' @param sampling.ratio The bound on the sampling ratio.
#' @param center The center of symmetry.
#' @param delta The tolerance parameter.
#' 
#' @return Xcdf.weighted A weighted version of Xcdf that maximizes hat{mu}
#' 
#' @export hajek.constrained.symmetric
hajek.constrained.symmetric = function(Xcdf, xvals, sampling.ratio, center, delta) {
  
  K = length(Xcdf)
  X.density = Xcdf - c(0, Xcdf[-K])
  
  center.bucket = which.min(abs(xvals - center))[1]
  R.min = min(center.bucket - 1, K - center.bucket)
  
  # Computes F.hat(center + y) for y in the range covered by (1:R.min)
  integrate.up = outer(center.bucket + (1:R.min), 1:K, function(x, y) as.numeric(x >= y) * X.density[y])
  # Computes F.hat(center - y) for y in the range covered by (1:R.min)
  integrate.down = outer(center.bucket - (1:R.min), 1:K, function(x, y) as.numeric(x >= y) * X.density[y])
  
  # Symmetry implies that M * f is nearly 1 (to within tolerance 2*delta)
  M = integrate.up + integrate.down
  
  objective.in = c(xvals * X.density, 0)
  const.mat = rbind(
    cbind(M, 0),
    cbind(-M, 0),
    cbind(diag(1, K, K), -1),
    cbind(diag(-1, K, K), sampling.ratio),
    c(X.density, 0)
  )
  const.rhs = c(rep(1 - delta, R.min),
                rep(-1 - delta, R.min),
                rep(0, 2*K),
                1)
  lp.out = lpSolve::lp(direction="max",
                       objective.in= objective.in,
                       const.mat=const.mat,
                       const.dir=c(rep(">=", length(const.rhs) - 1), "=="),
                       const.rhs=const.rhs)
  
  if (lp.out$status != 0) return(NA)
  
  Xcdf.weighted = cumsum(lp.out$solution[-(K+1)] * X.density)
  Xcdf.weighted
}