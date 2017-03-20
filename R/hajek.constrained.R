#' Computational workhorse function. Maximizes the value of hat{mu} among all Hajek
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

#' Computational workhorse function. Maximizes the value of hat{mu} among all Hajek
#' ratio estimators that satisfy the following properties:
#' - The ratio of the sampling weights gamma is bounded
#' - The resulting CDF is symmetric to within tolerance delta.
#' 
#' @param Fhat The raw, unweighted, empirical CDF.
#' @param xvals The points at which Fhat is specified.
#' @param center The center of symmetry.
#' @param delta The tolerance parameter.
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