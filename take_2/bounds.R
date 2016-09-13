library(Rglpk)
source("bounds_orig.R")

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

bounds.upper.internal = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = NULL, xmax = NULL, buckets = 100) {
	
	X = sort(X)
	
	n = length(X)
	
	if(is.null(xmin)) { xmin = min(X) }
	if(is.null(xmax)) { xmax = max(X) }
	
	if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
	
	xvals = seq(xmin, xmax, length.out = buckets + 1)

	orig.out = bounds(X, 1, sampling.ratio, solver="glpk")
	w = orig.out$weights_h
	Hhat.upper = approx(X, cumsum(w) / sum(w), xvals)$y
	Lhat.upper = get.logconc.upper(Hhat.upper)
	
	ret = data.frame(
		xvals=xvals,
		Fhat=ecdf(X)(xvals),
		Hhat.upper=Hhat.upper,
		Lhat.upper=Lhat.upper
	)
	
	return(ret)
}

bounds.upper.mean = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = min(X), xmax = max(X), buckets = 100) {
	
	ret = bounds.upper.internal(X, alpha, sampling.ratio, xmin, xmax, buckets)
	sum(ret$xvals * (c(ret$Lhat.upper[-1], 1) - ret$Lhat.upper))
	
}

bounds.logconc = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = min(X), xmax = max(X), buckets = 100) {

	up = bounds.upper.mean(X, alpha, sampling.ratio, xmin, xmax, buckets)
	low = -bounds.upper.mean(-X, alpha, sampling.ratio, -xmax, -xmin, buckets)
	
	c(low=low, up=up)
}

plot.ret = function(ret) {
	plot(ret$xvals, ret$Fhat)
	lines(ret$xvals, ret$Fhat.upper)
	lines(ret$xvals, ret$Hhat.upper, col = 4)
	lines(ret$xvals, ret$Lhat.upper, col = 2)
}

