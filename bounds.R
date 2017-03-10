library(kolmim)
library(lpSolve)
library(zoo) # for plotting only
library(splines) #for plotting only...

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

# Inspired by "A Modified Kolmogorov-Smirnov Test Sensitive to Tail Alternatives"
# By Mason & Schuenemeyer

get.mod.ks.upper = function(Fhat, n, alpha = 0.05) {
	
	if(alpha != 0.05) { stop("Only alpha = 0.05 implemented...") }
	
	delta = 2.8 #calibrated by hand for alpha = 0.05.... 
	Fhat.AD = pmax(Fhat - 2 * delta * sqrt(Fhat * (1 - Fhat) / n), 0) #control body
	Fhat.tail = pmax(Fhat - delta * (Fhat * (1 - Fhat)), 0) #control tail
	pmax(Fhat.AD, Fhat.tail)
}

get.ks.upper = function(Fhat, n, alpha = 0.05) {
	
	ks.thresh = uniroot(
		function(D) { pkolmim(D, n) - 1 + alpha},
		c(0.3/sqrt(n), 2/sqrt(n))
	)$root
	
	pmax(Fhat - ks.thresh, 0)
}

get.hajek.upper = function(Fhat.upper, sampling.ratio) {
	Q = 1/sampling.ratio
	sapply(Fhat.upper, function(u) {
		Q * u / (1 - u + Q * u)
	})
}

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
	
	lp.out = lp(direction="max",
				objective.in= objective.in,
			    const.mat=const.mat,
			    const.dir=c(rep(">=", length(const.rhs) - 1), "=="),
			    const.rhs=const.rhs)
			    
	Fhat.weighted = cumsum(lp.out$solution[-(K+1)] * fhat)
	Fhat.weighted
	
}

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

bounds.upper.internal = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = NULL, xmax = NULL, buckets = 1000) {
	
	n = length(X)
	
	if(is.null(xmin)) { xmin = min(X) }
	if(is.null(xmax)) { xmax = max(X) }
	
	if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
	
	xvals = seq(xmin, xmax, length.out = buckets + 1)

	Fhat = ecdf(X)(xvals)
	Fhat.upper = get.ks.upper(Fhat, n, alpha)
	Hhat.upper = get.hajek.upper(Fhat.upper, sampling.ratio)
	#Lhat.upper = get.logconc.upper(Hhat.upper)
	Lhat.upper = get.logconc.upper.scan(Fhat.upper, sampling.ratio)
	Fhat.weighted = hajek.constrained(Fhat, xvals, Lhat.upper, sampling.ratio)
	Fhat.AL = hajek.constrained(Fhat, xvals, 0 * Fhat, sampling.ratio)
	
	ret = data.frame(
		xvals=xvals,
		Fhat=Fhat,
		Fhat.upper=Fhat.upper,
		Hhat.upper=Hhat.upper,
		Lhat.upper=Lhat.upper,
		Fhat.weighted=Fhat.weighted, 
		Fhat.AL = Fhat.AL
	)
	
	return(ret)
}

bounds.upper.mean = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = min(X), xmax = max(X), buckets = 1000) {
	
	ret = bounds.upper.internal(X, alpha, sampling.ratio, xmin, xmax, buckets)
	K = length(ret$Fhat.weighted)
	c(sum(ret$xvals * (ret$Fhat.weighted - c(0, ret$Fhat.weighted[-K]))),
	 sum(ret$xvals * (ret$Fhat.AL - c(0, ret$Fhat.AL[-K]))))
	
}

bounds.logconc = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = min(X), xmax = max(X), buckets = 1000) {

	up = bounds.upper.mean(X, alpha, sampling.ratio, xmin, xmax, buckets)
	low = -bounds.upper.mean(-X, alpha, sampling.ratio, -xmax, -xmin, buckets)
	
	c(low=low, up=up)
}

plot.ret = function(ret) {
	plot(ret$xvals, ret$Fhat, lwd = 2, xlab="y", ylab="F(y)")
	lines(ret$xvals, ret$Fhat.upper, lwd = 2)
	lines(ret$xvals, ret$Lhat.upper, col = 6, lwd = 2)
	lines(ret$xvals, ret$Fhat.weighted, col = 2, lwd = 2)
	lines(ret$xvals, ret$Fhat.AL, col = 4, lwd = 2, lty = 1)
	legend("topleft", c("S.hat", "S.hat+", "L.hat", "F.hat (us)", "F.hat (AL)"), pch = c(1, NA, NA, NA, NA), lwd = c(NA, 2, 2, 2, 2), col = c(1, 1, 6, 2, 4))
}

zeropad = function(x) c(rep(0, 100), x, rep(0, 100))
unpad = function(x) (x[101:(length(x) - 100)])

plot.ret2 = function(ret, df=9) {
	
	K = length(ret$Fhat.weighted)
	
	fhat = ret$Fhat - c(0, ret$Fhat[-K])
	fhat.weighted = ret$Fhat.weighted - c(0, ret$Fhat.weighted[-K])
	fhat.AL = ret$Fhat.AL - c(0, ret$Fhat.AL[-K])
	
	fhat.smooth = unpad(predict(glm(zeropad(fhat) ~ ns(1:(length(fhat) + 200), df = df), family=quasi(link = "log", variance = "mu")), type = "response"))
	
	fhat.weighted.ratio = smooth.spline(na.approx(fhat.weighted / fhat), df = df)$y
	fhat.AL.ratio = na.approx(fhat.AL / fhat)
	
	plot(ret$xvals, fhat.smooth, lwd = 2, ylim = range(fhat.smooth, fhat.smooth * fhat.weighted.ratio, fhat.smooth * fhat.AL.ratio), xlab="y", ylab="f(y)")
	lines(ret$xvals, fhat.smooth * fhat.weighted.ratio, col = 2, lwd = 2)
	lines(ret$xvals, fhat.smooth * fhat.AL.ratio, col = 4, lwd = 2)
	legend("topleft", c("Uncorrected", "AL + Log-concave constr.", "Plain Aronow-Lee"), lwd = 2, col = c(1, 2, 4))
}

