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

get.berkjones.upper = function(Fhat, n, alpha = 0.05) {
	# This is an extremely inaccurate threshold, based on
	# crude scaling starting from an accurate calculation at n = 100.
	# The software in C is available at
	# http://www.stat.washington.edu/jaw/RESEARCH/SOFTWARE/software.list.html
	# TODO: Integrate exact C code into our package.
	# (This is the Berk-Jones test, based on Owen's paper)
	# (See also Jager/Wellner)
	delta.BJ = 0.53766 / sqrt(n)

	BJ.fun = function(A) {
		if(A == 0) return(0)
		if(A == 1) {A = 1 - 1/n}
		uniroot(
			function(B) (A*log(A/B) + (1 - A)*log((1 - A) / (1 - B)) - delta.BJ),
			interval = c(0, A)
		)$root
	}
	
	Fhat.upper = sapply(Fhat, Vectorize(BJ.fun))
	Fhat.upper
}

get.hajek.upper = function(Fhat.upper, sampling.ratio) {
	rho = 1/sampling.ratio
	sapply(Fhat.upper, function(u) {
		rho * u / (1 - u + rho * u)
	})
}


bounds.upper.internal = function(X, alpha = 0.05, sampling.ratio = 5,
	xmin = NULL, xmax = NULL, buckets = 100) {
	
	n = length(X)
	
	if(is.null(xmin)) { xmin = min(X) }
	if(is.null(xmax)) { xmax = max(X) }
	
	if(xmin > min(X) | xmax < max(X)) { stop ("support too short") }
	
	xvals = seq(xmin, xmax, length.out = buckets + 1)

	Fhat = ecdf(X)(xvals)
	#Fhat.upper = get.berkjones.upper(Fhat, n, alpha)
	Fhat.upper = get.mod.ks.upper(Fhat, n)
	Hhat.upper = get.hajek.upper(Fhat.upper, sampling.ratio)
	Lhat.upper = get.logconc.upper(Hhat.upper)
	
	ret = data.frame(
		xvals=xvals,
		Fhat=Fhat,
		Fhat.upper=Fhat.upper,
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

