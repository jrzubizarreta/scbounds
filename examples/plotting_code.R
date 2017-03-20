library(zoo) # for plotting only
library(splines) #for plotting only...


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
