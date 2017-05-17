rm(list = ls())

library(foreign)

detach("package:bounds", unload=TRUE)
library(bounds)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

d = read.dta("minorities/all08_11_indigenous.dta")
nrow(d)
head(d)
table(d$ethnic_student_08)

d2 = subset(d, ethnic_student_08=="Aymara")	# Mapuche
nrow(d2)

d3 = subset(d2, !is.na(mate_students_11))
d3 = subset(d2, mate_students_11!=0)
d3 = d3[order(d3$mate_students_11), ]
head(d3)

n = nrow(d3)
y = d3$mate_students_11

hist(y, breaks=21, col="grey")

cat("Number units:", n, "\n")

##############################################
# Test code
##############################################

out.al.9 = bounds(y, sampling.ratio = 9, constraint="none")
out.lc.9 = bounds(y, sampling.ratio = 9, constraint="logconcave")
out.sy.9 = bounds(y, sampling.ratio = 9, constraint="symmetric")

# out.al.5 = bounds(y, sampling.ratio = 5, constraint="none")
# out.lc.5 = bounds(y, sampling.ratio = 5, constraint="logconcave")
# out.sy.5 = bounds(y, sampling.ratio = 5, constraint="symmetric")
# 
# out.al.2 = bounds(y, sampling.ratio = 2, constraint="none")
# out.lc.2 = bounds(y, sampling.ratio = 2, constraint="logconcave")
# out.sy.2 = bounds(y, sampling.ratio = 2, constraint="symmetric")

xvals.all = out.al.9$upper.bound.internal$raw$xvals
good.idx = which(xvals.all >= min(y) & xvals.all <= max(y))

K = length(good.idx)
xvals = xvals.all[good.idx]

Xhat = out.al.9$upper.bound.internal$raw$Xhat[good.idx]
Fhat.al = out.al.9$upper.bound.internal$raw$Xhat.weighted[good.idx]
Fhat.lc = out.lc.9$upper.bound.internal$raw$Xhat.weighted[good.idx]
Fhat.sy = out.sy.9$upper.bound.internal$raw$Xhat.weighted[good.idx]

xhat = (Xhat - c(0, Xhat[-K])) / (xvals[2] - xvals[1])
fhat.al = (Fhat.al - c(0, Fhat.al[-K])) / (xvals[2] - xvals[1])
fhat.lc = (Fhat.lc - c(0, Fhat.lc[-K])) / (xvals[2] - xvals[1])
fhat.sy = (Fhat.sy - c(0, Fhat.sy[-K])) / (xvals[2] - xvals[1])

df = 9
xhat.smooth = unpad(predict(glm(zeropad(xhat) ~ ns(1:(length(xhat) + 200), df = df), family=quasi(link = "log", variance = "mu")), type = "response"))

#fhat.lc.ratio = smooth.spline(na.approx(fhat.lc / fhat), df = df)$y

my.na.approx = function(xxx) { 
  aaa = approx((1:length(xxx))[!is.na(xxx)], xxx[!is.na(xxx)], 1:length(xxx))$y
  max.idx = max(which(!is.na(aaa)))
  aaa[1:length(xxx) > max.idx] = aaa[max.idx]
  min.idx = min(which(!is.na(aaa)))
  aaa[1:length(xxx) < min.idx] = aaa[min.idx]
  aaa
}

fhat.al.ratio = my.na.approx(fhat.al / xhat)
fhat.lc.ratio = smooth.spline(my.na.approx(fhat.lc / xhat), df = df)$y
fhat.sy.ratio = smooth.spline(my.na.approx(fhat.sy / xhat), df = df)$y

rng = range(xhat.smooth, xhat.smooth * fhat.al.ratio, xhat.smooth * fhat.lc.ratio, xhat.smooth * fhat.sy.ratio, na.rm = TRUE)

pdf("minorities_DEN.pdf")
plot(xvals, xhat.smooth, lwd = 2, ylim = rng, xlab="y", ylab="f(y)", type = "l")
lines(xvals, xhat.smooth * fhat.lc.ratio, col = 2, lwd = 2)
lines(xvals, xhat.smooth * fhat.al.ratio, col = 4, lwd = 2)
lines(xvals, xhat.smooth * fhat.sy.ratio, col = 5, lwd = 2)
legend("topleft", c("Uncorrected", "Symmetry const.", "Log-concave constr.", "Plain Aronow-Lee"), lwd = 2, col = c(1, 5, 2, 4))
dev.off()