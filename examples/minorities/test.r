##############################################
# Load libraries and functions
##############################################

library(foreign)
library(Rcplex)
library(Rglpk)
library(SDMTools)
library(weights)

source("bounds.r")

##############################################
# Read and explore data
##############################################

d = read.dta("all08_11_indigenous.dta")
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

# Aronow and Lee	  
alpha = 0.1
beta = 0.9
out_al = bounds(y, alpha, beta, solver="glpk")
out_al$mu_l
out_al$mu_h

# Do we get the same back from shapebounds when we have no real constraints on the shape?
grid_length = length(y)
target = "normal"
#delta = 0.025
delta = beta / ((n-1)*alpha + beta)
cat("delta =", delta, "\n")
# If delta_shape = 1 then we are within 1 of ANY CDF so there are no effective shape constraints
check_out = bounds(y, alpha, beta, target_shape="normal", delta_shape=1, theta_shape=c(500, 100), grid_length_shape=length(y), solver="glpk")
check_out$mu_l
check_out$mu_h

# Shape constraints
delta = 0.1
out_shape_0.1 = bounds(y, alpha, beta, target_shape="normal", delta_shape=delta, theta_shape=c(500, 100), grid_length_shape=.5*n, solver="glpk")
out_shape_0.1$mu_l
out_shape_0.1$mu_h

# Plot 1
par(mfrow=c(2, 1))
wtd.hist(y, weight=out_al$weights_l, col=rgb(0.1, 0.1, 0.1, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="Weighted histograms and bounds under constrained probabilities", font.main = 1, freq=T, ylim=c(0, 550), xlab="Test scores")
abline(v=out_al$mu_l, col=rgb(0.1, 0.1, 0.1, 0.5), lwd=3)
wtd.hist(y, weight=out_al$weights_h, col=rgb(0.7, 0.7, 0.7, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="", add=T, freq=T, ylim=c(0, 550), xlab="")
abline(v=out_al$mu_h, col=rgb(0.7, 0.7, 0.7, 0.5), lwd=3)
legend("topright", c(expression(mu[L]), expression(mu[H])), lwd=c(3, 3), lty=c(1, 1), col=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), bty="n")
wtd.hist(y, weight=out_shape_0.1$weights_l, col=rgb(0.1, 0.1, 0.1, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="Weighted histograms and bounds under shape constraints", font.main = 1, freq=T, ylim=c(0, 550), xlab="Test scores")
abline(v=out_shape_0.1$mu_l, col=rgb(0.1, 0.1, 0.1, 0.5), lwd=3)
wtd.hist(y, weight=out_shape_0.1$weights_h, col=rgb(0.7, 0.7, 0.7, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="", add=T, freq=T, ylim=c(0, 550), xlab="")
abline(v=out_shape_0.1$mu_h, col=rgb(0.7, 0.7, 0.7, 0.5), lwd=3)
legend("topright", c(expression(mu[L]), expression(mu[H])), lwd=c(3, 3), lty=c(1, 1), col=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), bty="n")
# 
((out_al$mu_h-out_shape_0.1$mu_l)-(out_shape_0.1$mu_h-out_shape_0.1$mu_l))/(out_al$mu_h-out_shape_0.1$mu_l)

# Proportional change bounds	
delta_prop = 2
aux2 = bounds(y, alpha, beta, class=NULL, delta_prop, target_shape=NULL, theta_shape=NULL, delta_shape=NULL, grid_length_shape=NULL, solver="glpk")
aux2$mu_h
aux2$mu_l
aux2$weights_h
aux2$weights_l

a = out_al$mu_h-out_al$mu_l
b = aux2$mu_h-aux2$mu_l
(a-b)/a

# Plot 2
par(mfrow=c(2, 1))
wtd.hist(y, weight=out_al$weights_l, col=rgb(0.1, 0.1, 0.1, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="Weighted histograms and bounds under constrained probabilities", freq=T, ylim=c(0, 550), xlab="Test scores")
abline(v=out_al$mu_l, col=rgb(0.1, 0.1, 0.1, 0.5), lwd=3)
wtd.hist(y, weight=out_al$weights_h, col=rgb(0.7, 0.7, 0.7, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="", add=T, freq=T, ylim=c(0, 550), xlab="")
abline(v=out_al$mu_h, col=rgb(0.7, 0.7, 0.7, 0.5), lwd=3)
legend("topright", c(expression(mu[L]), expression(mu[H])), lwd=c(3, 3), lty=c(1, 1), col=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), bty="n")
wtd.hist(y, weight=aux2$weights_l, col=rgb(0.1, 0.1, 0.1, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="Weighted histograms and worst-case bounds under change constraints", freq=T, ylim=c(0, 550), xlab="Test scores")
abline(v=aux2$mu_l, col=rgb(0.1, 0.1, 0.1, 0.5), lwd=3)
wtd.hist(y, weight=aux2$weights_h, col=rgb(0.7, 0.7, 0.7, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="", add=T, freq=T, ylim=c(0, 550), xlab="")
abline(v=aux2$mu_h, col=rgb(0.7, 0.7, 0.7, 0.5), lwd=3)
legend("topright", c(expression(mu[L]), expression(mu[H])), lwd=c(3, 3), lty=c(1, 1), col=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), bty="n")

##############################################
# Run grid search
############################################## 

grid_length = .25*n
target = "normal"

delta = 0.1

# Make grid
grid_res = 51
gridH = matrix(NA, nrow=grid_res, ncol=grid_res)
gridL = matrix(NA, nrow=grid_res, ncol=grid_res)
mu_list = seq(out_al$mu_l, out_al$mu_h, length.out=grid_res)
mu_list

sigma_bounds = function(y, mu, w_min = 0.1, w_max=0.9, solver="glpk") {
	ystar = (y-mu)^2
	aux = bounds(ystar, w_min, w_max, class=NULL, delta_prop=NULL, target_shape=NULL, theta_shape=NULL, delta_shape=NULL, grid_length_shape=NULL, solver="glpk")
	list(sigma_h = sqrt(aux$mu_h), sigma_l = sqrt(aux$mu_l), weights_h = aux$weights_h, weights_l = aux$weights_l)
}
sd_bounds = sapply(mu_list, function(mu0) {
	lims = sigma_bounds(y, mu0, alpha, beta)
	c(lims$sigma_l, lims$sigma_h)
})
rng = range(sd_bounds)
rng

# Determine which sigma will be tested on grid
sigma_list = seq(rng[1], rng[2], length.out=grid_res)

for (i in 1:length(mu_list)) {
	mu0 = mu_list[i]
	for (j in 1:length(sigma_list)) {		
		sigma0 = sigma_list[j]
        cat("\n", "i = ", i, ", j = ", j, "; mu = ", mu0, ", sigma = ", sigma0, "\n", sep="")

		# Solve if in the range of possible sigmas
		if (sigma0 >= sd_bounds[1, i] && sigma0 <= sd_bounds[2, i]) {

			out_shape = bounds(y, alpha, beta, target_shape=target, delta_shape=delta, theta_shape=c(mu0, sigma0), grid_length_shape=grid_length, solver="glpk")
			
			# feasible
			if (is.null(out_shape$mu_h)) {
				gridH[i, j] = -1
			} else {
				gridH[i, j] = out_shape$mu_h
			}
			# feasible
			if (is.null(out_shape$mu_l)) {
			    gridL[i, j] = -1
			} else {
			    gridL[i, j] = out_shape$mu_l
			}
			cat("mu_l = ", out_shape$mu_l, ", mu_h = ", out_shape$mu_h, "\n", sep="")
				            
		}
	}
}

# Upper bound
library(reshape2)
mlt = melt(gridH, varnames=c("mu", "sigma"), na.rm=TRUE) 
mlt$mu = mu_list[mlt$mu]
mlt$sigma = sigma_list[mlt$sigma]
mlt = subset(mlt, value!=-1)
head(mlt)
plot(value~mu, data=mlt)
abline(0, 1)
# Lower bound
mlt = melt(gridL, varnames=c("mu", "sigma"), na.rm=TRUE) 
mlt$mu = mu_list[mlt$mu]
mlt$sigma = sigma_list[mlt$sigma]
mlt = subset(mlt, value!=-1)
head(mlt)
plot(value~mu, data=mlt)
abline(0, 1)

# Compare ranges
rng = c(min(gridL, na.rm=TRUE), max(gridH, na.rm=TRUE))
df = data.frame(method=c("A & L", "Shape"), low=c(out_al$mu_l, rng[1]), high=c(out_al$mu_h, rng[2]))
df$width = df$high-df$low
df$per = 100*df$width/max(df$width)
print(df)

coor_l = mlt[which.min(abs(mlt$value-rng[1])), ][, 1:2]
coor_h = mlt[which.min(abs(mlt$value-rng[2])), ][, 1:2]

# Map of the possible bounds depending on theta (given A & L constraints)
par(mfrow=c(1, 2))
img = as.numeric(gridL)
brks = c(-1, seq(min(img[img>-1], na.rm=TRUE), max(img, na.rm=TRUE), length.out=12))
summary(img)
image(mu_list, sigma_list, gridL, breaks=brks, col=c(gray.colors(12)), xlab=expression(mu[0]), ylab=expression(sigma[0]))
title(main=expression(paste("Possible values of ", mu[L])))
lines(sd_bounds[1,] ~ mu_list, ylim=range(sd_bounds), type="l")
lines(sd_bounds[2,] ~ mu_list, type="l")
text(coor_l[1], coor_l[2], expression(paste("            min ", mu[L])), col="black")
abline(h = rng, col="red")
summary(img[img>-1])
img = as.numeric(gridH)
brks = c(-1, seq(min(img[img>-1], na.rm=TRUE), max(img,na.rm=TRUE), length.out=12))
summary(img)
image(mu_list, sigma_list, gridH, breaks=brks, col=c(gray.colors(12)), xlab=expression(mu[0]), ylab=expression(sigma[0]))
title(main=expression(paste("Possible values of ", mu[H])))
lines(sd_bounds[1,] ~ mu_list, ylim=range(sd_bounds), type="l")
lines(sd_bounds[2,] ~ mu_list, type="l")
text(coor_h[1]-23, coor_h[2], expression(paste("max ", mu[H])), col="black")

abline(h = rng, col="red")
summary(img[img>-1])
# Key to chart:
# white is out of bounds
# black is infeasible
# color is maximum mu for that given theta

##############################################
# Notes
############################################## 

par(mfrow=c(2, 1))
wtd.hist(y, weight=out_al$weights_l, col=rgb(0.1, 0.1, 0.1, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="Weighted histograms and worst-case bounds under bounded probabilities", freq=T, ylim=c(0, 550), xlab="Test scores")
abline(v=out_al$mu_l, col=rgb(0.1, 0.1, 0.1, 0.5), lwd=3)
wtd.hist(y, weight=out_al$weights_h, col=rgb(0.7, 0.7, 0.7, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="", add=T, freq=T, ylim=c(0, 550), xlab="")
abline(v=out_al$mu_h, col=rgb(0.7, 0.7, 0.7, 0.5), lwd=3)
legend("topright", c(expression(mu[L]), expression(mu[H])), lwd=c(3, 3), lty=c(1, 1), col=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), bty="n")
wtd.hist(y, weight=out_shape_0.2$weights_l, col=rgb(0.1, 0.1, 0.1, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="Weighted histograms and worst-case bounds under shape constraints", freq=T, ylim=c(0, 550), xlab="Test scores")
abline(v=out_shape_0.2$mu_l, col=rgb(0.1, 0.1, 0.1, 0.5), lwd=3)
wtd.hist(y, weight=out_shape_0.2$weights_h, col=rgb(0.7, 0.7, 0.7, 0.5), breaks=c(seq(min(y), 799, by=29), max(y)), main="", add=T, freq=T, ylim=c(0, 550), xlab="")
abline(v=out_shape_0.2$mu_h, col=rgb(0.7, 0.7, 0.7, 0.5), lwd=3)
legend("topright", c(expression(mu[L]), expression(mu[H])), lwd=c(3, 3), lty=c(1, 1), col=c(rgb(0.1, 0.1, 0.1, 0.5), rgb(0.7, 0.7, 0.7, 0.5)), bty="n")

par(mfrow=c(2, 1))
plot(density(y), ylim=c(0, 0.01), main="Weighted densities of y: Aronow and Lee's weights", xlab="y")
lines(density(y, weights=out_al$weights_h/sum(out_al$weights_h)), col="red")
lines(density(y, weights=out_al$weights_l/sum(out_al$weights_l)), col="blue")
abline(v=mean(y), col="black")
abline(v=weighted.mean(y, out_al$weights_l), col="blue")
abline(v=weighted.mean(y, out_al$weights_h), col="red")	
legend("topright", c("True", "High", "Low"), lty=c(1, 1, 1), col=c("black", "red", "blue"), bty="n")
plot(density(y), ylim=c(0, 0.01), main="Weighted densities of y: proposed weights", xlab="y")
lines(density(y, weights=aux2$weights_h/sum(aux2$weights_h)), col="red")
lines(density(y, weights=aux2$weights_l/sum(aux2$weights_l)), col="blue")
abline(v=mean(y), col="black")
abline(v=weighted.mean(y, aux2$weights_l), col="blue")
abline(v=weighted.mean(y, aux2$weights_h), col="red")
legend("topright", c("True", "High", "Low"), lty=c(1, 1, 1), col=c("black", "red", "blue"), bty="n")
		
par(mfrow=c(2, 1))
plot(density(y), ylim=c(0, 0.01), main="Weighted densities of y: Aronow and Lee's weights", xlab="y")
lines(density(y, weights=out_al$weights_h/sum(out_al$weights_h)), col="red")
lines(density(y, weights=out_al$weights_l/sum(out_al$weights_l)), col="blue")
abline(v=mean(y), col="black")
abline(v=weighted.mean(y, out_al$weights_l), col="blue")
abline(v=weighted.mean(y, out_al$weights_h), col="red")	
legend("topright", c("True", "High", "Low"), lty=c(1, 1, 1), col=c("black", "red", "blue"), bty="n")
plot(density(y), ylim=c(0, 0.01), main="Weighted densities of y: proposed weights", xlab="y")
lines(density(y, weights=out_shape_0.1$weights_h/sum(out_shape_0.1$weights_h)), col="red")
lines(density(y, weights=out_shape_0.1$weights_l/sum(out_shape_0.1$weights_l)), col="blue")
abline(v=mean(y), col="black")
abline(v=weighted.mean(y, out_shape_0.1$weights_l), col="blue")
abline(v=weighted.mean(y, out_shape_0.1$weights_h), col="red")
legend("topright", c("True", "High", "Low"), lty=c(1, 1, 1), col=c("black", "red", "blue"), bty="n")
	
