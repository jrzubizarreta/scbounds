rm(list = ls())
library(foreign)
source("bounds.R")
source("bounds_orig.R")

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

sampling.ratio = 9

out_al = bounds(y, 1, sampling.ratio, solver="glpk")
c("low"=out_al$mu_l, "high"=out_al$mu_h)
bounds.logconc(y, sampling.ratio = 9)

ret = bounds.upper.internal(y, sampling.ratio = 9)

pdf("minorities_CDF.pdf")
plot.ret(ret)
dev.off()

pdf("minorities_DEN.pdf")
plot.ret2(ret)
dev.off()
