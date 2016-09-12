n = 10000

x = (1:n)/(n + 1)

res = replicate(10000, {
U = sort(runif(n))

DAD = max(abs(x - U) / sqrt(x * (1 - x))) * sqrt(n)
Dtail = max(abs(x - U) / (x * (1 - x)))
Dtot = max(DAD, Dtail)
c(DAD=DAD, Dtail=Dtail, Dtot=Dtot)
})