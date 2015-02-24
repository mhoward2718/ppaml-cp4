setwd("~/galois/ppaml/challenge-problems/cp4/")

xy <- read.csv("problem-1-data.csv")
x <- xy[,2:6]
y <- xy[,7]


theData <- list(
    x=x
  , y=y
  , K=ncol(x)
  , N=nrow(x)
)

library(rstan)
sm <- stan_model("cp4_1.stan")

# For MAP estimate
map <- optimizing(sm, data=theData)

# To sample from posterior
fit <- sampling(sm, data=theData)

samp <- extract(fit)
write.table(samp$w,"cp4_1_w_posterior.csv",row.names=F, col.names=F, sep=',')