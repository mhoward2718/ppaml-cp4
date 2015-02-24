# Problem 1
# Tom Dietterich, June 2014
#
# The coefficients (including the intercept) are generated from a normal-inverse-wishart
# with mean drawn from norm(0,2) and the identity matrix for the matrix parameter
#
# The data are then generated with additional zero-mean gaussian observation noise
# the standard deviation is drawn from an inverse gamma with shape=2 and rate=1
#
library(MASS)
d <- 4                # dimensions of the input
dplus <- d+1           # for the bias term
n <- 500                # sample size
# prior mean
prior.mu <- rnorm(dplus, 0, 2)      # mean = 0, standard deviation 2
# create an identity matrix
prior.S <- matrix(0, dplus, dplus)
for (i in 1:dplus) prior.S[i,i]=1
# Draw a wishart 
prior.Sigma <- rWishart(1, dplus, prior.S)
# force it to be a matrix
dim(prior.Sigma) <- c(dplus,dplus)
# compute precision matrix
prior.SigmaInv <- solve(prior.Sigma)
# draw the weight vector
w <- mvrnorm(1, prior.mu, prior.SigmaInv)

# noise variance (precision) prior; 
tau <- rgamma(1, 0.5, 2) # shape=0.5, rate=2

## the features will be sampled from [-1,1]
x <- matrix(data = runif(n*(d+1), -1, 1), nrow = n, ncol = dplus)
x[,1] <- -1                # first column corresponds to the intercept term
eps <- rnorm(n, 0, 1/tau)  # additive noise
y <- x %*% w + eps         # observations

# Compare the MLE. Hard coded for d=4
df <- data.frame(x,y)
# omitting the intercept term since it will be supplied by lm
model1 <- lm(y ~ X2 + X3 + X4 + X5, df) 

# write out the data
write.csv(df, "problem-1-data.csv")

