# Problem 3: Discrete Time Hidden Markov Model
# 
d <- 5  # number of states
steps <- 20 # number of steps

# initial distribution
p0 <- rep(0,d)
p0[1] <- 1.0   # initial state is 1 w.p. 1.0

# transitions
stay <- 1/3  # probability of staying
transitions <- matrix(0.0, d, d)
# each state will stay w.p. stay and move to one of the next two states w.p. (1-stay)/2
for (state in 1:d) {
  transitions[state, state] <- stay
  transitions[state, (state %% d) + 1] <- (1 - stay)/2
  transitions[state, ((state+1) %% d) +1] <- (1 - stay)/2
}
# emissions
# each state emits its true identity with this probability
faithful <- 0.6
unfaithful <- (1-faithful)/(d-1)
# otherwise each state emits an output chosen uniformly at random
observations <- matrix(unfaithful, d, d)
for (state in 1:d) {
  observations[state, state] <- faithful
}

# generate the observation sequence
true.state <- matrix(0, nrow=d, ncol=steps)
outputs <- rep(0, steps)
true.state[,1] <- rmultinom(1, rep(1, d), p0)
for (i in 2:steps) {
  true.state[,i] <- rmultinom(1, rep(1, d), true.state[,i-1] %*% transitions)
}
for (i in 1:steps) {
  outputs[i] <- which.max( rmultinom(1, rep(1, d), true.state[,i] %*% observations) )
}


# save to file
# I have manually edited these files to add meaningful row and column names
write.csv(outputs, "problem-3-outputs.csv")
write.csv(true.state, "problem-3-true-state.csv")

