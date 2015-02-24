# Generate a synthetic QMR-DT network
#
# QMR-DT 534 disease 4040 findings 40740 dependency arcs
#
# 4040 = 534 * 7.57
# Used for scaling the number of findings based on the number of diseases
total.findings.per.disease = 7.57

# 76.29 findings caused by each disease on average
caused.findings.per.disease = 76.29

# 
# diseases are marginally independent
# causal independence (noise-OR)
# P(finding-positive | diesease-positive)
# P(finding=positive | D) = 1 - prod_{i \in parents(finding)} [1 - P(finding=positive|disease_i = prsent)
# largest leak probability in QMR-DT was 0.153. 
# smallest was 5.8e-8
# I decided not to include leak probabilities
#
# disease incidence. A multinomial distribution with exponentially decaying probabilities
# mean # of inbound arcs is 10.08
# 
# weight on each disease link beta(shape1=1, shape2=3)
# 

num.diseases <- 20
num.findings <- floor(num.diseases*total.findings.per.disease)
edges <- matrix(0, nrow = num.diseases, ncol=num.findings)
for (idisease in 1:num.diseases) {
  # how many findings should this disease have?
  # we will model this as a binomial distribution with prob = 0.25 plus 1
  # and size (number of trials) = (76.29-1)/0.25
  size <- round((caused.findings.per.disease - 1) / 0.25)
  # least one finding per disease.  
  n.findings <- rbinom(1, size, 0.25) + 1
  # now choose the findings uniformly at random without replacement
  findings <- sample(1:num.findings, n.findings)
  # link weights
  weights <- rbeta(n.findings, 1, 3)
  edges[idisease, findings] <- weights
}
# how many findings are unused
unused <- c()
for (ifinding in 1:num.findings) {
  if (sum(edges[,ifinding]>0)==0) {unused <- c(unused, ifinding)}
}

# prior probabilities of disease
# distributed iid beta(1, 20)
disease.priors <- rbeta(num.diseases, 1, 20)

# now generate a patient
gen.diseases <- function () {
   true.diseases <- rep(0, num.diseases)
   for (idisease in 1:num.diseases) {
     true.diseases[idisease] <- rbinom(1, 1, disease.priors[idisease]);
   }
   return(true.diseases)
}

# repeat until we get a sick patient
gen.case <- function () {
   diseases <- c()
   tries <- 1
   while (sum(diseases <- gen.diseases())==0) {tries <- tries + 1}
   return(diseases)
}

noisy.or <- function (probs) {
  result <- 1.0
  for (p in probs) result <- result * (1 - p)
  return(1 - result)
}

# now generate the findings for the diseases of this patient
gen.findings <- function(diseases) {
  findings <- rep(0, num.findings)
  for (i in 1:num.findings) {
     findings[i] <- rbinom(1, 1, noisy.or(edges[diseases == 1, i]))
  }
  return(findings)
}

# Treatment costs: exponentially distributed
# log(cost(d)) ~ lognormal(0,1)
treatment.costs <- exp(rlnorm(num.diseases,0, 1))

# Observation costs: gamma distributed
# observation.cost(f) ~ 1 + gamma(4, 1)
observation.costs <- 1 + rgamma(num.findings, 4, 1)

# Probability of measuring finding f propto exp(-observation.ccsts[f])
z <- sum(exp(-observation.costs))
observation.probability <- exp(-observation.costs) / z


# generate the cases
num.cases <- 4
diseases <- matrix(0, num.diseases, num.cases)
findings <- matrix(0, num.findings, num.cases)
partial.findings <- matrix(NA, num.findings, num.cases)
for (icase in 1:num.cases) {
  # generate the true diseases of this patient
  diseases[,icase] <- gen.case()
  # generate the true findings
  findings[,icase] <- gen.findings(diseases[,icase])
  # Sample 20% of the findings to get the observed findings
  s <- sample(1:num.findings, size = round(0.2 * num.findings), prob = observation.probability)
  partial.findings[s, icase] <- findings[s, icase]
}

# Save the true network parameters
# note: In the release, I have manually edited these files to improve
# readability and to name the rows and columns.   
write.csv(disease.priors, "problem-2-disease-priors.csv")
write.csv(edges, "problem-2-edges.csv")
write.csv(observation.costs, "problem-2-observation-costs.csv")
write.csv(treatment.costs, "problem-2-treatment-costs.csv")
write.csv(diseases, "problem-2-cases-ground-truth.csv")
write.csv(findings, "problem-2-cases-findings.csv")
write.csv(partial.findings, "problem-2-cases-partial-findings.csv")
