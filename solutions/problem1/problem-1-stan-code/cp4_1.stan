data {
  int<lower=0> N;
  int<lower=1> K;
  matrix[N, K] x;
  vector[N] y;
}

parameters {
  vector[K] w_mean;
  cov_matrix[K] w_prec;
  vector[K] w;
  real<lower=0> noise_sd;
}

model {
  // Just the identity matrix
  matrix[K,K] idK;
  idK <- diag_matrix(rep_vector(1.0,K));

  w_mean ~ normal(0,2);
  w_prec ~ wishart(K, idK);
  //w ~ multi_normal_prec(w_mean, w_prec);
  w ~ normal(w_mean,1);
  // The second parameter is the scale 0.5, 
  // corresponding to a rate parameter of 2.
  noise_sd ~ inv_gamma(0.5, 0.5); 

  y ~ normal(x * w, noise_sd);
}