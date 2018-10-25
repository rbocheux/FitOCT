data {
  // GP
  int<lower=0>       Nn;          // Nb. of control points (CP) for GP
  real<lower=0>      lambda_rate;// Scale of CP S.D.
}
parameters {
  real<lower=0>       lambda;  // Lasso rate
  vector<lower=0>[Nn] u;       // Intermediate variable
  vector[Nn]          yGP;     // GP controle values
}
model {
  // Hierarchical Prior for Lasso (adapted from Mallick2014 and ARD)
  lambda ~ gamma(2,lambda_rate);
  for (n in 1:Nn)
    u[n] ~ gamma(2,lambda);
  for (n in 1:Nn)
    yGP[n] ~ normal(0, u[n]);
}
