// Model standard deviation of heteroscedastic residuals
// by exponential: uy(x) ~ Normal(mean = 0, sd = a*exp(-x/b) )

data {
  int<lower=1> N; // size of vectors
  vector[N]    x; // positions
  vector[N]    y; // values
}
parameters {
  vector<lower=0, upper=1e4>[2] theta; // upper limit: homoscedasticity
}
transformed parameters {
  vector[N] sig;
  for (n in 1:N)
    sig[n] = theta[1]*exp(-x[n]/theta[2]);
}
model {
  for (n in 1:N)
    y[n] ~ normal(0., sig[n]);
}
