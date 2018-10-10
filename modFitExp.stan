// Monoexponential fit

functions{
  // Exp. decay model
  vector phys_mod(real[] x, vector p) {
    int N = size(x);
    vector[N] m;
    for (n in 1:N)
      m[n] = p[1] + p[2] * exp(-x[n]/p[3]);
    return m;
  }

}
data {
  // Calibration dataset
  int<lower=1>       N;
  real               x[N];
  vector[N]          y;
  vector<lower=0>[N] uy;

  // Decay model
  int<lower=1>       Np; // Nb params in decay: should be 3 !!!
}
parameters {
  vector<lower=0>[Np] theta;   // Decay parameters
}
transformed parameters {
  vector[N]     m = phys_mod(x,theta);
  vector[N]     resid;

  for (n in 1:N)
    resid[n] = (y[n] - m[n]); // Rediduals
}
model {

  // Likelihood
  for (n in 1:N)
    0 ~ normal(resid[n], uy[n]);

}
generated quantities{
  // Birge Ratio
  real br = quad_form(inverse(diag_matrix(uy .* uy)),resid)
          / (N-Np);
}

