data {
  // GP
  int<lower=0>           Nn;          // Nb. of control points (CP) for GP
  real<lower=0>          lambda_rate;// Weight of EN penalty
  real<lower=0, upper=1> beta;        // Mixture coed of Lasso in EN
}
transformed data{
  real<lower=0> lambda = 1./lambda_rate;
}
parameters {
  vector[Nn]          yGP;     // GP controle values
}
model {
  // Ridge Regression
  target += -0.5*lambda*(1-beta) * dot_self(yGP);
  //Lasso
  target += -lambda* beta * sum(fabs(yGP));
}
