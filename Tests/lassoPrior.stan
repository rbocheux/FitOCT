data {
  // GP
  int<lower=0>       Nn;          // Nb. of control points (CP) for GP
  real<lower=0>      lambda_scale;// Scale of CP S.D.
}
parameters {
  vector[Nn]          yGP;     // GP controle values
}
model {
  for (n in 1:Nn)
    target += - lambda_scale * fabs(yGP[n]);
  target += -lambda_scale * dot_self(yGP);
}
