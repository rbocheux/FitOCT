rm(list = ls()); gc() # Clean environment

# Libraries ####
libs =c('parallel','rstan','inlmisc','purrr')

for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

# Options ####
options(mc.cores = parallel::detectCores(), width = 90)
rstan_options(auto_write = TRUE)
set.seed(1234) # Initialise la graine du RNG


# RUN ####
model ="
data {
  real<lower=0> lambda_scale;
}
parameters {
  real<lower=0> lambda;
}
model {
  lambda ~ exponential(1./lambda_scale);
}

"

modHP = (purrr::quietly(stan_model)(model_code = model))$result

stanData = list(
  lambda_scale= 10
)

# Bayesian Lasso

pars = c('lambda')

fit = sampling(object  = modHP,
               data    = stanData,
               pars    = pars,
               control = list(adapt_delta=0.99,max_treedepth=12),
               warmup  = 500, iter = 50000, chains = 4,
               verbose=FALSE)

print(fit)
stan_hist(fit,pars=pars)
stan_ac(fit,pars=pars)
stan_trace(fit,pars=pars)
