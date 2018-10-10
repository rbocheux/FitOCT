rm(list = ls()); gc() # Clean environment

# Libraries ####
libs =c('parallel','rstan','inlmisc')

for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

# Options ####
options(mc.cores = parallel::detectCores(), width = 90)
rstan_options(auto_write = TRUE)
set.seed(1234) # Initialise la graine du RNG

# Color schemes ####
cols    = inlmisc::GetTolColors(8)
col_tr  = inlmisc::GetTolColors(8,alpha=0.1)
col_tr2 = inlmisc::GetTolColors(8,alpha=0.4) # For legends

# Get stan models ####
modHet   = stan_model(file = 'modHetero.stan')
modExp   = stan_model(file = 'modFitExp.stan')
modExpGP = stan_model(file = 'modFitExpGP.stan')

# Interface functions to stan models ####
estimateNoise <- function(x, y, tag, df = 15,
                          model = modHet, sample = FALSE,
                          nb_chains = 4, nb_warmup = 500,
                          nb_iter = nb_warmup + 500) {

  # Smoothing
  ySpl = smooth.spline(x,y,df=df)$y
  resSpl = y-ySpl

  # Fit smoothing residuals by exponential variance model

  stanData = list(N =length(x), x=x, y=resSpl)
  init = list(theta = c(max(resSpl),mean(x)))

  if (sample) {
    # Generate MCMC sample
    parOpt = c('theta')
    pars   = c(parOpt)
    fit = sampling(model,
                   data = stanData,
                   pars = pars,
                   init = function() {init},
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)

    pdf(file = paste0('Results/controleResiduals_',tag,'.pdf'))
    print(rstan::traceplot(fit, inc_warmup=TRUE, pars = parOpt))
    pairs(fit, pars = c(parOpt,'lp__'), gap=0)
    dev.off()

    sink(file = paste0('Results/controleResiduals_',tag,'.txt'))
    print(fit,pars=parOpt)
    sink()
    # Estimate data uncertainty
    S = as.matrix(fit)
    theta = colMeans(S)[1:2]

  } else {
    # Optimize
    fit = optimizing(model,
                     data = stanData,
                     init = init,
                     as_vector = FALSE,
                     verbose   = TRUE)

    sink(file = paste0('Results/controleResiduals_',tag,'.txt'))
    cat('theta :',fit$par$theta)
    sink()

    # Estimate data uncertainty
    theta = fit$par$theta

  }
  sig = theta[1]*exp(-x/theta[2])

  return(list(uy = sig, ySmooth = ySpl))
}

fitMonoExp <- function(x, y, uy, tag,
                       model = modExp, sample = FALSE,
                       nb_chains = 4, nb_warmup = 500,
                       nb_iter = nb_warmup + 500) {

  stanData = list(
    N =length(x),
    x=x, y=y, uy=uy,
    Np=3
  )

  init = function() {
    list(theta  = c(min(y),max(y)-min(y),mean(x)))
  }

  if(sample) {
    parOpt = c('theta')
    pars   = c(parOpt,'resid','m','br')

    fit = sampling(model,
                   data = stanData,
                   pars = pars,
                   init = init,
                   control = list(adapt_delta=0.99,
                                  max_treedepth=12),
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)


    sink(file = paste0('Results/controleExp_',tag,'.txt'))
    print(fit,pars=c(parOpt,'br'))
    sink()

    parP = c('theta')
    S = as.matrix(fit,pars=c(parP,'lp__'))

    pdf(file = paste0('Results/controleExp_',tag,'.pdf'))
    print(rstan::traceplot(fit, inc_warmup=TRUE, pars = parOpt))
    pairs(fit, pars = c(parP,'lp__'), gap=0)
    par(mfrow=c(1,1),pty='m')
    rgumlib::SAPlot(S)
    dev.off()

    # Estimate decay params
    theta   = extract(fit,'theta')[[1]]
    theta0  = colMeans(theta)

  } else {

    fit = optimizing(model,
                     data = stanData,
                     init = init,
                     as_vector = FALSE,
                     verbose   = FALSE)

    sink(file = paste0('Results/controleExp_',tag,'.txt'))
    cat('theta  :',fit$par$theta,'\n')
    cat('br     :',fit$par$br,'\n')
    sink()

    # Estimate decay params
    theta0  = fit$par$theta

  }
  return(list(fit = fit, best.theta = theta0))
}

fitExpGP <- function(x, y, uy, tag,
                     Nn = 10, # Nb of control points
                     theta0 = NULL, ru_theta    = 0.1, # Theta prior
                     model = modExpGP,
                     sample = FALSE, # Flag MCMC / optimize
                     iter = 10000,   # Max. iterations of optimizer
                     nb_chains = 4, nb_warmup = 500,
                     nb_iter = nb_warmup + 500) {

  # Grid of GP control points
  dx  = diff(range(x))/(Nn+1)
  xGP = seq(min(x)+dx/2,max(x)-dx/2,length.out = Nn)

  # Initial monoexp params
  if(is.null(theta0))
    theta0 = c(min(y),max(y)-min(y),mean(x))

  stanData = list(
    N =length(x),
    x=x, y=y, uy=uy,
    Np = 3,
    Nn = Nn,
    xGP = xGP,
    alpha_scale = 0.1,
    rho_scale   = 0.1,
    theta0      = theta0,
    ru_theta    = ru_theta
  )

  init = list(
    theta  = theta0,
    yGP    = 0.01*rnorm(Nn),
    lambda = 20.0,
    sigma  = 1.0
  )


  if(sample) {
    parOpt = c('theta','yGP','lambda','sigma')
    pars   = c(parOpt,'resid','br','m','dL')

    fit = sampling(model,
                   data = stanData,
                   pars = pars,
                   init = function() {init},
                   control = list(adapt_delta=0.99,max_treedepth=12),
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)


    sink(file = paste0('Results/controleGP_',tag,'.txt'))
    print(fit,pars=c(parOpt,'br'))
    sink()

    parP = c('theta','lambda','sigma')
    S = as.matrix(fit,pars=c(parP,'lp__'))

    pdf(file = paste0('Results/controleGP_',tag,'.pdf'))
    print(rstan::traceplot(fit, inc_warmup=TRUE, pars = parOpt))
    pairs(fit, pars = c(parP,'lp__'), gap=0)
    pairs(fit, pars = c('yGP','lp__'), gap=0)
    par(mfrow=c(1,1),pty='m')
    rgumlib::SAPlot(S)
    dev.off()

  } else {

    fit = optimizing(model,
                     data = stanData,
                     init = init,
                     as_vector = FALSE,
                     verbose   = FALSE,
                     algorithm = 'LBFGS',
                     hessian = TRUE,
                     draws = 500,
                     iter = iter,
                     refresh = 500)
    if(is.null(fit$theta_tilde))
      # No Hessian-based sampling done => try to refine solution
      fit = optimizing(model,
                       data = stanData,
                       init = fit$par,
                       as_vector = FALSE,
                       verbose   = FALSE,
                       algorithm = 'LBFGS',
                       hessian = TRUE,
                       draws = 500,
                       iter = iter,
                       refresh = 500)

    sink(file = paste0('Results/controleGP_',tag,'.txt'))
    cat('theta  :',fit$par$theta,'\n')
    cat('yGP    :',fit$par$yGP,'\n')
    cat('lambda :',fit$par$lambda,'\n')
    cat('sigma  :',fit$par$sigma,'\n')
    cat('br     :',fit$par$br,'\n')
    sink()

  }

  return(list(fit = fit, sample = sample, xGP = xGP))
}

# RUN ####

setwd("~/Bureau/Collabs/Romain")
dataDirs = c("DataWl","Data1","DataSynth")
for (dataDir in dataDirs) {

  dataSets = list.dirs(path=dataDir,full.names = FALSE)[-1]

  for(dataSet in dataSets) {

    tag = paste0('GP_',dataDir,'_',dataSet)
    cat(tag,'------------------------','\n')

    # Get Data ####
    D = read.csv(paste0(dataDir,'/',dataSet,'/Courbe.csv'))
    x=D[,1]; y=D[,2]
    sel = x > 20 & x<=500 # Exclude aberrant points
    x = x[sel]; y = y[sel]

    # Estimate data uncertainty
    fits = estimateNoise(x, y, tag)
    uy   = fits$uy      # Used by next stages
    ySpl = fits$ySmooth # Used by plotMonoExp

    # Inference of exponential decay parameters
    fitm   = fitMonoExp(x, y, uy, tag)
    theta0 = fitm$best.theta   # Used by next stage
    source ("./plotMonoExp.R")

    # Inference of modulated decay parameters
    fitGP = fitExpGP(x, y, uy, tag,
                     sample = TRUE,
                     theta0 = theta0)
    source ("./plotExpGP.R")

  }
}
