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
# quietly avoids heavy useless compilation warnings
modHet        = (purrr::quietly(stan_model)(file = 'modHetero.stan'))$result
modExp        = (purrr::quietly(stan_model)(file = 'modFitExp.stan'))$result
modExpGP      = (purrr::quietly(stan_model)(file = 'modFitExpGP.stan'))$result
modExpGPLasso = (purrr::quietly(stan_model)(file = 'modFitExpGPLasso.stan'))$result

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
    ndf = stanData$N - stanData$Np
    CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf))/ndf
    cat('CI95 for br: ',CI95)
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
    thetaCor= cor(theta)

  } else {

    fit = optimizing(model,
                     data = stanData,
                     init = init,
                     as_vector = FALSE,
                     verbose   = FALSE,
                     hessian   = TRUE )

    sink(file = paste0('Results/controleExp_',tag,'.txt'))
    cat('theta  :',fit$par$theta,'\n')
    cat('br     :',fit$par$br,'\n')
    ndf = stanData$N - stanData$Np
    CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf))/ndf
    cat('CI95 for br: ',CI95)
    sink()

    # Estimate decay params
    theta0  = fit$par$theta
    thetaCor= cov2cor(solve(-fit$hessian))

  }
  return(list(fit = fit, best.theta = theta0, cor.theta = thetaCor))
}

fitExpGP <- function(x, y, uy, tag,
                     Nn = 10, # Nb of control points
                     theta0 = NULL, cor_theta = NULL, ru_theta    = 0.1, # Theta prior
                     lambda_scale = 0.1, # Scale of ctrl points prior
                     model = modExpGP,
                     sample = FALSE, # Flag MCMC / optimize
                     iter = 10000,   # Max. iterations of optimizer
                     prior_PD = 0,   # Flag to sample from prior only
                     nb_chains = 4, nb_warmup = 500,
                     nb_iter = nb_warmup + 500) {

  # Adjust tag
  tag1 = tag
  if(prior_PD != 0)
    tag1=paste0('priPD_',tag)

  # Grid of GP control points
  dx  = diff(range(x))/(Nn+1)
  xGP = seq(min(x)+dx/2,max(x)-dx/2,length.out = Nn)

  # Initial monoexp params
  if(is.null(theta0))
    theta0 = c(min(y),max(y)-min(y),mean(x))

  if(is.null(cor_theta))
    cor_theta = diag(1,length(theta0))

  stanData = list(
    N =length(x),
    x=x, y=y, uy=uy,
    Np = 3,
    Nn = Nn,
    xGP = xGP,
    alpha_scale = 0.1,
    rho_scale   = 0.1,
    theta0      = theta0,
    cor_theta   = cor_theta,
    ru_theta    = ru_theta,
    prior_PD    = prior_PD,
    lambda_scale= lambda_scale
  )

  init = list(
    theta  = theta0,
    yGP    = 0.01*rnorm(Nn),
    lambda = 10.0,
    sigma  = 1.0
  )


  if(sample) {
    parOpt = c('theta','yGP','lambda','sigma')
    pars   = parOpt
    if(prior_PD == 0)
      pars   = c(pars,'resid','br','m','dL')

    fit = sampling(model,
                   data = stanData,
                   pars = pars,
                   init = function() {init},
                   control = list(adapt_delta=0.99,max_treedepth=12),
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)


    sink(file = paste0('Results/controleGP_',tag1,'.txt'))
    pars   = parOpt
    if(prior_PD == 0)
      pars=c(parOpt,'br')
    print(fit,pars=pars)
    if(prior_PD == 0) {
      ndf = stanData$N - stanData$Np - stanData$Nn - 2
      CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf))/ndf
      cat('CI95 for br: ',CI95)
    }
    sink()

    parP = c('theta','lambda','sigma')
    S = as.matrix(fit,pars=c(parP,'lp__'))

    pdf(file = paste0('Results/controleGP_',tag1,'.pdf'))
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

    sink(file = paste0('Results/controleGP_',tag1,'.txt'))
    cat('theta  :',fit$par$theta,'\n')
    cat('yGP    :',fit$par$yGP,'\n')
    cat('lambda :',fit$par$lambda,'\n')
    cat('sigma  :',fit$par$sigma,'\n')
    if(prior_PD == 0) {
      cat('br     :',fit$par$br,'\n')
      ndf = stanData$N - stanData$Np - stanData$Nn - 2
      CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf))/ndf
      cat('CI95 for br: ',CI95)
    }
    sink()

  }

  return(list(fit = fit, sample = sample, xGP = xGP, prior_PD = prior_PD))
}

fitExpGPLasso <- function(x, y, uy, tag,
                     Nn = 10, # Nb of control points
                     theta0 = NULL, cor_theta = NULL, ru_theta    = 0.1, # Theta prior
                     lambda_scale = 0.1, # Scale of ctrl points prior
                     model = modExpGPLasso,
                     sample = FALSE, # Flag MCMC / optimize
                     iter = 10000,   # Max. iterations of optimizer
                     prior_PD = 0,   # Flag to sample from prior only
                     nb_chains = 4, nb_warmup = 500,
                     nb_iter = nb_warmup + 500) {

  # Adjust tag
  tag1 = tag
  if(prior_PD != 0)
    tag1=paste0('priPD_',tag1)

  # Grid of GP control points
  dx  = diff(range(x))/(Nn+1)
  xGP = seq(min(x)+dx/2,max(x)-dx/2,length.out = Nn)

  # Initial monoexp params
  if(is.null(theta0))
    theta0 = c(min(y),max(y)-min(y),mean(x))

  if(is.null(cor_theta))
    cor_theta = diag(1,length(theta0))

  stanData = list(
    N =length(x),
    x=x, y=y, uy=uy,
    Np = 3,
    Nn = Nn,
    xGP = xGP,
    alpha_scale = 0.1,
    rho_scale   = 0.1,
    theta0      = theta0,
    cor_theta   = cor_theta,
    ru_theta    = ru_theta,
    prior_PD    = prior_PD,
    lambda_scale= lambda_scale
  )

  init = list(
    theta  = theta0,
    yGP    = 0.01*rnorm(Nn),
    sigma  = 1.0
  )


  if(sample) {
    parOpt = c('theta','yGP','sigma')
    pars   = parOpt
    if(prior_PD == 0)
      pars   = c(pars,'resid','br','m','dL')

    fit = sampling(model,
                   data = stanData,
                   pars = pars,
                   init = function() {init},
                   control = list(adapt_delta=0.99,max_treedepth=12),
                   iter = nb_iter, chains = nb_chains,
                   warmup = nb_warmup, verbose=FALSE)


    sink(file = paste0('Results/controleGP_',tag1,'.txt'))
    pars   = parOpt
    if(prior_PD == 0)
      pars=c(parOpt,'br')
    print(fit,pars=pars)
    if(prior_PD == 0) {
      ndf = stanData$N - stanData$Np - stanData$Nn - 2
      CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf))/ndf
      cat('CI95 for br: ',CI95)
    }
    sink()

    parP = c('theta','sigma')
    S = as.matrix(fit,pars=c(parP,'lp__'))

    pdf(file = paste0('Results/controleGP_',tag1,'.pdf'))
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

    sink(file = paste0('Results/controleGP_',tag1,'.txt'))
    cat('theta  :',fit$par$theta,'\n')
    cat('yGP    :',fit$par$yGP,'\n')
    cat('sigma  :',fit$par$sigma,'\n')
    if(prior_PD == 0) {
      cat('br     :',fit$par$br,'\n')
      ndf = stanData$N - stanData$Np - stanData$Nn - 2
      CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf))/ndf
      cat('CI95 for br: ',CI95)
    }
    sink()

  }

  return(list(fit = fit, sample = sample, xGP = xGP, prior_PD = prior_PD))
}
# Misc. functions ####
plotPriPos <- function(pri,pos,tag,xlim=range(c(pri,pos))) {
  d = density(pri)
  d$y = d$y/max(d$y)
  plot(d, type = 'l', col = cols[4],
       main = tag,
       xlab = '', xlim = xlim,
       ylab = 'Norm. density', ylim = c(0,1.1))
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[4],border=NA)
  d = density(pos)
  d$y = d$y/max(d$y)
  lines(d$x,d$y,col=cols[6])
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[6],border=NA)
}

# RUN ####
dataDirs = c("DataWl","Data1","DataSynth")[3]
for (dataDir in dataDirs) {

  dataSets = list.dirs(path=dataDir,full.names = FALSE)[-1]

  for(dataSet in dataSets) {

    tag = paste0('GP_',dataDir,'_',dataSet)
    cat(tag,'------------------------','\n')

    # Get Data
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
    cor.theta = fitm$cor.theta
    source ("./plotMonoExp.R")

    # Inference of modulated decay parameters
    lasso = TRUE # Use lasso prior
    if(lasso)
      tag=paste0('lasso_',tag)
    # - Prior (Predictive) Distribution
    fitGP_pri = fitExpGPLasso(x, y, uy, tag,
                              sample = TRUE,
                              theta0 = theta0,
                              cor_theta = cor.theta,
                              ru_theta = 0.05,
                              prior_PD = 1,
                              lambda_scale=0.05)
    fitOut = fitGP_pri; source ("./plotExpGP.R")

    # - Posterior Distribution
    fitGP = fitExpGPLasso(x, y, uy, tag,
                          sample = TRUE,
                          theta0 = theta0,
                          cor_theta = cor.theta,
                          ru_theta = 0.05,
                          prior_PD = 0,
                          lambda_scale=0.05)
    fitOut = fitGP    ; source ("./plotExpGP.R")

    # - Compare prior/posterior marginal pdfs
    theta_pri   = extract(fitGP_pri$fit,'theta')[[1]]
    yGP_pri     = extract(fitGP_pri$fit,'yGP')[[1]]
    if(!lasso) lambda_pri  = extract(fitGP_pri$fit,'lambda')[[1]]
    sigma_pri   = extract(fitGP_pri$fit,'sigma')[[1]]
    theta_pos   = extract(fitGP$fit,'theta')[[1]]
    yGP_pos     = extract(fitGP$fit,'yGP')[[1]]
    if(!lasso) lambda_pos  = extract(fitGP$fit,'lambda')[[1]]
    sigma_pos   = extract(fitGP$fit,'sigma')[[1]]
    png(filename = paste0('Results/pripos_',tag,'.png'),
        width=2000, height=2000)
    par(mfrow=c(4,4),pty='s',yaxs='i',
        mar=c(1.5,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
        lwd=6, cex=3.5)
    for(i in 1:ncol(theta_pri))
      plotPriPos(theta_pri[,i],theta_pos[,i],paste0('theta_',i))
    if(!lasso) plotPriPos(lambda_pri,lambda_pos,'lambda')
    plotPriPos(sigma_pri,sigma_pos,'sigma')
    for(i in 1:ncol(yGP_pri))
      plotPriPos(yGP_pri[,i],yGP_pos[,i],paste0('yGP_',i),0.5*c(-1,1))
    dev.off()


  }
}
