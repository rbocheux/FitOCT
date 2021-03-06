rm(list = ls()); gc() # Clean environment

# Libraries ####
libs =c('parallel','rstan','inlmisc','FitOCTLib')

for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

# Options ####
options(mc.cores = parallel::detectCores(), width = 90)
rstan_options(auto_write = TRUE)
set.seed(1234) # Initialise la graine du RNG

# Graphical parameters ####
gPars = list(
  cols    = inlmisc::GetColors(8),
  col_tr  = inlmisc::GetColors(8,alpha=0.1), # Light, for spaghetti curves
  col_tr2 = inlmisc::GetColors(8,alpha=0.4), # Darker, for legends...
  pty     = 's',
  mar     = c(3,3,1.6,.2),
  mgp     = c(2,.75,0),
  tcl     = -0.5,
  lwd     = 6,
  cex     = 3.5,
  modScale = 1,
  xlabel  = 'stromal depth (µm)',
  plot_title = '',
  graphTable = TRUE
)

# Control parameters ####

### Default values / Set values here
ctrlPars = list(
  depthSel    = NULL, # Otherwise c(xmin,xmax)
  dataType    = 2,    # Intensity
  subSample   = 1,
  smooth_df   = 15,
  method      = c('sample','optim','vb')[1],
  nb_warmup   = 500,
  nb_sample   = 1000,
  modRange    = 0.5,
  ru_theta    = 0.05,
  lambda_rate = 0.1,
  gridType    = 'internal',
  Nn          = 10,
  rho_scale   = 0.1,
  priPost     = TRUE, # Compare prior and posterior pdf ?
  priorType   = 'abc'
)

# Override parameters with control file
ctrlFile = 'ctrlParams.yaml'
if (file.exists(ctrlFile)) {
  ## Get from file
  lPars = rlist::list.load(ctrlFile)
  ## Expose parameters
  for (n in names(lPars))
    ctrlPars[[n]]= lPars[[n]]
}
cat('Configuration Parameters\n')
cat('------------------------\n')
str(ctrlPars,give.head=FALSE, give.length=FALSE)

# Expose parameters
for (n in names(ctrlPars))
  assign(n,rlist::list.extract(ctrlPars,n))

# RUN ####
dataDirs = c("DataWl","Data1","DataSynth")[2]
for (dataDir in dataDirs) {

  dataSets = list.dirs(path=dataDir,full.names = FALSE)[-1]

  for(dataSet in dataSets[4]) {

    tag = paste0(dataDir,'_',dataSet)
    cat(tag,'------------------------','\n')

    ### Get ans select Data
    D = read.csv(paste0(dataDir,'/',dataSet,'/Courbe.csv'))
    C = FitOCTLib::selX(D[,1],D[,2],depthSel,subSample)
    x = C$x; y = C$y

    ### Estimate data uncertainty
    fits = FitOCTLib::estimateNoise(x, y, df = smooth_df)
    uy   = fits$uy      # Used by next stages
    ySpl = fits$ySmooth # Used by plotMonoExp
    source ("./plotNoise.R")

    ### MAP Inference of exponential decay parameters
    fitm   = FitOCTLib::fitMonoExp(x, y, uy, dataType = dataType)
    theta0    = fitm$best.theta   # Used by next stage
    cor.theta = fitm$cor.theta
    source ("./plotMonoExp.R")

    if(is.null(br$alert)) break # Monoexp fit OK, no need to try fitExpGP

    ### Prior
    priExp = FitOCTLib::estimateExpPrior(
      x, uy, dataType, priorType,
      out = fitm, ru_theta = ru_theta,
      eps = 1e-3
    )

    # - Posterior Distribution
    fitGP = FitOCTLib::fitExpGP(
      x, y, uy,
      dataType      = dataType,
      Nn            = Nn,
      gridType      = gridType,
      method        = method,
      theta0        = priExp$theta0,
      Sigma0        = priExp$Sigma0,
      lambda_rate   = lambda_rate,
      rho_scale     = ifelse(rho_scale==0, 1./Nn, rho_scale),
      nb_warmup     = nb_warmup,
      nb_iter       = nb_warmup + nb_sample,
      prior_PD      = 0,
      open_progress = FALSE
    )
    fitOut = fitGP
    source ("./plotExpGP.R")

    # - Compare prior/posterior marginal pdfs
    if(priPost & fitGP$method != 'optim')
       source("./priPost.R")

  }
}
