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

# Color schemes ####
cols    = inlmisc::GetTolColors(8)
col_tr  = inlmisc::GetTolColors(8,alpha=0.1)
col_tr2 = inlmisc::GetTolColors(8,alpha=0.4) # For legends

# Misc. functions
expDecayModel <- function(x,c) {
  return( c[1] + c[2] * exp(-2*x/c[3]) )
}


# RUN ####
dataDirs = c("DataWl","Data1","DataSynth")[2]
for (dataDir in dataDirs) {

  dataSets = list.dirs(path=dataDir,full.names = FALSE)[-1]

  for(dataSet in dataSets[4]) {

    tag = paste0(dataDir,'_',dataSet)
    cat(tag,'------------------------','\n')

    ### Get Data
    D = read.csv(paste0(dataDir,'/',dataSet,'/Courbe.csv'))
    x=D[,1]; y=D[,2]
    # sel = x > 20 & x<=500 # Exclude aberrant points
    # x = x[sel]; y = y[sel]


    ### Estimate data uncertainty
    fits = FitOCTLib::estimateNoise(x, y)
    uy   = fits$uy      # Used by next stages
    ySpl = fits$ySmooth # Used by plotMonoExp


    ### Inference of exponential decay parameters
    fitm   = FitOCTLib::fitMonoExp(x, y, uy)
    theta0    = fitm$best.theta   # Used by next stage
    cor.theta = fitm$cor.theta
    source ("./plotMonoExp.R")

    ### Inference of modulated decay parameters
    Nn = 10
    lambda_rate = 0.1
    ru_theta    = 0.05
    priPost     = TRUE
    method      = c('sample','optim','vb')[1]

    # - Posterior Distribution
    fitGP = FitOCTLib::fitExpGP(x, y, uy, Nn,
                     method = method,
                     theta0 = theta0,
                     cor_theta = cor.theta,
                     ru_theta = ru_theta,
                     lambda_rate=lambda_rate,
                     nb_warmup = 100, nb_iter = 200,
                     open_progress = FALSE)
    fitOut = fitGP
    source ("./plotExpGP.R")

    # - Compare prior/posterior marginal pdfs
    if(priPost & fitGP$method != 'optim')
       source("./priPost.R")

  }
}
