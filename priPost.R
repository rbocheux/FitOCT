# - Prior (Predictive) Distribution
fitGP_pri = FitOCTLib::fitExpGP(
  x, y, uy,
  dataType      = dataType,
  Nn            = Nn,
  gridType      = gridType,
  method        = 'sample',
  theta0        = theta0,
  cor_theta     = cor.theta,
  ru_theta      = ru_theta,
  lambda_rate   = lambda_rate,
  rho_scale     = ifelse(rho_scale==0, 1./Nn, rho_scale),
  nb_warmup     = nb_warmup,
  nb_iter       = nb_warmup + nb_sample,
  prior_PD = 1,
  open_progress = FALSE
)
fitOut = fitGP_pri
source ("./plotExpGP.R")

png(filename = paste0('Results/',tag,'_pripos.png'),
    width=2400, height=2400)
FitOCTLib::plotPriPostAll(fitGP_pri$fit, fitGP$fit, gPars = gPars)
dev.off()
