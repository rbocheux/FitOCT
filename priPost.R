# - Prior (Predictive) Distribution
fitGP_pri = FitOCTLib::fitExpGP(x, y, uy,
                                method = 'sample',
                                theta0 = theta0,
                                cor_theta = cor.theta,
                                ru_theta = ru_theta,
                                lambda_rate=lambda_rate,
                                prior_PD = 1,
                                open_progress = FALSE)
fitOut = fitGP_pri
source ("./plotExpGP.R")

png(filename = paste0('Results/pripos_',tag,'.png'),
    width=2400, height=2400)
FitOCTLib::plotPriPostAll(fitGP_pri$fit, fitGP$fit, gPars = gPars)
dev.off()
