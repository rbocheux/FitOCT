# Ctrl output
sink(file=paste0('Results/',tag,'_ctrl.txt'),append = TRUE)
opt = theta0
cat('\n MonoExp decay parameters:\n')
for(i in 1:length(opt))
  cat(paste0('b_',i,' : '),signif(opt[i],3),'\n')
cat('\n\n')

# Probability Interval for Birge's ratio
FitOCTLib::printBr(fitm$fit)
sink()

# Plots
fit     = fitm$fit
resid   = fit$par$resid
mod     = fit$par$m

tagOut  = paste0('Results/',tag,'_MonoExp')
png(filename = paste0(tagOut,'_results.png'),
    width=2000, height=1000)

FitOCTLib::plotMonoExp(x, y, uy, ySpl, mod, resid, gPars)

dev.off()
