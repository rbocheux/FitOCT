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

theta_pri   = extract(fitGP_pri$fit,'theta')[[1]]
yGP_pri     = extract(fitGP_pri$fit,'yGP')[[1]]
lambda_pri  = extract(fitGP_pri$fit,'lambda')[[1]]
sigma_pri   = extract(fitGP_pri$fit,'sigma')[[1]]
theta_pos   = extract(fitGP$fit,'theta')[[1]]
yGP_pos     = extract(fitGP$fit,'yGP')[[1]]
lambda_pos  = extract(fitGP$fit,'lambda')[[1]]
sigma_pos   = extract(fitGP$fit,'sigma')[[1]]

png(filename = paste0('Results/pripos_',tag,'.png'),
    width=2400, height=2400)
par(mfrow=c(4,4),pty='s',yaxs='i',
    mar=c(2,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
    lwd=6, cex=3.5)
for(i in 1:ncol(theta_pri))
  plotPriPos(theta_pri[,i],theta_pos[,i],paste0('theta_',i))
plotPriPos(lambda_pri,lambda_pos,'lambda')
plotPriPos(sigma_pri,sigma_pos,'sigma')
for(i in 1:ncol(yGP_pri))
  plotPriPos(yGP_pri[,i],yGP_pos[,i],paste0('yGP_',i),0.5*c(-1,1))
dev.off()
