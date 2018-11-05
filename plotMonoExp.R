# Ctrl output
sink(file=paste0('Results/',tag,'_ctrl.txt'))
cat('\n Noise fit parameters:\n')
a = fits$theta
for(i in 1:length(a))
  cat(paste0('a_',i,' : '),signif(a[i],3),'\n')
cat('\n\n')

opt = theta0
cat('\n MonoExp decay parameters:\n')
for(i in 1:length(opt))
  cat(paste0('b_',i,' : '),signif(opt[i],3),'\n')
br = fitm$fit$par$br
N  = length(x)
Np = 3
ndf = N - Np
CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf)) / ndf
if(prod(CI95-br) >= 0)
  cat('\n!!! WARNING !!! \n')
cat('br       :',signif(br,2),'\n')
cat('CI95(br) :',paste0(signif(CI95,2),collapse='-'))
cat('\n\n')
sink()

# Plots
fit    = fitm$fit

theta   = fit$par$theta
resid   = fit$par$resid
mod     = fit$par$m
tagOut  = paste0('Results/',tag,'_MonoExp')

png(filename = paste0(tagOut,'_results.png'),
    width=2000, height=1000)
par(mfrow=c(1,2),pty='s',
    mar=c(3,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
    lwd=6, cex=3.5)

# Fit
plot(x,y,pch=20,cex=0.5,col=cols[6],
     main='Data fit',
     xlab='depth (a.u.)',
     ylab='mean OCT signal (a.u.)')
lines(x,mod,col=cols[7])
legend('topright', bty='n',
       title = '(a) ', title.adj = 1,
       legend=c('data','best fit'),
       pch=c(20,NA),lty=c(-1,1),
       col=c(cols[6],cols[7])
)
grid(lwd=3);box()

# Residus
ylim=1.2*max(abs(resid))*c(-1,1)
res = resid
plot(x,res,type='n',
     ylim=ylim, main='Residuals',
     xlab='depth (a.u.)',
     ylab='residuals (a.u.)')
grid(lwd=3); abline(h=0)
polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
points(x,res,pch=20,cex=0.75,col=cols[6])
lines(x, ySpl-mod, col=cols[7])
legend('topright', bty='n',
       title = '(b) ', title.adj = 1,
       legend=c('mean resid.','data 95% uncert.','best fit - smooth'),
       pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,50,6),
       col=c(cols[6],col_tr2[4],cols[7])
)
box()

dev.off()
