# Plots
fit    = fitm$fit

theta   = fit$par$theta
resid   = fit$par$resid
mod     = fit$par$m
tagOut  = fitm$tag

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
