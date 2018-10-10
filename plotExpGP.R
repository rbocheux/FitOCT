nMC = 100  # Nb of plotted spaghetti lines

fit    = fitGP$fit
sample = fitGP$sample
xGP    = fitGP$xGP

png(filename = paste0('Results/results_',tag,'.png'),
    width=2000, height=2000)
par(mfrow=c(2,2),pty='s',
    mar=c(3,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
    lwd=6, cex=3.5)

if(sample) {

  theta   = extract(fit,'theta')[[1]]
  resid   = extract(fit,'resid')[[1]]
  mod     = extract(fit,'m')[[1]]
  dL      = extract(fit,'dL')[[1]]
  yGP     = extract(fit,'yGP')[[1]]
  lp       = extract(fit,'lp__')[[1]]

  iMC = sample.int(nrow(lp),nMC)

  map = which.max(lp)
  y_map = mod[map,]

  # Fit
  plot(x,y,pch=20,cex=0.5,col=cols[6],
       main='Data fit',
       xlab='depth (a.u.)',
       ylab='mean OCT signal (a.u.)')
  if(nMC >0)
    for (i in 1:nMC)
      lines(x, mod[iMC[i],], col=col_tr[4])
  # lines(x,y_map,col=cols[7])
  lines(x,theta[map,1]+theta[map,2]*exp(-x/theta[map,3]),col=cols[7])

  legend('topright', bty='n',
         title = '(a) ', title.adj = 1,
         legend=c('data','expo. best fit','post. sample'),
         pch=c(20,NA,NA),lty=c(-1,1,1),
         col=c(cols[6],cols[7], col_tr2[4])
  )
  grid(lwd=3);box()

  # Residuals
  ylim=1.2*max(abs(resid))*c(-1,1)
  res = colMeans(resid)
  plot(x,res,type='n',
       ylim=ylim, main='Residuals',
       xlab='depth (a.u.)',
       ylab='residuals (a.u.)')
  grid(lwd=3); abline(h=0)
  polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
  points(x,res,pch=20,cex=0.75,col=cols[6])
  legend('topright', bty='n',
         title = '(b) ', title.adj = 1,
         legend=c('mean resid.','data 95% uncert.'),
         pch=c(20,NA),lty=c(-1,1),lwd=c(1,50),
         col=c(cols[6],col_tr2[4])
  )
  box()

  # Local deviations
  plot(x, dL[map,], type = 'n',
       ylim = c(-0.3,0.3),
       col  = cols[4],
       main='Deviation from mean depth',
       xlab = 'depth (a.u.)',
       ylab = 'relative deviation')
  abline(h=0); grid(lwd=3)
  if(nMC >0)
    for (i in 1:nMC)
      lines(x, dL[iMC[i],], col=col_tr[4])
  Q = t(apply(yGP,2,
              function(x)
                quantile(x,probs = c(0.025,0.25,0.75,0.975))
  )
  )
  segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])        # 95 %
  segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=24) # 50 %

  legend('topright', bty='n',
         title = '(c) ', title.adj = 1,
         legend=c('50% CI','95% CI','post. sample'),
         pch=NA ,lty=c(1,1,1),lwd=c(24,6,6),
         col=c(cols[6],cols[7],col_tr2[4])
  )
  box()


} else {

  theta = fit$par$theta
  mod   = fit$par$m
  resid = fit$par$resid
  dL    = fit$par$dL
  yGP   = fit$par$yGP

  plot(x,y,pch=20,cex=0.5,col=cols[6],
       main='Data fit',
       xlab='depth (a.u.)',
       ylab='mean OCT signal (a.u.)')
  lines(x,theta[1]+theta[2]*exp(-x/theta[3]),col=cols[7])
  lines(x,mod, col=cols[4])

  legend('topright', bty='n',
         title = '(a) ', title.adj = 1,
         legend=c('data','expo. best fit','best fit'),
         pch=c(20,NA,NA),lty=c(-1,1,1),
         col=c(cols[6],cols[7], col_tr2[4])
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
  legend('topright', bty='n',
         title = '(b) ', title.adj = 1,
         legend=c('resid.','data 95% uncert.'),
         pch=c(20,NA),lty=c(-1,1),lwd=c(1,50),
         col=c(cols[6],col_tr2[4])
  )
  box()

  # Local deviations
  plot(x, dL, type = 'l',
       ylim = c(-0.3,0.3),
       col  = cols[4],
       main='Deviation from mean depth',
       xlab = 'depth (a.u.)',
       ylab = 'relative deviation')

  abline(h=0); grid(lwd=3)

  if(!is.null(fit$theta_tilde)) {
    S = fit$theta_tilde
    iMC = sample.int(nrow(S),nMC)

    c = which(grepl(pattern = 'dL\\[', x=colnames(S)))
    for (i in 1:nMC)
      lines(x, S[iMC[i],c], col=col_tr[4])

    c = which(grepl(pattern = 'yGP\\[', x=colnames(S)))
    yGP = S[iMC,c]
    Q = t(apply(yGP,2,
                function(x)
                  quantile(x,probs = c(0.025,0.25,0.75,0.975))
      )
    )
    segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])        # 95 %
    segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=24) # 50 %
    legend('topright', bty='n',
           title = '(c) ', title.adj = 1,
           legend=c('50% CI','95% CI','post. sample'),
           pch=NA ,lty=c(1,1,1),lwd=c(24,6,6),
           col=c(cols[6],cols[7],col_tr2[4])
    )

  } else {
    points(xGP,yGP,pch=19,col=cols[7])
    segments(xGP,yGP,xGP,0*yGP,col=cols[7])
    legend('topright', bty='n',
           title = '(c) ', title.adj = 1,
           legend=c('ctrl points','modulation'),
           pch=c(19,NA) ,lty=c(1,1),lwd=c(-1,6),
           col=c(cols[7],cols[4])
    )
  }


  box()

}

# Plot true modulation for synthetic signals
fName = paste0(dataDir,'/',dataSet,'/Modulation.csv')
if(file.exists(fName)) {
  M = read.csv(fName)
  lines(M[,1],M[,2],lty=2)
}

dev.off()


