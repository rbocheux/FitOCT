library(RandomFields)
RFoptions(spConform=FALSE)
library(inlmisc)
col_tr = inlmisc::GetTolColors(8,alpha=0.1)

set.seed(NULL)
n <- 10
xdat = seq(0.05,0.95,length.out = n)
ydat = rnorm(n)

np <- 100
xp <- seq(0,1,length.out = np)

par(mfrow=c(5,5),pty='s',
    mar=c(3,3,1.6,.2),
    mgp=c(2,.75,0),tcl=-0.5)

nRun=100
for (scale in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
  for (scv in c(0.05,0.1,0.5,1,2)) {
    set.seed(12345)
    plot(xdat,ydat,pch=20,col=2,ylim=c(-3,3),
         main=paste(scale,'/',scv))
    cond <- RFsimulate(model = RMgauss(var=scv*sd(ydat),scale=scale),
                       x=xp,
                       given = list(x=xdat),
                       data=list(y=ydat),
                       n=nRun )
    matlines(xp,cond,col=col_tr[4])
    grid(); box()
  }
}
