# Comparison of Hierarchical prior (Bayesian Lasso)
# and Elastic Net prior

rm(list = ls()); gc() # Clean environment

# Libraries ####
libs =c('parallel','rstan','inlmisc','purrr')

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

# Misc. functions ####
plotPri <- function(pri,tag,xlim=range(pri)) {
  d = density(pri)
  d$y = log10(d$y/max(d$y))
  plot(d, type = 'l', col = cols[4],
       main = tag,
       xlab = '', xlim = xlim,
       ylab = 'Norm. density', ylim = c(-4,0))
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[4],border=NA)
}
plotPri2 <- function(pri,pri2,tag,xlim=range(pri),legend=FALSE) {
  d = density(pri)
  d$y = log10(d$y/max(d$y))
  d2 = density(pri2)
  d2$y = log10(d2$y/max(d2$y))
  plot(d, type = 'l', col = cols[4],
       main = tag,
       xlab = '', xlim = xlim,
       ylab = 'Norm. density', ylim = c(-4,0))
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[4],border=NA)
  lines(d2, col=cols[6])
  polygon(d2$x,d2$y,rev(d2$w),0*d2$y,col=col_tr2[6],border=NA)
  if(legend)
    legend('topright',
           legend = c('HP','EN'),
           col = c(cols[4],cols[6]),
           lty = 1, cex=0.75, bty = 'n',
           pch=NA)
}
# RUN ####

Nn = 4

# Bayesian Lasso ####
modHP = (purrr::quietly(stan_model)(file = 'hierPrior.stan'))$result

stanData = list(
  Nn = Nn,
  lambda_rate= 0.05
)

init = list(
  yGP    = 0.01*rnorm(Nn),
  lambda = 20
)

pars = c('lambda','yGP')

fitHP = sampling(object  = modHP,
                 data    = stanData,
                 pars    = pars,
                 init    = function() {init},
                 control = list(adapt_delta=0.99,max_treedepth=13),
                 warmup  = 500, iter = 50000, chains = 4,
                 verbose=FALSE)

print(fitHP)
#
# yGP_pri     = extract(fitHP,'yGP')[[1]]
# png(filename = paste0('Results/pdfHP.png'),
#     width=2500, height=600)
# par(mfrow=c(1,4),pty='s',yaxs='i',
#     mar=c(1.5,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
#     lwd=6, cex=3.5)
# for(i in 1:ncol(yGP_pri))
#   plotPri(yGP_pri[,i],paste0('yGP_',i),1.5*c(-1,1))
# dev.off()

# Elastic Net ####

modEN = (purrr::quietly(stan_model)(file = 'elasticNetPrior.stan'))$result

stanData = list(
  Nn = Nn,
  lambda_rate= 0.1,
  beta = 0.5
)

pars = c('yGP')

fitEN = sampling(object  = modEN,
                 data    = stanData,
                 pars    = pars,
                 control = list(adapt_delta=0.99,max_treedepth=12),
                 warmup  = 500, iter = 50000, chains = 4,
                 verbose = FALSE
)

print(fitEN)

# yGP_pri     = extract(fitEN,'yGP')[[1]]
# png(filename = paste0('Results/pdfEN.png'),
#     width=2500, height=600)
# par(mfrow=c(1,4),pty='s',yaxs='i',
#     mar=c(1.5,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
#     lwd=6, cex=3.5)
# for(i in 1:ncol(yGP_pri))
#   plotPri(yGP_pri[,i],paste0('yGP_',i),1*c(-1,1))
# dev.off()


yGP_HP = extract(fitHP,'yGP')[[1]]
yGP_EN = extract(fitEN,'yGP')[[1]]
png(filename = paste0('Results/pdfHPvsEN.png'),
    width=2500, height=600)
par(mfrow=c(1,4),pty='s',yaxs='i',
    mar=c(2,3,1.6,.2),mgp=c(2,.75,0),tcl=-0.5,
    lwd=6, cex=3.5)
for(i in 1:ncol(yGP_HP))
  plotPri2(yGP_HP[,i],yGP_EN[,i],paste0('yGP_',i),1.5*c(-1,1),
           legend = ifelse(i==1,TRUE,FALSE))
dev.off()
