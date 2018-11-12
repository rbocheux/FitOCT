# Libraries ####
libs =c('shiny','parallel','rstan','knitr',
        'inlmisc','shinycssloaders','DT',
        'RandomFields','stringr','devtools')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

lib ='FitOCTLib'
if(!require(lib,character.only = TRUE))
  devtools::install_github("ppernot/FitOCTlib")
library(lib,character.only = TRUE)

# Options ####
Sys.setlocale(category = "LC_NUMERIC", locale="C")
options(mc.cores = parallel::detectCores(),
        width = 60,
        warn  = 0)
rstan::rstan_options(auto_write = TRUE)
RandomFields::RFoptions(spConform=FALSE)

# set.seed(1234) # Initialise la graine du RNG

# Color schemes ####
cols    = inlmisc::GetColors(8)
# Transparents for spaghetti
col_tr  = inlmisc::GetColors(8,alpha=0.1)
# Darker for legends or fillings
col_tr2 = inlmisc::GetColors(8,alpha=0.4)

# Graphical parameters and functions ####
pty = 's'
mar = c(3,3,1.5,.5)
mgp = c(2,.75,0)
tcl = -0.5
lwd = 2
cex = 1

plotPriPos <- function(pri,pos,tag,xlim=range(c(pri,pos))) {
  # Plot overlapped densities for 2 samples
  d = density(pri)
  d$y = d$y/max(d$y)
  plot(d, type = 'l', col = cols[4],
       main = tag,
       xlab = '', xlim = xlim,
       ylab = 'Norm. density', ylim = c(0,1.1), yaxs = 'i')
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[4],border=NA)
  d = density(pos)
  d$y = d$y/max(d$y)
  lines(d$x,d$y,col=cols[6])
  polygon(d$x,d$y,rev(d$w),0*d$y,col=col_tr2[6],border=NA)
}

expDecayModel <- function(x,c) {
  return( c[1] + c[2] * exp(-2*x/c[3]) )
}

# Global variables ####
Inputs = reactiveValues(
  x           = NULL,
  y           = NULL,
  xSel        = NULL,
  outSmooth   = NULL,
  outMonoExp  = NULL,
  outExpGP    = NULL,
  outPriExpGP = NULL,
  fitOut      = NULL
)

# MAIN ####
function(input, output, session) {

  # Directory to store stan logs
  session_dir = file.path(tempdir(),stringr::str_sub(session$token, 1, 8))
  dir.create(session_dir, showWarnings = FALSE)
  session$onSessionEnded(function() {unlink(session_dir, TRUE)})

  observeEvent(
    input$dataFile,
    {
      dat = try(
        read.csv(input$dataFile[['datapath']]),
        silent = TRUE
      )
      if(class(dat) == 'try-error')
        return(NULL)

      # Update depts range selector
      rangeX = range(dat[,1])
      updateSliderInput(
        session,
        inputId = 'depthSel',
        min     = round(rangeX[1]),
        max     = round(rangeX[2]),
        value   = round(rangeX),
        step    = 1
      )

      # Store data and empty buffers
      Inputs$x          <<- dat[,1]
      Inputs$y          <<- dat[,2]
      Inputs$outSmooth  <<- NULL
      Inputs$outMonoExp <<- NULL
      Inputs$outExpGP   <<- NULL
      Inputs$fitOut     <<- NULL
      Inputs$xSel       <<- 1:length(Inputs$x)
    }
  )

  selX <- function(x,y) {
    # Apply selectors to inputs

    xSel = which(x >= input$depthSel[1] &
                 x <= input$depthSel[2]  )
    x = x[xSel]; y = y[xSel]

    if(input$subSample != 1) {
      xSel = seq(1,length(x),by=input$subSample)
      x = x[xSel]; y = y[xSel]
    }
    return(list(x=x,y=y))
  }

  output$plotNoise   <- renderPlot({
    if(is.null(Inputs$x))
      return(NULL)

    C = selX(x = Inputs$x, y = Inputs$y)
    x = C$x; y = C$y

    if(!is.finite(IQR(x)))
      return(NULL) # Protect smooth.spline from zero tol

    out = FitOCTLib::estimateNoise(
      x  = x,
      y  = y,
      df = input$smooth_df
      )

    Inputs$outSmooth  <<- out
    Inputs$outMonoExp <<- NULL
    Inputs$outExpGP   <<- NULL

    ySmooth = out$ySmooth
    resid   = y - ySmooth
    uy      = out$uy

    par(mfrow=c(1,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    # Smooth Fit
    plot(x,y,pch=20,cex=0.5,col=cols[6],
         main='Data fit',
         xlab='depth (a.u.)',
         ylab='mean OCT signal (a.u.)')
    grid()
    lines(x,ySmooth,col=cols[7])
    legend('topright', bty='n',
           title = '', title.adj = 1,
           legend=c('data','smoother'),
           pch=c(20,NA),lty=c(-1,1),
           col=c(cols[6],cols[7])
    )
    box()

    # Residuals
    ylim=1.2*max(abs(resid))*c(-1,1)
    res = resid
    plot(x,res,type='n',
         ylim=ylim, main='Residuals',
         xlab='depth (a.u.)',
         ylab='residuals (a.u.)')
    grid()
    abline(h=0)
    polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
    points(x,res,pch=20,cex=0.75,col=cols[6])
    legend('topright', bty='n',
           title = '', title.adj = 1,
           legend=c('resid.','data 95% uncert.'),
           pch=c(20,NA),lty=c(-1,1),lwd=c(1,10),
           col=c(cols[6],col_tr2[4])
    )
    box()

  })

  output$plotMonoExp <- renderPlot({
    if(is.null(Inputs$outSmooth))
      return(NULL)

    C = selX(x = Inputs$x, y = Inputs$y)
    x = C$x; y = C$y

    out     = Inputs$outSmooth
    ySmooth = out$ySmooth
    uy      = out$uy

    out = FitOCTLib::fitMonoExp(
      x  = x,
      y  = y,
      uy = uy
    )

    Inputs$outMonoExp <<- out
    Inputs$outExpGP   <<- NULL

    fit    = out$fit
    theta  = fit$par$theta
    resid  = fit$par$resid
    mod    = fit$par$m

    par(mfrow=c(1,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    # Fit
    plot(x,y,pch=20,cex=0.5,col=cols[6],
         main='Data fit',
         xlab='depth (a.u.)',
         ylab='mean OCT signal (a.u.)')
    grid()
    lines(x,mod,col=cols[7])
    legend('topright', bty='n',
           title = '', title.adj = 1,
           legend=c('data','best fit'),
           pch=c(20,NA),lty=c(-1,1),
           col=c(cols[6],cols[7])
    )
    box()

    # Residus
    ylim=1.2*max(abs(resid))*c(-1,1)
    res = resid
    plot(x,res,type='n',
         ylim=ylim, main='Residuals',
         xlab='depth (a.u.)',
         ylab='residuals (a.u.)')
    grid()
    abline(h=0)
    polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
    points(x,res,pch=20,cex=0.75,col=cols[6])
    lines(x, ySmooth-mod, col=cols[7])
    legend('topright', bty='n',
           title = '(b) ', title.adj = 1,
           legend=c('mean resid.','data 95% uncert.','best fit - smooth'),
           pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,10,2),
           col=c(cols[6],col_tr2[4],cols[7])
    )
    box()

  })

  output$resMonoExp  <- renderPrint({
    if (is.null(Inputs$outMonoExp))
      return(NULL)

    out  = Inputs$outSmooth
    a    = out$theta

    outm = Inputs$outMonoExp
    fit  = outm$fit

    if(outm$method == 'optim') {
      for(i in 1:length(a))
        cat(paste0('a_',i,' : '),signif(a[i],3),'\n')

      opt = unlist(fit$par[['theta']],use.names = TRUE)
      for(i in 1:length(opt))
        cat(paste0('b_',i,' : '),signif(opt[i],3),'\n')

      br = fit$par$br
    } else {
      br = mean(extract(fit,pars='br')[[1]])
    }

    # Probability Interval for Birge's ratio
    N  = length(Inputs$x)
    Np = 3
    ndf = N - Np
    CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf)) / ndf

    if(prod(CI95-br) >= 0)
      cat('\n!!! WARNING !!! \n')
    cat('br       :',signif(br,2),'\n')
    cat('CI95(br) :',paste0(signif(CI95,2),collapse='-'))

  })

  runExpGP <-function(prior_PD = 0) {
    if(is.null(Inputs$outMonoExp))
      return(NULL)

    outm  = Inputs$outMonoExp

    isolate({
      C = selX(x = Inputs$x, y = Inputs$y)
      x = C$x; y = C$y
    })

    log_file = file.path(session_dir, "stan.log")
    dummy = suppressWarnings(unlink(log_file))
    sink(log_file)
    out <- FitOCTLib::fitExpGP(
      x         = x,
      y         = y,
      uy        = Inputs$outSmooth[['uy']],
      method    = ifelse(prior_PD == 0,
                         input$method,
                         'sample'),
      theta0    = outm$best.theta,
      cor_theta = outm$cor.theta,
      ru_theta  = input$ru_theta,
      nb_warmup = input$nb_warmup,
      prior_PD  =  prior_PD,
      nb_iter   = input$nb_warmup + input$nb_sample,
      Nn        = input$Nn,
      rho_scale = ifelse(input$rho_scale==0,
                         1. / input$Nn,
                         input$rho_scale),
      lambda_rate = input$lambda_rate,
      gridType  = input$gridType,
      open_progress = FALSE
    )
    sink()

    return(out)
  }

  doExpGP <- eventReactive(
    input$runExpGP,
    {
      out = runExpGP()
      Inputs$outExpGP <<- out
      return(out)
    }
  )

  doPriExpGP <- eventReactive(
    input$runPriExpGP,
    {
      out = runExpGP(prior_PD = 1)
      Inputs$outPriExpGP <<- out
      return(out)
    }
  )

  do_progress = function(file) {
    if (!file.exists(file))
      return(NULL)
    r = readLines(file, warn = FALSE)
    if (length(r) == 0)
      return(NULL)
    r = unlist(stringr::str_extract_all(r, "Chain \\d+.*"))
    r = r[length(r)]
    frac_s = stringr::str_match(r, "(\\d+)%")
    if (nrow(frac_s) == 0) return(NULL)
    frac = as.numeric(frac_s[1,2])
    chain = as.integer(stringr::str_match(r, "Chain (\\d+)")[1,2])
    complete = floor(((chain - 1)*100 + frac)/4)
    print(complete)
    return(complete)
  }

  file_info = reactiveFileReader(
    intervalMillis = 1000,
    session       = session,
    filePath      = file.path(session_dir, "stan.log"),
    readFunc      = readLines
  )

  output$outExpGP <- renderPrint({
    file_info()
    #cat(file_info(),'%')
  })
  #outputOptions(output, "outExpGP",suspendWhenHidden = FALSE)

  output$plotExpGP   <- renderPlot({
    if (is.null(out <- doExpGP()))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    C = selX(x = Inputs$x, y = Inputs$y)
    x = C$x; y = C$y

    outS    = Inputs$outSmooth
    ySmooth = outS$ySmooth
    uy      = outS$uy

    fit      = out$fit
    method   = out$method
    xGP      = out$xGP
    prior_PD = out$prior_PD
    lasso    = out$lasso

    modScale = input$modRange

    nMC = 100

    par(mfrow=c(2,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    if(method == 'sample') {

      theta   = extract(fit,'theta')[[1]]
      yGP     = extract(fit,'yGP')[[1]]
      sigma   = mean(extract(fit,'sigma')[[1]])
      if(prior_PD == 0) {
        resid   = extract(fit,'resid')[[1]]
        mod     = extract(fit,'m')[[1]]
        dL      = extract(fit,'dL')[[1]]
        lp      = extract(fit,'lp__')[[1]]
        map     = which.max(lp)
        y_map   = mod[map,]
      }

      iMC = sample.int(nrow(theta),nMC)

      # Fit
      plot(x,y,pch=20,cex=0.5,col=cols[6],
           main='Data fit',
           xlab='depth (a.u.)',
           ylab='mean OCT signal (a.u.)')
      grid()
      if(prior_PD == 0) {
        if(nMC >0) {
          for (i in 1:nMC)
            lines(x, mod[iMC[i],], col=col_tr[4])
        }
        # Calculate AVerage Exponential Decay
        mExp = x*0
        for (i in 1:nrow(theta))
           mExp = mExp + expDecayModel(x,theta[i,1:3])
        mExp = mExp/nrow(theta)
        lines(x,mExp,col=cols[7])
      } else {
        if(nMC >0)
          for (i in 1:nMC)
            lines(x,expDecayModel(x,theta[i,1:3]),col=col_tr[7])
      }

      legend('topright', bty='n',
             legend=c('data','mean exp. fit','post. sample'),
             pch=c(20,NA,NA),lty=c(-1,1,1),
             col=c(cols[6],cols[7], col_tr2[4])
      )
      box()

      if(prior_PD == 0) {
        # Residuals
        res = colMeans(resid)
        ylim=1.2*max(abs(res))*c(-1,1)
        plot(x,res,type='n',
             ylim=ylim, main='Residuals',
             xlab='depth (a.u.)',
             ylab='residuals (a.u.)')
        grid(lwd=3); abline(h=0)
        polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
        points(x,res,pch=20,cex=0.75,col=cols[6])
        lines(x, ySmooth-y_map, col=cols[7])
        legend('topright', bty='n',
               legend=c('mean resid.','data 95% uncert.','smooth - fit'),
               pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
               col=c(cols[6],col_tr2[4],cols[7])
        )
        box()
      }

      # Local deviations
      if(prior_PD == 0) {
        plot(x, dL[map,], type = 'n',
             ylim = modScale * c(-1,1),
             col  = cols[4],
             main='Deviation from mean depth',
             xlab = 'depth (a.u.)',
             ylab = 'relative deviation')
        abline(h=0)
        grid()
        if(nMC >0)
          for (i in 1:nMC)
            lines(x, dL[iMC[i],], col=col_tr[4])

      } else {
        plot(x, x, type = 'n',
             ylim = modScale*c(-1,1),
             col  = cols[4],
             main='Deviation from mean depth',
             xlab = 'depth (a.u.)',
             ylab = 'relative deviation')
        abline(h=0)
        grid()
      }

      Q = t(apply(yGP,2,
                  function(x)
                    quantile(x,probs = c(0.025,0.25,0.75,0.975))
      )
      )
      segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])       # 95 %
      segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=6) # 50 %

      legend('topright', bty='n',
             legend=c('50% CI','95% CI','post. sample'),
             pch=NA ,lty=c(1,1,1),lwd=c(6,2,2),
             col=c(cols[6],cols[7],col_tr2[4])
      )
      box()

    } else {

      theta = fit$par$theta
      mod   = fit$par$m
      resid = fit$par$resid
      dL    = fit$par$dL
      yGP   = fit$par$yGP
      sigma = fit$par$sigma

      plot(x,y,pch=20,cex=0.5,col=cols[6],
           main='Data fit',
           xlab='depth (a.u.)',
           ylab='mean OCT signal (a.u.)')
      grid()
      lines(x,expDecayModel(x,theta),col=cols[7])
      lines(x,mod, col=cols[4])

      legend('topright', bty='n',
             legend=c('data','expo. best fit','best fit'),
             pch=c(20,NA,NA),lty=c(-1,1,1),
             col=c(cols[6],cols[7], col_tr2[4])
      )
      box()

      # Residus
      ylim=1.2*max(abs(resid))*c(-1,1)
      res = resid
      plot(x,res,type='n',
           ylim=ylim, main='Residuals',
           xlab='depth (a.u.)',
           ylab='residuals (a.u.)')
      grid()
      abline(h=0)
      polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
      points(x,res,pch=20,cex=0.75,col=cols[6])
      lines(x, ySmooth - mod, col=cols[7])
      legend('topright', bty='n',
             legend=c('mean resid.','data 95% uncert.','smooth - fit'),
             pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
             col=c(cols[6],col_tr2[4],cols[7])
      )
      box()

      # Local deviations
      plot(x, dL, type = 'l',
           ylim = modScale*c(-1,1),
           col  = cols[4],
           main='Deviation from mean depth',
           xlab = 'depth (a.u.)',
           ylab = 'relative deviation')
      grid()
      abline(h=0)

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
        segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])       # 95 %
        segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=6) # 50 %
        legend('topright', bty='n',
               legend=c('50% CI','95% CI','post. sample'),
               pch=NA ,lty=c(1,1,1),lwd=c(6,2,2),
               col=c(cols[6],cols[7],col_tr2[4])
        )

      } else {
        points(xGP,yGP,pch=19,col=cols[7])
        segments(xGP,yGP,xGP,0*yGP,col=cols[7])
        legend('topright', bty='n',
               legend=c('ctrl points','modulation'),
               pch=c(19,NA) ,lty=c(1,1),lwd=c(-1,2),
               col=c(cols[7],cols[4])
        )
      }
      box()
    }

    if(prior_PD == 0) {
      # Plot true modulation for synthetic signals
      fName = paste0('./Modulation.csv')
      if(file.exists(fName)) {
        M = read.csv(fName)
        lines(M[,1],M[,2],lty=2)
      }
    }

    resw = res/(uy*sigma)
    xlim = range(resw)
    qqnorm(resw, xlim=xlim,ylim=xlim,
           main='Norm. Q-Q plot of wghtd resid.',
           col=cols[6], pch=20)
    abline(a=0,b=1,lty=2)
    grid();box()

  })

  output$plotPriExpGP   <- renderPlot({
    if (is.null(out <- doPriExpGP()))
      return(NULL)

    if(is.null(Inputs$outPriExpGP))
      return(NULL)

    C = selX(x = Inputs$x, y = Inputs$y)
    x = C$x; y = C$y

    outS    = Inputs$outSmooth
    ySmooth = outS$ySmooth
    uy      = outS$uy

    fit      = out$fit
    method   = out$method
    xGP      = out$xGP
    prior_PD = out$prior_PD
    lasso    = out$lasso

    modScale = input$modRange

    nMC = 100

    par(mfrow=c(2,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    if(method == 'sample') {

      theta   = extract(fit,'theta')[[1]]
      yGP     = extract(fit,'yGP')[[1]]
      if(prior_PD == 0) {
        resid   = extract(fit,'resid')[[1]]
        mod     = extract(fit,'m')[[1]]
        dL      = extract(fit,'dL')[[1]]
        lp      = extract(fit,'lp__')[[1]]
        map     = which.max(lp)
        y_map   = mod[map,]
      }

      iMC = sample.int(nrow(theta),nMC)

      # Fit
      plot(x,y,pch=20,cex=0.5,col=cols[6],
           main='Data fit',
           xlab='depth (a.u.)',
           ylab='mean OCT signal (a.u.)')
      grid()
      if(prior_PD == 0) {
        if(nMC >0) {
          for (i in 1:nMC)
            lines(x, mod[iMC[i],], col=col_tr[4])
        }
        # Calculate AVerage Exponential Decay
        mExp = x*0
        for (i in 1:nrow(theta))
          mExp = mExp + expDecayModel(x,theta[i,1:3])
        mExp = mExp/nrow(theta)
        lines(x,mExp,col=cols[7])
      } else {
        if(nMC >0)
          for (i in 1:nMC)
            lines(x,expDecayModel(x,theta[i,1:3]),col=col_tr[7])
      }

      legend('topright', bty='n',
             legend=c('data','mean exp. fit','post. sample'),
             pch=c(20,NA,NA),lty=c(-1,1,1),
             col=c(cols[6],cols[7], col_tr2[4])
      )
      box()

      if(prior_PD == 0) {
        # Residuals
        res = colMeans(resid)
        ylim=1.2*max(abs(res))*c(-1,1)
        plot(x,res,type='n',
             ylim=ylim, main='Residuals',
             xlab='depth (a.u.)',
             ylab='residuals (a.u.)')
        grid(lwd=3); abline(h=0)
        polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
        points(x,res,pch=20,cex=0.75,col=cols[6])
        lines(x, ySmooth-y_map, col=cols[7])
        legend('topright', bty='n',
               legend=c('mean resid.','data 95% uncert.','smooth - fit'),
               pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
               col=c(cols[6],col_tr2[4],cols[7])
        )
        box()
      }

      # Local deviations
      if(prior_PD == 0) {
        plot(x, dL[map,], type = 'n',
             ylim = modScale * c(-1,1),
             col  = cols[4],
             main='Deviation from mean depth',
             xlab = 'depth (a.u.)',
             ylab = 'relative deviation')
        abline(h=0)
        grid()
        if(nMC >0)
          for (i in 1:nMC)
            lines(x, dL[iMC[i],], col=col_tr[4])

      } else {
        plot(x, x, type = 'n',
             ylim = modScale*c(-1,1),
             col  = cols[4],
             main='Deviation from mean depth',
             xlab = 'depth (a.u.)',
             ylab = 'relative deviation')
        abline(h=0)
        grid()
      }

      Q = t(apply(yGP,2,
                  function(x)
                    quantile(x,probs = c(0.025,0.25,0.75,0.975))
      )
      )
      segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])       # 95 %
      segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=6) # 50 %

      legend('topright', bty='n',
             legend=c('50% CI','95% CI','post. sample'),
             pch=NA ,lty=c(1,1,1),lwd=c(6,2,2),
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
      grid()
      lines(x,expDecayModel(x,theta),col=cols[7])
      lines(x,mod, col=cols[4])

      legend('topright', bty='n',
             legend=c('data','expo. best fit','best fit'),
             pch=c(20,NA,NA),lty=c(-1,1,1),
             col=c(cols[6],cols[7], col_tr2[4])
      )
      box()

      # Residus
      ylim=1.2*max(abs(resid))*c(-1,1)
      res = resid
      plot(x,res,type='n',
           ylim=ylim, main='Residuals',
           xlab='depth (a.u.)',
           ylab='residuals (a.u.)')
      grid()
      abline(h=0)
      polygon(c(x,rev(x)),c(-2*uy,rev(2*uy)),col=col_tr2[4],border = NA)
      points(x,res,pch=20,cex=0.75,col=cols[6])
      lines(x, ySmooth - mod, col=cols[7])
      legend('topright', bty='n',
             legend=c('mean resid.','data 95% uncert.','smooth - fit'),
             pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
             col=c(cols[6],col_tr2[4],cols[7])
      )
      box()

      # Local deviations
      plot(x, dL, type = 'l',
           ylim = modScale*c(-1,1),
           col  = cols[4],
           main='Deviation from mean depth',
           xlab = 'depth (a.u.)',
           ylab = 'relative deviation')
      grid()
      abline(h=0)

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
        segments(xGP,Q[,1],xGP,Q[,4],col=cols[7])       # 95 %
        segments(xGP,Q[,2],xGP,Q[,3],col=cols[6],lwd=6) # 50 %
        legend('topright', bty='n',
               legend=c('50% CI','95% CI','post. sample'),
               pch=NA ,lty=c(1,1,1),lwd=c(6,2,2),
               col=c(cols[6],cols[7],col_tr2[4])
        )

      } else {
        points(xGP,yGP,pch=19,col=cols[7])
        segments(xGP,yGP,xGP,0*yGP,col=cols[7])
        legend('topright', bty='n',
               legend=c('ctrl points','modulation'),
               pch=c(19,NA) ,lty=c(1,1),lwd=c(-1,2),
               col=c(cols[7],cols[4])
        )
      }
      box()
    }

    if(prior_PD == 0) {
      # Plot true modulation for synthetic signals
      fName = paste0('./Modulation.csv')
      if(file.exists(fName)) {
        M = read.csv(fName)
        lines(M[,1],M[,2],lty=2)
      }
    }

  })

  output$summaryPriOut  <- DT::renderDataTable({
    if (is.null(out <- try(doPriExpGP(),silent=TRUE)))
      return(NULL)

    if(is.null(Inputs$outPriExpGP))
      return(NULL)

    fit  = out$fit

    pars = c('theta','yGP','lambda','sigma')

    sum  = rstan::summary(fit,pars=pars,
                          use_cache=FALSE,
                          c(0.025,0.5,0.975))$summary[,-2]
    `%>%` <- DT::`%>%`
    DT::datatable(data = signif(sum,digits = 3),
                  options = list(
                    ordering    = FALSE,
                    searching   = FALSE,
                    paging      = FALSE,
                    info        = FALSE,
                    pageLength  = 16,
                    deferRender = TRUE,
                    scrollY     = FALSE,
                    scrollX     = TRUE,
                    stateSave   = FALSE
                  ) ) %>%
      DT::formatStyle(columns = "Rhat",
                      color = DT::styleInterval(1.1, c("blue", "red"))) %>%
      DT::formatRound(columns = "n_eff", digits = 0)

  })

  output$tracesPriExpGP <- renderPlot({
    if (is.null(out <- doPriExpGP()))
      return(NULL)

    if(is.null(Inputs$outPriExpGP))
      return(NULL)

    fit  = out$fit
    pars = c('theta','yGP','lambda','sigma')
    print(rstan::traceplot(fit, inc_warmup=TRUE, pars = pars))

  })

  nzCtrlPts   <- function(p=0.95){
    # Count non-zero (at p% level) ctrl points
    if (is.null(out <- doExpGP()))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    fit      = out$fit
    method   = out$method

    # Get a sample or return NULL
    if(method == 'sample') {
      yGP     = extract(fit,'yGP')[[1]]

    } else {

      if(!is.null(fit$theta_tilde)) {
        S = fit$theta_tilde
        c = which(grepl(pattern = 'yGP\\[', x=colnames(S)))
        yGP = S[,c]

      } else {
        return(NULL)

      }
    }

    # Estimate p% CI on yGP
    Q = t(
      apply(
        yGP,2,
        function(x)
          quantile(x,probs = 0.5 + 0.5*p*c(-1,1))
      )
    )

    # Count zeros out of p% CI
    nz = sum( apply(Q,1,prod) > 0)

    return(nz)
  }

  output$resExpGP    <- renderPrint({
    if (is.null(out <- doExpGP()))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    fit = out$fit

    if(out$method == 'optim') {
      br = fit$par$br
    } else {
      br = mean(extract(fit,pars='br')[[1]])
    }

    # Probability Interval for Birge's ratio
    N  = length(Inputs$x)
    Np = 3
    Nn = length(out$xGP)
    ndf = N - (Np + Nn + 2)
    CI95 = c(qchisq(0.025,df=ndf),qchisq(0.975,df=ndf)) / ndf

    cat('raw br\n')
    cat('--------------------------------------\n')
    if(prod(CI95-br) >= 0)
      cat('!!! WARNING !!! \n')
    cat('br       :',signif(br,2),' ( ndf =',ndf,')\n')
    cat('CI95(br) :',paste0(signif(CI95,2),collapse='-'),'\n')

    # Correct br by counting only non-zero ctrl points
    cat('\n')
    cat('correction from inactive ctrl points\n')
    cat('------------------------------------\n')

    nz = nzCtrlPts(0.5)
    if(is.null(nz)) {
      # Let the user decide by himself
      for(n in rev(0:Nn)) {
        nf = N - (Np + n + 2)
        CI95 = c(qchisq(0.025,df=nf),qchisq(0.975,df=nf)) / nf
        br1  = br * ndf / nf
        if(prod(CI95-br1) < 0)
          break
      }
      cat('Fit OK if there are less than \n',
          n+1,' active ctrl points\n')

    } else {
      cat('Estimated',nz,'active ctrl points\n')
      nf = N - (Np + nz + 2)
      CI95 = c(qchisq(0.025,df=nf),qchisq(0.975,df=nf)) / nf
      br1  = br * ndf / nf
      if(prod(CI95-br1) >= 0)
        cat('!!! WARNING !!! \n')
      cat('br       :',signif(br1,2),' ( ndf =',nf,')\n')
      cat('CI95(br) :',paste0(signif(CI95,2),collapse='-'),'\n')
    }

  })

  output$summaryOut  <- DT::renderDataTable({
    if (is.null(out <- try(doExpGP(),silent=TRUE)))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    fit  = out$fit

    if (out$method == 'optim') {

      pars = c('theta','yGP','lambda','sigma')
      opt = list()
      for (par in pars)
        opt[[par]] = fit$par[[par]]
      opt = unlist(opt,use.names = TRUE)

      if(!is.null(fit$hessian)) {
        H = fit$hessian
        tags = colnames(H)
        tags = gsub('\\.','',tags)
        colnames(H) = rownames(H) = tags
        se = list()
        for (par in names(opt))
          se[[par]]  = sqrt(-1/H[par,par])
        se = unlist(se)
      }
      sum = data.frame(opt = opt, sd = se)

      # br = fit$par$br

      DT::datatable(data = signif(sum,digits = 3),
                    options = list(
                      ordering    = FALSE,
                      searching   = FALSE,
                      paging      = FALSE,
                      info        = FALSE,
                      pageLength  = 16,
                      deferRender = TRUE,
                      scrollY     = FALSE,
                      scrollX     = TRUE,
                      # scrollCollapse = FALSE,
                      stateSave   = FALSE
                    ) )

    } else {

      pars = c('theta','yGP','lambda','sigma','br')

      sum  = rstan::summary(fit,pars=pars,
                            use_cache=FALSE,
                            c(0.025,0.5,0.975))$summary[,-2]
      # br = mean(extract(fit,pars='br')[[1]])

      `%>%` <- DT::`%>%`
      DT::datatable(data = signif(sum,digits = 3),
                    options = list(
                      ordering    = FALSE,
                      searching   = FALSE,
                      paging      = FALSE,
                      info        = FALSE,
                      pageLength  = 16,
                      deferRender = TRUE,
                      scrollY     = FALSE,
                      scrollX     = TRUE,
                      # scrollCollapse = FALSE,
                      stateSave   = FALSE
                    ) ) %>%
        DT::formatStyle(columns = "Rhat",
                        color = DT::styleInterval(1.1, c("blue", "red"))) %>%
        DT::formatRound(columns = "n_eff", digits = 0)
    }


  })

  output$tracesExpGP <- renderPlot({
    if (is.null(out <- doExpGP()))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    if (out$method == 'optim')
      return(NULL)

    fit  = out$fit
    pars = c('theta','yGP','lambda','sigma','br')
    print(rstan::traceplot(fit, inc_warmup=TRUE, pars = pars))

  })

  output$priPostExpGP   <- renderPlot({
    if (is.null(fitGP <- doExpGP()) ||
        is.null(fitGP_pri <- doPriExpGP()))
      return(NULL)
    if(fitGP$method != 'sample')
      return(NULL)

    theta_pri   = rstan::extract(fitGP_pri$fit,'theta')[[1]]
    yGP_pri     = rstan::extract(fitGP_pri$fit,'yGP')[[1]]
    lambda_pri  = rstan::extract(fitGP_pri$fit,'lambda')[[1]]
    sigma_pri   = rstan::extract(fitGP_pri$fit,'sigma')[[1]]

    theta_pos   = rstan::extract(fitGP$fit,'theta')[[1]]
    yGP_pos     = rstan::extract(fitGP$fit,'yGP')[[1]]
    lambda_pos  = rstan::extract(fitGP$fit,'lambda')[[1]]
    sigma_pos   = rstan::extract(fitGP$fit,'sigma')[[1]]

    nPar = ncol(theta_pri) + ncol(yGP_pri) + 2
    ncol = 4
    nrow = floor(nPar / ncol)
    if(ncol*nrow < nPar) nrow = nrow + 1
    par(mfrow=c(nrow,ncol),pty=pty,mar=mar,
        mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    for(i in 1:ncol(theta_pri))
      plotPriPos(theta_pri[,i],theta_pos[,i],paste0('theta_',i))
    plotPriPos(lambda_pri,lambda_pos,'lambda')
    plotPriPos(sigma_pri,sigma_pos,'sigma')
    for(i in 1:ncol(yGP_pri))
      plotPriPos(yGP_pri[,i],yGP_pos[,i],paste0('yGP_',i),0.5*c(-1,1))

  })

  output$plotGP      <- renderPlot({
    if(is.null(out <- Inputs$outMonoExp))
      return(NULL)

    # Nb control points
    n = input$NnTest

    # Grid
    dx  = 1/(n+1)
    if(input$gridTypeTest == 'internal')
      xdat = seq(dx/2,1-dx/2,length.out = n)
    else
      xdat = seq(0.0,1.0,length.out = n)

    # Output values and reference
    xp <- (Inputs$x-min(Inputs$x)) / (max(Inputs$x)-min(Inputs$x))
    yref = 0.3*scale(out$fit$par$resid)
    yref = smooth.spline(xp,yref,df=input$smooth_df)$y
    ydat = yref[1+round(xdat*(length(xp)-1))]

    # Simulated GP
    nRun=100

    rho = input$rho_scaleTest
    if(rho == 0)
      rho = 1/n

    cond = RandomFields::RFsimulate(
      model = RMgauss(var=input$alpha_scaleTest*sd(ydat),
                      scale=rho),
      x=xp,
      given = list(x=xdat),
      data  = list(y=ydat),
      n=nRun,
    )

    par(mfrow=c(1,1),pty=pty,mar=mar,mgp=mgp,
        tcl=tcl,lwd=1.5*lwd, cex=1.5*cex)
    plot(xdat,ydat,pch=20,col=cols[7],
         xlim = c(0,1),
         ylim = 1.2*c(-1,1),
         xlab = 'depth (a.u.)',
         ylab = 'relative deviation')
    grid()
    lines(xp,yref,col=1,lty=2)
    matlines(xp,cond,col=col_tr[4])
    lines(xp,rowMeans(cond),col=cols[6])
    legend('topright',
           legend = c('ctrl. points','GP pred.',
                      'mean GP pred.','ref. data'),
           col = c(cols[7],col_tr2[4],cols[6],1),
           pch = c(20,NA,NA,NA),
           lty = c(NA,1,1,2),
           bty='n'
           )
    box()

  })
}


