# Libraries ####
libs =c('shiny','parallel','rstan','FitOCTLib',
        'inlmisc','shinycssloaders','DT')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

# Options ####
Sys.setlocale(category = "LC_NUMERIC", locale="C")
options(mc.cores = parallel::detectCores(),
        width = 20,
        warn  = 0)
rstan_options(auto_write = TRUE)
# set.seed(1234) # Initialise la graine du RNG

# Color schemes ####
cols    = inlmisc::GetTolColors(8)
# Transparents for spaghetti
col_tr  = inlmisc::GetTolColors(8,alpha=0.1)
# Darker for legends or fillings
col_tr2 = inlmisc::GetTolColors(8,alpha=0.4)

# Graphical parameters ####
pty = 's'
mar = c(3,3,1.5,.5)
mgp = c(2,.75,0)
tcl = -0.5
lwd = 2
cex = 1

# Global variables ####
Inputs = reactiveValues(
  x          = NULL,
  y          = NULL,
  outSmooth  = NULL,
  outMonoExp = NULL,
  outExpGP   = NULL
)

# MAIN ####
function(input, output, session) {

  observeEvent(
    input$dataFile,
    {
      dat = try(
        read.csv(input$dataFile[['datapath']]),
        silent = TRUE
      )
      if(class(dat) == 'try-error')
        return(NULL)

      # Store data and empty buffers
      Inputs$x          <<- dat[,1]
      Inputs$y          <<- dat[,2]
      Inputs$outSmooth  <<- NULL
      Inputs$outMonoExp <<- NULL
      Inputs$outExpGP   <<- NULL
    }
  )

  output$plotNoise   <- renderPlot({

    if(is.null(Inputs$x))
      return(NULL)

    x = Inputs$x
    y = Inputs$y

    out = FitOCTLib::estimateNoise(x=x,y=y,df=input$smooth_df)

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

    x = Inputs$x
    y = Inputs$y

    out = Inputs$outSmooth
    ySmooth = out$ySmooth
    uy      = out$uy

    out = FitOCTLib::fitMonoExp(x=x,y=y,uy=uy)

    Inputs$outMonoExp <<- out
    Inputs$outExpGP   <<- NULL

    fit    = out$fit
    theta  = fit$par$theta
    resid  = fit$par$resid
    mod    = fit$par$m

    par(mfrow=c(1,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    # Smooth Fit
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

    out = Inputs$outMonoExp
    fit = out$fit

    if(out$method == 'optim') {
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
      cat('!!! WARNING !!! \n')
    cat('br       :',signif(br,2),'\n')
    cat('CI95(br) :',paste0(signif(CI95,2),collapse='-'))

  })

  runExpGP <-function() {
    if(is.null(Inputs$outMonoExp))
      return(NULL)

    outm  = Inputs$outMonoExp
    fit   = outm$fit
    theta = fit$par$theta

    out = FitOCTLib::fitExpGP(x      = Inputs$x,
                              y      = Inputs$y,
                              uy     = Inputs$outSmooth[['uy']],
                              method = input$method,
                              theta0 = theta,
                              nb_warmup = 100, nb_iter = 200)

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

  output$plotExpGP   <- renderPlot({
    if (is.null(out <- doExpGP()))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    x = Inputs$x
    y = Inputs$y

    outS = Inputs$outSmooth
    ySmooth = outS$ySmooth
    uy      = outS$uy

    fit      = out$fit
    method   = out$method
    xGP      = out$xGP
    prior_PD = out$prior_PD
    lasso    = out$lasso

    nMC = 100

    par(mfrow=c(2,2),pty=pty,mar=mar,mgp=mgp,tcl=tcl,lwd=lwd, cex=cex)

    if(method == 'sample') {

      theta   = extract(fit,'theta')[[1]]
      yGP     = extract(fit,'yGP')[[1]]
      if(prior_PD == 0) {
        resid   = extract(fit,'resid')[[1]]
        mod     = extract(fit,'m')[[1]]
        dL      = extract(fit,'dL')[[1]]
        lp       = extract(fit,'lp__')[[1]]
        map = which.max(lp)
        y_map = mod[map,]
      }

      iMC = sample.int(nrow(theta),nMC)

      # Fit
      plot(x,y,pch=20,cex=0.5,col=cols[6],
           main='Data fit',
           xlab='depth (a.u.)',
           ylab='mean OCT signal (a.u.)')
      grid()
      if(prior_PD == 0) {
        if(nMC >0)
          for (i in 1:nMC)
            lines(x, mod[iMC[i],], col=col_tr[4])
        lines(x,theta[map,1]+theta[map,2]*exp(-x/theta[map,3]),col=cols[7])
      } else {
        if(nMC >0)
          for (i in 1:nMC)
            lines(x,theta[i,1]+theta[i,2]*exp(-x/theta[i,3]),col=col_tr[7])
      }

      legend('topright', bty='n',
             legend=c('data','expo. best fit','post. sample'),
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
        lines(x, y_map - ySmooth, col=cols[7])
        legend('topright', bty='n',
               legend=c('mean resid.','data 95% uncert.','best.fit - smooth'),
               pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
               col=c(cols[6],col_tr2[4],cols[7])
        )
        box()
      }

      # Local deviations
      if(prior_PD == 0) {
        plot(x, dL[map,], type = 'n',
             ylim = 0.5 * c(-1,1),
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
             ylim = 0.5*c(-1,1),
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
      lines(x,theta[1]+theta[2]*exp(-x/theta[3]),col=cols[7])
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
      lines(x, mod - ySmooth, col=cols[7])
      legend('topright', bty='n',
             legend=c('mean resid.','data 95% uncert.','best.fit - smooth'),
             pch=c(20,NA,NA),lty=c(-1,1,1),lwd=c(1,6,1),
             col=c(cols[6],col_tr2[4],cols[7])
      )
      box()

      # Local deviations
      plot(x, dL, type = 'l',
           ylim = 0.5*c(-1,1),
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

    if(prod(CI95-br) >= 0)
      cat('!!! WARNING !!! \n')
    cat('br       :',signif(br,2),'\n')
    cat('CI95(br) :',paste0(signif(CI95,2),collapse='-'))

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


}


