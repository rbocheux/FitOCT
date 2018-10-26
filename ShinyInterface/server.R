# Libraries ####
libs =c('shiny','parallel','rstan','FitOCTLib','inlmisc')
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

x = 1:200
y = 100 + 100*exp(-x/50)
y = y + rnorm(length(y),0,0.2)
# y = y + rnorm(length(y),0,0.1*sqrt(y-100+5))

function(input, output, session) {

    output$plotNoise <- renderPlot({

      out = FitOCTLib::estimateNoise(x=x,y=y,df=input$set_df)
      ySmooth = out$ySmooth
      resid   = y - ySmooth
      uy      = out$uy

      par(mfrow=c(1,2),pty='s',
          mar=c(3,3,1.5,.2),mgp=c(2,.75,0),tcl=-0.5,
          lwd=2, cex=1)

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

      out = FitOCTLib::estimateNoise(x=x,y=y,df=input$set_df)
      ySmooth = out$ySmooth
      uy      = out$uy

      out = FitOCTLib::fitMonoExp(x=x,y=y,uy=uy)
      fit    = out$fit
      theta  = fit$par$theta
      resid  = fit$par$resid
      mod    = fit$par$m

      par(mfrow=c(1,2),pty='s',
          mar=c(3,3,1.5,.2),mgp=c(2,.75,0),tcl=-0.5,
          lwd=2, cex=1)

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

    output$summary <- renderPrint({
        summary(cars)
    })

    output$table <- DT::renderDataTable({
        DT::datatable(cars)
    })
}


