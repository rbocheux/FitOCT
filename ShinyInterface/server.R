# Libraries ####
libs =c('shiny','parallel','rstan','knitr',
        'inlmisc','shinycssloaders','DT',
        'RandomFields','stringr','devtools',
        'rmarkdown','rlist')
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

# Graphical parameters ####
gPars = list(
  cols    = inlmisc::GetColors(8),
  # Transparents for spaghetti
  col_tr  = inlmisc::GetColors(8,alpha=0.1),
  # Darker for legends or fillings
  col_tr2 = inlmisc::GetColors(8,alpha=0.4),
  pty = 's',
  mar = c(3,3,1.5,.5),
  mgp = c(2,.75,0),
  tcl = -0.5,
  lwd = 2,
  cex = 1
)

# Tables output parameters ###
DTopts = list(
  ordering    = FALSE,
  searching   = FALSE,
  paging      = FALSE,
  info        = FALSE,
  pageLength  = 16,
  deferRender = TRUE,
  scrollY     = FALSE,
  scrollX     = TRUE,
  stateSave   = FALSE
)

# Misc functions ####
summaryNoise    <- function(out){
  fit = out$fit
  if (out$method == 'optim') {

    pars = c('theta')
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

    DT::datatable(
      data = signif(sum,digits = 3),
      options = DTopts )

  } else {
    pars = c('theta')

    sum  = rstan::summary(fit,pars=pars,
                          use_cache=FALSE,
                          c(0.025,0.5,0.975))$summary[,-2]

    `%>%` <- DT::`%>%`
    DT::datatable(
      data = signif(sum,digits = 3),
      options = DTopts
    ) %>%
      DT::formatStyle(
        columns = "Rhat",
        color = DT::styleInterval(1.1, c("blue", "red"))
      ) %>%
      DT::formatRound(
        columns = "n_eff",
        digits = 0
      )
  }
}
summaryMonoExp  <- function(out){
  fit = out$fit
  if (out$method == 'optim') {

    pars = c('theta')
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
                  options = DTopts )

  } else {
    pars = c('theta','br')

    sum  = rstan::summary(fit,pars=pars,
                          use_cache=FALSE,
                          c(0.025,0.5,0.975))$summary[,-2]

    `%>%` <- DT::`%>%`
    DT::datatable(
      data = signif(sum,digits = 3),
      options = DTopts
    ) %>%
      DT::formatStyle(
        columns = "Rhat",
        color = DT::styleInterval(1.1, c("blue", "red"))
      ) %>%
      DT::formatRound(
        columns = "n_eff",
        digits = 0
      )
  }
}
summaryExpGP    <- function(out){
  fit = out$fit
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
                    stateSave   = FALSE
                  ) ) %>%
    DT::formatStyle(
      columns = "Rhat",
      color = DT::styleInterval(1.1, c("blue", "red"))) %>%
    DT::formatRound(columns = "n_eff", digits = 0)
  }
}
summaryPriExpGP <- function(out) {
  fit  = out$fit

  pars = c('theta','yGP','lambda','sigma')

  sum  = rstan::summary(fit,pars=pars,
                        use_cache=FALSE,
                        c(0.025,0.5,0.975))$summary[,-2]
  `%>%` <- DT::`%>%`
  DT::datatable(
    data = signif(sum,digits = 3),
    options = DTopts
  ) %>%
  DT::formatStyle(
    columns = "Rhat",
    color = DT::styleInterval(1.1, c("blue", "red"))
  ) %>%
  DT::formatRound(
    columns = "n_eff",
    digits = 0
  )

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
  session_dir = file.path(tempdir(),
                          stringr::str_sub(session$token, 1, 8))
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

  selX <- function(x,y,depthSel,subSample) {
    # Apply selectors to inputs

    xSel = which(x >= depthSel[1] &
                 x <= depthSel[2]  )
    x = x[xSel]; y = y[xSel]

    if(subSample != 1) {
      xSel = seq(1,length(x),by=subSample)
      x = x[xSel]; y = y[xSel]
    }
    return(list(x=x,y=y))
  }

  # Noise estimation ####
  output$plotNoise   <- renderPlot({
    if(is.null(Inputs$x))
      return(NULL)

    C = selX(Inputs$x,Inputs$y,input$depthSel,input$subSample)

    if(!is.finite(IQR(C$x)))
      return(NULL) # Protect smooth.spline from zero tol

    out = FitOCTLib::estimateNoise(
      x=C$x, y=C$y, df=input$smooth_df
    )

    Inputs$outSmooth  <<- out
    Inputs$outMonoExp <<- NULL
    Inputs$outExpGP   <<- NULL

    FitOCTLib::plotNoise(
      x=C$x, y=C$y, uy=out$uy, ySmooth=out$ySmooth, gPars=gPars
    )

  })

  # Mono-exponential fit ####
  output$plotMonoExp <- renderPlot({
    if(is.null(Inputs$outSmooth))
      return(NULL)

    C = selX(Inputs$x,Inputs$y,input$depthSel,input$subSample)

    out  = Inputs$outSmooth

    outm = FitOCTLib::fitMonoExp(
      x=C$x, y=C$y, uy=out$uy
    )

    Inputs$outMonoExp <<- outm
    Inputs$outExpGP   <<- NULL

    FitOCTLib::plotMonoExp(
      x=C$x, y=C$y, uy=out$uy, ySmooth=out$ySmooth,
      mod=outm$fit$par$m, resid=outm$fit$par$resid,
      gPars=gPars
    )

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

    }
    cat('\n')
    FitOCTLib::printBr(fit)

  })

  # Modulated exponential fit ####
  runExpGP <-function(prior_PD = 0) {
    if(is.null(Inputs$outMonoExp))
      return(NULL)

    outm  = Inputs$outMonoExp

    isolate({
      C = selX(Inputs$x,Inputs$y,input$depthSel,input$subSample)
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

    C = selX(Inputs$x,Inputs$y,input$depthSel,input$subSample)
    outS = Inputs$outSmooth
    FitOCTLib::plotExpGP(
      C$x, C$y, outS$uy, outS$ySmooth,
      out = out, modScale = input$modRange,
      gPars = gPars
    )
  })

  output$resExpGP    <- renderPrint({
    if (is.null(out <- doExpGP()))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    # Probability Interval for Birge's ratio
    FitOCTLib::printBr(out$fit)

  })

  output$summaryOut  <- DT::renderDataTable({
    if (is.null(out <- try(doExpGP(),silent=TRUE)))
      return(NULL)

    if(is.null(Inputs$outExpGP))
      return(NULL)

    summaryExpGP(out)

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

  # Prior pdf ####
  doPriExpGP <- eventReactive(
    input$runPriExpGP,
    {
      out = runExpGP(prior_PD = 1)
      Inputs$outPriExpGP <<- out
      return(out)
    }
  )

  output$plotPriExpGP   <- renderPlot({
    if (is.null(out <- doPriExpGP()))
      return(NULL)

    if(is.null(Inputs$outPriExpGP))
      return(NULL)

    C = selX(Inputs$x,Inputs$y,input$depthSel,input$subSample)
    outS = Inputs$outSmooth
    FitOCTLib::plotExpGP(
      C$x, C$y, outS$uy, outS$ySmooth,
      out = out, modScale = input$modRange,
      gPars = gPars
    )
  })

  output$summaryPriOut  <- DT::renderDataTable({
    if (is.null(out <- try(doPriExpGP(),silent=TRUE)))
      return(NULL)

    if(is.null(Inputs$outPriExpGP))
      return(NULL)

    summaryPriExpGP(out)

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

  output$priPostExpGP   <- renderPlot({
    if (is.null(fitGP <- doExpGP()) ||
        is.null(fitGP_pri <- doPriExpGP()))
      return(NULL)
    if(fitGP$method != 'sample')
      return(NULL)

    FitOCTLib::plotPriPostAll(fitGP_pri$fit,fitGP$fit,
                              gPars = gPars)

  })

  # GP-Design ####
  doApplyGPDesign <- observeEvent(
    input$applyGPDesign,
    {
      updateRadioButtons(session, 'gridType',
                         selected = input$gridTypeTest)
      updateNumericInput(session, 'Nn',
                         value   = input$NnTest)
      updateSliderInput(session, 'rho_scale',
                        value   = input$rho_scaleTest)
    }
  )

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
    C = selX(Inputs$x,Inputs$y,input$depthSel,input$subSample)
    xp <- (C$x-min(C$x)) / (max(C$x)-min(C$x))
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

    # Extract graphical params
    for (n in names(gPars))
      assign(n,list.extract(gPars,n))

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

  # Save ####
  listCtrlParams <- function() {
    list(
      depthSel    = input$depthSel,
      subSample   = input$subSample,
      smooth_df   = input$smooth_df,
      method      = input$method,
      nb_warmup   = input$nb_warmup,
      nb_sample   = input$nb_sample,
      modRange    = input$modRange,
      ru_theta    = input$ru_theta,
      lambda_rate = input$lambda_rate,
      gridType    = input$gridType,
      Nn          = input$Nn,
      rho_scale   = input$rho_scale
    )
  }

  output$report = downloadHandler(
    filename = "fitOCTReport.html",
    content = function(file) {

      parList = listCtrlParams()

      src <- normalizePath('reportTemplate.Rmd')
      owd <- setwd(tempdir())
      on.exit(expr = setwd(owd))
      file.copy(src, 'reportTemplate.Rmd', overwrite = TRUE)
      rmarkdown::render('reportTemplate.Rmd', output_file = file)

    }
  )

  output$params = downloadHandler(
    filename = "save_ctrlParams.yaml",
    content = function(file) {
      rlist::list.save(listCtrlParams(),file)
    }
  )

}


