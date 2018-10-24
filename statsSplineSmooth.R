# Assess the effect of parameter df in smooth.spline

is=0
dataDirs = c("DataWl","Data1","DataSynth")[3]
for (dataDir in dataDirs) {
  dataSets = list.dirs(path=dataDir,full.names = FALSE)[-1]
  for(dataSet in dataSets) {
    is=is+1
    tag = paste0('GP1_',dataDir,'_',dataSet); cat(tag,'---------','\n')

    # Get Data ####
    D = read.csv(paste0(dataDir,'/',dataSet,'/Courbe.csv'))
    x=D[,1]; y=D[,2]
    sel = x > 20 & x<=500
    x = x[sel]; y = y[sel]

    stat = c()
    i=0
    dfs = 2:20
    for(df in dfs) {
      ySpl = smooth.spline(x,y,df=df)$y
      resSpl = y-ySpl
      i=i+1; stat[i] = sd(resSpl)
    }

    if(is==1)
      plot(dfs,stat/stat[1],type='l',ylim=c(0,1),
           xlab = 'df', ylab= 'sd(residuals)')
    else
      lines(dfs,stat/stat[1],col=is)

    # l = 10
    # x1 = dfs[l:max(dfs)]
    # y1 = stat[l:max(dfs)]
    # abline(reg=lm(y1~x1),col=4)

  }
}
