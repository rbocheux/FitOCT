# Generate synthetic datasets for FitOCT

x  = 20:500
a  = 1000
b  = 2000
l0 = 150
s  = 0.5       # Noise factor

# Mono-exponential ####
y0 = a + b*exp(-x/l0)
y = y0 + rnorm(length(x), 0,s*sqrt(y0-a+1))

plot(x,y)
lines(x,y0,col=2)

write.csv(cbind(x,y),file = 'DataSynth/monoExp/Courbe.csv',
          row.names = FALSE)


# Sinc-Modulated exponential ####
m  = 10 * sin(x/50) / x
y1 =  a + b*exp(-x/(l0*(1+m)))
y = y1 + rnorm(length(x), 0, s*sqrt(y0-a+1))

plot(x,y,pch=20)
lines(x,y1,col=4,lwd=2)
lines(x,y0,col=2,lwd=2)

write.csv(cbind(x,y),file = 'DataSynth/sincExp/Courbe.csv',
          row.names = FALSE)
write.csv(cbind(x,m),file = 'DataSynth/sincExp/Modulation.csv',
          row.names = FALSE)

# Sinc-Modulated exponential ####
m  = 10 * sin(x/25) / x
y1 =  a + b*exp(-x/(l0*(1+m)))
y = y1 + rnorm(length(x), 0, s*sqrt(y0-a+1))

plot(x,y,pch=20)
lines(x,y1,col=4,lwd=2)
lines(x,y0,col=2,lwd=2)

write.csv(cbind(x,y),file = 'DataSynth/sincExp1/Courbe.csv',
          row.names = FALSE)
write.csv(cbind(x,m),file = 'DataSynth/sincExp1/Modulation.csv',
          row.names = FALSE)

# Sinc-Modulated exponential ####
m  = 10 * sin(x/75) / x
y1 =  a + b*exp(-x/(l0*(1+m)))
y = y1 + rnorm(length(x), 0, s*sqrt(y0-a+1))

plot(x,y,pch=20)
lines(x,y1,col=4,lwd=2)
lines(x,y0,col=2,lwd=2)

write.csv(cbind(x,y),file = 'DataSynth/sincExp2/Courbe.csv',
          row.names = FALSE)
write.csv(cbind(x,m),file = 'DataSynth/sincExp2/Modulation.csv',
          row.names = FALSE)

# Sinc-Modulated exponential: Shot ####
m  = 1 * sin((x-250)/20) / (x-250+0.1) ; plot(x,m)
y1 =  a + b*exp(-x/(l0*(1+m)))
y = y1 + rnorm(length(x), 0, s*sqrt(y0-a+1))

plot(x,y,pch=20)
lines(x,y1,col=4,lwd=2)
lines(x,y0,col=2,lwd=2)

write.csv(cbind(x,y),file = 'DataSynth/sincExp3/Courbe.csv',
          row.names = FALSE)
write.csv(cbind(x,m),file = 'DataSynth/sincExp3/Modulation.csv',
          row.names = FALSE)
