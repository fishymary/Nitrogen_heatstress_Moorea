# -----------------------------------------------------------------------
# Nutrient pollution alters coral bleaching across the seascape
# Donovan et al. 
# 5 - model outputs from 4
# -----------------------------------------------------------------------

# This script manipulates JAGS outputs to produce figures and data summaries

# Initialization ----------------------------------------------------------
rm(list=ls())
library(rjags)
# library(dplyr)
library(plotrix)
library(wesanderson)
library(boot) # for inv.logit()

# data --------------------------------------------------------------------
poc.hpc.ap14 <- readRDS('hpc_out/poc_avgTotN_060519_q.Rdata')
acr.hpc.ap14 <- readRDS('hpc_out/acr_avgTotN_060519_q.Rdata')
moorea <- read.csv('data/moorea_withavgnuts.csv')

# pocillopora data inputs
temp <- moorea[moorea$Taxa=="Pocillopora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
# temp <- temp[sample(nrow(temp),100),]
poc.y <- temp$Percent_bleached
poc.y <- poc.y/100

# split the data into discrete and continuous components
poc.y.d <- ifelse(poc.y == 1 | poc.y == 0, poc.y, NA)
poc.y.discrete <- ifelse(is.na(poc.y.d), 0, 1)
poc.y.d <- poc.y.d[!is.na(poc.y.d)]
poc.n.discrete <- length(poc.y.d)

which.cont <- which(poc.y < 1 & poc.y > 0)
poc.y.c <- ifelse(poc.y < 1 & poc.y > 0, poc.y, NA)
poc.y.c <- poc.y.c[!is.na(poc.y.c)]
poc.n.cont <- length(poc.y.c)

# acropora data inputs
temp <- moorea[moorea$Taxa=="Acropora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
acr.y <- temp$Percent_bleached
acr.y <- acr.y/100

# split the data into discrete and continuous components
acr.y.d <- ifelse(acr.y == 1 | acr.y == 0, acr.y, NA)
acr.y.discrete <- ifelse(is.na(acr.y.d), 0, 1)
acr.y.d <- acr.y.d[!is.na(acr.y.d)]
acr.n.discrete <- length(acr.y.d)

which.cont <- which(acr.y < 1 & acr.y > 0)
acr.y.c <- ifelse(acr.y < 1 & acr.y > 0, acr.y, NA)
acr.y.c <- acr.y.c[!is.na(acr.y.c)]
acr.n.cont <- length(acr.y.c)

# convergence and diagnostics ---------------------------------------------
# convergence
param_check <-
  c(
    'beta_colony[1]',
    'beta_colony[2]',
    'beta_point[1]',
    'beta_point[2]',
    'beta_point[3]',
    'B0_habitat[1]',
    'B0_habitat[2]',
    'beta_colony_p[1]',
    'beta_colony_p[2]',
    'beta_point_p[1]',
    'beta_point_p[2]',
    'beta_point_p[3]',
    'B0_habitat_p[1]',
    'B0_habitat_p[2]'
  )

gel_check <- cbind(gelman.diag(mcmc.list(
  poc.hpc.ap14[[1]][, param_check], poc.hpc.ap14[[2]][, param_check], poc.hpc.ap14[[3]][, param_check]))[[1]],
  gelman.diag(mcmc.list(
    acr.hpc.ap14[[1]][, param_check], acr.hpc.ap14[[2]][, param_check], acr.hpc.ap14[[3]][, param_check]))[[1]])
gel_check

write.csv(gel_check, 'outputs/TableS1.csv', row.names = T)

# fit
# pocillopora prevalence
poc.ynew.d <- rep(NA,poc.n.discrete)
for(i in 1:poc.n.discrete) {
  poc.ynew.d[i] <- mean(c(poc.hpc.ap14[[1]][, paste0('y.new.d[', i, ']')], 
                          poc.hpc.ap14[[2]][, paste0('y.new.d[', i, ']')], 
                          poc.hpc.ap14[[3]][, paste0('y.new.d[', i, ']')]))
}
hist(poc.ynew.d)

poc.d.check <- ifelse(poc.ynew.d >= mean(poc.y.d),0,1)
sum(poc.d.check)/length(poc.ynew.d)

poc.ynew.p <- rep(NA,length(poc.y.c))
for(i in 1:length(poc.y.c)){
  poc.ynew.p[i] <- mean(c(poc.hpc.ap14[[1]][,paste0('y.new.p[',i,']')],
                          poc.hpc.ap14[[2]][,paste0('y.new.p[',i,']')],
                          poc.hpc.ap14[[3]][,paste0('y.new.p[',i,']')]))
}
hist(poc.ynew.p)

poc.p.check <- ifelse(poc.ynew.p >= mean(poc.y.c),0,1)
sum(poc.p.check)/length(poc.ynew.p)
plot(poc.y.c,poc.ynew.p)

# acropora prevalence
acr.ynew.d <- rep(NA,acr.n.discrete)
for(i in 1:acr.n.discrete){
  acr.ynew.d[i] <- mean(c(acr.hpc.ap14[[1]][,paste0('y.new.d[',i,']')],
                          acr.hpc.ap14[[2]][,paste0('y.new.d[',i,']')],
                          acr.hpc.ap14[[3]][,paste0('y.new.d[',i,']')]))
}
hist(acr.ynew.d)

acr.d.check <- ifelse(acr.ynew.d >= mean(acr.y.d),0,1)
sum(acr.d.check)/length(acr.ynew.d)

acr.ynew.p <- rep(NA,length(acr.y))
for(i in 1:length(acr.y)){
  acr.ynew.p[i] <- mean(c(acr.hpc.ap14[[1]][,paste0('y.new.p[',i,']')],
                          acr.hpc.ap14[[2]][,paste0('y.new.p[',i,']')],
                          acr.hpc.ap14[[3]][,paste0('y.new.p[',i,']')]))
}
hist(acr.ynew.p)

acr.p.check <- ifelse(acr.ynew.p >= mean(acr.y.discrete),0,1)
sum(acr.p.check)/length(acr.ynew.p)
plot(acr.y.discrete,acr.ynew.p)

# pocillopora severity
poc.ynew.s <- rep(NA,poc.n.cont)
for(i in 1:poc.n.cont){
  poc.ynew.s[i] <- mean(c(poc.hpc.ap14[[1]][,paste0('y.new.s[',i,']')],
                          poc.hpc.ap14[[2]][,paste0('y.new.s[',i,']')],
                          poc.hpc.ap14[[3]][,paste0('y.new.s[',i,']')]))
}
hist(poc.ynew.s)

poc.s.check <- ifelse(poc.ynew.s >= mean(poc.y.c),0,1)
sum(poc.s.check)/length(poc.ynew.s)

# acropora severity
acr.ynew.s <- rep(NA,acr.n.cont)
for(i in 1:acr.n.cont){
  acr.ynew.s[i] <- mean(c(acr.hpc.ap14[[1]][,paste0('y.new.s[',i,']')],
                          acr.hpc.ap14[[2]][,paste0('y.new.s[',i,']')],
                          acr.hpc.ap14[[3]][,paste0('y.new.s[',i,']')]))
}
hist(acr.ynew.s)

acr.s.check <- ifelse(acr.ynew.s >= mean(acr.y.c),0,1)
sum(acr.s.check)/length(acr.ynew.s)

# R2
varF <- sd(poc.y.d)^2
varE <- sd(poc.y.d-c(poc.ynew.d))^2
varF/(varF+varE)

varF <- sd(poc.y.c)^2
varE <- sd(poc.y.c-c(poc.ynew.s))^2
varF/(varF+varE)

varF <- sd(acr.y.d)^2
varE <- sd(acr.y.d-c(acr.ynew.d))^2
varF/(varF+varE)

varF <- sd(acr.y.discrete)^2
varE <- sd(acr.y.discrete-c(acr.ynew.p))^2
varF/(varF+varE)


varF <- sd(acr.y.c)^2
varE <- sd(acr.y.c-c(acr.ynew.s))^2
varF/(varF+varE)

# Figure 4 - parameter summary --------------------------------------------
poc.avgTotN.s <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
for(i in 1:2){ 
  poc.avgTotN.s[i,2] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.5)
  poc.avgTotN.s[i,3] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.975)
  poc.avgTotN.s[i,4] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.025)
  poc.avgTotN.s[i,5] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.90)
  poc.avgTotN.s[i,6] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.10)
}
for(i in 1:3){ 
  poc.avgTotN.s[i+2,2] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.50)
  poc.avgTotN.s[i+2,3] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.975)
  poc.avgTotN.s[i+2,4] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.025)
  poc.avgTotN.s[i+2,5] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.90)
  poc.avgTotN.s[i+2,6] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.10)
}
poc.avgTotN.s

poc.avgTotN.p <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat'),
                            mean=NA,up=NA,down=NA,up80=NA,down80=NA)
for(i in 1:2){ 
  poc.avgTotN.p[i,2] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.50)
  poc.avgTotN.p[i,3] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.975)
  poc.avgTotN.p[i,4] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.025)
  poc.avgTotN.p[i,5] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  poc.avgTotN.p[i,6] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){ 
  poc.avgTotN.p[i+2,2] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.50)
  poc.avgTotN.p[i+2,3] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.975)
  poc.avgTotN.p[i+2,4] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.025)
  poc.avgTotN.p[i+2,5] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  poc.avgTotN.p[i+2,6] <- quantile(c(poc.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
poc.avgTotN.p

acr.avgTotN.s <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
for(i in 1:2){ 
  acr.avgTotN.s[i,2] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.5)
  acr.avgTotN.s[i,3] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.975)
  acr.avgTotN.s[i,4] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.025)
  acr.avgTotN.s[i,5] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.90)
  acr.avgTotN.s[i,6] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony[',i,']')]),0.10)
}
for(i in 1:3){ 
  acr.avgTotN.s[i+2,2] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.50)
  acr.avgTotN.s[i+2,3] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.975)
  acr.avgTotN.s[i+2,4] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.025)
  acr.avgTotN.s[i+2,5] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.90)
  acr.avgTotN.s[i+2,6] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point[',i,']')]),0.10)
}
acr.avgTotN.s

acr.avgTotN.p <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat'),
                            mean=NA,up=NA,down=NA,up80=NA,down80=NA)
for(i in 1:2){ 
  acr.avgTotN.p[i,2] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.50)
  acr.avgTotN.p[i,3] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.975)
  acr.avgTotN.p[i,4] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.025)
  acr.avgTotN.p[i,5] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  acr.avgTotN.p[i,6] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.hpc.ap14[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){ 
  acr.avgTotN.p[i+2,2] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.50)
  acr.avgTotN.p[i+2,3] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.975)
  acr.avgTotN.p[i+2,4] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.025)
  acr.avgTotN.p[i+2,5] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  acr.avgTotN.p[i+2,6] <- quantile(c(acr.hpc.ap14[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.hpc.ap14[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
acr.avgTotN.p

png(file='outputs/Figure4.png',height=1500,width=3500,res=300)
par(mfrow=c(1,2),mar=c(3.5,1.5,2,1),mgp=c(2,1,0),oma=c(1,7,0,0))
# acropora
plot(acr.avgTotN.p$mean, c(5,4,3,2,1), xlim=c(min(acr.avgTotN.p$down),max(acr.avgTotN.p$up))
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.5,cex.lab=1.5,bty='l')
text(-7.8,5.5, expression("A)"~italic(Acropora)),cex=1.5)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.5)
abline(v=0,lty=2)
#prevelence
plotCI(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgTotN.p$up, li=acr.avgTotN.p$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgTotN.p$up80, li=acr.avgTotN.p$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# severity
plotCI(acr.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=acr.avgTotN.s$up, li=acr.avgTotN.s$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(acr.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=acr.avgTotN.s$up80, li=acr.avgTotN.s$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(acr.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9),pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# pocillopora
plot(poc.avgTotN.p$mean, c(5,4,3,2,1), xlim=c(min(poc.avgTotN.p$down),max(poc.avgTotN.p$up))
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n'
     ,cex.axis=1.5,cex.lab=1.5,bty='l')
text(-3.5,5.5, expression("B)"~italic(Pocillopora)),cex=1.5)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.3)
abline(v=0,lty=2)
#prevelence
plotCI(poc.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=poc.avgTotN.p$up, li=poc.avgTotN.p$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(poc.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=poc.avgTotN.p$up80, li=poc.avgTotN.p$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(poc.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# severity
plotCI(poc.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=poc.avgTotN.s$up, li=poc.avgTotN.s$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(poc.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=poc.avgTotN.s$up80, li=poc.avgTotN.s$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(poc.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9),pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

legend('topright',legend=c('Prevalence','Severity'),pch=19,col=c(wes_palette("Darjeeling1")[5],wes_palette("Darjeeling1")[1]),bty='n',cex=1.5)

mtext('Beta',outer=T,side=1,cex=1.5,line=-1.1)
mtext('Colony size',outer=T,side=2,cex=1.5,line=-0.9,at=0.81,las=2)
mtext('Depth',outer=T,side=2,cex=1.5,line=-0.9,at=0.66,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.515,las=2)
mtext('Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.37,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.25,las=2)
mtext('x Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.19,las=2)
dev.off()

# Figure 5 - interaction --------------------------------------------------

ind.minmax <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.minmax$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interA.max.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interA.max.1[1]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interB.max.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interB.max.1[1]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interC.max.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interC.max.1[1]'])))),
                     inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interA.max.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interA.max.1[2]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interB.max.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interB.max.1[2]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interC.max.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interC.max.1[2]'])))))

ind.minmax$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interA.max.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interA.max.1[1]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interB.max.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interB.max.1[1]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interC.max.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interC.max.1[1]']), 0.975))),
                   inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interA.max.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interA.max.1[2]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interB.max.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interB.max.1[2]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interC.max.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interC.max.1[2]']), 0.975))))

ind.minmax$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.max.1[1]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.max.1[1]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.max.1[1]']), 0.025))),
                     inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.max.1[2]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.max.1[2]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.max.1[2]']), 0.025))))

ind.minmax$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.max.1[1]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.max.1[1]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.max.1[1]']), 0.9))),
                     inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.max.1[2]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.max.1[2]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.max.1[2]']), 0.9))))

ind.minmax$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interA.max.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interA.max.1[1]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interB.max.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interB.max.1[1]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interC.max.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interC.max.1[1]']), 0.1))),
                       inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.max.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interA.max.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interA.max.1[2]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.max.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interB.max.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interB.max.1[2]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.max.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interC.max.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interC.max.1[2]']), 0.1))))

ind.minmax2 <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.minmax2$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interA.max.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interA.max.2[1]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interB.max.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interB.max.2[1]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interC.max.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interC.max.2[1]'])))),
                      inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interA.max.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interA.max.2[2]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interB.max.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interB.max.2[2]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interC.max.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interC.max.2[2]'])))))

ind.minmax2$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interA.max.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interA.max.2[1]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interB.max.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interB.max.2[1]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interC.max.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interC.max.2[1]']),0.975))),
                    inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interA.max.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interA.max.2[2]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interB.max.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interB.max.2[2]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interC.max.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interC.max.2[2]']),0.975))))

ind.minmax2$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.max.2[1]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.max.2[1]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.max.2[1]']),0.025))),
                      inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.max.2[2]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.max.2[2]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.max.2[2]']),0.025))))

ind.minmax2$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.max.2[1]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.max.2[1]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.max.2[1]']),0.9))),
                      inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.max.2[2]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.max.2[2]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.max.2[2]']),0.9))))

ind.minmax2$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interA.max.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interA.max.2[1]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interB.max.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interB.max.2[1]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interC.max.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interC.max.2[1]']),0.1))),
                        inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.max.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interA.max.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interA.max.2[2]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interB.max.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interB.max.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interB.max.2[2]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interC.max.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interC.max.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interC.max.2[2]']),0.1))))

png(file='outputs/Figure5.png',height=2000,width=2700,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0))
xseq <- c(1,2,3,5,6,7)
plot(xseq,ind.minmax$mean,ylim=c(0,1),xlim=c(0.5,7.5),type='n',ylab='Predicted Bleaching Severity',xlab='Nutrients',xaxt='n',cex.axis=1.3,cex.lab=1.4)
axis(1,at=xseq,labels=c('Low','Medium','High','Low','Medium','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
plotCI(xseq-0.1,ind.minmax$mean,ui=ind.minmax$up,li=ind.minmax$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,ind.minmax$mean,ui=ind.minmax$up80,li=ind.minmax$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,ind.minmax$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

plotCI(xseq+0.1,ind.minmax2$mean,ui=ind.minmax2$up,li=ind.minmax2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,ind.minmax2$mean,ui=ind.minmax2$up80,li=ind.minmax2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,ind.minmax2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

mtext('Low Heat Stress',side=3,at=2,cex=1.4)
mtext('High Heat Stress',side=3,at=6,cex=1.4)
legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)
dev.off()


ind.q <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.q$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interA.q.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interA.q.1[1]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interB.q.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interB.q.1[1]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interC.q.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interC.q.1[1]'])))),
                     inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interA.q.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interA.q.1[2]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interB.q.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interB.q.1[2]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interC.q.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interC.q.1[2]'])))))

ind.q$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interA.q.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interA.q.1[1]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interB.q.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interB.q.1[1]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interC.q.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interC.q.1[1]']), 0.975))),
                   inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interA.q.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interA.q.1[2]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interB.q.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interB.q.1[2]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interC.q.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interC.q.1[2]']), 0.975))))

ind.q$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.q.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.q.1[1]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.q.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.q.1[1]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.q.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.q.1[1]']), 0.025))),
                     inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.q.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.q.1[2]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.q.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.q.1[2]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.q.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.q.1[2]']), 0.025))))

ind.q$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.q.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.q.1[1]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.q.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.q.1[1]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.q.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.q.1[1]']), 0.9))),
                     inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA.q.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA.q.1[2]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB.q.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB.q.1[2]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC.q.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC.q.1[2]']), 0.9))))

ind.q$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interA.q.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interA.q.1[1]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interB.q.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interB.q.1[1]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interC.q.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interC.q.1[1]']), 0.1))),
                       inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA.q.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interA.q.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interA.q.1[2]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interB.q.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interB.q.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interB.q.1[2]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interC.q.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interC.q.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interC.q.1[2]']), 0.1))))

ind.q2 <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.q2$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interA.q.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interA.q.2[1]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interB.q.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interB.q.2[1]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interC.q.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interC.q.2[1]'])))),
                      inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interA.q.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interA.q.2[2]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interB.q.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interB.q.2[2]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interC.q.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interC.q.2[2]'])))))

ind.q2$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interA.q.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interA.q.2[1]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interB.q.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interB.q.2[1]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interC.q.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interC.q.2[1]']),0.975))),
                    inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interA.q.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interA.q.2[2]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interB.q.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interB.q.2[2]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interC.q.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interC.q.2[2]']),0.975))))

ind.q2$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.q.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.q.2[1]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.q.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.q.2[1]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.q.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.q.2[1]']),0.025))),
                      inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.q.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.q.2[2]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.q.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.q.2[2]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.q.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.q.2[2]']),0.025))))

ind.q2$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.q.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.q.2[1]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.q.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.q.2[1]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.q.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.q.2[1]']),0.9))),
                      inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interA.q.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interA.q.2[2]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interB.q.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interB.q.2[2]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interC.q.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interC.q.2[2]']),0.9))))

ind.q2$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interA.q.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interA.q.2[1]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interB.q.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interB.q.2[1]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interC.q.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interC.q.2[1]']),0.1))),
                        inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA.q.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interA.q.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interA.q.2[2]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interB.q.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interB.q.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interB.q.2[2]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interC.q.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interC.q.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interC.q.2[2]']),0.1))))


png(file='outputs/Figure5_z.png',height=1800,width=3400,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0),oma=c(0,0,2,1))

# min
xseq <- c(1,2,3)
go <- ind.minmax[ind.minmax$temp=='low',]
go2 <- ind.minmax2[ind.minmax2$temp=='low',]
plot(xseq,go$mean,ylim=c(0,1),xlim=c(0.5,15.5),type='n',ylab='Predicted Bleaching Severity',xlab='Nutrients',xaxt='n',cex.axis=1.3,cex.lab=1.4)
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# low quantile
xseq <- c(5,6,7)
go <- ind.q[ind.q$temp=='low',]
go2 <- ind.q2[ind.q2$temp=='low',]
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=8,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# up quantile
xseq <- c(9,10,11)
go <- ind.q[ind.q$temp=='high',]
go2 <- ind.q2[ind.q2$temp=='high',]
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=12,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# max temp
xseq <- c(13,14,15)
go <- ind.minmax[ind.minmax$temp=='high',]
go2 <- ind.minmax2[ind.minmax2$temp=='high',]
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
# abline(v=12,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

mtext('Minimum',side=3,at=2,cex=1.4)
mtext('Low Quartile',side=3,at=6,cex=1.4)
mtext('Up Quartile',side=3,at=10,cex=1.4)
mtext('Maximum',side=3,at=14,cex=1.4)
mtext('Heat Stress',side=3,outer=T,cex=1.4)

legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)

dev.off()

###### prevalence
ind.minmax_p <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.minmax_p$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[1]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[1]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[1]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[1]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[1]'])))),
                     inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[2]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[2]'])),
                                 mean(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[2]'], 
                                        poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[2]'], 
                                        poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[2]'])))))

ind.minmax_p$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[1]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[1]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[1]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[1]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[1]']), 0.975))),
                   inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[2]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[2]']), 0.975),
                               quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[2]'], 
                                          poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[2]'], 
                                          poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[2]']), 0.975))))

ind.minmax_p$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[1]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[1]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[1]']), 0.025))),
                     inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[2]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[2]']), 0.025),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[2]']), 0.025))))

ind.minmax_p$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[1]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[1]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[1]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[1]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[1]']), 0.9))),
                     inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[2]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[2]']), 0.9),
                                 quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[2]'], 
                                            poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[2]'], 
                                            poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[2]']), 0.9))))

ind.minmax_p$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[1]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[1]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[1]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[1]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[1]']), 0.1))),
                       inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.max.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interA_p.max.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interA_p.max.1[2]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.max.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interB_p.max.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interB_p.max.1[2]']), 0.1),
                                   quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.max.1[2]'], 
                                              poc.hpc.ap14[[2]][, 'pred.interC_p.max.1[2]'], 
                                              poc.hpc.ap14[[3]][, 'pred.interC_p.max.1[2]']), 0.1))))

ind.minmax2_p <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.minmax2_p$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interA_p.max.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interA_p.max.2[1]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interB_p.max.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interB_p.max.2[1]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[1]'],
                                         poc.hpc.ap14[[2]][,'pred.interC_p.max.2[1]'],
                                         poc.hpc.ap14[[3]][,'pred.interC_p.max.2[1]'])))),
                      inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interA_p.max.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interA_p.max.2[2]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interB_p.max.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interB_p.max.2[2]'])),
                                  mean(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[2]'],
                                         poc.hpc.ap14[[2]][,'pred.interC_p.max.2[2]'],
                                         poc.hpc.ap14[[3]][,'pred.interC_p.max.2[2]'])))))

ind.minmax2_p$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interA_p.max.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interA_p.max.2[1]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interB_p.max.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interB_p.max.2[1]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[1]'],
                                           poc.hpc.ap14[[2]][,'pred.interC_p.max.2[1]'],
                                           poc.hpc.ap14[[3]][,'pred.interC_p.max.2[1]']),0.975))),
                    inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interA_p.max.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interA_p.max.2[2]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interB_p.max.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interB_p.max.2[2]']),0.975),
                                quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[2]'],
                                           poc.hpc.ap14[[2]][,'pred.interC_p.max.2[2]'],
                                           poc.hpc.ap14[[3]][,'pred.interC_p.max.2[2]']),0.975))))

ind.minmax2_p$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interA_p.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interA_p.max.2[1]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interB_p.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interB_p.max.2[1]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interC_p.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interC_p.max.2[1]']),0.025))),
                      inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interA_p.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interA_p.max.2[2]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interB_p.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interB_p.max.2[2]']),0.025),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interC_p.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interC_p.max.2[2]']),0.025))))

ind.minmax2_p$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interA_p.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interA_p.max.2[1]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interB_p.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interB_p.max.2[1]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[1]'],
                                             poc.hpc.ap14[[2]][,'pred.interC_p.max.2[1]'],
                                             poc.hpc.ap14[[3]][,'pred.interC_p.max.2[1]']),0.9))),
                      inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interA_p.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interA_p.max.2[2]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interB_p.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interB_p.max.2[2]']),0.9),
                                  quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[2]'],
                                             poc.hpc.ap14[[2]][,'pred.interC_p.max.2[2]'],
                                             poc.hpc.ap14[[3]][,'pred.interC_p.max.2[2]']),0.9))))

ind.minmax2_p$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interA_p.max.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interA_p.max.2[1]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interB_p.max.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interB_p.max.2[1]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[1]'],
                                               poc.hpc.ap14[[2]][,'pred.interC_p.max.2[1]'],
                                               poc.hpc.ap14[[3]][,'pred.interC_p.max.2[1]']),0.1))),
                        inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.max.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interA_p.max.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interA_p.max.2[2]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.max.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interB_p.max.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interB_p.max.2[2]']),0.1),
                                    quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.max.2[2]'],
                                               poc.hpc.ap14[[2]][,'pred.interC_p.max.2[2]'],
                                               poc.hpc.ap14[[3]][,'pred.interC_p.max.2[2]']),0.1))))


ind.q_p <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.q_p$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[1]'], 
                                   poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[1]'], 
                                   poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[1]'])),
                            mean(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[1]'], 
                                   poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[1]'], 
                                   poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[1]'])),
                            mean(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[1]'], 
                                   poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[1]'], 
                                   poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[1]'])))),
                inv.logit(c(mean(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[2]'], 
                                   poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[2]'], 
                                   poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[2]'])),
                            mean(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[2]'], 
                                   poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[2]'], 
                                   poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[2]'])),
                            mean(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[2]'], 
                                   poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[2]'], 
                                   poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[2]'])))))

ind.q_p$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[1]'], 
                                     poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[1]'], 
                                     poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[1]']), 0.975),
                          quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[1]'], 
                                     poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[1]'], 
                                     poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[1]']), 0.975),
                          quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[1]'], 
                                     poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[1]'], 
                                     poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[1]']), 0.975))),
              inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[2]'], 
                                     poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[2]'], 
                                     poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[2]']), 0.975),
                          quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[2]'], 
                                     poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[2]'], 
                                     poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[2]']), 0.975),
                          quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[2]'], 
                                     poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[2]'], 
                                     poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[2]']), 0.975))))

ind.q_p$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[1]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[1]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[1]']), 0.025),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[1]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[1]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[1]']), 0.025),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[1]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[1]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[1]']), 0.025))),
                inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[2]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[2]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[2]']), 0.025),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[2]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[2]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[2]']), 0.025),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[2]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[2]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[2]']), 0.025))))

ind.q_p$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[1]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[1]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[1]']), 0.9),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[1]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[1]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[1]']), 0.9),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[1]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[1]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[1]']), 0.9))),
                inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[2]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[2]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[2]']), 0.9),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[2]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[2]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[2]']), 0.9),
                            quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[2]'], 
                                       poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[2]'], 
                                       poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[2]']), 0.9))))

ind.q_p$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[1]'], 
                                         poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[1]'], 
                                         poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[1]']), 0.1),
                              quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[1]'], 
                                         poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[1]'], 
                                         poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[1]']), 0.1),
                              quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[1]'], 
                                         poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[1]'], 
                                         poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[1]']), 0.1))),
                  inv.logit(c(quantile(c(poc.hpc.ap14[[1]][, 'pred.interA_p.q.1[2]'], 
                                         poc.hpc.ap14[[2]][, 'pred.interA_p.q.1[2]'], 
                                         poc.hpc.ap14[[3]][, 'pred.interA_p.q.1[2]']), 0.1),
                              quantile(c(poc.hpc.ap14[[1]][, 'pred.interB_p.q.1[2]'], 
                                         poc.hpc.ap14[[2]][, 'pred.interB_p.q.1[2]'], 
                                         poc.hpc.ap14[[3]][, 'pred.interB_p.q.1[2]']), 0.1),
                              quantile(c(poc.hpc.ap14[[1]][, 'pred.interC_p.q.1[2]'], 
                                         poc.hpc.ap14[[2]][, 'pred.interC_p.q.1[2]'], 
                                         poc.hpc.ap14[[3]][, 'pred.interC_p.q.1[2]']), 0.1))))

ind.q2_p <-
  data.frame(
    temp = c(rep('low', 3), rep('high', 3)),
    nuts = c('low', 'medium', 'high', 'low', 'medium', 'high'),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )

ind.q2_p$mean <- c(inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[1]'],
                                    poc.hpc.ap14[[2]][,'pred.interA_p.q.2[1]'],
                                    poc.hpc.ap14[[3]][,'pred.interA_p.q.2[1]'])),
                             mean(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[1]'],
                                    poc.hpc.ap14[[2]][,'pred.interB_p.q.2[1]'],
                                    poc.hpc.ap14[[3]][,'pred.interB_p.q.2[1]'])),
                             mean(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[1]'],
                                    poc.hpc.ap14[[2]][,'pred.interC_p.q.2[1]'],
                                    poc.hpc.ap14[[3]][,'pred.interC_p.q.2[1]'])))),
                 inv.logit(c(mean(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[2]'],
                                    poc.hpc.ap14[[2]][,'pred.interA_p.q.2[2]'],
                                    poc.hpc.ap14[[3]][,'pred.interA_p.q.2[2]'])),
                             mean(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[2]'],
                                    poc.hpc.ap14[[2]][,'pred.interB_p.q.2[2]'],
                                    poc.hpc.ap14[[3]][,'pred.interB_p.q.2[2]'])),
                             mean(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[2]'],
                                    poc.hpc.ap14[[2]][,'pred.interC_p.q.2[2]'],
                                    poc.hpc.ap14[[3]][,'pred.interC_p.q.2[2]'])))))

ind.q2_p$up <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[1]'],
                                      poc.hpc.ap14[[2]][,'pred.interA_p.q.2[1]'],
                                      poc.hpc.ap14[[3]][,'pred.interA_p.q.2[1]']),0.975),
                           quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[1]'],
                                      poc.hpc.ap14[[2]][,'pred.interB_p.q.2[1]'],
                                      poc.hpc.ap14[[3]][,'pred.interB_p.q.2[1]']),0.975),
                           quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[1]'],
                                      poc.hpc.ap14[[2]][,'pred.interC_p.q.2[1]'],
                                      poc.hpc.ap14[[3]][,'pred.interC_p.q.2[1]']),0.975))),
               inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[2]'],
                                      poc.hpc.ap14[[2]][,'pred.interA_p.q.2[2]'],
                                      poc.hpc.ap14[[3]][,'pred.interA_p.q.2[2]']),0.975),
                           quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[2]'],
                                      poc.hpc.ap14[[2]][,'pred.interB_p.q.2[2]'],
                                      poc.hpc.ap14[[3]][,'pred.interB_p.q.2[2]']),0.975),
                           quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[2]'],
                                      poc.hpc.ap14[[2]][,'pred.interC_p.q.2[2]'],
                                      poc.hpc.ap14[[3]][,'pred.interC_p.q.2[2]']),0.975))))

ind.q2_p$down <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[1]'],
                                        poc.hpc.ap14[[2]][,'pred.interA_p.q.2[1]'],
                                        poc.hpc.ap14[[3]][,'pred.interA_p.q.2[1]']),0.025),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[1]'],
                                        poc.hpc.ap14[[2]][,'pred.interB_p.q.2[1]'],
                                        poc.hpc.ap14[[3]][,'pred.interB_p.q.2[1]']),0.025),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[1]'],
                                        poc.hpc.ap14[[2]][,'pred.interC_p.q.2[1]'],
                                        poc.hpc.ap14[[3]][,'pred.interC_p.q.2[1]']),0.025))),
                 inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[2]'],
                                        poc.hpc.ap14[[2]][,'pred.interA_p.q.2[2]'],
                                        poc.hpc.ap14[[3]][,'pred.interA_p.q.2[2]']),0.025),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[2]'],
                                        poc.hpc.ap14[[2]][,'pred.interB_p.q.2[2]'],
                                        poc.hpc.ap14[[3]][,'pred.interB_p.q.2[2]']),0.025),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[2]'],
                                        poc.hpc.ap14[[2]][,'pred.interC_p.q.2[2]'],
                                        poc.hpc.ap14[[3]][,'pred.interC_p.q.2[2]']),0.025))))

ind.q2_p$up80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[1]'],
                                        poc.hpc.ap14[[2]][,'pred.interA_p.q.2[1]'],
                                        poc.hpc.ap14[[3]][,'pred.interA_p.q.2[1]']),0.9),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[1]'],
                                        poc.hpc.ap14[[2]][,'pred.interB_p.q.2[1]'],
                                        poc.hpc.ap14[[3]][,'pred.interB_p.q.2[1]']),0.9),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[1]'],
                                        poc.hpc.ap14[[2]][,'pred.interC_p.q.2[1]'],
                                        poc.hpc.ap14[[3]][,'pred.interC_p.q.2[1]']),0.9))),
                 inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[2]'],
                                        poc.hpc.ap14[[2]][,'pred.interA_p.q.2[2]'],
                                        poc.hpc.ap14[[3]][,'pred.interA_p.q.2[2]']),0.9),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[2]'],
                                        poc.hpc.ap14[[2]][,'pred.interB_p.q.2[2]'],
                                        poc.hpc.ap14[[3]][,'pred.interB_p.q.2[2]']),0.9),
                             quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[2]'],
                                        poc.hpc.ap14[[2]][,'pred.interC_p.q.2[2]'],
                                        poc.hpc.ap14[[3]][,'pred.interC_p.q.2[2]']),0.9))))

ind.q2_p$down80 <- c(inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[1]'],
                                          poc.hpc.ap14[[2]][,'pred.interA_p.q.2[1]'],
                                          poc.hpc.ap14[[3]][,'pred.interA_p.q.2[1]']),0.1),
                               quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[1]'],
                                          poc.hpc.ap14[[2]][,'pred.interB_p.q.2[1]'],
                                          poc.hpc.ap14[[3]][,'pred.interB_p.q.2[1]']),0.1),
                               quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[1]'],
                                          poc.hpc.ap14[[2]][,'pred.interC_p.q.2[1]'],
                                          poc.hpc.ap14[[3]][,'pred.interC_p.q.2[1]']),0.1))),
                   inv.logit(c(quantile(c(poc.hpc.ap14[[1]][,'pred.interA_p.q.2[2]'],
                                          poc.hpc.ap14[[2]][,'pred.interA_p.q.2[2]'],
                                          poc.hpc.ap14[[3]][,'pred.interA_p.q.2[2]']),0.1),
                               quantile(c(poc.hpc.ap14[[1]][,'pred.interB_p.q.2[2]'],
                                          poc.hpc.ap14[[2]][,'pred.interB_p.q.2[2]'],
                                          poc.hpc.ap14[[3]][,'pred.interB_p.q.2[2]']),0.1),
                               quantile(c(poc.hpc.ap14[[1]][,'pred.interC_p.q.2[2]'],
                                          poc.hpc.ap14[[2]][,'pred.interC_p.q.2[2]'],
                                          poc.hpc.ap14[[3]][,'pred.interC_p.q.2[2]']),0.1))))


png(file='outputs/Figure5_poc_prevalance.png',height=1800,width=3400,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0),oma=c(0,0,2,1))

# min
xseq <- c(1,2,3)
go <- ind.minmax_p[ind.minmax_p$temp=='low',]
go2 <- ind.minmax2_p[ind.minmax2_p$temp=='low',]
plot(xseq,go$mean,ylim=c(0,1),xlim=c(0.5,15.5),type='n',ylab='Predicted Bleaching Prevalence',xlab='Nutrients',xaxt='n',cex.axis=1.3,cex.lab=1.4)
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# low quantile
xseq <- c(5,6,7)
go <- ind.q_p[ind.q_p$temp=='low',]
go2 <- ind.q2_p[ind.q2_p$temp=='low',]
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=8,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# up quantile
xseq <- c(9,10,11)
go <- ind.q_p[ind.q_p$temp=='high',]
go2 <- ind.q2_p[ind.q2_p$temp=='high',]
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
abline(v=12,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# max temp
xseq <- c(13,14,15)
go <- ind.minmax_p[ind.minmax_p$temp=='high',]
go2 <- ind.minmax2_p[ind.minmax2_p$temp=='high',]
axis(1,at=xseq,labels=c('Low','Med','High'),cex.axis=1.3,cex.lab=1.4)
# abline(v=12,lty=2)
plotCI(xseq-0.1,go$mean,ui=go$up,li=go$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,go$mean,ui=go$up80,li=go$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,go$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
plotCI(xseq+0.1,go2$mean,ui=go2$up,li=go2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,go2$mean,ui=go2$up80,li=go2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,go2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

mtext('Minimum',side=3,at=2,cex=1.4)
mtext('Low Quartile',side=3,at=6,cex=1.4)
mtext('Up Quartile',side=3,at=10,cex=1.4)
mtext('Maximum',side=3,at=14,cex=1.4)
mtext('Heat Stress',side=3,outer=T,cex=1.4)

legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)

dev.off()


# continuous plot ---------------------------------------------------------

poc.hpc.fine <- readRDS('hpc_out/poc_avgTotN_fine_070119.Rdata')

# "a","B0_p","beta_colony_p","B0","beta_colony"

param_go <- 'a'
a.post <- mcmc.list(
  poc.hpc.fine[[1]][, param_go], poc.hpc.fine[[2]][, param_go], poc.hpc.fine[[3]][, param_go])

dim(a.post[[1]])

param_go <- grep("B0_p",colnames(poc.hpc.fine[[1]]))
B0_p.post <- mcmc.list(
  poc.hpc.fine[[1]][, param_go], poc.hpc.fine[[2]][, param_go], poc.hpc.fine[[3]][, param_go])

param_go <- grep("B0",colnames(poc.hpc.fine[[1]]))
B0.post <- mcmc.list(
  poc.hpc.fine[[1]][, param_go], poc.hpc.fine[[2]][, param_go], poc.hpc.fine[[3]][, param_go])

param_go <- grep("B0",colnames(poc.hpc.fine[[1]]))
B0.post <- mcmc.list(
  poc.hpc.fine[[1]][, param_go], poc.hpc.fine[[2]][, param_go], poc.hpc.fine[[3]][, param_go])


library(ggmcmc)
poc.hpc.fine.gs <- ggs(poc.hpc.fine)
a.post <- subset(poc.hpc.fine.gs, Parameter=='a')$value
n.stored <- length(a.post)
n <- length(poc.y)
P.discrete <- array(dim=c(n.stored, n))
for (i in 1:n){
  P.discrete[, i] <- inv.logit(a.post)
}
pdd <- reshape2::melt(P.discrete, varnames = c("iteration", "site"), 
            value.name = "Pr.discrete")

temp <- moorea[moorea$Taxa=="Pocillopora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
x <- temp$Depth_m
pdd$x <- x[pdd$site]
ggplot(pdd) +
  geom_line(aes(x=x, y=Pr.discrete, group=iteration), 
            alpha=0.05, color="grey")
