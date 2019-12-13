# -----------------------------------------------------------------------
# Nitrogen pollution interacts with heat stress to increase coral bleaching across the seascape
# Donovan et al. 
# 5 - model outputs from 4
# -----------------------------------------------------------------------

# This script manipulates JAGS outputs to produce figures and data summaries

# Initialization ----------------------------------------------------------
rm(list=ls())
library(rjags)
library(plotrix)
library(wesanderson)
library(boot) # for inv.logit()
library(dplyr)

# data --------------------------------------------------------------------
poc.binom <- readRDS('data/hpc_out/binom_poc.Rdata')
poc.binom.sum <- summary(poc.binom)

poc.beta <- readRDS('data/hpc_out/beta_poc.Rdata')
poc.beta.sum <- summary(poc.beta)

acr.binom <- readRDS('data/hpc_out/binom_acr.Rdata')
acr.binom.sum <- summary(acr.binom)

acr.beta <- readRDS('data/hpc_out/beta_acr.Rdata')
acr.beta.sum <- summary(acr.beta)

moorea <- read.csv('data/moorea_2016_bleaching.csv')

# convergence and diagnostics ---------------------------------------------

param_check <-
  c('beta_colony_p[1]',
    'beta_colony_p[2]',
    'beta_point_p[1]',
    'beta_point_p[2]',
    'beta_point_p[3]',
    'B0_habitat_p[1]',
    'B0_habitat_p[2]',
    'mu_coast_p[1]',
    'mu_coast_p[2]',
    'mu_coast_p[3]'
  )
mod_check <- 
  c('pval.mean',
    'pval.sd',
    'R2',
    'pval.cv',
    'pval.fit',
    'pval.min',
    'pval.max'
  )

# convergence check
gel_check <- cbind(
  gelman.diag(mcmc.list(
    poc.binom[[1]][, param_check], poc.binom[[2]][, param_check], poc.binom[[3]][, param_check]))[[1]][,1],
  gelman.diag(mcmc.list(
    poc.beta[[1]][, param_check], poc.beta[[2]][, param_check], poc.beta[[3]][, param_check]))[[1]][,1],
  gelman.diag(mcmc.list(
    acr.binom[[1]][, param_check], acr.binom[[2]][, param_check], acr.binom[[3]][, param_check]))[[1]][,1],
  gelman.diag(mcmc.list(
    acr.beta[[1]][, param_check], acr.beta[[2]][, param_check], acr.beta[[3]][, param_check]))[[1]][,1]
)
gel_check
temp <- data.frame(gel_check,var=rownames(gel_check))
colnames(temp) <- c('poc.binom','poc.beta','acr.binom','acr.beta','var')
write.csv(temp, 'outputs/TableS1.csv', row.names = T)

# model fit
round(poc.binom.sum$statistics[mod_check,'Mean'],2)
round(poc.beta.sum$statistics[mod_check,'Mean'],2)
round(acr.binom.sum$statistics[mod_check,'Mean'],2)
round(acr.beta.sum$statistics[mod_check,'Mean'],2)

# poc severity - observed versus predicted
grepgo <- grep('y.new',colnames(poc.beta[[1]]))
poc.beta.ynew <- poc.beta.sum$statistics[grepgo,'Mean']
poc.beta.ynew <- poc.beta.sum$quantiles[grepgo,'50%']

temp <- moorea[moorea$Taxa=="Pocillopora",]; temp <- temp[!is.na(temp$Percent_bleached),]; temp <- temp[!is.na(temp$avgTotN),]
temp <- temp[temp$Percent_bleached > 0,]
y <- temp$Percent_bleached; y <- y/100; y <- ifelse(y == 1, 0.99, y)

plot(y,poc.beta.ynew,xlim=c(0,1),ylim=c(0,1)); abline(a=0,b=1)

d.check <- ifelse(poc.beta.ynew >= mean(y),0,1)
sum(d.check)/length(d.check)

hist(y,breaks=20,freq=F)
lines(density(poc.beta.ynew),col='blue')

# prev
temp <- moorea[moorea$Taxa=="Pocillopora",]; temp <- temp[!is.na(temp$Percent_bleached),]; temp <- temp[!is.na(temp$avgTotN),]
y <- temp$Percent_bleached; y <- y/100; y_p <- ifelse(y > 0, 1, 0)
grepgo <- grep('y.new',colnames(poc.binom[[1]]))
poc.binom.ynew <- poc.binom.sum$quantiles[grepgo,'50%']
plot(y_p,poc.binom.ynew,xlim=c(0,1),ylim=c(0,1)); abline(a=0,b=1)
d.check <- ifelse(poc.binom.ynew >= mean(y_p),0,1)
sum(d.check)/length(d.check)


# acr severity - observed versus predicted
grepgo <- grep('y.new',colnames(acr.beta[[1]]))
acr.beta.ynew <- acr.beta.sum$statistics[grepgo,'Mean']
acr.beta.ynew <- acr.beta.sum$quantiles[grepgo,'50%']

temp <- moorea[moorea$Taxa=="Acropora",]; temp <- temp[!is.na(temp$Percent_bleached),]; temp <- temp[!is.na(temp$avgTotN),]
temp <- temp[temp$Percent_bleached > 0,]
y <- temp$Percent_bleached; y <- y/100; y <- ifelse(y == 1, 0.99, y)

plot(y,acr.beta.ynew,xlim=c(0,1),ylim=c(0,1)); abline(a=0,b=1)

d.check <- ifelse(acr.beta.ynew >= mean(y),0,1)
sum(d.check)/length(d.check)

hist(y,breaks=20,freq=F)
lines(density(acr.beta.ynew),col='blue')

varF <- (sd(y))^2
varE <- (sd(y-acr.beta.ynew))^2
varF/(varF+varE)


# parameter summaries -----------------------------------------------------

poc.binom.sum$quantiles[param_check,]
poc.beta.sum$quantiles[param_check,]
acr.binom.sum$quantiles[param_check,]
acr.beta.sum$quantiles[param_check,]

write.csv(cbind(
  inv.logit(poc.binom.sum$quantiles[param_check[6:7],]),
  inv.logit(poc.beta.sum$quantiles[param_check[6:7],]),
  inv.logit(acr.binom.sum$quantiles[param_check[6:7],]),
  inv.logit(acr.beta.sum$quantiles[param_check[6:7],])),
          'outputs/TableS3.csv')

poc.avgTotN.s <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(poc.beta[[1]]))
poc.avgTotN.s$mean <- poc.beta.sum$quantiles[grepgo,'50%']
poc.avgTotN.s$up <- poc.beta.sum$quantiles[grepgo,'97.5%']
poc.avgTotN.s$down <- poc.beta.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  poc.avgTotN.s[i,5] <- quantile(c(poc.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  poc.avgTotN.s[i,6] <- quantile(c(poc.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  poc.avgTotN.s[i+2,5] <- quantile(c(poc.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.beta[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  poc.avgTotN.s[i+2,6] <- quantile(c(poc.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.beta[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
poc.avgTotN.s

poc.avgTotN.p <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(poc.binom[[1]]))
poc.avgTotN.p$mean <- poc.binom.sum$quantiles[grepgo,'50%']
poc.avgTotN.p$up <- poc.binom.sum$quantiles[grepgo,'97.5%']
poc.avgTotN.p$down <- poc.binom.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  poc.avgTotN.p[i,5] <- quantile(c(poc.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  poc.avgTotN.p[i,6] <- quantile(c(poc.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  poc.avgTotN.p[i+2,5] <- quantile(c(poc.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.binom[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  poc.avgTotN.p[i+2,6] <- quantile(c(poc.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.binom[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
poc.avgTotN.p


acr.avgTotN.s <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(acr.beta[[1]]))
acr.avgTotN.s$mean <- acr.beta.sum$quantiles[grepgo,'50%']
acr.avgTotN.s$up <- acr.beta.sum$quantiles[grepgo,'97.5%']
acr.avgTotN.s$down <- acr.beta.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  acr.avgTotN.s[i,5] <- quantile(c(acr.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  acr.avgTotN.s[i,6] <- quantile(c(acr.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  acr.avgTotN.s[i+2,5] <- quantile(c(acr.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.beta[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  acr.avgTotN.s[i+2,6] <- quantile(c(acr.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.beta[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
acr.avgTotN.s

acr.avgTotN.p <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(acr.binom[[1]]))
acr.avgTotN.p$mean <- acr.binom.sum$quantiles[grepgo,'50%']
acr.avgTotN.p$up <- acr.binom.sum$quantiles[grepgo,'97.5%']
acr.avgTotN.p$down <- acr.binom.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  acr.avgTotN.p[i,5] <- quantile(c(acr.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  acr.avgTotN.p[i,6] <- quantile(c(acr.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  acr.avgTotN.p[i+2,5] <- quantile(c(acr.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.binom[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  acr.avgTotN.p[i+2,6] <- quantile(c(acr.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.binom[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
acr.avgTotN.p


# png(file='outputs/Figure4.png',height=1500,width=3500,res=300)
# par(mfrow=c(1,2),mar=c(3.5,1.5,2,1),mgp=c(2,1,0),oma=c(1,7,0,0))
# # acropora
# plot(acr.avgTotN.p$mean, c(5,4,3,2,1), xlim=c(-0.92,1.3)
#      , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.5,cex.lab=1.5,bty='l')
# text(-0.99,5.6, expression("A)"~italic(Acropora)),cex=1.5,pos=4,xpd=T)
# axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.5)
# abline(v=0,lty=2)
# #prevelence
# plotCI(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgTotN.p$up, li=acr.avgTotN.p$down
#        , err='x',add=T,lwd=1,pch=NA,sfrac=0)
# plotCI(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgTotN.p$up80, li=acr.avgTotN.p$down80
#        , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
# points(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# 
# # prevalence
# plotCI(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgTotN.p$up, li=acr.avgTotN.p$down
#        , err='x',add=T,lwd=1,pch=NA,sfrac=0)
# plotCI(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1), ui=acr.avgTotN.p$up80, li=acr.avgTotN.p$down80
#        , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
# points(acr.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# 
# # severity
# plotCI(acr.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=acr.avgTotN.s$up, li=acr.avgTotN.s$down
#        , err='x',add=T,lwd=1,pch=NA,sfrac=0)
# plotCI(acr.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=acr.avgTotN.s$up80, li=acr.avgTotN.s$down80
#        , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
# points(acr.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9),pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
# 
# # pocillopora
# plot(poc.avgTotN.p$mean, c(5,4,3,2,1), xlim=c(-0.92,1.3)
#      , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n'
#      ,cex.axis=1.5,cex.lab=1.5,bty='l')
# text(-0.99,5.6, expression("B)"~italic(Pocillopora)),cex=1.5,pos=4,xpd=T)
# axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.3)
# abline(v=0,lty=2)
# #prevelence
# plotCI(poc.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=poc.avgTotN.p$up, li=poc.avgTotN.p$down
#        , err='x',add=T,lwd=1,pch=NA,sfrac=0)
# plotCI(poc.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=poc.avgTotN.p$up80, li=poc.avgTotN.p$down80
#        , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
# points(poc.avgTotN.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# # severity
# plotCI(poc.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=poc.avgTotN.s$up, li=poc.avgTotN.s$down
#        , err='x',add=T,lwd=1,pch=NA,sfrac=0)
# plotCI(poc.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=poc.avgTotN.s$up80, li=poc.avgTotN.s$down80
#        , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
# points(poc.avgTotN.s$mean, c(4.9,3.9,2.9,1.9,0.9),pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)
# 
# legend('topright',legend=c('Prevalence','Severity'),pch=19,col=c(wes_palette("Darjeeling1")[5],wes_palette("Darjeeling1")[1]),bty='n',cex=1.5)
# 
# mtext('Beta',outer=T,side=1,cex=1.5,line=-1.1)
# mtext('Colony size',outer=T,side=2,cex=1.5,line=-0.9,at=0.81,las=2)
# mtext('Depth',outer=T,side=2,cex=1.5,line=-0.9,at=0.66,las=2)
# mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.515,las=2)
# mtext('Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.37,las=2)
# mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.25,las=2)
# mtext('x Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.19,las=2)
# dev.off()


png(file='outputs/Figure4_final.png',height=1500,width=3500,res=300)
par(mfrow=c(1,2),mar=c(3.5,1.5,2,1),mgp=c(2,1,0),oma=c(1,7,0,0))
# severity
plot(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), xlim=c(-0.5,0.5)
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.5,cex.lab=1.7,bty='l')
text(-0.52,5.6,'A) Severity',cex=1.5,pos=4,xpd=T)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.8)
abline(v=0,lty=2)

plotCI(acr.avgTotN.s$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.s$up[1:5], li=acr.avgTotN.s$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[5]
       ))
plotCI(acr.avgTotN.s$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.s$up80[1:5], li=acr.avgTotN.s$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(
         wes_palette("Darjeeling1")[5]
       ))
points(acr.avgTotN.s$mean[1:5], c(5,4,3,2,1),pch=21,bg=c(
  wes_palette("Darjeeling1")[5]
),cex=1.8)

plotCI(poc.avgTotN.s$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.s$up[1:5], li=poc.avgTotN.s$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
plotCI(poc.avgTotN.s$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.s$up80[1:5], li=poc.avgTotN.s$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
points(poc.avgTotN.s$mean[1:5], c(5,4,3,2,1)+0.2,pch=21,bg=c(
  wes_palette("Darjeeling1")[1]
),cex=1.8)
# 

#prevelence
plot(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), xlim=c(-0.92,1.3)
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.5,cex.lab=1.7,bty='l')
text(-0.95,5.6,'B) Prevalence',cex=1.5,pos=4,xpd=T)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.8)
abline(v=0,lty=2)

plotCI(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.p$up[1:5], li=acr.avgTotN.p$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[5]
       ))
plotCI(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.p$up80[1:5], li=acr.avgTotN.p$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(wes_palette("Darjeeling1")[5]))
points(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1),pch=21,bg=c(
  wes_palette("Darjeeling1")[5]
),cex=1.8)

plotCI(poc.avgTotN.p$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.p$up[1:5], li=poc.avgTotN.p$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
plotCI(poc.avgTotN.p$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.p$up80[1:5], li=poc.avgTotN.p$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
points(poc.avgTotN.p$mean[1:5], c(5,4,3,2,1)+0.2,pch=21,bg=c(
  wes_palette("Darjeeling1")[1]
),cex=1.8)
# 
# legend('topright',legend=c(expression(~italic(Pocillopora)),expression(~italic(Acropora))),pch=19,col=c(wes_palette("Darjeeling1")[1],wes_palette("Darjeeling1")[5]),bty='n',cex=1.5)

mtext('Effect',outer=T,side=1,cex=1.5,line=-1.1)
mtext('Colony size',outer=T,side=2,cex=1.5,line=-0.9,at=0.81,las=2)
mtext('Depth',outer=T,side=2,cex=1.5,line=-0.9,at=0.66,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.515,las=2)
mtext('Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.37,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.25,las=2)
mtext('x Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.19,las=2)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend(0.3,-0.85,legend='Acropora',pch=21,lty=1,cex=1.5,bty='n',text.font=3,pt.bg=wes_palette("Darjeeling1")[5],col=wes_palette("Darjeeling1")[5],lwd=2)
legend(0.7,-0.85,legend='Pocillopora',pch=21,lty=1,cex=1.5,bty='n',text.font=3,pt.bg=wes_palette("Darjeeling1")[1],col=wes_palette("Darjeeling1")[1],lwd=2)

dev.off()



setEPS()
postscript('outputs/Figure4_final.eps',height=4.3,width=8.6)
# png(file='outputs/Figure4_final.png',height=1500,width=3500,res=300)
par(mfrow=c(1,2),mar=c(3.5,1.5,2,1),mgp=c(1.7,.7,0),oma=c(1,5,0,0))
# severity
plot(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), xlim=c(-0.5,0.5)
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.2,cex.lab=1.7,bty='l')
text(-0.55,5.6,'A) Severity',cex=1.2,pos=4,xpd=T)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.2,cex.lab=1.2)
abline(v=0,lty=2)

plotCI(acr.avgTotN.s$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.s$up[1:5], li=acr.avgTotN.s$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[5]
       ))
plotCI(acr.avgTotN.s$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.s$up80[1:5], li=acr.avgTotN.s$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(
         wes_palette("Darjeeling1")[5]
       ))
points(acr.avgTotN.s$mean[1:5], c(5,4,3,2,1),pch=21,bg=c(
  wes_palette("Darjeeling1")[5]
),cex=1.5)

plotCI(poc.avgTotN.s$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.s$up[1:5], li=poc.avgTotN.s$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
plotCI(poc.avgTotN.s$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.s$up80[1:5], li=poc.avgTotN.s$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
points(poc.avgTotN.s$mean[1:5], c(5,4,3,2,1)+0.2,pch=21,bg=c(
  wes_palette("Darjeeling1")[1]
),cex=1.5)
# 

#prevelence
plot(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), xlim=c(-0.92,1.3)
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.2,cex.lab=1.7,bty='l')
text(-0.99,5.6,'B) Prevalence',cex=1.2,pos=4,xpd=T)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.2,cex.lab=1.2)
abline(v=0,lty=2)

plotCI(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.p$up[1:5], li=acr.avgTotN.p$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[5]
       ))
plotCI(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1), ui=acr.avgTotN.p$up80[1:5], li=acr.avgTotN.p$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(wes_palette("Darjeeling1")[5]))
points(acr.avgTotN.p$mean[1:5], c(5,4,3,2,1),pch=21,bg=c(
  wes_palette("Darjeeling1")[5]
),cex=1.5)

plotCI(poc.avgTotN.p$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.p$up[1:5], li=poc.avgTotN.p$down[1:5]
       , err='x',add=T,lwd=2,pch=NA,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
plotCI(poc.avgTotN.p$mean[1:5], c(5,4,3,2,1)+0.2, ui=poc.avgTotN.p$up80[1:5], li=poc.avgTotN.p$down80[1:5]
       , err='x',add=T,pch=NA,lwd=4,sfrac=0,col=c(
         wes_palette("Darjeeling1")[1]
       ))
points(poc.avgTotN.p$mean[1:5], c(5,4,3,2,1)+0.2,pch=21,bg=c(
  wes_palette("Darjeeling1")[1]
),cex=1.5)
# 
# legend('topright',legend=c(expression(~italic(Pocillopora)),expression(~italic(Acropora))),pch=19,col=c(wes_palette("Darjeeling1")[1],wes_palette("Darjeeling1")[5]),bty='n',cex=1.5)

mtext('Effect',outer=T,side=1,cex=1.2,line=-1.8)
mtext('Colony size',outer=T,side=2,cex=1.2,line=-0.9,at=0.8,las=2)
mtext('Depth',outer=T,side=2,cex=1.2,line=-0.9,at=0.66,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.2,line=-0.9,at=0.52,las=2)
mtext('Heat Stress',outer=T,side=2,cex=1.2,line=-0.9,at=0.38,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.2,line=-0.9,at=0.26,las=2)
mtext('x Heat Stress',outer=T,side=2,cex=1.2,line=-0.9,at=0.21,las=2)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend(-1,-0.85,legend='Acropora',pch=21,lty=1,cex=1.2,bty='n',text.font=3,pt.bg=wes_palette("Darjeeling1")[5],col=wes_palette("Darjeeling1")[5],lwd=2)
legend(-0.6,-0.85,legend='Pocillopora',pch=21,lty=1,cex=1.2,bty='n',text.font=3,pt.bg=wes_palette("Darjeeling1")[1],col=wes_palette("Darjeeling1")[1],lwd=2)

dev.off()

# interaction -------------------------------------------------------------

#### Acropora - data set up
temp <- moorea[moorea$Taxa=="Acropora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
temp <- temp[temp$Percent_bleached > 0,]
temp$point <- as.numeric(as.factor(as.character(temp$Point)))
site_preds <- temp %>% group_by(point) %>% summarise(cumtemp=unique(cumstress),totN=unique(avgTotN),habitat=unique(Habitat),coast=unique(Island_shore))
grepgo <- grep('pred_X_1',colnames(acr.beta[[1]]))
acr_site_pred <- data.frame(site_preds,
                             fringe_median=inv.logit(acr.beta.sum$quantiles[grepgo,'50%']),
                             fringe_up=inv.logit(acr.beta.sum$quantiles[grepgo,'97.5%']),
                             fringe_down=inv.logit(acr.beta.sum$quantiles[grepgo,'2.5%']))
grepgo <- grep('pred_X_2',colnames(acr.beta[[1]]))
acr_site_pred <- data.frame(acr_site_pred,
                             back_median=inv.logit(acr.beta.sum$quantiles[grepgo,'50%']),
                             back_up=inv.logit(acr.beta.sum$quantiles[grepgo,'97.5%']),
                             back_down=inv.logit(acr.beta.sum$quantiles[grepgo,'2.5%']))
head(acr_site_pred)
acr_site_pred$cumtemp_s <- scale(acr_site_pred$cumtemp)[,1]
acr_site_pred$totN_s <- scale(acr_site_pred$totN)[,1]

acr_site_pred %>% distinct(habitat,cumtemp) %>% arrange(cumtemp) # inspect for values to use for preds

#### Pocillopora - data set up
temp <- moorea[moorea$Taxa=="Pocillopora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
temp <- temp[temp$Percent_bleached > 0,]
temp$point <- as.numeric(as.factor(as.character(temp$Point)))
site_preds <- temp %>% group_by(point) %>% summarise(cumtemp=unique(cumstress),totN=unique(avgTotN),habitat=unique(Habitat),coast=unique(Island_shore))

grepgo <- grep('pred_X_1',colnames(poc.beta[[1]]))
poc_site_pred <- data.frame(site_preds,
                             fringe_median=inv.logit(poc.beta.sum$quantiles[grepgo,'50%']),
                             fringe_up=inv.logit(poc.beta.sum$quantiles[grepgo,'97.5%']),
                             fringe_down=inv.logit(poc.beta.sum$quantiles[grepgo,'2.5%']))
grepgo <- grep('pred_X_2',colnames(poc.beta[[1]]))
poc_site_pred <- data.frame(poc_site_pred,
                            back_median=inv.logit(poc.beta.sum$quantiles[grepgo,'50%']),
                             back_up=inv.logit(poc.beta.sum$quantiles[grepgo,'97.5%']),
                             back_down=inv.logit(poc.beta.sum$quantiles[grepgo,'2.5%']))
head(poc_site_pred)
poc_site_pred$cumtemp_s <- scale(poc_site_pred$cumtemp)[,1]
poc_site_pred$totN_s <- scale(poc_site_pred$totN)[,1]

poc_site_pred %>% distinct(habitat,cumtemp) %>% arrange(cumtemp) # inspect for values to use for preds


#### Acropora severity - low heat stress
acr_min_back <- acr_site_pred %>% 
  filter(habitat=='Lagoon') %>% 
  filter(cumtemp < 2)
plot(acr_min_back$totN_s,acr_min_back$back_median,ylim=c(0.1,0.83),xlim=c(-2.5,3.5))
plotCI(acr_min_back$totN_s,acr_min_back$back_median,ui=acr_min_back$back_up,li=acr_min_back$back_down,add=T)

#### Acropora severity - moderate heat stress
acr_mod_back <- acr_site_pred %>% 
  filter(habitat=='Lagoon') %>% 
  filter(cumtemp > 2.1 & cumtemp < 2.3)
plot(acr_mod_back$totN_s,acr_mod_back$back_median,ylim=c(0.1,0.83),xlim=c(-2.5,3.5))
plotCI(acr_mod_back$totN_s,acr_mod_back$back_median,ui=acr_mod_back$back_up,li=acr_mod_back$back_down,add=T)

#### Acropora severity - high heat stress
acr_high_back <- acr_site_pred %>% 
  filter(habitat=='Lagoon') %>% 
  filter(cumtemp > 2.7)
plot(acr_high_back$totN_s,acr_high_back$back_median,ylim=c(0.1,0.83),xlim=c(-2.5,3.5))
plotCI(acr_high_back$totN_s,acr_high_back$back_median,ui=acr_high_back$back_up,li=acr_high_back$back_down,add=T)

#### Pocillopora severity - low heat stress
poc_min_back <- poc_site_pred %>% 
  filter(habitat=='Lagoon') %>% 
  filter(cumtemp < 2)
plot(poc_min_back$totN_s,poc_min_back$back_median,ylim=c(0.1,0.83),xlim=c(-2.5,3.5))
plotCI(poc_min_back$totN_s,poc_min_back$back_median,ui=poc_min_back$back_up,li=poc_min_back$back_down,add=T)

#### Pocillopora severity - moderate heat stress
poc_mod_back <- poc_site_pred %>% 
  filter(habitat=='Lagoon') %>% 
  filter(cumtemp > 2.09 & cumtemp < 2.14)
plot(poc_mod_back$totN_s,poc_mod_back$back_median,ylim=c(0.1,0.83),xlim=c(-2.5,3.5))
plotCI(poc_mod_back$totN_s,poc_mod_back$back_median,ui=poc_mod_back$back_up,li=poc_mod_back$back_down,add=T)

#### Pocillopora severity - high heat stress
poc_high_back <- poc_site_pred %>% 
  filter(habitat=='Lagoon') %>% 
  filter(cumtemp > 2.8)
plot(poc_high_back$totN_s,poc_high_back$back_median,ylim=c(0.1,0.83),xlim=c(-2.5,3.5))
plotCI(poc_high_back$totN_s,poc_high_back$back_median,ui=poc_high_back$back_up,li=poc_high_back$back_down,add=T)


####### combine
# png(file='outputs/Figure5.png',height=2300,width=3100,res=300)
# par(mgp=c(2.2,.9,0),oma=c(4,4,5,0),mfrow=c(2,3),mar=c(1.5,1,1,1))
# 
# plot(acr_min_back$totN_s,acr_min_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
#      xlab='',ylab='',cex.axis=1.7,bty='l',xaxt='n')
# axis(1,at=c(-2,-1,0,1,2,3),labels=c('','','','','',''))
# plotCI(acr_min_back$totN_s,acr_min_back$back_median,ui=acr_min_back$back_up,li=acr_min_back$back_down,add=T
#        ,pch=19,cex=1.3,sfrac=0,lwd=1.1)
# 
# plot(acr_mod_back$totN_s,acr_mod_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
#      xlab='',ylab='',cex.axis=1.7,bty='l',xaxt='n',yaxt='n')
# axis(1,at=c(-2,-1,0,1,2,3),labels=c('','','','','',''))
# axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
# plotCI(acr_mod_back$totN_s,acr_mod_back$back_median,ui=acr_mod_back$back_up,li=acr_mod_back$back_down,add=T
#        ,pch=19,cex=1.3,sfrac=0,lwd=1.1)
# 
# plot(acr_high_back$totN_s,acr_high_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
#      xlab='',ylab='',cex.axis=1.7,bty='l',xaxt='n',yaxt='n')
# axis(1,at=c(-2,-1,0,1,2,3),labels=c('','','','','',''))
# axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
# plotCI(acr_high_back$totN_s,acr_high_back$back_median,ui=acr_high_back$back_up,li=acr_high_back$back_down,add=T
#        ,pch=19,cex=1.3,sfrac=0,lwd=1.1)
# 
# plot(poc_min_back$totN_s,poc_min_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
#      xlab='',ylab='',cex.axis=1.7,bty='l')
# plotCI(poc_min_back$totN_s,poc_min_back$back_median,ui=poc_min_back$back_up,li=poc_min_back$back_down,add=T
#        ,pch=19,cex=1.3,sfrac=0,lwd=1.19)
# 
# plot(poc_mod_back$totN_s,poc_mod_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
#      xlab='',ylab='',cex.axis=1.7,bty='l',yaxt='n')
# axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
# plotCI(poc_mod_back$totN_s,poc_mod_back$back_median,ui=poc_mod_back$back_up,li=poc_mod_back$back_down,add=T
#        ,pch=19,cex=1.3,sfrac=0,lwd=1.1)
# 
# plot(poc_high_back$totN_s,poc_high_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
#      xlab='',ylab='',cex.axis=1.7,bty='l',yaxt='n')
# axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
# plotCI(poc_high_back$totN_s,poc_high_back$back_median,ui=poc_high_back$back_up,li=poc_high_back$back_down,add=T
#        ,pch=19,cex=1.3,sfrac=0,lwd=1.1)
# 
# mtext('Total N',side=1,cex=1.5,outer=T,line=1.5)
# mtext('Bleaching Severity',side=2,cex=1.5,outer=T,line=1.5)
# mtext('Low',side=3,cex=1.5,outer=T,at=0.15,line=1)
# mtext('Moderate',side=3,cex=1.5,outer=T,line=1)
# mtext('High',side=3,cex=1.5,outer=T,at=0.85,line=1)
# mtext('Heat Stress',side=3,cex=1.5,outer=T,line=3)
# mtext('Pocillopora',side=3,cex=1.5,outer=T,at=0.09,line=-26,font=3)
# mtext('Acropora',side=3,cex=1.5,outer=T,at=0.09,line=-1,font=3)
# 
# dev.off()


png(file='outputs/Figure5_final.png',height=1200,width=3100,res=300)
par(mgp=c(2.2,.9,0),oma=c(4,4,4,0),mfrow=c(1,3),mar=c(1.5,1,1,1))

plot(acr_min_back$totN_s,acr_min_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
     xlab='',ylab='',cex.axis=1.7,bty='l',xaxt='n')
axis(1,at=c(-2,-1,0,1,2,3),labels=c(0.35,0.41,0.48,0.54,0.60,0.66),cex.axis=1.7)
plotCI(acr_min_back$totN_s,acr_min_back$back_median,ui=acr_min_back$back_up,li=acr_min_back$back_down,add=T
        ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[5])
points(acr_min_back$totN_s,acr_min_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[5])

plotCI(poc_min_back$totN_s,poc_min_back$back_median,ui=poc_min_back$back_up,li=poc_min_back$back_down,add=T
        ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[1])
points(poc_min_back$totN_s,poc_min_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[1])

plot(acr_mod_back$totN_s,acr_mod_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
     xlab='',ylab='',cex.axis=1.7,bty='l',yaxt='n',xaxt='n')
axis(1,at=c(-2,-1,0,1,2,3),labels=c(0.35,0.41,0.48,0.54,0.60,0.66),cex.axis=1.7)
axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
plotCI(acr_mod_back$totN_s,acr_mod_back$back_median,ui=acr_mod_back$back_up,li=acr_mod_back$back_down,add=T
        ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[5])
points(acr_mod_back$totN_s,acr_mod_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[5])

plotCI(poc_mod_back$totN_s,poc_mod_back$back_median,ui=poc_mod_back$back_up,li=poc_mod_back$back_down,add=T
        ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[1])
points(poc_mod_back$totN_s,poc_mod_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[1])

plot(acr_high_back$totN_s,acr_high_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
     xlab='',ylab='',cex.axis=1.7,bty='l',yaxt='n',xaxt='n')
axis(1,at=c(-2,-1,0,1,2,3),labels=c(0.35,0.41,0.48,0.54,0.60,0.66),cex.axis=1.7)
axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
plotCI(acr_high_back$totN_s,acr_high_back$back_median,ui=acr_high_back$back_up,li=acr_high_back$back_down,add=T
        ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[5])
points(acr_high_back$totN_s,acr_high_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[5])

plotCI(poc_high_back$totN_s,poc_high_back$back_median,ui=poc_high_back$back_up,li=poc_high_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[1])
points(poc_high_back$totN_s,poc_high_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[1])

mtext('Total N',side=1,cex=1.2,outer=T,line=1.5)
mtext('Bleaching Severity',side=2,cex=1.2,outer=T,line=1.5)
mtext('Low',side=3,cex=1.2,outer=T,at=0.15,line=0)
mtext('Moderate',side=3,cex=1.2,outer=T,line=0)
mtext('High',side=3,cex=1.2,outer=T,at=0.85,line=0)
mtext('Heat Stress',side=3,cex=1.2,outer=T,line=2)
# mtext('Pocillopora',side=3,cex=1.5,outer=T,at=0.09,line=-26,font=3)
# mtext('Acropora',side=3,cex=1.5,outer=T,at=0.09,line=-1,font=3)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend(0.35,-0.85,legend='Acropora',pch=19,lty=1,lwd=1.1,cex=1.5,bty='n',text.font=3,col=wes_palette("Darjeeling1")[5])
legend(0.65,-0.85,legend='Pocillopora',pch=19,lty=1,lwd=1.1,cex=1.5,bty='n',text.font=3,col=wes_palette("Darjeeling1")[1])

dev.off()


setEPS()
postscript('outputs/Figure5_final.eps',height=4.3,width=8.6)
par(mgp=c(2.2,.9,0),oma=c(4,4,4,0),mfrow=c(1,3),mar=c(1.5,1,1,.7))

plot(acr_min_back$totN_s,acr_min_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
     xlab='',ylab='',cex.axis=1.7,bty='l',xaxt='n')
axis(1,at=c(-2,-1,0,1,2,3),labels=c(0.35,0.41,0.48,0.54,0.60,0.66),cex.axis=1.7)
plotCI(acr_min_back$totN_s,acr_min_back$back_median,ui=acr_min_back$back_up,li=acr_min_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[5])
points(acr_min_back$totN_s,acr_min_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[5])

plotCI(poc_min_back$totN_s,poc_min_back$back_median,ui=poc_min_back$back_up,li=poc_min_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[1])
points(poc_min_back$totN_s,poc_min_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[1])

plot(acr_mod_back$totN_s,acr_mod_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
     xlab='',ylab='',cex.axis=1.7,bty='l',yaxt='n',xaxt='n')
axis(1,at=c(-2,-1,0,1,2,3),labels=c(0.35,0.41,0.48,0.54,0.60,0.66),cex.axis=1.7)
axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
plotCI(acr_mod_back$totN_s,acr_mod_back$back_median,ui=acr_mod_back$back_up,li=acr_mod_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[5])
points(acr_mod_back$totN_s,acr_mod_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[5])

plotCI(poc_mod_back$totN_s,poc_mod_back$back_median,ui=poc_mod_back$back_up,li=poc_mod_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[1])
points(poc_mod_back$totN_s,poc_mod_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[1])

plot(acr_high_back$totN_s,acr_high_back$back_median,ylim=c(0.24,0.8),xlim=c(-2.5,3.5),
     xlab='',ylab='',cex.axis=1.7,bty='l',yaxt='n',xaxt='n')
axis(1,at=c(-2,-1,0,1,2,3),labels=c(0.35,0.41,0.48,0.54,0.60,0.66),cex.axis=1.7)
axis(2,at=c(0.3,0.4,0.5,0.6,0.7,0.8),labels=c('','','','','',''))
plotCI(acr_high_back$totN_s,acr_high_back$back_median,ui=acr_high_back$back_up,li=acr_high_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[5])
points(acr_high_back$totN_s,acr_high_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[5])

plotCI(poc_high_back$totN_s,poc_high_back$back_median,ui=poc_high_back$back_up,li=poc_high_back$back_down,add=T
       ,pch=19,cex=1.3,sfrac=0,lwd=1.3,col=wes_palette("Darjeeling1")[1])
points(poc_high_back$totN_s,poc_high_back$back_median,cex=1.3,pch=21,bg=wes_palette("Darjeeling1")[1])

mtext('Total N',side=1,cex=1.2,outer=T,line=1.5)
mtext('Bleaching Severity',side=2,cex=1.2,outer=T,line=1.5)
mtext('Low',side=3,cex=1.2,outer=T,at=0.15,line=0)
mtext('Moderate',side=3,cex=1.2,outer=T,line=0)
mtext('High',side=3,cex=1.2,outer=T,at=0.85,line=0)
mtext('Heat Stress',side=3,cex=1.2,outer=T,line=2)
# mtext('Pocillopora',side=3,cex=1.5,outer=T,at=0.09,line=-26,font=3)
# mtext('Acropora',side=3,cex=1.5,outer=T,at=0.09,line=-1,font=3)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
# legend(0.25,-0.85,legend='Acropora',pch=19,lty=1,lwd=1.1,cex=1.7,bty='n',text.font=3,col=wes_palette("Darjeeling1")[5])
# legend(0.65,-0.85,legend='Pocillopora',pch=19,lty=1,lwd=1.1,cex=1.7,bty='n',text.font=3,col=wes_palette("Darjeeling1")[1])
legend(-1,-0.85,legend='Acropora',pch=21,lty=1,cex=1.7,bty='n',text.font=3,pt.bg=wes_palette("Darjeeling1")[5],col=wes_palette("Darjeeling1")[5],lwd=2)
legend(-0.6,-0.85,legend='Pocillopora',pch=21,lty=1,cex=1.7,bty='n',text.font=3,pt.bg=wes_palette("Darjeeling1")[1],col=wes_palette("Darjeeling1")[1],lwd=2)
dev.off()






