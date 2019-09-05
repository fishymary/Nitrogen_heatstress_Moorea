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
poc.binom <- readRDS('hpc_out/binom_poc_avgTotN_081519.Rdata')
poc.binom.sum <- summary(poc.binom)

poc.beta <- readRDS('hpc_out/beta_poc_avgTotN_081519.Rdata')
poc.beta.sum <- summary(poc.beta)

acr.binom <- readRDS('hpc_out/binom_acr_avgTotN_081519.Rdata')
acr.binom.sum <- summary(acr.binom)

acr.beta <- readRDS('hpc_out/beta_acr_avgTotN_081519.Rdata')
acr.beta.sum <- summary(acr.beta)

moorea <- read.csv('data/moorea_withavgnuts.csv')

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


png(file='outputs/Figure4_sep.png',height=1500,width=3500,res=300)
par(mfrow=c(1,2),mar=c(3.5,1.5,2,1),mgp=c(2,1,0),oma=c(1,7,0,0))
# acropora
plot(acr.avgTotN.p$mean, c(5,4,3,2,1), xlim=c(min(acr.avgTotN.p$down),max(acr.avgTotN.p$up))
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.5,cex.lab=1.5,bty='l')
text(-4.1,5.5, expression("A)"~italic(Acropora)),cex=1.5)
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
plot(poc.avgTotN.p$mean, c(5,4,3,2,1), xlim=c(min(poc.avgTotN.s$down),3)
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n'
     ,cex.axis=1.5,cex.lab=1.5,bty='l')
text(-1.7,5.5, expression("B)"~italic(Pocillopora)),cex=1.5)
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


# interaction -------------------------------------------------------------

ind.minmax <-
  data.frame(
    temp = c(rep('min', 3), rep('lower', 3), rep('upper', 3), rep('max', 3)),
    nuts = rep(c('low', 'medium', 'high'),4),
    'mean' = rep(NA, 6),
    'up' = rep(NA, 6),
    'down' = rep(NA, 6),
    'up80' = rep(NA, 6),
    'down80' = rep(NA, 6)
  )


###### Pocillopora prevalence
grepgo <- grep('pred.inter.cat.1',colnames(poc.binom[[1]]))
poc.binom.sum$quantiles[grepgo,]
poc.binom.inter <- ind.minmax
poc.binom.inter$mean <- inv.logit(poc.binom.sum$statistics[grepgo,'Mean'])
poc.binom.inter$up <- inv.logit(poc.binom.sum$quantiles[grepgo,'97.5%'])
poc.binom.inter$down <- inv.logit(poc.binom.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  poc.binom.inter$down80[i] <- inv.logit(quantile(c(poc.binom[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                    poc.binom[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                    poc.binom[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.10))
}
for(i in 1:12){
  poc.binom.inter$up80[i] <- inv.logit(quantile(c(poc.binom[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                  poc.binom[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                  poc.binom[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.90))
}

grepgo <- grep('pred.inter.cat.2',colnames(poc.binom[[1]]))
poc.binom.sum$quantiles[grepgo,]
poc.binom.inter2 <- ind.minmax
poc.binom.inter2$mean <- inv.logit(poc.binom.sum$statistics[grepgo,'Mean'])
poc.binom.inter2$up <- inv.logit(poc.binom.sum$quantiles[grepgo,'97.5%'])
poc.binom.inter2$down <- inv.logit(poc.binom.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  poc.binom.inter2$down80[i] <- inv.logit(quantile(c(poc.binom[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                     poc.binom[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                     poc.binom[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.10))
}
for(i in 1:12){
  poc.binom.inter2$up80[i] <- inv.logit(quantile(c(poc.binom[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                   poc.binom[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                   poc.binom[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.90))
}


png(file='outputs/Figure5_poc_prev.png',height=1800,width=3400,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0),oma=c(0,0,2,1))

xseq <- c(1,2,3,5,6,7,9,10,11,13,14,15)
plot(xseq,poc.binom.inter$mean,ylim=c(0,1),xlim=c(0.5,15.5),type='n',ylab='Predicted Bleaching Severity',xlab='Nutrients',xaxt='n',cex.axis=1.3,cex.lab=1.4)
plotCI(xseq-0.1,poc.binom.inter$mean,ui=poc.binom.inter$up,li=poc.binom.inter$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,poc.binom.inter$mean,ui=poc.binom.inter$up80,li=poc.binom.inter$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,poc.binom.inter$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

plotCI(xseq+0.1,poc.binom.inter2$mean,ui=poc.binom.inter2$up,li=poc.binom.inter2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,poc.binom.inter2$mean,ui=poc.binom.inter2$up80,li=poc.binom.inter2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,poc.binom.inter2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)


axis(1,at=xseq,labels=rep(c('Low','Med','High'),4),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
abline(v=8,lty=2)
abline(v=12,lty=2)

mtext('Minimum',side=3,at=2,cex=1.4)
mtext('Low Quartile',side=3,at=6,cex=1.4)
mtext('Up Quartile',side=3,at=10,cex=1.4)
mtext('Maximum',side=3,at=14,cex=1.4)
mtext('Heat Stress',side=3,outer=T,cex=1.4)

legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)

dev.off()


###### Acropora prevalence

grepgo <- grep('pred.inter.cat.1',colnames(acr.binom[[1]]))
acr.binom.sum$quantiles[grepgo,]
acr.binom.inter <- ind.minmax
acr.binom.inter$mean <- inv.logit(acr.binom.sum$statistics[grepgo,'Mean'])
acr.binom.inter$up <- inv.logit(acr.binom.sum$quantiles[grepgo,'97.5%'])
acr.binom.inter$down <- inv.logit(acr.binom.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  acr.binom.inter$down80[i] <- inv.logit(quantile(c(acr.binom[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                    acr.binom[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                    acr.binom[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.10))
}
for(i in 1:12){
  acr.binom.inter$up80[i] <- inv.logit(quantile(c(acr.binom[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                  acr.binom[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                  acr.binom[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.90))
}

grepgo <- grep('pred.inter.cat.2',colnames(acr.binom[[1]]))
acr.binom.sum$quantiles[grepgo,]
acr.binom.inter2 <- ind.minmax
acr.binom.inter2$mean <- inv.logit(acr.binom.sum$statistics[grepgo,'Mean'])
acr.binom.inter2$up <- inv.logit(acr.binom.sum$quantiles[grepgo,'97.5%'])
acr.binom.inter2$down <- inv.logit(acr.binom.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  acr.binom.inter2$down80[i] <- inv.logit(quantile(c(acr.binom[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                     acr.binom[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                     acr.binom[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.10))
}
for(i in 1:12){
  acr.binom.inter2$up80[i] <- inv.logit(quantile(c(acr.binom[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                   acr.binom[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                   acr.binom[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.90))
}


png(file='outputs/Figure5_acr_prev.png',height=1800,width=3400,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0),oma=c(0,0,2,1))

xseq <- c(1,2,3,5,6,7,9,10,11,13,14,15)
plot(xseq,acr.binom.inter$mean,ylim=c(0,1),xlim=c(0.5,15.5),type='n',ylab='Predicted Bleaching Severity',xlab='Nutrients',xaxt='n',cex.axis=1.3,cex.lab=1.4)
plotCI(xseq-0.1,acr.binom.inter$mean,ui=acr.binom.inter$up,li=acr.binom.inter$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,acr.binom.inter$mean,ui=acr.binom.inter$up80,li=acr.binom.inter$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,acr.binom.inter$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

plotCI(xseq+0.1,acr.binom.inter2$mean,ui=acr.binom.inter2$up,li=acr.binom.inter2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,acr.binom.inter2$mean,ui=acr.binom.inter2$up80,li=acr.binom.inter2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,acr.binom.inter2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)


axis(1,at=xseq,labels=rep(c('Low','Med','High'),4),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
abline(v=8,lty=2)
abline(v=12,lty=2)

mtext('Minimum',side=3,at=2,cex=1.4)
mtext('Low Quartile',side=3,at=6,cex=1.4)
mtext('Up Quartile',side=3,at=10,cex=1.4)
mtext('Maximum',side=3,at=14,cex=1.4)
mtext('Heat Stress',side=3,outer=T,cex=1.4)

legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)

dev.off()


###### Pocillopora severity
grepgo <- grep('pred.inter.cat.1',colnames(poc.beta[[1]]))
poc.beta.sum$quantiles[grepgo,]
poc.beta.inter <- ind.minmax
poc.beta.inter$mean <- inv.logit(poc.beta.sum$statistics[grepgo,'Mean'])
poc.beta.inter$up <- inv.logit(poc.beta.sum$quantiles[grepgo,'97.5%'])
poc.beta.inter$down <- inv.logit(poc.beta.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  poc.beta.inter$down80[i] <- inv.logit(quantile(c(poc.beta[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                   poc.beta[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                   poc.beta[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.10))
}
for(i in 1:12){
  poc.beta.inter$up80[i] <- inv.logit(quantile(c(poc.beta[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                 poc.beta[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                 poc.beta[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.90))
}

grepgo <- grep('pred.inter.cat.2',colnames(poc.beta[[1]]))
poc.beta.sum$quantiles[grepgo,]
poc.beta.inter2 <- ind.minmax
poc.beta.inter2$mean <- inv.logit(poc.beta.sum$statistics[grepgo,'Mean'])
poc.beta.inter2$up <- inv.logit(poc.beta.sum$quantiles[grepgo,'97.5%'])
poc.beta.inter2$down <- inv.logit(poc.beta.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  poc.beta.inter2$down80[i] <- inv.logit(quantile(c(poc.beta[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                    poc.beta[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                    poc.beta[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.10))
}
for(i in 1:12){
  poc.beta.inter2$up80[i] <- inv.logit(quantile(c(poc.beta[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                  poc.beta[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                  poc.beta[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.90))
}


png(file='outputs/Figure5_poc_sev.png',height=1800,width=3400,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0),oma=c(0,0,2,1))

xseq <- c(1,2,3,5,6,7,9,10,11,13,14,15)
plot(xseq,poc.beta.inter$mean,ylim=c(0,1),xlim=c(0.5,15.5),type='n',ylab='Predicted Bleaching Severity',xlab='Nitrogen',xaxt='n',cex.axis=1.3,cex.lab=1.4)
plotCI(xseq-0.1,poc.beta.inter$mean,ui=poc.beta.inter$up,li=poc.beta.inter$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,poc.beta.inter$mean,ui=poc.beta.inter$up80,li=poc.beta.inter$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,poc.beta.inter$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

plotCI(xseq+0.1,poc.beta.inter2$mean,ui=poc.beta.inter2$up,li=poc.beta.inter2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,poc.beta.inter2$mean,ui=poc.beta.inter2$up80,li=poc.beta.inter2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,poc.beta.inter2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)


axis(1,at=xseq,labels=rep(c('Low','Med','High'),4),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
abline(v=8,lty=2)
abline(v=12,lty=2)

mtext('Minimum',side=3,at=2,cex=1.4)
mtext('Low Quartile',side=3,at=6,cex=1.4)
mtext('Up Quartile',side=3,at=10,cex=1.4)
mtext('Maximum',side=3,at=14,cex=1.4)
mtext('Heat Stress',side=3,outer=T,cex=1.4)

legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)

dev.off()


###### Acropora severity

grepgo <- grep('pred.inter.cat.1',colnames(acr.beta[[1]]))
acr.beta.sum$quantiles[grepgo,]
acr.beta.inter <- ind.minmax
acr.beta.inter$mean <- inv.logit(acr.beta.sum$statistics[grepgo,'Mean'])
acr.beta.inter$up <- inv.logit(acr.beta.sum$quantiles[grepgo,'97.5%'])
acr.beta.inter$down <- inv.logit(acr.beta.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  acr.beta.inter$down80[i] <- inv.logit(quantile(c(acr.beta[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                   acr.beta[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                   acr.beta[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.10))
}
for(i in 1:12){
  acr.beta.inter$up80[i] <- inv.logit(quantile(c(acr.beta[[1]][,paste0('pred.inter.cat.1[',i,']')],
                                                 acr.beta[[2]][,paste0('pred.inter.cat.1[',i,']')],
                                                 acr.beta[[3]][,paste0('pred.inter.cat.1[',i,']')]),0.90))
}

grepgo <- grep('pred.inter.cat.2',colnames(acr.beta[[1]]))
acr.beta.sum$quantiles[grepgo,]
acr.beta.inter2 <- ind.minmax
acr.beta.inter2$mean <- inv.logit(acr.beta.sum$statistics[grepgo,'Mean'])
acr.beta.inter2$up <- inv.logit(acr.beta.sum$quantiles[grepgo,'97.5%'])
acr.beta.inter2$down <- inv.logit(acr.beta.sum$quantiles[grepgo,'2.5%'])
for(i in 1:12){
  acr.beta.inter2$down80[i] <- inv.logit(quantile(c(acr.beta[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                    acr.beta[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                    acr.beta[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.10))
}
for(i in 1:12){
  acr.beta.inter2$up80[i] <- inv.logit(quantile(c(acr.beta[[1]][,paste0('pred.inter.cat.2[',i,']')],
                                                  acr.beta[[2]][,paste0('pred.inter.cat.2[',i,']')],
                                                  acr.beta[[3]][,paste0('pred.inter.cat.2[',i,']')]),0.90))
}


png(file='outputs/Figure5_acr_sev.png',height=1800,width=3400,res=300)
par(mar=c(4,4,2,1),mgp=c(2.2,0.9,0),oma=c(0,0,2,1))

xseq <- c(1,2,3,5,6,7,9,10,11,13,14,15)
plot(xseq,acr.beta.inter$mean,ylim=c(0,1),xlim=c(0.5,15.5),type='n',ylab='Predicted Bleaching Severity',xlab='Nutrients',xaxt='n',cex.axis=1.3,cex.lab=1.4)
plotCI(xseq-0.1,acr.beta.inter$mean,ui=acr.beta.inter$up,li=acr.beta.inter$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq-0.1,acr.beta.inter$mean,ui=acr.beta.inter$up80,li=acr.beta.inter$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq-0.1,acr.beta.inter$mean,pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

plotCI(xseq+0.1,acr.beta.inter2$mean,ui=acr.beta.inter2$up,li=acr.beta.inter2$down,add=T,pch=NA,lwd=1,sfrac=0)
plotCI(xseq+0.1,acr.beta.inter2$mean,ui=acr.beta.inter2$up80,li=acr.beta.inter2$down80,add=T,pch=NA,lwd=4,sfrac=0)
points(xseq+0.1,acr.beta.inter2$mean,pch=21,bg='white',col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)


axis(1,at=xseq,labels=rep(c('Low','Med','High'),4),cex.axis=1.3,cex.lab=1.4)
abline(v=4,lty=2)
abline(v=8,lty=2)
abline(v=12,lty=2)

mtext('Minimum',side=3,at=2,cex=1.4)
mtext('Low Quartile',side=3,at=6,cex=1.4)
mtext('Up Quartile',side=3,at=10,cex=1.4)
mtext('Maximum',side=3,at=14,cex=1.4)
mtext('Heat Stress',side=3,outer=T,cex=1.4)

legend('topleft',legend=c('Fringing Reef','Backreef'),pch=c(19,21),col=wes_palette("Darjeeling1")[1],bty='n',cex=1.4)

dev.off()



# dN15 --------------------------------------------------------------------
poc.dN15.binom <- readRDS('hpc_out/binom_poc_avgdN15_081519.Rdata')
poc.dN15.binom.sum <- summary(poc.dN15.binom)

poc.dN15.beta <- readRDS('hpc_out/beta_poc_avgdN15_081519.Rdata')
poc.dN15.beta.sum <- summary(poc.dN15.beta)

acr.dN15.binom <- readRDS('hpc_out/binom_acr_avgdN15_081519.Rdata')
acr.dN15.binom.sum <- summary(acr.dN15.binom)

acr.dN15.beta <- readRDS('hpc_out/beta_acr_avgdN15_081519.Rdata')
acr.dN15.beta.sum <- summary(acr.dN15.beta)


poc.avgdN15.s <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(poc.dN15.beta[[1]]))
poc.avgdN15.s$mean <- poc.dN15.beta.sum$quantiles[grepgo,'50%']
poc.avgdN15.s$up <- poc.dN15.beta.sum$quantiles[grepgo,'97.5%']
poc.avgdN15.s$down <- poc.dN15.beta.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  poc.avgdN15.s[i,5] <- quantile(c(poc.dN15.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  poc.avgdN15.s[i,6] <- quantile(c(poc.dN15.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  poc.avgdN15.s[i+2,5] <- quantile(c(poc.dN15.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.beta[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  poc.avgdN15.s[i+2,6] <- quantile(c(poc.dN15.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.beta[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
poc.avgdN15.s

poc.avgdN15.p <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(poc.dN15.binom[[1]]))
poc.avgdN15.p$mean <- poc.dN15.binom.sum$quantiles[grepgo,'50%']
poc.avgdN15.p$up <- poc.dN15.binom.sum$quantiles[grepgo,'97.5%']
poc.avgdN15.p$down <- poc.dN15.binom.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  poc.avgdN15.p[i,5] <- quantile(c(poc.dN15.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  poc.avgdN15.p[i,6] <- quantile(c(poc.dN15.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   poc.dN15.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  poc.avgdN15.p[i+2,5] <- quantile(c(poc.dN15.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.binom[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  poc.avgdN15.p[i+2,6] <- quantile(c(poc.dN15.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     poc.dN15.binom[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
poc.avgdN15.p


acr.avgdN15.s <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(acr.dN15.beta[[1]]))
acr.avgdN15.s$mean <- acr.dN15.beta.sum$quantiles[grepgo,'50%']
acr.avgdN15.s$up <- acr.dN15.beta.sum$quantiles[grepgo,'97.5%']
acr.avgdN15.s$down <- acr.dN15.beta.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  acr.avgdN15.s[i,5] <- quantile(c(acr.dN15.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  acr.avgdN15.s[i,6] <- quantile(c(acr.dN15.beta[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.beta[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.beta[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  acr.avgdN15.s[i+2,5] <- quantile(c(acr.dN15.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.beta[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  acr.avgdN15.s[i+2,6] <- quantile(c(acr.dN15.beta[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.beta[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.beta[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
acr.avgdN15.s

acr.avgdN15.p <- data.frame(beta=c('Colony Size','Depth','Nutrients','Cum Heat','Nutrients x Heat')
                            ,mean=NA,up=NA,down=NA,up80=NA,down80=NA)
grepgo <- grep('beta',colnames(acr.dN15.binom[[1]]))
acr.avgdN15.p$mean <- acr.dN15.binom.sum$quantiles[grepgo,'50%']
acr.avgdN15.p$up <- acr.dN15.binom.sum$quantiles[grepgo,'97.5%']
acr.avgdN15.p$down <- acr.dN15.binom.sum$quantiles[grepgo,'2.5%']
for(i in 1:2){
  acr.avgdN15.p[i,5] <- quantile(c(acr.dN15.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.90)
  acr.avgdN15.p[i,6] <- quantile(c(acr.dN15.binom[[1]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.binom[[2]][,paste0('beta_colony_p[',i,']')],
                                   acr.dN15.binom[[3]][,paste0('beta_colony_p[',i,']')]),0.10)
}
for(i in 1:3){
  acr.avgdN15.p[i+2,5] <- quantile(c(acr.dN15.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.binom[[3]][,paste0('beta_point_p[',i,']')]),0.90)
  acr.avgdN15.p[i+2,6] <- quantile(c(acr.dN15.binom[[1]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.binom[[2]][,paste0('beta_point_p[',i,']')],
                                     acr.dN15.binom[[3]][,paste0('beta_point_p[',i,']')]),0.10)
}
acr.avgdN15.p


png(file='outputs/Figure4_sep_dN15.png',height=1500,width=3500,res=300)
par(mfrow=c(1,2),mar=c(3.5,1.5,2,1),mgp=c(2,1,0),oma=c(1,7,0,0))
# acropora
plot(acr.avgdN15.p$mean, c(5,4,3,2,1), xlim=c(min(acr.avgdN15.p$down),max(acr.avgdN15.p$up))
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n',cex.axis=1.5,cex.lab=1.5,bty='l')
text(-2.6,5.5, expression("A)"~italic(Acropora)),cex=1.5)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.5)
abline(v=0,lty=2)
#prevelence
plotCI(acr.avgdN15.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgdN15.p$up, li=acr.avgdN15.p$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(acr.avgdN15.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=acr.avgdN15.p$up80, li=acr.avgdN15.p$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(acr.avgdN15.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# severity
plotCI(acr.avgdN15.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=acr.avgdN15.s$up, li=acr.avgdN15.s$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(acr.avgdN15.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=acr.avgdN15.s$up80, li=acr.avgdN15.s$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(acr.avgdN15.s$mean, c(4.9,3.9,2.9,1.9,0.9),pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

# pocillopora
plot(poc.avgdN15.p$mean, c(5,4,3,2,1), xlim=c(-2.5,2.5)
     , ylim=c(0.7,5.6),type='n',ylab='',xlab='',yaxt='n'
     ,cex.axis=1.5,cex.lab=1.5,bty='l')
text(-1.5,5.5, expression("B)"~italic(Pocillopora)),cex=1.5)
axis(2,at=c(5,4,3,2,1),labels = rep("",5),las=1,cex.axis=1.3,cex.lab=1.3)
abline(v=0,lty=2)
#prevelence
plotCI(poc.avgdN15.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=poc.avgdN15.p$up, li=poc.avgdN15.p$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(poc.avgdN15.p$mean, c(5.1,4.1,3.1,2.1,1.1), ui=poc.avgdN15.p$up80, li=poc.avgdN15.p$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(poc.avgdN15.p$mean, c(5.1,4.1,3.1,2.1,1.1),pch=19,col=rep(wes_palette("Darjeeling1")[5],5),cex=1.5)
# severity
plotCI(poc.avgdN15.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=poc.avgdN15.s$up, li=poc.avgdN15.s$down
       , err='x',add=T,lwd=1,pch=NA,sfrac=0)
plotCI(poc.avgdN15.s$mean, c(4.9,3.9,2.9,1.9,0.9), ui=poc.avgdN15.s$up80, li=poc.avgdN15.s$down80
       , err='x',add=T,pch=NA,lwd=3,sfrac=0,col='black')
points(poc.avgdN15.s$mean, c(4.9,3.9,2.9,1.9,0.9),pch=19,col=rep(wes_palette("Darjeeling1")[1],5),cex=1.5)

legend('topright',legend=c('Prevalence','Severity'),pch=19,col=c(wes_palette("Darjeeling1")[5],wes_palette("Darjeeling1")[1]),bty='n',cex=1.5)

mtext('Beta',outer=T,side=1,cex=1.5,line=-1.1)
mtext('Colony size',outer=T,side=2,cex=1.5,line=-0.9,at=0.81,las=2)
mtext('Depth',outer=T,side=2,cex=1.5,line=-0.9,at=0.66,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.515,las=2)
mtext('Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.37,las=2)
mtext('Nitrogen',outer=T,side=2,cex=1.5,line=-0.9,at=0.25,las=2)
mtext('x Heat Stress',outer=T,side=2,cex=1.5,line=-0.9,at=0.19,las=2)
dev.off()

