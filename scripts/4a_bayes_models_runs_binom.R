# -----------------------------------------------------------------------
# Nutrient pollution alters coral bleaching across the seascape
# Donovan et al. 
# 4 - Bayesian hierarchical model for testing interaction between nutrients and temperature on bleaching
# binomial and beta models seperated
# -----------------------------------------------------------------------

# This script runs a Bayesian hierarchical models for testing interaction between nutrients and temperature on bleaching
# Prevalence is modeled seperately from severity

# Initialization ----------------------------------------------------------
library(rjags)
library(dplyr)
library(parallel)

# data --------------------------------------------------------------------
moorea <- read.csv('data/moorea_withavgnuts.csv')

# prevalence --------------------------------------------------------------

cat("model{
    # Colony - prevalence 
    for(i in 1:n){
      y[i] ~ dbern(mu[i])
      logit(mu[i]) <- B0_p[point[i]] + inprod(beta_colony_p,X_colony[i,])
      y.new[i] ~ dbern(mu[i]) 
    }
    
    for(b in 1:B){
      beta_colony_p[b] ~ dnorm(0,0.01)
    }
    
    # Point - prevalence
    for(p in 1:P){
      B0_p[p] ~ dnorm(g_p[p],tau_point_p)
      g_p[p] <- B0_habitat_p[habitat[p]] + inprod(beta_point_p,X_point[p,])
      
      B0.new_p[p] ~ dnorm(g_p[p],tau_point_p) # posterior predictive distribution of B0 for plotting
    }
    for(c in 1:C){
      beta_point_p[c] ~ dnorm(0,0.01)
    }
    sigma_point_p ~ dunif(0,50)
    tau_point_p <- pow(sigma_point_p,-2) 
    
    # habitat - prevelence 
    for(k in 1:K){
     B0_habitat_p[k] ~ dnorm(mu_coast_p[coast[k]],tau_habitat_p)
    }
    sigma_habitat_p ~ dunif(0,50)
    tau_habitat_p <- pow(sigma_habitat_p,-2) 
    
    # coast - prevelence
    for(m in 1:M){
      mu_coast_p[m] ~ dnorm(mu_overall_p,tau_coast_p)
    }
    sigma_coast_p ~ dunif(0,50)
    tau_coast_p <- pow(sigma_coast_p,-2) 
    
    # overall - prevelence
    mu_overall_p ~ dnorm(0,0.01)
    
    
    # Posterior predictive checks
    y.mean <- mean(y)
    ynew.mean <- mean(y.new)
    pval.mean <- step(y.mean-ynew.mean)
    
    sd.y <- sd(y)
    sd.y.new <- sd(y.new)
    pval.sd <- step(sd.y.new-sd.y)

    cv.y <- sd(y)/mean(y)
    cv.y.new <- sd(y.new)/mean(y.new)
    pval.cv <- step(cv.y.new - cv.y)
    
    for(i in 1:n){
      sq[i] <- (y[i] - mu[i])^2
      sq.new[i] <- (y.new[i] - mu[i])^2
    }
    fit <- sum(sq)
    fit.new <- sum(sq.new)
    pval.fit <- step(fit.new - fit)
    
    min.y <- min(y)
    min.y.new <- min(y.new)
    pval.min <- step(min.y.new - min.y)
    
    max.y <- max(y)
    max.y.new <- max(y.new)
    pval.max <- step(max.y.new - max.y)
     
    # R-squared
    varF <- sd(y)^2
    varE <- sd(y - y.new)^2
    R2 <- varF/(varF+varE)

   
    # predictions
    for(z in 1:50){
        pred.interA_p.1[z] <- B0_habitat_p[1] + inprod(beta_point_p,interA[z,])
        pred.interB_p.1[z] <- B0_habitat_p[1] + inprod(beta_point_p,interB[z,])
        pred.interC_p.1[z] <- B0_habitat_p[1] + inprod(beta_point_p,interC[z,])
        
        pred.interA_p.2[z] <- B0_habitat_p[2] + inprod(beta_point_p,interA[z,])
        pred.interB_p.2[z] <- B0_habitat_p[2] + inprod(beta_point_p,interB[z,])
        pred.interC_p.2[z] <- B0_habitat_p[2] + inprod(beta_point_p,interC[z,])
    }

   for(z in 1:12){
        pred.inter.cat.1[z] <- B0_habitat_p[1] + inprod(beta_point_p,interAll[z,])
        pred.inter.cat.2[z] <- B0_habitat_p[2] + inprod(beta_point_p,interAll[z,])
    }
    
    }  ", file="binom_hierarchical.jags")

# Pocillopora & TotN ------------------------------------------------------
temp <- moorea[moorea$Taxa=="Pocillopora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
# temp <- temp[sample(nrow(temp),100),]
y <- temp$Percent_bleached
y <- y/100
y_p <- ifelse(y > 0, 1, 0)

# point
temp$point <- as.numeric(as.factor(as.character(temp$Point)))
point <- temp$point

# colony level predictors
x_colony <- as.matrix(data.frame(temp$Colony_size_class,temp$Depth_m))
for(i in 1:ncol(x_colony)) x_colony[,i] <- scale(x_colony[,i])[,1]
B <- ncol(x_colony)

# point level preds
site_preds <- temp %>% group_by(point) %>% summarise(cumtemp=unique(cumstress),totN=unique(avgTotN),habitat=unique(Habitat),coast=unique(Island_shore))
# summary(site_preds %>% group_by(point) %>% summarise(length(point)))
site_preds$cumtempXtotN <- site_preds$cumtemp*site_preds$totN

X_point <- as.matrix(data.frame(site_preds$totN,site_preds$cumtemp,site_preds$totN*site_preds$cumtemp))
for(i in 1:ncol(X_point)) X_point[,i] <- scale(X_point[,i])[,1]

# new data for predictions
totN_new <- as.matrix(data.frame('totN'=seq(from=min(X_point[,1]),to=max(X_point[,1]),length=50),
                                 'cumtemp'=rep(mean(X_point[,2]),length=50),
                                 'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

cumtemp_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                    'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                                    'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

inter_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                  'cumtemp'=rep(mean(X_point[,2]),length=50),
                                  'cumtempXtotN'=seq(from=min(X_point[,3]),to=max(X_point[,3]),length=50)))

# interaction as 3 categories defined by lower, middle, upper quartiles
interA <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25),50))*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))

interB <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),50)*seq(from=min(X_point[,2])
                                                                                       ,to=max(X_point[,2]),length=50)))

interC <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),50)*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))


interA.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interB.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interC.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interA.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interB.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interC.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

inter_min <- as.matrix(data.frame(
  # min heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(min(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(min(X_point[,2]),3)
))
inter_low <- as.matrix(data.frame(
  # lower q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.25),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.25),3)
))
inter_up <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.75),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.75),3)
))
inter_max <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(max(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(max(X_point[,2]),3)
))
interAll <- rbind(inter_min,inter_low,inter_up,inter_max)


jd <- list(B=B,
           y=y_p,
           n=length(y),
           X_colony=x_colony,
           point=point,
           P=length(unique(temp$point)),
           habitat=as.numeric(as.factor(as.character(site_preds$habitat))),
           X_point=X_point,
           C=ncol(X_point),
           K=2,
           coast=as.numeric(as.factor(as.character(site_preds$coast))),
           M=3,
           interA=interA,interB=interB,interC=interC,interAll=interAll
)

nXcol <- ncol(x_colony)
nXpoint <- ncol(X_point)

initFunc <- function(){return(list(
  beta_colony_p=rnorm(nXcol,0,1),
  beta_point_p=rnorm(nXpoint,-1,1),
  sigma_point_p=runif(1,0,20),
  sigma_habitat_p=runif(1,0,20),
  sigma_coast_p=runif(1,0,1),
  mu_overall_p=rnorm(1,0,1)
))}

n.adapt <- 500; n.update <- 5000; n.iter <- 10000
# n.adapt <- 100; n.update <- 500; n.iter <- 1000

# run chains in parallel
cl <- makeCluster(3) # this determines the number of chains, must be less than the number of cores on computer
clusterExport(cl, c('jd','n.adapt','n.update','n.iter','initFunc','nXcol','nXpoint'))

out <- clusterEvalQ(cl,{
  library(rjags)
  tmod <- jags.model("binom_hierarchical.jags", data= jd, n.chains=1, n.adapt=n.adapt,inits = initFunc())
  update(tmod, n.iter = n.update)
  zmCore <- coda.samples(tmod, c(
                                 "beta_colony_p","beta_point_p","B0_habitat_p",
                                 "mu_coast_p","pval.mean","y.new",'pval.sd','R2',
                                 'pred.interA_p.1','pred.interB_p.1','pred.interC_p.1',
                                 'pred.interA_p.2','pred.interB_p.2','pred.interC_p.2',
                                  'pred.interA_p.max.1','pred.interB_p.max.1','pred.interC_p.max.1',
                                 'pred.interA_p.max.2','pred.interB_p.max.2','pred.interC_p.max.2',
                                 'pred.interA_p.q.1','pred.interB_p.q.1','pred.interC_p.q.1',
                                 'pred.interA_p.q.2','pred.interB_p.q.2','pred.interC_p.q.2',
                                 'pred.inter.cat.1',
                                 'pred.inter.cat.2',
                                 'pval.cv','pval.fit','pval.min','pval.max'),n.iter=n.iter, n.thin=1)
  return(as.mcmc(zmCore))
})
zmPrp <- mcmc.list(out)
stopCluster(cl)

saveRDS(zmPrp, 'binom_poc_avgTotN_081519.Rdata')


# Acropora & TotN ------------------------------------------------------
temp <- moorea[moorea$Taxa=="Acropora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgTotN),]
# temp <- temp[sample(nrow(temp),100),]
y <- temp$Percent_bleached
y <- y/100
y_p <- ifelse(y > 0, 1, 0)

# point
temp$point <- as.numeric(as.factor(as.character(temp$Point)))
point <- temp$point

# colony level predictors
x_colony <- as.matrix(data.frame(temp$Colony_size_class,temp$Depth_m))
for(i in 1:ncol(x_colony)) x_colony[,i] <- scale(x_colony[,i])[,1]
B <- ncol(x_colony)

# point level preds
site_preds <- temp %>% group_by(point) %>% summarise(cumtemp=unique(cumstress),totN=unique(avgTotN),habitat=unique(Habitat),coast=unique(Island_shore))
# summary(site_preds %>% group_by(point) %>% summarise(length(point)))
site_preds$cumtempXtotN <- site_preds$cumtemp*site_preds$totN

X_point <- as.matrix(data.frame(site_preds$totN,site_preds$cumtemp,site_preds$totN*site_preds$cumtemp))
for(i in 1:ncol(X_point)) X_point[,i] <- scale(X_point[,i])[,1]

# new data for predictions
totN_new <- as.matrix(data.frame('totN'=seq(from=min(X_point[,1]),to=max(X_point[,1]),length=50),
                                 'cumtemp'=rep(mean(X_point[,2]),length=50),
                                 'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

cumtemp_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                    'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                                    'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

inter_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                  'cumtemp'=rep(mean(X_point[,2]),length=50),
                                  'cumtempXtotN'=seq(from=min(X_point[,3]),to=max(X_point[,3]),length=50)))

# interaction as 3 categories defined by lower, middle, upper quartiles
interA <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25),50))*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))

interB <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),50)*seq(from=min(X_point[,2])
                                                                                       ,to=max(X_point[,2]),length=50)))

interC <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),50)*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))


interA.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interB.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interC.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interA.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interB.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interC.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

inter_min <- as.matrix(data.frame(
  # min heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(min(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(min(X_point[,2]),3)
))
inter_low <- as.matrix(data.frame(
  # lower q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.25),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.25),3)
))
inter_up <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.75),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.75),3)
))
inter_max <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(max(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(max(X_point[,2]),3)
))
interAll <- rbind(inter_min,inter_low,inter_up,inter_max)


jd <- list(B=B,
           y=y_p,
           n=length(y),
           X_colony=x_colony,
           point=point,
           P=length(unique(temp$point)),
           habitat=as.numeric(as.factor(as.character(site_preds$habitat))),
           X_point=X_point,
           C=ncol(X_point),
           K=2,
           coast=as.numeric(as.factor(as.character(site_preds$coast))),
           M=3,
           interA=interA,interB=interB,interC=interC,interAll=interAll
)

nXcol <- ncol(x_colony)
nXpoint <- ncol(X_point)

initFunc <- function(){return(list(
  beta_colony_p=rnorm(nXcol,0,1),
  beta_point_p=rnorm(nXpoint,-1,1),
  sigma_point_p=runif(1,0,20),
  sigma_habitat_p=runif(1,0,20),
  sigma_coast_p=runif(1,0,1),
  mu_overall_p=rnorm(1,0,1)
))}

n.adapt <- 500; n.update <- 5000; n.iter <- 10000
# n.adapt <- 100; n.update <- 500; n.iter <- 1000

# run chains in parallel
cl <- makeCluster(3) # this determines the number of chains, must be less than the number of cores on computer
clusterExport(cl, c('jd','n.adapt','n.update','n.iter','initFunc','nXcol','nXpoint'))

out <- clusterEvalQ(cl,{
  library(rjags)
  tmod <- jags.model("binom_hierarchical.jags", data= jd, n.chains=1, n.adapt=n.adapt,inits = initFunc())
  update(tmod, n.iter = n.update)
  zmCore <- coda.samples(tmod, c(
    "beta_colony_p","beta_point_p","B0_habitat_p",
    "mu_coast_p","pval.mean","y.new",'pval.sd','R2',
    'pred.interA_p.1','pred.interB_p.1','pred.interC_p.1',
    'pred.interA_p.2','pred.interB_p.2','pred.interC_p.2',
    'pred.interA_p.max.1','pred.interB_p.max.1','pred.interC_p.max.1',
    'pred.interA_p.max.2','pred.interB_p.max.2','pred.interC_p.max.2',
    'pred.interA_p.q.1','pred.interB_p.q.1','pred.interC_p.q.1',
    'pred.interA_p.q.2','pred.interB_p.q.2','pred.interC_p.q.2',
    'pred.inter.cat.1',
    'pred.inter.cat.2',
    'pval.cv','pval.fit','pval.min','pval.max'),n.iter=n.iter, n.thin=1)
  return(as.mcmc(zmCore))
})
zmPrp <- mcmc.list(out)
stopCluster(cl)

saveRDS(zmPrp, 'binom_acr_avgTotN_081519.Rdata')



# Pocillopora & avgdN15 ------------------------------------------------------
temp <- moorea[moorea$Taxa=="Pocillopora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgdN15),]
# temp <- temp[sample(nrow(temp),100),]
y <- temp$Percent_bleached
y <- y/100
y_p <- ifelse(y > 0, 1, 0)

# point
temp$point <- as.numeric(as.factor(as.character(temp$Point)))
point <- temp$point

# colony level predictors
x_colony <- as.matrix(data.frame(temp$Colony_size_class,temp$Depth_m))
for(i in 1:ncol(x_colony)) x_colony[,i] <- scale(x_colony[,i])[,1]
B <- ncol(x_colony)

# point level preds
site_preds <- temp %>% group_by(point) %>% summarise(cumtemp=unique(cumstress),totN=unique(avgdN15),habitat=unique(Habitat),coast=unique(Island_shore))
# summary(site_preds %>% group_by(point) %>% summarise(length(point)))
site_preds$cumtempXtotN <- site_preds$cumtemp*site_preds$totN

X_point <- as.matrix(data.frame(site_preds$totN,site_preds$cumtemp,site_preds$totN*site_preds$cumtemp))
for(i in 1:ncol(X_point)) X_point[,i] <- scale(X_point[,i])[,1]

# new data for predictions
totN_new <- as.matrix(data.frame('totN'=seq(from=min(X_point[,1]),to=max(X_point[,1]),length=50),
                                 'cumtemp'=rep(mean(X_point[,2]),length=50),
                                 'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

cumtemp_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                    'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                                    'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

inter_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                  'cumtemp'=rep(mean(X_point[,2]),length=50),
                                  'cumtempXtotN'=seq(from=min(X_point[,3]),to=max(X_point[,3]),length=50)))

# interaction as 3 categories defined by lower, middle, upper quartiles
interA <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25),50))*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))

interB <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),50)*seq(from=min(X_point[,2])
                                                                                       ,to=max(X_point[,2]),length=50)))

interC <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),50)*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))


interA.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interB.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interC.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interA.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interB.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interC.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

inter_min <- as.matrix(data.frame(
  # min heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(min(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(min(X_point[,2]),3)
))
inter_low <- as.matrix(data.frame(
  # lower q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.25),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.25),3)
))
inter_up <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.75),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.75),3)
))
inter_max <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(max(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(max(X_point[,2]),3)
))
interAll <- rbind(inter_min,inter_low,inter_up,inter_max)


jd <- list(B=B,
           y=y_p,
           n=length(y),
           X_colony=x_colony,
           point=point,
           P=length(unique(temp$point)),
           habitat=as.numeric(as.factor(as.character(site_preds$habitat))),
           X_point=X_point,
           C=ncol(X_point),
           K=2,
           coast=as.numeric(as.factor(as.character(site_preds$coast))),
           M=3,
           interA=interA,interB=interB,interC=interC,interAll=interAll
)

nXcol <- ncol(x_colony)
nXpoint <- ncol(X_point)

initFunc <- function(){return(list(
  beta_colony_p=rnorm(nXcol,0,1),
  beta_point_p=rnorm(nXpoint,-1,1),
  sigma_point_p=runif(1,0,20),
  sigma_habitat_p=runif(1,0,20),
  sigma_coast_p=runif(1,0,1),
  mu_overall_p=rnorm(1,0,1)
))}

n.adapt <- 500; n.update <- 5000; n.iter <- 10000
# n.adapt <- 100; n.update <- 500; n.iter <- 1000

# run chains in parallel
cl <- makeCluster(3) # this determines the number of chains, must be less than the number of cores on computer
clusterExport(cl, c('jd','n.adapt','n.update','n.iter','initFunc','nXcol','nXpoint'))

out <- clusterEvalQ(cl,{
  library(rjags)
  tmod <- jags.model("binom_hierarchical.jags", data= jd, n.chains=1, n.adapt=n.adapt,inits = initFunc())
  update(tmod, n.iter = n.update)
  zmCore <- coda.samples(tmod, c(
    "beta_colony_p","beta_point_p","B0_habitat_p",
    "mu_coast_p","pval.mean","y.new",'pval.sd','R2',
    'pred.interA_p.1','pred.interB_p.1','pred.interC_p.1',
    'pred.interA_p.2','pred.interB_p.2','pred.interC_p.2',
    'pred.interA_p.max.1','pred.interB_p.max.1','pred.interC_p.max.1',
    'pred.interA_p.max.2','pred.interB_p.max.2','pred.interC_p.max.2',
    'pred.interA_p.q.1','pred.interB_p.q.1','pred.interC_p.q.1',
    'pred.interA_p.q.2','pred.interB_p.q.2','pred.interC_p.q.2',
    'pred.inter.cat.1',
    'pred.inter.cat.2',
    'pval.cv','pval.fit','pval.min','pval.max'),n.iter=n.iter, n.thin=1)
  return(as.mcmc(zmCore))
})
zmPrp <- mcmc.list(out)
stopCluster(cl)

saveRDS(zmPrp, 'binom_poc_avgdN15_081519.Rdata')


# Acropora & avgdN15 ------------------------------------------------------
temp <- moorea[moorea$Taxa=="Acropora",]
temp <- temp[!is.na(temp$Percent_bleached),]
temp <- temp[!is.na(temp$avgdN15),]
# temp <- temp[sample(nrow(temp),100),]
y <- temp$Percent_bleached
y <- y/100
y_p <- ifelse(y > 0, 1, 0)

# point
temp$point <- as.numeric(as.factor(as.character(temp$Point)))
point <- temp$point

# colony level predictors
x_colony <- as.matrix(data.frame(temp$Colony_size_class,temp$Depth_m))
for(i in 1:ncol(x_colony)) x_colony[,i] <- scale(x_colony[,i])[,1]
B <- ncol(x_colony)

# point level preds
site_preds <- temp %>% group_by(point) %>% summarise(cumtemp=unique(cumstress),totN=unique(avgdN15),habitat=unique(Habitat),coast=unique(Island_shore))
# summary(site_preds %>% group_by(point) %>% summarise(length(point)))
site_preds$cumtempXtotN <- site_preds$cumtemp*site_preds$totN

X_point <- as.matrix(data.frame(site_preds$totN,site_preds$cumtemp,site_preds$totN*site_preds$cumtemp))
for(i in 1:ncol(X_point)) X_point[,i] <- scale(X_point[,i])[,1]

# new data for predictions
totN_new <- as.matrix(data.frame('totN'=seq(from=min(X_point[,1]),to=max(X_point[,1]),length=50),
                                 'cumtemp'=rep(mean(X_point[,2]),length=50),
                                 'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

cumtemp_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                    'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                                    'cumtempXtotN'=rep(mean(X_point[,3]),length=50)))

inter_new <- as.matrix(data.frame('totN'=rep(mean(X_point[,1]),50),
                                  'cumtemp'=rep(mean(X_point[,2]),length=50),
                                  'cumtempXtotN'=seq(from=min(X_point[,3]),to=max(X_point[,3]),length=50)))

# interaction as 3 categories defined by lower, middle, upper quartiles
interA <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25),50))*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))

interB <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),50)*seq(from=min(X_point[,2])
                                                                                       ,to=max(X_point[,2]),length=50)))

interC <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),50),
                               'cumtemp'=seq(from=min(X_point[,2]),to=max(X_point[,2]),length=50),
                               'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),50)*seq(from=min(X_point[,2])
                                                                                        ,to=max(X_point[,2]),length=50)))


interA.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interB.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interC.max <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                   'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                   'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(min(X_point[,2]),max(X_point[,2]))))

interA.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.25)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.25)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interB.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.5)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.5)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

interC.q <- as.matrix(data.frame('totN'=rep(quantile(X_point[,1],c(0.75)),2),
                                 'cumtemp'=c(min(X_point[,2]),max(X_point[,2])),
                                 'cumtempXtotN'=rep(quantile(X_point[,1],c(0.75)),2)*c(quantile(X_point[,2],0.25),quantile(X_point[,2],0.75))))

inter_min <- as.matrix(data.frame(
  # min heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(min(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(min(X_point[,2]),3)
))
inter_low <- as.matrix(data.frame(
  # lower q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.25),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.25),3)
))
inter_up <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(quantile(X_point[,2],0.75),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(quantile(X_point[,2],0.75),3)
))
inter_max <- as.matrix(data.frame(
  # upper q heat stress, low-med-high nuts
  'totN' = quantile(X_point[,1],c(0.25,0.5,0.75)),
  'cumtemp' = rep(max(X_point[,2]),3),
  'cumtempXtotN' = quantile(X_point[,1],c(0.25,0.5,0.75))*rep(max(X_point[,2]),3)
))
interAll <- rbind(inter_min,inter_low,inter_up,inter_max)


jd <- list(B=B,
           y=y_p,
           n=length(y),
           X_colony=x_colony,
           point=point,
           P=length(unique(temp$point)),
           habitat=as.numeric(as.factor(as.character(site_preds$habitat))),
           X_point=X_point,
           C=ncol(X_point),
           K=2,
           coast=as.numeric(as.factor(as.character(site_preds$coast))),
           M=3,
           interA=interA,interB=interB,interC=interC,interAll=interAll
)

nXcol <- ncol(x_colony)
nXpoint <- ncol(X_point)

initFunc <- function(){return(list(
  beta_colony_p=rnorm(nXcol,0,1),
  beta_point_p=rnorm(nXpoint,-1,1),
  sigma_point_p=runif(1,0,20),
  sigma_habitat_p=runif(1,0,20),
  sigma_coast_p=runif(1,0,1),
  mu_overall_p=rnorm(1,0,1)
))}

n.adapt <- 500; n.update <- 5000; n.iter <- 10000
# n.adapt <- 100; n.update <- 500; n.iter <- 1000

# run chains in parallel
cl <- makeCluster(3) # this determines the number of chains, must be less than the number of cores on computer
clusterExport(cl, c('jd','n.adapt','n.update','n.iter','initFunc','nXcol','nXpoint'))

out <- clusterEvalQ(cl,{
  library(rjags)
  tmod <- jags.model("binom_hierarchical.jags", data= jd, n.chains=1, n.adapt=n.adapt,inits = initFunc())
  update(tmod, n.iter = n.update)
  zmCore <- coda.samples(tmod, c(
    "beta_colony_p","beta_point_p","B0_habitat_p",
    "mu_coast_p","pval.mean","y.new",'pval.sd','R2',
    'pred.interA_p.1','pred.interB_p.1','pred.interC_p.1',
    'pred.interA_p.2','pred.interB_p.2','pred.interC_p.2',
    'pred.interA_p.max.1','pred.interB_p.max.1','pred.interC_p.max.1',
    'pred.interA_p.max.2','pred.interB_p.max.2','pred.interC_p.max.2',
    'pred.interA_p.q.1','pred.interB_p.q.1','pred.interC_p.q.1',
    'pred.interA_p.q.2','pred.interB_p.q.2','pred.interC_p.q.2',
    'pred.inter.cat.1',
    'pred.inter.cat.2',
    'pval.cv','pval.fit','pval.min','pval.max'),n.iter=n.iter, n.thin=1)
  return(as.mcmc(zmCore))
})
zmPrp <- mcmc.list(out)
stopCluster(cl)

saveRDS(zmPrp, 'binom_acr_avgdN15_081519.Rdata')





