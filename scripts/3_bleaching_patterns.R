# -----------------------------------------------------------------------
# Nutrient pollution alters coral bleaching across the seascape
# Donovan et al. 
# 3 - Bleaching patterns summary
# -----------------------------------------------------------------------

# This script summarizes bleaching patterns and creates rasters for mapping in Figure 3.

# Initialization ----------------------------------------------------------
library(dplyr)
library(rjags)
library(boot)

# data --------------------------------------------------------------------
# bleaching
moorea <- read.csv('data/moorea_withavgnuts.csv')

# # lidar
# m.lidar <- read.table('data/bathy_LIDAR_Riegl_820_05m.xyz')
# head(m.lidar)
# colnames(m.lidar) <- c('longitude','latitude','depth')
# 
# # summarize ---------------------------------------------------------------
poc_prev <- moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached)) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanPrev'=mean(y)) %>% ungroup()
hist(poc_prev$meanPrev)
summary(poc_prev$meanPrev)
moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached)) %>% group_by(Island_shore,Point,Latitude,Longitude) %>% summarise('meanPrev'=mean(y)) %>% ungroup() %>% group_by(Island_shore) %>% summarise('mean'=mean(meanPrev),'max'=max(meanPrev),'min'=min(meanPrev)) %>% ungroup()

acr_prev <- moorea %>% filter(Taxa=="Acropora" & !is.na(Percent_bleached)) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanPrev'=mean(y)) %>% ungroup()
hist(acr_prev$meanPrev)

poc_sev <- moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached) & y==1) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup()
summary(poc_sev$meanSev)
moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached) & y==1) %>% group_by(Island_shore,Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup() %>% group_by(Island_shore) %>% summarise('median'=median(meanSev),'max'=max(meanSev),'min'=min(meanSev)) %>% ungroup()

acr_sev <- moorea %>% filter(Taxa=="Acropora" & !is.na(Percent_bleached) & y==1) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup()
hist(acr_sev$meanSev)
summary(acr_sev$meanSev)
moorea %>% filter(Taxa=="Acropora" & !is.na(Percent_bleached) & y==1) %>% group_by(Island_shore,Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup() %>% group_by(Island_shore) %>% summarise('median'=median(meanSev),'max'=max(meanSev),'min'=min(meanSev)) %>% ungroup()
# 
# 
# # create mask -------------------------------------------------------------
# # convert LiDAR data to a regular grid for drawing contours
# pts <- m.lidar
# colnames(pts) <- c('x','y','z')
# coordinates(pts) <- ~x+y # create a SpatialPointsDataFrame									   
# rast <- raster(ext=extent(pts), resolution=.00005) # create an empty raster object to the extent of the points
# rasOut <- rasterize(pts, rast, pts$z, fun = mean) # rasterize to regularly grid 
# plot(rasOut)
# summary(rasOut)
# 
# # make a raster with the values we want being 1 and the rest being NA
# ras.use <- rasOut # save a copy
# ras.use[ras.use > 0] <- 1 # all depths in the water are 1
# ras.use[ras.use < 0] <- NA # everything else NA
# plot(ras.use)
# class(ras.use) #RasterLayer
# # pol.use <- rasterToPolygons(ras.use,dissolve=T,fun=function(x){x==1}) #convert to polygons- super slow 
# # writeOGR(obj=pol.use, dsn="data", layer="pol_mask", driver="ESRI Shapefile") #write out output so don't have to run above again
# pol.go <- readOGR('data/pol_mask.shp') # read back in
# class(pol.go) # check class preserved in read out/in
# plot(pol.go)
# 
# # create spatial points
# xy <- SpatialPointsDataFrame(matrix(c(poc_prev$Longitude,poc_prev$Latitude),ncol=2),data=poc_prev,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# points(xy)
# 
# # get the min/max range for lat/long to make an empty grid 
# x.range <- as.numeric(c(min(poc_prev$Longitude)-0.01, max(poc_prev$Longitude)+0.01))  # min/max longitude of the interpolation area
# y.range <- as.numeric(c(min(poc_prev$Latitude)-0.01, max(poc_prev$Latitude)+0.01))  # min/max latitude of the interpolation area  
# # from the range, exapnd the coordinates to make a regular grid
# grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.0005), 
#                    y = seq(from = y.range[1], to = y.range[2], by = 0.0005))  # expand 
# coordinates(grd) <- ~x + y
# sp::gridded(grd) <- TRUE
# proj4string(grd) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
# # quick plot to see what we have
# plot(grd, cex = 1.5, col = "grey")
# points(xy, pch = 1, col = "red", cex = 1)
# 
# #clip to extent of LiDAR
# grd.mask.p <- crop(grd,pol.go) # slow
# plot(grd.mask.p)
# # convert to correct data type
# grd.mask.r <- as.data.frame(grd.mask.p)
# grd.mask.r <- SpatialPointsDataFrame(grd.mask.r,data=data.frame(z=rep(1,length(grd.mask.r$x)),t=rep(2,length(grd.mask.r$x))))
# proj4string(grd.mask.r) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
# grd.mask.rs <- SpatialPixels(grd.mask.r)
# plot(grd.mask.rs)
# 
# # interpolate -------------------------------------------------------------
# 
# xy <- SpatialPointsDataFrame(matrix(c(poc_sev$Longitude,poc_sev$Latitude),ncol=2),data=poc_sev,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# interp_poc_s <- idw(formula=poc_sev$meanSev~1,locations = xy, newdata = grd.mask.rs)
# spplot(interp_poc_s['var1.pred']) # quick plot
# 
# xy <- SpatialPointsDataFrame(matrix(c(poc_prev$Longitude,poc_prev$Latitude),ncol=2),data=poc_prev,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# interp_poc_p <- idw(formula=poc_prev$meanPrev~1,locations = xy, newdata = grd.mask.rs)
# spplot(interp_poc_p['var1.pred']) # quick plot
# 
# xy <- SpatialPointsDataFrame(matrix(c(acr_sev$Longitude,acr_sev$Latitude),ncol=2),data=acr_sev,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# interp_acr_s <- idw(formula=acr_sev$meanSev~1,locations = xy, newdata = grd.mask.rs)
# spplot(interp_acr_s['var1.pred']) # quick plot
# 
# xy <- SpatialPointsDataFrame(matrix(c(acr_prev$Longitude,acr_prev$Latitude),ncol=2),data=acr_prev,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# interp_acr_p <- idw(formula=acr_prev$meanPrev~1,locations = xy, newdata = grd.mask.rs)
# spplot(interp_acr_p['var1.pred']) # quick plot
# 
# 
# # export ------------------------------------------------------------------
# interp_poc_p_r <- raster(interp_poc_p)
# interp_poc_s_r <- raster(interp_poc_s)
# interp_acr_p_r <- raster(interp_acr_p)
# interp_acr_s_r <- raster(interp_acr_s)
# 
# writeRaster(interp_poc_p_r, filename='data/Pocillopora_prevalence.tif', format='GTiff')
# writeRaster(interp_poc_s_r, filename='data/Pocillopora_severity.tif', format='GTiff')
# writeRaster(interp_acr_p_r, filename='data/Acropora_prevalence.tif', format='GTiff')
# writeRaster(interp_acr_s_r, filename='data/Acropora_severity.tif', format='GTiff')
# 

# bayesian null model -----------------------------------------------------

cat("model{
  # priors
    B0_p ~ dnorm(0, 0.001)
    t0 ~ dnorm(0, 0.01) 
    tau <- exp(t0)
    
    # likelihood for alpha
    for (i in 1:n){
      logit(alpha[i]) <- a[point[i]]
      y.discrete[i] ~ dbern(alpha[i])
    }
    
    # likelihood for gamma
    for (i in 1:n.discrete){
      y.d[i] ~ dbern(mu[i])
      logit(mu[i]) <- B0_p
    }
    
    # likelihood for mu and tau
    for (i in 1:n.cont){
      y.c[i] ~ dbeta(p[i], q[i])
      p[i] <- mu2[i] * tau
      q[i] <- (1 - mu2[i]) * tau
      logit(mu2[i]) <- B0[point.c[i]]
    }
    
    # Point - prevelence
    for(p in 1:P){
      a[p] ~ dnorm(mu_p,tau_mu_p)
      B0[p] ~ dnorm(mu_s,tau_mu_s)
    }
    
    sigma_mu_p ~ dunif(0,50)
    tau_mu_p <- pow(sigma_mu_s,-2) 
    sigma_mu_s ~ dunif(0,50)
    tau_mu_s <- pow(sigma_mu_s,-2) 
    
    # overall 
    mu_p ~ dnorm(0,0.01)
    mu_s ~ dnorm(0,0.01)

    }  ", file="beinf_nox.txt")


# pocillopora -------------------------------------------------------------

# moorea.poc <- subset(moorea, Taxa=='Pocillopora')
# temp <- moorea.poc[!is.na(moorea.poc$Percent_bleached),]
# y <- temp$Percent_bleached
# y <- y/100
# 
# # point
# temp$point <- as.numeric(as.factor(as.character(temp$Point)))
# point <- temp$point
# 
# # split the data into discrete and continuous components
# y.d <- ifelse(y == 1 | y == 0, y, NA)
# y.discrete <- ifelse(is.na(y.d), 0, 1)
# y.d <- y.d[!is.na(y.d)]
# n.discrete <- length(y.d)
# point.d <- point[y.discrete == 1]
# 
# which.cont <- which(y < 1 & y > 0)
# y.c <- ifelse(y < 1 & y > 0, y, NA)
# y.c <- y.c[!is.na(y.c)]
# n.cont <- length(y.c)
# point.c <- point[which.cont]
# 
# 
# jd <- list(n=length(y),
#            y.discrete=temp$y,
#            n.discrete=n.discrete,
#            y.d=y.d, 
#            point=temp$point,
#            n.cont=n.cont,
#            y.c=y.c, 
#            point.c=point.c,
#            P=length(unique(temp$point))
# )
# 
# initFunc <- function(){return(list(
#   # a=rnorm(1,0,1),
#   t0=rnorm(1,0,1),
#   sigma_mu_p=runif(1,0,1),
#   mu_p=rnorm(1,0,1),
#   sigma_mu_s=runif(1,0,1),
#   mu_s=rnorm(1,0,1)
# ))}
# 
# n.adapt <- 500; n.update <- 1500; n.iter <- 10000
# 
# mod <- jags.model("beinf_nox.txt", data= jd, n.chains=3, n.adapt=n.adapt,inits = initFunc())
# update(mod, n.iter = n.update)
# poc.prev.nox <- coda.samples(mod, c("a"),n.iter=n.iter, n.thin=1)
# # gelman.diag(poc.prev.nox)# check for convergence
# # summary(poc.prev.nox)
# 
# poc.b0.p <- rep(NA,length(unique(moorea.poc$Point)))
# grepgo <- grep('a',colnames(poc.prev.nox[[1]])); poc.b0.p <- summary(poc.prev.nox)$quantiles[grepgo,'50%']
# hist(poc.b0.p)
# hist(inv.logit(poc.b0.p))
# 
# poc.sev.nox <- coda.samples(mod, c("B0"),n.iter=n.iter, n.thin=1)
# # gelman.diag(poc.prev.nox)# check for convergence
# # summary(poc.prev.nox)
# 
# poc.b0.s <- rep(NA,length(unique(moorea.poc$Point)))
# grepgo <- grep('B0',colnames(poc.sev.nox[[1]])); poc.b0.s <- summary(poc.sev.nox)$quantiles[grepgo,'50%']
# hist(poc.b0.s)
# hist(inv.logit(poc.b0.s))
# 
# str(temp)
# temp <- temp %>% group_by(point) %>% summarise('meanPrev'=mean(y),'meanSev'=mean(Percent_bleached/100))
# plot(temp$meanPrev,inv.logit(poc.b0.p)); abline(a=0,b=1)
# plot(temp$meanSev,inv.logit(poc.b0.s)); abline(a=0,b=1)
# 
# 
# cat("model{
#     # likelihood for alpha
#     for (i in 1:n){
#       logit(alpha[i]) <- a[point[i]]
#       y.discrete[i] ~ dbern(alpha[i])
#     }
#     
#     # Point - prevelence
#     for(p in 1:P){
#     a[p] ~ dnorm(mu,tau)
#     }
#     
#     sigma ~ dunif(0,50)
#     tau <- pow(sigma,-2) 
#     mu ~ dnorm(0,0.01)
#     
#     }  ", file="simp_logistic.txt")
# 
# moorea.poc <- subset(moorea, Taxa=='Pocillopora')
# temp <- moorea.poc[!is.na(moorea.poc$Percent_bleached),]
# y <- temp$Percent_bleached
# y <- y/100
# 
# # point
# temp$point <- as.numeric(as.factor(as.character(temp$Point)))
# point <- temp$point
# 
# # split the data into discrete and continuous components
# y.d <- ifelse(y == 1 | y == 0, y, NA)
# y.discrete <- ifelse(is.na(y.d), 0, 1)
# 
# jd <- list(n=length(y),
#            y.discrete=temp$y,
#            point=temp$point,
#            P=length(unique(temp$point))
# )
# 
# initFunc <- function(){return(list(
#   sigma=runif(1,0,1),
#   mu=rnorm(1,0,1)
# ))}
# 
# n.adapt <- 500; n.update <- 1500; n.iter <- 10000
# 
# mod <- jags.model("simp_logistic.txt", data= jd, n.chains=3, n.adapt=n.adapt,inits = initFunc())
# update(mod, n.iter = n.update)
# poc.prev.simplog <- coda.samples(mod, c("a"),n.iter=n.iter, n.thin=1)
# 
# poc.a <- rep(NA,length(unique(moorea.poc$Point)))
# grepgo <- grep('a',colnames(poc.prev.simplog[[1]])); poc.a <- summary(poc.prev.simplog)$quantiles[grepgo,'50%']
# hist(poc.a)
# hist(inv.logit(poc.a))
# 
# temp <- temp %>% group_by(point) %>% summarise('meanPrev'=mean(y),'meanSev'=mean(Percent_bleached/100))
# plot(temp$meanPrev,inv.logit(poc.a)); abline(a=0,b=1)
# 
# plot(inv.logit(poc.b0.p),inv.logit(poc.a)); abline(a=0,b=1)
# 
# library(lme4)
# moorea.poc <- subset(moorea, Taxa=='Pocillopora')
# temp <- moorea.poc[!is.na(moorea.poc$Percent_bleached),]
# y <- temp$Percent_bleached
# y <- y/100
# 
# # point
# temp$point <- as.numeric(as.factor(as.character(temp$Point)))
# point <- temp$point
# 
# # split the data into discrete and continuous components
# y.d <- ifelse(y == 1 | y == 0, y, NA)
# y.discrete <- ifelse(is.na(y.d), 0, 1)
# 
# poc.prev.glm <- glmer(y.discrete~1+(1|temp$point),family='binomial')
# summary(poc.prev.glm)
# plot(coef(poc.prev.glm)$`temp$point`[,1],poc.a)
# 
# oldrez <- read.csv('/Users/mary/github/moorea_bleach_bayes/outputs/temp.csv')
# plot(coef(poc.prev.glm)$`temp$point`[,1],oldrez$x)
# plot(poc.a,oldrez$x)
