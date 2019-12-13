
# -----------------------------------------------------------------------
# Nitrogen pollution interacts with thermal stress to increase coral bleaching across the seascape
# Donovan et al. 
# 1 - Temperature Patterns
# -----------------------------------------------------------------------

# This script analyzes temperature patterns across Moorea LTER sites, plots patterns across time and space, and calculates cumulative heat stress for use in modeling.

# Initialization ----------------------------------------------------------
rm(list=ls())
library(lubridate)
library(dplyr)

# Data --------------------------------------------------------------------
# data from Moorea Coral Reef LTER core time series, bottom mounted termistors:
# http://mcrlter.msi.ucsb.edu/cgi-bin/showDataset.cgi?docid=knb-lter-mcr.1035 # accessed December 3, 2018
# Leichter, J, K. Seydel and C. Gotschalk of Moorea Coral Reef LTER. 2018. MCR LTER: Coral Reef: Benthic Water Temperature, ongoing since 2005. knb-lter-mcr.1035.11

# temperature <- read.csv(file.path(getwd(),'data',paste('MCR_LTER','01','_BottomMountThermistors_20181024.csv',sep='')))
# for(i in c('02','03','04','05','06')){
#   temperature <- rbind(temperature,read.csv(file.path(getwd(),'data',paste('MCR_LTER',i,'_BottomMountThermistors_20181024.csv',sep=''))))
# }
# head(temperature)
# saveRDS(temperature,file='temperature_combined.Rdata')
temperature <- readRDS(file='data/temperature_combined.Rdata')

# Figure 1 - Temperature trends -------------------------------------------
# format time and date
temperature$time_use <- ymd_hms(temperature$time_local)
temperature$day <- format(temperature$time_use, '%Y-%m-%d')

# check metadata
temperature %>% filter(reef_type_code=='BAK') %>% distinct(site,reef_type_code,sensor_type,sensor_depth_m)
temperature %>% filter(reef_type_code=='FRI') %>% distinct(site,reef_type_code,sensor_type,sensor_depth_m)

# median temperature by day (and filter out forereef)
# temperature.day <- temperature %>% group_by(site,reef_type_code,sensor_type,sensor_depth_m,day) %>% filter(reef_type_code=='BAK' | reef_type_code=='FRI') %>% summarise(temp_c = median(temperature_c)) %>% ungroup()
# temperature.day$day <- ymd(temperature.day$day)
# saveRDS(temperature.day, file='data/temperature_day.Rdata')
temperature.day <- readRDS('data/temperature_day.Rdata')
temperature.day$day <- as.Date(temperature.day$day, '%Y-%m-%d')

samples <- temperature.day %>% group_by(site,reef_type_code,sensor_type,sensor_depth_m) %>% filter(reef_type_code=='BAK'| reef_type_code=='FRI') %>% summarise('start'=min(day),'end'=max(day), n=length(day)) %>% ungroup() %>% print(n=46)

# filter to include only time series we are interested in
samples <- filter(samples,  
                  # site== 'LTER01' & reef_type_code=='BAK' & sensor_depth_m==1 | ##missing data during sampling
                  site== 'LTER01' & reef_type_code=='FRI' & sensor_depth_m==1 |
                    site== 'LTER02' & reef_type_code=='BAK' & sensor_depth_m==1.5 |
                    site== 'LTER02' & reef_type_code=='FRI' & sensor_depth_m==0.900 |
                    site== 'LTER03' & reef_type_code=='BAK' & sensor_depth_m==1.8 |
                    site== 'LTER03' & reef_type_code=='FRI' & sensor_depth_m==1 |
                    site== 'LTER04' & reef_type_code=='BAK' & sensor_depth_m==1.8 |
                    site== 'LTER04' & reef_type_code=='FRI' & sensor_depth_m==0.800 |
                    # site== 'LTER05' & reef_type_code=='BAK' & sensor_depth_m==1.8 | ##missing data during sampling
                    site== 'LTER05' & reef_type_code=='FRI' & sensor_depth_m==0.600 |
                    site== 'LTER06' & reef_type_code=='BAK' & sensor_depth_m==2 |
                    site== 'LTER06' & reef_type_code=='FRI' & sensor_depth_m==1.2
)
samples %>% distinct(site, reef_type_code,sensor_depth_m)

# create an index for joining with data
samples$ind <- paste(samples$site, samples$reef_type_code, samples$sensor_depth_m, sep='_')
temperature.day$ind <- paste(temperature.day$site, temperature.day$reef_type_code, temperature.day$sensor_depth_m, sep='_')
temperature.day.sub <- temperature.day[temperature.day$ind %in% samples$ind,] # select only data corresponding to metadata in samples
temperature.day.sub$ind <- as.factor(temperature.day.sub$ind)
temperature.day.sub %>% distinct(site, reef_type_code,sensor_depth_m)

# subset date range so consistent across sites
temperature.day.sub <- filter(temperature.day.sub, day >= '2009-12-31' & day <= '2016-12-31') 

# check metadata again and plot
temperature.day.sub %>% group_by(site,reef_type_code,sensor_type,sensor_depth_m) %>% filter(reef_type_code=='BAK'| reef_type_code=='FRI') %>% summarise('start'=min(day),'end'=max(day), n=length(day)) %>% print(n=46)

# par(mfrow=c(5,2),mar=c(2,2,1,1),mgp=c(1,.5,0))
# for(i in 1:10){
#   plot(temperature.day.sub$day[temperature.day.sub$ind==levels(temperature.day.sub$ind)[i]],(temperature.day.sub$temp_c[temperature.day.sub$ind==levels(temperature.day.sub$ind)[i]]),ylab="",xlab="",ylim=c(24,30),main=levels(temperature.day.sub$ind)[i])
#   abline(v=ymd('2016-05-08'),col='blue',lwd=2)
# }

# note missing data in some cases, need to fill in NA rows for missing dates
## create a dataframe with the full date range
time.seq <- data.frame(day=unique(temperature.day.sub$day),ind=seq(1:length(unique(temperature.day.sub$day))))
time.seq$year <- year(time.seq$day)
time.seq$week.num <- week(time.seq$day)
time.seq.week <- time.seq %>% group_by(week.num,year) %>% summarise(day=min(day)) %>% ungroup()

## for each site/hab join with full date range and fill in NAs
alldays <- temperature.day.sub %>% dplyr::select(site,reef_type_code,sensor_depth_m,temp_c,day) %>% filter(site==samples$site[1] & reef_type_code==samples$reef_type_code[1] & sensor_depth_m==samples$sensor_depth_m[1])
alldays <- left_join(time.seq[c('day')],alldays,by='day')
alldays$site[is.na(alldays$site)] <- samples$site[1]
alldays$reef_type_code[is.na(alldays$reef_type_code)] <- samples$reef_type_code[1]
alldays$sensor_depth_m[is.na(alldays$sensor_depth_m)] <- samples$sensor_depth_m[1]
alldays <- alldays[with(alldays, order(day)),]
for(i in 2:10){
  temp <- temperature.day.sub %>% dplyr::select(site,reef_type_code,sensor_depth_m,temp_c,day) %>% filter(site==samples$site[i] & reef_type_code==samples$reef_type_code[i] & sensor_depth_m==samples$sensor_depth_m[i])
  temp <- left_join(time.seq[c('day')],temp,by='day')
  temp$site[is.na(temp$site)] <- samples$site[i]
  temp$reef_type_code[is.na(temp$reef_type_code)] <- samples$reef_type_code[i]
  temp$sensor_depth_m[is.na(temp$sensor_depth_m)] <- samples$sensor_depth_m[i]
  temp <- temp[with(temp, order(day)),]
  alldays <- rbind(alldays, temp)
}
str(alldays)

# check metadata again, should have same date range and number of rows
alldays %>% group_by(site,reef_type_code,sensor_depth_m) %>%  summarise('start'=min(day),'end'=max(day), n=length(day)) %>% ungroup() %>% print(n=46)

alldays$year <- year(alldays$day)
alldays$week.num <- week(alldays$day)

temperature.week <- alldays %>% group_by(site,reef_type_code,sensor_depth_m,year,week.num) %>% summarise(temp_c = mean(temp_c)) %>% ungroup()
temperature.week <- temperature.week[with(temperature.week, order(site,reef_type_code,sensor_depth_m,year,week.num)),]


# calculate maximum monthly temperature
# The threshold for thermal stress was defined as the median of maximum monthly summer temperatures [austral summer = Feb–Apr] 
# and accumulated heat stress (in °C-weeks) was calculated for all temperatures exceeding the MMM value.

temp <- temperature.day
temp$Month <- format(temperature.day$day, '%m')
temp$Year <- format(temperature.day$day, '%Y')
temp <- temp %>% group_by(Month,Year,site,reef_type_code) %>% summarise(m.max=max(temp_c)) %>% ungroup()
str(temp)
temp$Month <- as.numeric(temp$Month)
temp <- subset(temp, Month > 1 & Month < 5) # subset for austral summer
mma <- temp %>% summarise(mmm=median(m.max)) ; mma <- mma$mmm
mma

mma_ref <- 29

# calculate accumulated heat stress
temperature.week$hotspot <- temperature.week$temp_c - mma_ref
temperature.week$hotspot[temperature.week$hotspot < 0] <- 0

temperature.week$cumstress <- NA

# reset factor levels
temperature.week$reef_type_code <- as.character(temperature.week$reef_type_code)
temperature.week$reef_type_code <- as.factor(temperature.week$reef_type_code)

temperature.week$site <- as.character(temperature.week$site)
temperature.week$site <- as.factor(temperature.week$site)

samples$reef_type_code <- as.character(samples$reef_type_code)
samples$reef_type_code <- as.factor(samples$reef_type_code)

# 12 week running sum 
out <- temperature.week %>% filter(site==samples$site[1] & reef_type_code==samples$reef_type_code[1])
for(i in 13:nrow(out)){
  out$cumstress[i] <- sum(out$hotspot[(i-12):i],na.rm=T)
}

for(k in c(2:10)){
  temp <- temperature.week %>% filter(site==samples$site[k] & reef_type_code==samples$reef_type_code[k])
  for(i in 13:nrow(temp)){
    temp$cumstress[i] <- sum(temp$hotspot[(i-12):i],na.rm=F)
  }
  out <- rbind(out, temp)
}
str(out)
out <- left_join(out, time.seq.week, by=c('year','week.num'))


# plot
png(file='outputs/Figure1.png',height=3600,width=2400,res=300)
par(mgp=c(2,.7,0),oma=c(0,0,0,2),mfrow=c(2,1),mar=c(3.5,4,1,1))
# # temp <- subset(temperature.day, reef_type_code=='FRI')
# # temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
# # temp$Month <- format(temp$day, '%m')
# # temp$Year <- format(temp$day, '%Y')
# # temp$dayz <- format(temp$day, '%d')
# plot(mean_temp~day,data=temp[c(temp$Year==2016 &temp$reef_type_code=='FRI'),],type='n',
#      ylab=expression('Temperature '~degree~C),xlab='',cex.lab=1.5,cex.axis=1.5,ylim=c(26,30))
temp <- subset(temperature.day, reef_type_code=='FRI')
temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
temp$Month <- format(temp$day, '%m')
temp$Year <- format(temp$day, '%Y')
temp$dayz <- format(temp$day, '%d')
temp <- temp %>% group_by(Month,dayz) %>% summarise(mean_mean=mean(mean_temp), se=sd(mean_temp)/sqrt(length(mean_temp))) %>% ungroup()
temp$date <- as.Date(paste(temp$Month,temp$dayz,'2016',sep='-'),format='%m-%d-%Y')
plot(temp$date,temp$mean_mean,type='n',
     ylab=expression('Temperature '~degree~C),xlab='',cex.lab=1.8,cex.axis=1.5,ylim=c(26,30),bty='l')

text(as.Date('2016-01-05'),29.9,'A',cex=2.2)

abline(v=as.Date('2016-05-08'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-09'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-10'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-11'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-12'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-13'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-14'),col=rgb(190,190,190,200,max=255),lwd=1.8)

# temp <- subset(temperature.day, reef_type_code=='FRI')
# temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
# temp$Month <- format(temp$day, '%m')
# temp$Year <- format(temp$day, '%Y')
# temp$dayz <- format(temp$day, '%d')
# temp <- temp %>% group_by(Month,dayz) %>% summarise(mean_mean=mean(mean_temp), se=sd(mean_temp)/sqrt(length(mean_temp))) %>% ungroup()
# str(temp)
# temp$date <- as.Date(paste(temp$Month,temp$dayz,'2016',sep='-'),format='%m-%d-%Y')
points(temp$date,temp$mean_mean,type='l')
points(temp$date,temp$mean_mean+1.96*temp$se,type='l',lty=2)
points(temp$date,temp$mean_mean-1.96*temp$se,type='l',lty=2)

temp <- subset(temperature.day, reef_type_code=='FRI')
temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
temp$Month <- format(temp$day, '%m')
temp$Year <- format(temp$day, '%Y')
temp$dayz <- format(temp$day, '%d')
points(mean_temp~day,data=temp[c(temp$Year==2016 &temp$reef_type_code=='FRI'),],type='l',lwd=2,col='blue')
abline(h=mma_ref,lwd=3,col='red',lty=2)


temp <- temperature.day %>% group_by(site,reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
temp$Month <- format(temp$day, '%m')
temp$Year <- format(temp$day, '%Y')
temp$dayz <- format(temp$day, '%d')
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER02'),],type='l',lwd=2,col=rgb(0,255,0,100,max=255))
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER03'),],type='l',lwd=2,col=rgb(255,165,0,100,max=255))
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER04'),],type='l',lwd=2,col=rgb(255,0,255,100,max=255))
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER06'),],type='l',lwd=2,col=rgb(160,32,240,100,max=255))

st <- 255
pal <- c(rgb(190,125,216,st,max=255),
         rgb(141,177,71,st,max=255),
         rgb(96,163,222,st,max=255),
         rgb(213,144,71,st,max=255),
         rgb(86,185,135,st,max=255),
         rgb(229,112,134,st,max=255))

temp <- subset(out, year==2016)
plot(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER06'),],ylab='Cumulative Heat Stress',xlab='',type='n',cex.lab=1.8,cex.axis=1.5,bty='l')
text(as.Date('2016-01-05'),5,'B',cex=2.2)

abline(v=as.Date('2016-05-08'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-09'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-10'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-11'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-12'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-13'),col=rgb(190,190,190,200,max=255),lwd=1.8)
abline(v=as.Date('2016-05-14'),col=rgb(190,190,190,200,max=255),lwd=1.8)
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER01'),],type='l',lwd=4,col=pal[1])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER02'),],type='l',lwd=4,col=pal[2])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER03'),],type='l',lwd=4,col=pal[3])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER04'),],type='l',lwd=4,col=pal[4])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER05'),],type='l',lwd=4,col=pal[5])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER06'),],type='l',lwd=4,col=pal[6])
legend('topright',legend=c('LTER 1','LTER 2','LTER 3','LTER 4','LTER 5','LTER 6'),col=pal,lty=1,lwd=4,cex=1.5,bty='n')
dev.off()

# eps
setEPS()
postscript('outputs/Figure1.eps',height=12,width=7)
par(mgp=c(2,.7,0),oma=c(0,0,0,2),mfrow=c(2,1),mar=c(3.5,4,1,1))
# # temp <- subset(temperature.day, reef_type_code=='FRI')
# # temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
# # temp$Month <- format(temp$day, '%m')
# # temp$Year <- format(temp$day, '%Y')
# # temp$dayz <- format(temp$day, '%d')
# plot(mean_temp~day,data=temp[c(temp$Year==2016 &temp$reef_type_code=='FRI'),],type='n',
#      ylab=expression('Temperature '~degree~C),xlab='',cex.lab=1.5,cex.axis=1.5,ylim=c(26,30))
temp <- subset(temperature.day, reef_type_code=='FRI')
temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
temp$Month <- format(temp$day, '%m')
temp$Year <- format(temp$day, '%Y')
temp$dayz <- format(temp$day, '%d')
temp <- temp %>% group_by(Month,dayz) %>% summarise(mean_mean=mean(mean_temp), se=sd(mean_temp)/sqrt(length(mean_temp))) %>% ungroup()
temp$date <- as.Date(paste(temp$Month,temp$dayz,'2016',sep='-'),format='%m-%d-%Y')
plot(temp$date,temp$mean_mean,type='n',
     ylab=expression('Temperature '~degree~C),xlab='',cex.lab=1.7,cex.axis=1.5,ylim=c(26,30),bty='l')

text(as.Date('2016-01-05'),29.9,'A',cex=1.7)

abline(v=as.Date('2016-05-08'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-09'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-10'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-11'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-12'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-13'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-14'),col=rgb(190,190,190,255,max=255),lwd=1.8)

# temp <- subset(temperature.day, reef_type_code=='FRI')
# temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
# temp$Month <- format(temp$day, '%m')
# temp$Year <- format(temp$day, '%Y')
# temp$dayz <- format(temp$day, '%d')
# temp <- temp %>% group_by(Month,dayz) %>% summarise(mean_mean=mean(mean_temp), se=sd(mean_temp)/sqrt(length(mean_temp))) %>% ungroup()
# str(temp)
# temp$date <- as.Date(paste(temp$Month,temp$dayz,'2016',sep='-'),format='%m-%d-%Y')
points(temp$date,temp$mean_mean,type='l')
points(temp$date,temp$mean_mean+1.96*temp$se,type='l',lty=2)
points(temp$date,temp$mean_mean-1.96*temp$se,type='l',lty=2)

temp <- subset(temperature.day, reef_type_code=='FRI')
temp <- temp %>% group_by(reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
temp$Month <- format(temp$day, '%m')
temp$Year <- format(temp$day, '%Y')
temp$dayz <- format(temp$day, '%d')
points(mean_temp~day,data=temp[c(temp$Year==2016 &temp$reef_type_code=='FRI'),],type='l',lwd=2,col='blue')
abline(h=mma_ref,lwd=3,col='red',lty=2)


temp <- temperature.day %>% group_by(site,reef_type_code,day) %>% summarise('mean_temp'=mean(temp_c)) %>% ungroup()
temp$Month <- format(temp$day, '%m')
temp$Year <- format(temp$day, '%Y')
temp$dayz <- format(temp$day, '%d')
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER02'),],type='l',lwd=2,col=rgb(0,255,0,100,max=255))
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER03'),],type='l',lwd=2,col=rgb(255,165,0,100,max=255))
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER04'),],type='l',lwd=2,col=rgb(255,0,255,100,max=255))
# points(mean_temp~day,data=temp[c(temp$Year==2016 & temp$reef_type_code=='BAK' & temp$site=='LTER06'),],type='l',lwd=2,col=rgb(160,32,240,100,max=255))

st <- 255
pal <- c(rgb(190,125,216,st,max=255),
         rgb(141,177,71,st,max=255),
         rgb(96,163,222,st,max=255),
         rgb(213,144,71,st,max=255),
         rgb(86,185,135,st,max=255),
         rgb(229,112,134,st,max=255))

temp <- subset(out, year==2016)
plot(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER06'),],ylab='Cumulative Heat Stress',xlab='',type='n',cex.lab=1.7,cex.axis=1.5,bty='l')
text(as.Date('2016-01-05'),5,'B',cex=1.7)

abline(v=as.Date('2016-05-08'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-09'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-10'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-11'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-12'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-13'),col=rgb(190,190,190,255,max=255),lwd=1.8)
abline(v=as.Date('2016-05-14'),col=rgb(190,190,190,255,max=255),lwd=1.8)
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER01'),],type='l',lwd=4,col=pal[1])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER02'),],type='l',lwd=4,col=pal[2])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER03'),],type='l',lwd=4,col=pal[3])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER04'),],type='l',lwd=4,col=pal[4])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER05'),],type='l',lwd=4,col=pal[5])
points(cumstress~day,data=temp[c(temp$reef_type_code=='FRI' & temp$site=='LTER06'),],type='l',lwd=4,col=pal[6])
legend('topright',legend=c('LTER 1','LTER 2','LTER 3','LTER 4','LTER 5','LTER 6'),col=pal,lty=1,lwd=4,cex=1.3,bty='n')
dev.off()


# summarise temperature trends for results
out %>% filter(year == 2016 & cumstress > 0 & day < as.Date('2016-07-30')) %>% summarise(n_days=length(unique(day)),start=min(day),end=max(day)) %>% mutate(extent = end - start)
out %>% filter(year == 2016 & cumstress > 0) %>% group_by(site) %>% summarise(min(day)) %>% ungroup()
out %>% filter(year == 2016 & cumstress > 0 & day < as.Date('2016-07-30')) %>% group_by(site) %>% summarise(max(day)) %>% ungroup()

temp <- temperature.day
temp$year <- format(temp$day, '%Y')
temp %>% filter(year==2016) %>% group_by(reef_type_code,day) %>% summarise(mean_temp=mean(temp_c)) %>% filter(mean_temp > 29 & day < as.Date('2016-07-30')) %>% summarise(n_days=length(unique(day)),start=min(day),end=max(day)) %>% mutate(extent = end - start)

# join cumstress with data
out.may13 <- out[out$day==ymd('2016-05-13'),]
out.may13$LTER_site <- NA
out.may13$LTER_site[out.may13$site=='LTER01'] <- 'LTER_1'
out.may13$LTER_site[out.may13$site=='LTER02'] <- 'LTER_2'
out.may13$LTER_site[out.may13$site=='LTER03'] <- 'LTER_3'
out.may13$LTER_site[out.may13$site=='LTER04'] <- 'LTER_4'
out.may13$LTER_site[out.may13$site=='LTER05'] <- 'LTER_5'
out.may13$LTER_site[out.may13$site=='LTER06'] <- 'LTER_6'
out.may13$Habitat <- NA
out.may13$Habitat[out.may13$reef_type_code=='BAK'] <- 'Lagoon'
out.may13$Habitat[out.may13$reef_type_code=='FRI'] <- 'Fringe'

out.may13 <- rbind(out.may13,out.may13[2,])
out.may13[11,1] <- 'LTER01'
out.may13[11,10] <- 'LTER_1'
out.may13 <- rbind(out.may13,out.may13[9,])
out.may13[12,1] <- 'LTER05'
out.may13[12,10] <- 'LTER_5'


# join cum stress with data and export ------------------------------------

moorea <- read.csv('data/Poc_Acrop_bleaching_all_data.csv')

moorea <- left_join(moorea,out.may13[c('LTER_site','Habitat','cumstress')],by=c('LTER_site','Habitat'))
str(moorea)

write.csv(moorea,'data/Poc_Acrop_bleaching_all_data_wCumHeat.csv')

# supplemental figures ----------------------------------------------------

# Figure S5
temp <- temperature.day
temp$Month <- format(temperature.day$day, '%m')
temp$Year <- format(temperature.day$day, '%Y')
temp <- temp %>% group_by(Month,Year,site,reef_type_code) %>% summarise(m.max=max(temp_c)) %>% ungroup()
str(temp)
temp$Month <- as.numeric(temp$Month)
temp <- subset(temp, Month > 1 & Month < 5) # subset for austral summer
mma_site <- temp %>% group_by(site,reef_type_code) %>% summarise(mmm=median(m.max))
mma_site
mma_site$LTER_site <- 'LTER_1'
mma_site$LTER_site[mma_site$site=='LTER02'] <- 'LTER_2'
mma_site$LTER_site[mma_site$site=='LTER03'] <- 'LTER_3'
mma_site$LTER_site[mma_site$site=='LTER04'] <- 'LTER_4'
mma_site$LTER_site[mma_site$site=='LTER05'] <- 'LTER_5'
mma_site$LTER_site[mma_site$site=='LTER06'] <- 'LTER_6'
mma_site$Habitat <- as.character('Lagoon')
mma_site$Habitat[mma_site$reef_type_code=='FRI'] <- 'Fringe'

mma_site_adj <- mma_site
mma_site_adj$mmm[mma_site$LTER_site=='LTER_1' & mma_site$Habitat=='Lagoon'] <- mma_site$mmm[mma_site$LTER_site=='LTER_2' & mma_site$Habitat=='Lagoon']
mma_site_adj$mmm[mma_site$LTER_site=='LTER_5' & mma_site$Habitat=='Lagoon'] <- mma_site$mmm[mma_site$LTER_site=='LTER_6' & mma_site$Habitat=='Lagoon']

temp <- moorea %>% distinct(LTER_site,Habitat,avgTotN)
temp <- left_join(temp,mma_site_adj)
png(file='outputs/FigureS5.png',height=2000,width=2400,res=300)
plot(temp$avgTotN,temp$mmm,xlab='Total N',ylab='Mean Maximum Monthly Temperature',pch=19)
text(0.63,29.68,round(cor(temp$avgTotN,temp$mmm,use="complete.obs"),2),font=2)
dev.off()

# Figure S6
temp <- temperature.day
temp$Month <- format(temperature.day$day, '%m')
temp$Year <- format(temperature.day$day, '%Y')
temp <- temp %>% group_by(site,reef_type_code) %>% summarise(mean=mean(temp_c)) %>% ungroup()
str(temp)
# temp$Month <- as.numeric(temp$Month)
# temp <- subset(temp, Month > 1 & Month < 5) # subset for austral summer
# mma_site <- temp %>% group_by(site,reef_type_code) %>% summarise(mmm=median(m.max))
temp_mean_site <- temp
temp_mean_site$LTER_site <- 'LTER_1'
temp_mean_site$LTER_site[temp_mean_site$site=='LTER02'] <- 'LTER_2'
temp_mean_site$LTER_site[temp_mean_site$site=='LTER03'] <- 'LTER_3'
temp_mean_site$LTER_site[temp_mean_site$site=='LTER04'] <- 'LTER_4'
temp_mean_site$LTER_site[temp_mean_site$site=='LTER05'] <- 'LTER_5'
temp_mean_site$LTER_site[temp_mean_site$site=='LTER06'] <- 'LTER_6'
temp_mean_site$Habitat <- as.character('Lagoon')
temp_mean_site$Habitat[temp_mean_site$reef_type_code=='FRI'] <- 'Fringe'


temp <- temperature.day
temp$Month <- format(temperature.day$day, '%m')
temp$Year <- format(temperature.day$day, '%Y')
temp <- temp %>% group_by(site,reef_type_code) %>% summarise(sd=sd(temp_c)) %>% ungroup()
str(temp)
# temp$Month <- as.numeric(temp$Month)
# temp <- subset(temp, Month > 1 & Month < 5) # subset for austral summer
# mma_site <- temp %>% group_by(site,reef_type_code) %>% summarise(mmm=median(m.max))
temp_sd_site <- temp
temp_sd_site$LTER_site <- 'LTER_1'
temp_sd_site$LTER_site[temp_sd_site$site=='LTER02'] <- 'LTER_2'
temp_sd_site$LTER_site[temp_sd_site$site=='LTER03'] <- 'LTER_3'
temp_sd_site$LTER_site[temp_sd_site$site=='LTER04'] <- 'LTER_4'
temp_sd_site$LTER_site[temp_sd_site$site=='LTER05'] <- 'LTER_5'
temp_sd_site$LTER_site[temp_sd_site$site=='LTER06'] <- 'LTER_6'
temp_sd_site$Habitat <- as.character('Lagoon')
temp_sd_site$Habitat[temp_sd_site$reef_type_code=='FRI'] <- 'Fringe'

temp_mean_site$mean[temp_mean_site$LTER_site=='LTER_1' & temp_mean_site$Habitat=='Lagoon'] <- temp_mean_site$mean[temp_mean_site$LTER_site=='LTER_2' & temp_mean_site$Habitat=='Lagoon']
temp_mean_site$mean[temp_mean_site$LTER_site=='LTER_5' & temp_mean_site$Habitat=='Lagoon'] <- temp_mean_site$mean[temp_mean_site$LTER_site=='LTER_6' & temp_mean_site$Habitat=='Lagoon']
temp_sd_site$sd[temp_sd_site$LTER_site=='LTER_1' & temp_sd_site$Habitat=='Lagoon'] <- temp_sd_site$sd[temp_sd_site$LTER_site=='LTER_2' & temp_sd_site$Habitat=='Lagoon']
temp_sd_site$sd[temp_sd_site$LTER_site=='LTER_5' & temp_sd_site$Habitat=='Lagoon'] <- temp_sd_site$sd[temp_sd_site$LTER_site=='LTER_6' & temp_sd_site$Habitat=='Lagoon']

temp <- moorea %>% distinct(LTER_site,Habitat,avgTotN)
temp <- left_join(temp,temp_mean_site)
temp <- left_join(temp,temp_sd_site)

png(file='outputs/FigureS6.png',height=4000,width=2400,res=300)
par(mfrow=c(2,1),mar=c(4,4,2,1),oma=c(2,0,0,0),xpd=T)

plot(temp$avgTotN,temp$mean,xlab='Total N',ylab='Mean Temperature',pch=19)
text(0.63,27.92,round(cor(temp$avgTotN,temp$mean,use="complete.obs"),2),font=2)

plot(temp$avgTotN,temp$sd,xlab='Total N',ylab='SD Temperature',pch=19)
text(0.63,1.12,round(cor(temp$avgTotN,temp$sd,use="complete.obs"),2),font=2)
dev.off()

