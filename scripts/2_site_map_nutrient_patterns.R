# -----------------------------------------------------------------------
# Nutrient pollution alters coral bleaching across the seascape
# Donovan et al. 
# 2 - Site map and nutrient patterns
# -----------------------------------------------------------------------

# This script produces a site map with a continuous surface for nutrient data (using kriging).

# Initialization ----------------------------------------------------------
rm(list=ls())
library(raster)
library(rgdal)
library(dplyr)
library(sf)
library(kriging)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(rgeos)

  # Data --------------------------------------------------------------------
  # base layers
  dem <- raster('data/dtm_merged_5m.tif')
  dem_land <- raster("data/dem_land.tif")
  mo.shp <- readOGR('data/Moorea_outline/Isle_outline_gcs84.shp')
  mo.shp.dem <- spTransform(mo.shp,projection(dem_land)) # transform to same projection as dem
  mo.sf <- st_as_sf(mo.shp.dem) # save as an sf object
  
  # sample points
  lter.sites <- read.csv('data/LTER_Sites_LatLon.csv')
  moorea <- read.csv('data/Poc_Acrop_bleaching_all_data_wCumHeat.csv')
  
  # nutrient data
  may_nuts_all <- read.csv('data/UGA_May_isotope_data.csv')
  may_nuts_all$may_totN <- may_nuts_all$Total..N
  may_nuts_all$Point <- may_nuts_all$Sample.ID
  
  jan_nuts_all <- read.csv('data/UGA_Jan_isotope_data.csv')
  jan_nuts_all$jan_totN <- jan_nuts_all$Total..N
  
  all_nuts_all <- full_join(may_nuts_all[c('Point','may_totN')],jan_nuts_all[c('Point','jan_totN')],by='Point')
  all_nuts_all <- all_nuts_all[!all_nuts_all$Point=='*6',]
  all_nuts_all <- all_nuts_all[!all_nuts_all$Point=='*26',]
  all_nuts_all$avgTotN <- (all_nuts_all$may_totN + all_nuts_all$jan_totN)/2
  nrow(all_nuts_all[!is.na(all_nuts_all$avgTotN),])
  all_nuts_all <- all_nuts_all[!is.na(all_nuts_all$avgTotN),]
  
  nut_latlong <- read.csv('data/Turbinaria_sample_points.csv')
  nut_latlong$Point <- nut_latlong$Name
  all_nuts_all <- left_join(all_nuts_all,nut_latlong[,c('Point','Latitude','Longitude')],by='Point')
  
  # land --------------------------------------------------------------------
  # hill shade
  slp <- terrain(dem_land,opt='slope')
  asp <- terrain(dem_land,opt='aspect')
  hill <- hillShade(slp,asp)# compute hillshade 
  plot(hill)
  
  # transform rasters for ggplot
  dem.p  <-  rasterToPoints(dem_land)
  dem.df <-  data.frame(dem.p)
  colnames(dem.df) = c("x", "y", "alt")
  
  hill.p  <-  rasterToPoints(hill)
  hill.df <-  data.frame(hill.p)
  colnames(hill.df) = c("x", "y", "alt")
  
  # nutrients ---------------------------------------------------------------
  # convert nutrient data to spatial object and transform projection
  xy <- SpatialPointsDataFrame(matrix(c(all_nuts_all$Longitude,all_nuts_all$Latitude),ncol=2),data=all_nuts_all,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  nutsUTM <- spTransform(xy, projection(mo.shp.dem))
  nutsUTM.df <- as.data.frame(nutsUTM)
  
  # set up boundary for kriging
  # using outline from Tom, sent May 2 2019
  ashapem <- read.csv("data/nut_boundary2.csv", header = T)
  ashapem <- as.matrix(ashapem[,2:3])
  
  borderpolygon <- list(data.frame(ashapem[,1], ashapem[,2]))# the borderpolygon object is a list that will be part of the kriging function
  
  # old code for creating a reef outline
  # # create a convex hull around survey points
  # # cutz <- chull(matrix(c(all_nuts_all$Longitude,all_nuts_all$Latitude),ncol=2))
  # # cutz <- matrix(c(all_nuts_all$Longitude,all_nuts_all$Latitude),ncol=2)[c(cutz, cutz[1]), ]
  # 
  # # use lidar extent as boundary
  # # use outline of lidar extent created in ArcGIS
  # reef.outline <- readOGR('data/mask_shapes/lidar_extent_pygon_filled.shp')
  # reef.outline.pj <- spTransform(reef.outline,projection(dem_land)) # transform to same projection as dem
  # plot(reef.outline.pj)
  # 
  # # plot to make sure things are lining up
  # plot(reef.outline.pj)
  # plot(mo.shp.dem,add=T)
  # plot(nutsUTM,add=T)
  # 
  # # get the min/max range for lat/long to make a box
  # x.range <- as.numeric(c(min(all_nuts_all$Longitude)-0.01, max(all_nuts_all$Longitude)+0.01))  # min/max longitude of the interpolation area
  # y.range <- as.numeric(c(min(all_nuts_all$Latitude)-0.01, max(all_nuts_all$Latitude)+0.01))  # min/max latitude of the interpolation area  
  # 
  # square <- cbind(x.range[1],y.range[2],  # NW corner
  #                 x.range[2], y.range[2],  # NE corner
  #                 x.range[2],y.range[1],  # SE corner
  #                 x.range[1],y.range[1], # SW corner
  #                 x.range[1],y.range[2])  # NW corner again - close ploygon
  # id=1
  # polys <- SpatialPolygons(mapply(function(poly, id) 
  # {
  #   xy <- matrix(poly, ncol=2, byrow=TRUE)
  #   Polygons(list(Polygon(xy)), ID=id)
  # }, 
  # split(square, row(square)), id),
  # proj4string=CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84")))
  # plot(polys)
  # polys.utm <- spTransform(polys,projection(mo.shp.dem))
  # 
  # ashapesp.t.df <- slot(polys.utm, "polygons")[[1]]
  # ashapesp.t.df <- slot(ashapesp.t.df, 'Polygons')
  # ashapesp.t.df.go <- slot(ashapesp.t.df[[1]], 'coords')
  # 
  # borderpolygon <- list(data.frame(ashapesp.t.df.go[,1], ashapesp.t.df.go[,2]))
  
  # kriging
  krig1 <- kriging(x=nutsUTM.df$coords.x1, y=nutsUTM.df$coords.x2, response=nutsUTM.df$avgTotN, pixels=1000, polygons=borderpolygon) 
  str(krig1)
  
  krig2 <- krig1$map
  str(krig2)
  # krigxy <- SpatialPoints(coords=krig2[,c("x","y")],proj4string=CRS(projection(mo.shp.dem)))
  # krigxy$pred <- krig2$pred
  # # plot(krigxy)
  # krigxy.crop <- krigxy[reef.outline.pj,] # clip the output to the reef outline polygon
  # plot(krigxy.crop)
  krigLLdf <- as.data.frame(krig2)
  
  ggplot() +
    geom_raster(data=krigLLdf, aes(x=x, y=y, fill=pred)) +
    scale_fill_gradientn(colours = rev(rainbow(10))) +
    coord_sf()
  
  
  # sample points -----------------------------------------------------------
  # bsites <- moorea %>% filter(!is.na(Total..N)) %>% distinct(Point,Longitude,Latitude,Habitat,Island_shore) %>% ungroup() %>% rename(Shoreline = Island_shore)
  spdf <- SpatialPointsDataFrame(coords = all_nuts_all[c('Longitude','Latitude')], data = all_nuts_all, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  ptsUTM <- spTransform(spdf, projection(mo.shp.dem))
  ptsUTM.df <- as.data.frame(ptsUTM)
  str(ptsUTM.df)
  # ptsUTM.df$Habitat.2 <- as.character(ptsUTM.df$Habitat)
  # ptsUTM.df$Habitat.2[ptsUTM.df$Habitat=='Lagoon'] <- 'Backreef'
  # ptsUTM.df$Habitat.2[ptsUTM.df$Habitat=='Fringe'] <- 'Fringing Reef'
  # ptsUTM.df$Habitat <- ptsUTM.df$Habitat.2
  
  
  lsites <- subset(lter.sites, Habitat=='Back Reef'|Habitat=='Fringing Reef')
  spdf <- SpatialPointsDataFrame(coords = lsites[c('long','lat')], data = lsites, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  lterUTM <- spTransform(spdf, projection(mo.shp.dem))
  lterUTM.df <- as.data.frame(lterUTM)
  lterUTM.df.back <- subset(lterUTM.df, Habitat=='Back Reef')
  
  lterUTM.df.med <- lterUTM.df %>% group_by(Site) %>% summarise(lat=median(lat),long=median(long)) %>% ungroup()
  spdf <- SpatialPointsDataFrame(coords = lterUTM.df.med[c('long','lat')], data = lterUTM.df.med, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  lterUTM <- spTransform(spdf, projection(mo.shp.dem))
  lterUTM.df <- as.data.frame(lterUTM)
  
  # plot --------------------------------------------------------------------
  
  png(file='outputs/Figure2.png',height=2000,width=2400,res=300)
  ggplot() +
    theme_bw() +   
    # geom_raster(data=krigLLdf, aes(x=x, y=y, fill=pred),show.legend=FALSE) +
    # scale_fill_gradientn(colours = rev(rainbow(10))) +
    geom_point(data=krig2, aes(x=x, y=y, colour=pred), size=4) + 
    scale_colour_gradientn(name="% Nitrogen",colours = rev(rainbow(10))) + 
    theme(legend.position= c(0.9, 0.2),legend.text=element_text(size=12)) +
    
    new_scale("fill") +
    geom_raster(data=hill.df, aes(x=x,y=y,fill = alt),show.legend=FALSE) + scale_fill_gradientn(colors = grey(20:100/100)) +
    new_scale("fill") +
    # geom_raster(data=dem.df, aes(x=x,y=y,fill = log(alt+.45)), alpha=0.30,show.legend=FALSE) + scale_fill_gradientn(colors = terrain.colors(100)) +
    geom_raster(data=dem.df, aes(x=x,y=y,fill = log(alt+.45)), alpha=0.30,show.legend=FALSE) + scale_fill_gradientn(colors = grey(1:100/100)) +
    geom_sf(data = mo.sf, fill=NA) +
    geom_point(data=lterUTM.df, aes(long.1,lat.1),size=18,shape=22,stroke=1.15) + 
    new_scale("fill") +
    geom_point(data=ptsUTM.df, aes(Longitude.1,Latitude.1),size=1) +
    # geom_point(data=ptsUTM.df, aes(Longitude.1,Latitude.1,shape=Habitat),size=1.5) + scale_shape_manual(values=c(1,19)) +
    # geom_text_repel(data=lterUTM.df.back, aes(long.1,lat.1),label=lterUTM.df.back$Site) +
    coord_sf() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    # theme(panel.grid.major = element_blank(),
    #       plot.background=element_rect(fill='white')) + 
    theme(panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5),
          plot.background=element_rect(fill='white')) +
    # theme(legend.position = 'bottom',
    #       legend.justification='center',
    #       legend.key=element_rect(colour="white"),
    #       legend.text=element_text(size=16),
    #       legend.title=element_text(size=18)) +
    annotation_scale(location = "bl", width_hint = 0.25,text_cex=1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) +
    annotate("text", x = 198576.0, y = 8064920 + 1100, label = "LTER 1",size=4) +
    annotate("text", x = 202527.4, y = 8065298 + 1100, label = "LTER 2",size=4) +
    annotate("text", x = 206348.7 + 1500, y = 8061512, label = "LTER 3",size=4) +
    annotate("text", x = 205727.6 + 1500, y = 8058194, label = "LTER 4",size=4) +
    annotate("text", x = 195586.2 - 1500, y = 8054071, label = "LTER 5",size=4) +
    annotate("text", x = 190191.7 - 1500, y = 8060737, label = "LTER 6",size=4)
  dev.off()
  
  
