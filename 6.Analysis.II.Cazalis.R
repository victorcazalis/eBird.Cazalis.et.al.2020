library(DescTools) ; library(sf) ; library(sp) ; library(rgdal) ; library(rgeos) ; library(raster) ; library(spatialEco) ; library(exactextractr)
source_lines <- function(file, lines){source(textConnection(readLines(file)[lines]))} 



###########################################
### EXTRACT VALUES FOR BACKGROUND SITES ###
###########################################

### ADD EXTENT OF INTEREST 
limits<-st_transform(read_sf("D:/Victor/eBird/Data/10.Extent.of.interest/Limites.eBird.interHotspotBiomes.SundaMelan.shp"), st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
limits<-limits[limits$NAME==zone,]
limits<-limits[limits$FIRST_BIOM=="Tropical & Subtropical Moist Broadleaf Forests",]



### CREATE 1km buffers 
r <- raster(extent(limits), res=c(2000,2000), crs = st_crs(limits))    
r[] <- 1:length(r)

sp.r<-as(as(r, "SpatialGridDataFrame"), "SpatialPolygonsDataFrame")
pts.r<-SpatialPointsDataFrame(coords=data.frame(coordinates(r)[,1], coordinates(r)[,2]), data=as.data.frame(c(1:nrow(coordinates(r)))), proj4string=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
pts.r$zone<-over(pts.r, as_Spatial(limits))$NAME

bf1<-sp.r[is.na(pts.r$zone)==FALSE,] %>% st_as_sf() %>% st_centroid() %>% st_buffer(., 1000) 
st_crs(bf1)<-st_crs(limits)
bf<-st_transform(bf1, st_crs("+init=epsg:4238")) 

ggplot(bf[1:50,])+geom_sf() # Check they have the good shape and are adjacent




### RASTER EXTRACTS ###

### Altitude
alt<-raster("D:/Victor/PA GIS/GLOBEelevation.Soum/globe_v1_all/w001001.adf")
alt<-crop(alt, extent(bf), snap="out") # reduce range
bf$alt<-exact_extract(alt, bf, fun=function(value, cov_frac){median(value[cov_frac>0.01], na.rm=T)})

# # Plot to check extraction
# gridExtra::grid.arrange(ggplot()+
#   geom_raster(data=data.frame(alt=alt@data@values, x=coordinates(alt)[,1], y=coordinates(alt)[,2]), aes(x=x, y=y, fill=alt))+
#   geom_sf(data=bf, fill=NA, col="red"),
#   ggplot(bf)+geom_sf(aes(fill=alt)))


### Agricultural suitability
agri<-raster("D:/Victor/PA GIS/Agricultural suitability Zabel/overall_cropsuit_i_1981-2010.tif")
agri<-crop(agri, extent(bf), snap=  "out") 
bf$agri<-exact_extract(agri, bf, fun=function(value, cov_frac){median(value[cov_frac>0.01], na.rm=T)})

### Remoteness
remote<-raster("D:/Victor/PA GIS/Accessibility.to.cities.Weiss2018/accessibility_to_cities_2015_v1.0.tif")
remote<-crop(remote, extent(bf), snap=  "out") 
remote<-replace(remote, remote==(-9999), NA)
bf$remoteness<-exact_extract(remote, bf, fun=function(value, cov_frac){median(value[cov_frac>0.01], na.rm=T)})

### Canopy height
Canopy<-raster("D:/Victor/PA GIS/PA South America/CanopyHeight/Simard_Pinto_3DGlobalVeg_L3C.tif")
Canopy<-crop(Canopy, extent(bf), snap=  "out") 
bf$canopy<-exact_extract(Canopy, bf, fun=function(value, cov_frac){median(value[cov_frac>0.01], na.rm=T)})

### Human footprint
hfp<-raster("D:/Victor/PA GIS/Human Footprint/Footprint2.tif")
hfp<-crop(hfp, extent(bf), snap=  "out") 
bf$hfp<-exact_extract(hfp, bf, fun=function(value, cov_frac){median(value[cov_frac>0.01], na.rm=T)})

### Forest cover
glc<-raster("D:/Victor/PA GIS/CCI Land Cover/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif")
glcFor<-crop(glc, extent(bf), snap="out")
glcFor<-replace(glcFor, glcFor==210, NA) # Mettre les oceans en NA
glcFor2<-replace(glcFor, glcFor %in% c(50:90, 160,170), 1)
glcFor2<-replace(glcFor2, glcFor2 !=1, 0)

bf$For.prop<-exact_extract(glcFor2, bf, fun=function(value, cov_frac){mean(value[cov_frac>0.01], na.rm=T)})
bf$LUbin<-revalue(as.character(cut(bf$For.prop, breaks=c(-0.1,0.1,0.6,1.1))), c("(0.6,1.1]"="Forest", "(-0.1,0.1]"="Non-Forest", "(0.1,0.6]"="Intermed"))


### Hansen deforestation rates
if(zone %in% c("Tumbes", "Atlantic.Forest", "Mesoamerica", "Tropical.Andes")){hansen<-raster("D:/Victor/PA GIS/Hansen.hotspots.eBird/Data/Hansen.merged.Americas.tif")}
if(zone %in% c("Eastern.Afromontane")){hansen<-raster("D:/Victor/PA GIS/Hansen.hotspots.eBird/Data/Hansen.merged.Africa.tif")}
if(zone %in% c("Ghats.Lanka", "Indo-Burma", "Sundaland")){hansen<-raster("D:/Victor/PA GIS/Hansen.hotspots.eBird/Data/Hansen.merged.Asia.tif")}

bf$Hansen.loss<-exactextractr::exact_extract(hansen, bf, fun=function(value, cov_frac){mean(replace(value[cov_frac>0.01], value[cov_frac>0.01]>0, 1))})



### ADD PROTECTION VALUES ###
shp.PA<-st_transform(read_sf("D:/Victor/PA GIS/PA mondiales/WDPA_Nov2018-shapefile-polygons.fixed.30.11.18.shp"), st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
shp.PA<-st_intersection(st_buffer(shp.PA,0), st_buffer(limits,0))
shp.PA<-subset(shp.PA, shp.PA$STATUS %in% c("Designated", "Inscribed", "Established")) # Removes Adopted, Not Reported, Proposed (total = 1700 / 216 000)
shp.PA<-subset(shp.PA, shp.PA$DESIG_ENG != "UNESCO-MAB Biosphere Reserve")
shp.PA<-st_buffer(st_combine(shp.PA), 0)

bfSF<-st_transform(bf, st_crs(shp.PA))


# Binary
pts<-st_centroid(bfSF)
pts.inter<-st_intersection(pts, shp.PA)
bfSF$PA<-revalue(as.factor(bfSF$layer %in% pts.inter$layer), c("TRUE"="PA", "FALSE"="unPA"))



### Add lon and lat
CENTR<-st_centroid(bfSF)
bfSF$lon<-st_coordinates(CENTR)[,1]
bfSF$lat<-st_coordinates(CENTR)[,2]



### SAVE
st_write(bfSF, paste0("D:/Victor/eBird/Figures/Expand.plots/Area.characteristics/Shapefile.Analyse.II/Shape.Analyse.grid.habitat.Sept19.", zone, ".shp")) 

write.csv(data.frame(Layer=bfSF$layer, lon=bfSF$lon, lat=bfSF$lat, alt=bfSF$alt, agri=bfSF$agri, remoteness=bfSF$remoteness, canopy=bfSF$canopy, hfp=bfSF$hfp, For.prop=bfSF$For_prop, LUbin=bfSF$LUbin, PA=bfSF$PA
), paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/Tableau.Analyse.II.grid.habitat.", zone, ".csv"), row.names=FALSE)




####################
### ANALYSIS IIa ###
####################
grid<-read.csv(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/Tableau.Analyse.II.grid.habitat.", zone, ".csv"))


### Data management
grid<-subset(grid, grid$LUbin != "Intermed")
grid$LUbin<-factor(droplevels(grid$LUbin), c("Non-Forest", "Forest"))
grid$PA <- factor(grid$PA, c("unPA", "PA"))

### Forest model
M.For.2<-gam(LUbin ~ PA + s(remoteness) + s(alt) + s(agri) +te(lon,lat), data=grid, family="binomial")



#####################
### ANALYSIS IIa' ###
#####################
M.Hansen.2<-gam(log(0.001+Hansen.loss) ~ PA + s(remoteness) + s(alt) + s(agri) +te(lon,lat), data=grid, family="gaussian")

# I give here the example of how I calculate the percentage given in Extended Table 1 for habitat variables (predicting the difference if all protected sites haven't been protected)
gridPred<-grid[grid$PA=="PA",]
res$Pred.PA[res$zone==zone]<-exp(mean(predict(M.Hansen.2, newdata=gridPred, type="response"), na.rm=T))+0.001
gridPred$PA<-"unPA"
res$Pred.unPA[res$zone==zone]<-exp(mean(predict(M.Hansen.2, newdata=gridPred, type="response"), na.rm=T))+0.001




####################
### ANALYSIS IIb ###
####################
### Canopy model
gridFOR<-subset(grid, grid$LUbin == "Forest")
M.Can.2<-gam(canopy ~ PA + s(remoteness) + s(alt) + s(agri)  +te(lon,lat), data=gridFOR, family="gaussian")
# par(mfrow=c(2,2)) ; gam.check(M.Can.2) ; par(mfrow=c(1,1))
# visreg(M.Can.2)


# Model Wilderness
M.HFP.2<-gam((-1*hfp) ~ PA + s(remoteness) + s(alt) + s(agri)  +te(lon,lat), data=gridFOR, family="gaussian")


# Model Contiguity
M.Frag.2<-gam(For.prop ~ PA + s(remoteness) + s(alt) + s(agri)  +te(lon,lat), data=gridFOR, family="gaussian")




### Save results
res<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/Results.II.csv", sep=",")

res[res$zone==zone & res$Variable=="Forest", c("Coef", "SE", "P.value")]    <- c(coefficients(M.For.2)["PAPA"], summary(M.For.2)$se["PAPA"], summary(M.For.2)$p.pv["PAPA"])
res[res$zone==zone & res$Variable=="Hansen", c("Coef", "SE", "P.value")]    <- c(coefficients(M.Hansen.2)["PAPA"], summary(M.Hansen.2)$se["PAPA"], summary(M.Hansen.2)$p.pv["PAPA"])
res[res$zone==zone & res$Variable=="Canopy", c("Coef", "SE", "P.value")]    <- c(coefficients(M.Can.2)["PAPA"], summary(M.Can.2)$se["PAPA"], summary(M.Can.2)$p.pv["PAPA"])
res[res$zone==zone & res$Variable=="HFP", c("Coef", "SE", "P.value")]    <- c(coefficients(M.HFP.2)["PAPA"], summary(M.HFP.2)$se["PAPA"], summary(M.HFP.2)$p.pv["PAPA"])
res[res$zone==zone & res$Variable=="Fragmentation", c("Coef", "SE", "P.value")]    <- c(coefficients(M.Frag.2)["PAPA"], summary(M.Frag.2)$se["PAPA"], summary(M.Frag.2)$p.pv["PAPA"])

write.csv(res, "D:/Victor/eBird/Data/Expand.Analyses.tables/Results.II.csv", row.names=FALSE)

