`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_)) # Function I'll use in subset to exclude data (opposite to %in%)
source_lines <- function(file, lines){source(textConnection(readLines(file)[lines]))} # Function I'll use to source only part of a script


version="dec18" # Used across the analyses in file names

### Choose a zone (ie hotspot) between: "Atlantic.Forest", "Indo-Burma", "Mesoamerica", "Tropical.Andes", "Ghats.Lanka", "Caribbean.Islands", "Coastal.Africa", "Eastern.Afromontane", "Guinean.Africa", "Himalaya", "Madagascar", "Philippines", "Tumbes", "Sundaland
zone<-"Atlantic.Forest"




#########################
### EXPORT eBIRD DATA ###
#########################


## Go in the folder with all raw data from this zone
setwd(paste0("D:/Victor/eBird/Data/To export/", zone))

library(auk) ; library(plyr) ; library(readr)



### Loop to extract data from each file
for (i in 1:length(list.files())){

eb<-read_ebd(list.files()[i], rollup=T, unique=T) # Lit le fichier, utilise auk_rollup (supprime sp. et regroupe les sous-especes à l'espece) et auk_unique (supprime les listes dupliquées par les differents observateurs)



names(eb)<-replace(names(eb), names(eb)=="checklist_id", "checklist")

### Subset data
# Remove useless columns
eb[,c("global_unique_identifier", "last_edited_date", "subspecies_common_name", "breeding_bird_atlas_category", "effort_area_ha", "county", "subnational2_code", "iba_code", "bcr_code", "usfws_code", "atlas_block", "locality", "locality_id", "locality_type", "first_name", "last_name", "has_media", "reviewed", "x")]<-NULL
# Remove domestic species (Columba livia only)
eb<-subset(eb, eb$category != "domestic")
# Remove checklists that did not report all species
eb<-subset(eb, eb$all_species_reported==T)
# Remove disapproved observations
eb<-subset(eb, eb$approved==T)
# Select protocols
eb<-subset(eb, eb$protocol_type %in% c("Traveling","Stationary","Historical"))
# Remove data with no duration
eb<-subset(eb, is.na(eb$duration_minutes)==FALSE)
# Remove historical protocols with no distance reported (ie only stationary are allowed not to have a distance reported)
eb<-subset(eb, eb$protocol_type=="Stationary" | is.na(eb$effort_distance_km)==FALSE)


### Format date
eb$year<-as.numeric(format(eb$observation_date, "%Y"))
eb$month<-as.numeric(format(eb$observation_date, "%m"))
eb$day<-as.numeric(format(eb$observation_date, "%j"))




saveRDS(eb, file=paste("D:/Victor/eBird/Data/Exported/", zone,"/ebird.",version,".Export.", substr(list.files()[i],5,6), ".rds", sep=""))

cat(i)

}




################################
### MERGE DATA FROM THE ZONE ###
################################
setwd(paste0("D:/Victor/eBird/Data/Exported/", zone))


df<-readRDS(list.files()[1]) # Stocks the first file

if(length(list.files())>1){
  for(i in 2:length(list.files())){
    df<-rbind(df, readRDS(list.files()[i])) # Add each file, starting at i=2
  }
}

saveRDS(df, paste0("D:/Victor/eBird/Data/Exported/0.MERGED/Merged.data.", version, ".", zone, ".rds"))





###################
### FORMAT DATA ###
###################

library(plyr) ;library(dplyr) ; library(mgcv) ; library(gtools) ; library(sf) ; library(sp) ; library(rgdal) ; library(rgeos) ; library(SDMTools) ; library(shapefiles)



eb.raw<-readRDS(paste("Merged.data", version, zone, "rds", sep="."))



### Reduce dataset for good quality checklists
# Remove short or long checklists
eb<-subset(eb.raw, eb.raw$duration_minutes>30 & eb.raw$duration_minutes<600)
# Remove old data
eb<-subset(eb, eb$year>=2005)
# Remove long travelling distances
eb<-subset(eb, is.na(eb$effort_distance_km)==TRUE | eb$effort_distance_km<=5)
# Charge observer expertise data
if(zone %in% c("Eastern.Afromontane")){                                        obsqual<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obsKelling.Africa." , version, ".rds"))}
if(zone %in% c("Atlantic.Forest", "Mesoamerica", "Tropical.Andes", "Tumbes")){ obsqual<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obsKelling.America.", version, ".rds"))}
if(zone %in% c("Indo-Burma", "Ghats.Lanka", "Sundaland")){                     obsqual<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obsKelling.Asia."   , version, ".rds"))}



### Build observation tables (trick to keep X values)
eb$Ab<-as.numeric(replace(eb$observation_count, eb$observation_count=="X", 10000000000))
obs<-data.frame(round(tapply(eb$Ab, list(eb$checklist, eb$scientific_name), sum)))
obs1<-replace(obs, obs>100000000, "X")


### Build checklist table
chlist<-ddply(eb, ~eb$checklist, function(x){
  data.frame(duration=mean(x$duration_minutes), distance=mean(x$effort_distance_km), year=mean(x$year), day=mean(x$day), lat=mean(x$latitude), lon=mean(x$longitude), protocol=x$protocol_type[1], Country=x$country_code[1],
             Observers=paste(unique(unlist(strsplit(x$observer_id, ","))), collapse=";"))})
colnames(chlist)[1]<-"Liste"

# Add in chlist the expertise score of the most skilled observer
chlist$obsKelling<-chlist$N_obs<-NA
for(i in 1:nrow(chlist)){
  OBS<-unique(unlist(strsplit(as.character(chlist$Observers[i]), ";")))
  Qual<-obsqual$expertise[match(OBS, obsqual$observer)]
  chlist$obsKelling[i]<-max(Qual, na.rm=TRUE)
  chlist$N_obs[i]<-length(OBS)
}
chlist$obsKelling<-replace(chlist$obsKelling, chlist$obsKelling=="-Inf", NA)
chlist<-subset(chlist, is.na(chlist$obsKelling)==FALSE)

# Subset obs with checklists that are kept
obs<-obs[rownames(obs) %in% chlist$Liste,]
obs1<-obs1[rownames(obs1) %in% chlist$Liste,]
table(chlist$Liste == rownames(obs1)) # Should be 100% true




### Delete checklist too close the same day ###

## Load the function that enable to calculate distance between sites (could be redone with st_distance)
library(Imap)

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
} ## END OF FUNCTION


chlist$suppr<-NA # I'll write "S" in this column for checklists that should be removed
chlist$dates<-paste(chlist$day, chlist$year)
Dates<-levels(as.factor(chlist$dates))

for(i in 1:length(Dates)){  
  x<-subset(chlist, chlist$dates==Dates[i]) # Subset per date
  
  Dist<-round(GeoDistanceInMetresMatrix(x)/1000) ; diag(Dist)<-NA ; Dist[upper.tri(Dist)]<-NA # Distance matrix 
  
  if(min(Dist, na.rm=T)<2){ # if there are some distances <2km
    couples<-as.data.frame(which(Dist < 2, arr.ind=T)) # Show me the pairs closer to 2km
    
    for(j in 1:nrow(couples)){
      sub<-x[as.numeric(couples[j,]),] # For each pair, add an "S" for one of the checklist
      chlist$suppr[chlist$Liste == sub$Liste[sample(1:nrow(sub), 1)]] <- "S"
    }
  }
  if(round(i/100)==i/100){ cat(i, "  ")} # Print progress
}


chlist<-subset(chlist, is.na(chlist$suppr)==T) # Keep only the NA (i.e., the one without "S")
chlist$suppr<-chlist$dates<-NULL

### Check that the suppression of checklist the same day at less that 2km has worked
# ply<-ddply(chlist, .(paste(chlist$year, chlist$day)), function(x){
#   
#   if(nrow(x)==1){Dist<-9999999}else{ # Si il n'y a qu'une ligne : pas de matrice de distance
#     Dist<-round(GeoDistanceInMetresMatrix(x)/1000)
#     diag(Dist)<-NA}
#   
#   data.frame(
#     year=x$year[1],
#     day=x$day[1],
#     Nb_obs=nrow(x),
#     Dist_min=min(Dist, na.rm=T)
#   )
# })
# cat("Should be 100% FALSE")
# table(ply$Dist_min<2)


obs1<-subset(obs1, rownames(obs1) %in% chlist$Liste)
cat("Should be 100% TRUE")
table(rownames(obs1)==chlist$Liste)


### Save tables 
saveRDS(obs1, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script10.", version, ".", zone, ".rds"))
saveRDS(chlist, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script10.", version, ".", zone, ".rds"))








###############################
### EXTRACT GIS INFORMATION ###
###############################
library(DescTools) ; library(sf) ; library(sp) ; library(rgdal) ; library(rgeos) ; library(raster)


chlist<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script10.", version, ".", zone, ".rds"))
obs<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script10.", version, ".", zone, ".rds"))


### ADD EXTENT OF INTEREST
limits<-st_transform(read_sf("D:/Victor/eBird/Data/10.Extent.of.interest/Limites.eBird.shp"), st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
limits<-limits[limits$NAME==zone,]
limits<-limits[limits$FIRST_BIOM=="Tropical & Subtropical Moist Broadleaf Forests",]
limitsSP<-spTransform(readOGR("D:/Victor/eBird/Data/10.Extent.of.interest/Limites.eBird.shp"), CRSobj = CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
limitsSP<-limitsSP[limitsSP$NAME==zone,]
limitsSP<-limitsSP[limitsSP$FIRST_BIOM=="Tropical & Subtropical Moist Broadleaf Forests",]
limitsSP$id<-"1"

pts.full<-spTransform(SpatialPointsDataFrame(coords=data.frame(chlist$lon, chlist$lat), data=as.data.frame(chlist$Liste), proj4string=CRS("+init=epsg:4238")), CRSobj = CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
chlist$In.zone<-over(pts.full, limitsSP)$id

# Only keep sites within limitsSP
chlist<-subset(chlist, is.na(chlist$In.zone)==FALSE)
pts<-pts.full[pts.full$`chlist$Liste` %in% chlist$Liste,]



### ADD PROTECTED AREAS
shp.PA<-st_transform(read_sf("D:/Victor/PA GIS/PA mondiales/WDPA_Nov2018.shp"), st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
shp.PA<-st_intersection(st_buffer(shp.PA,0), st_buffer(limits,0))

# Remove PAs to exclude according to WDPA guidelines
shp.PA<-subset(shp.PA, shp.PA$STATUS %in% c("Designated", "Inscribed", "Established")) # Removes Adopted, Not Reported, Proposed (total = 1700 / 216 000)
shp.PA<-subset(shp.PA, shp.PA$DESIG_ENG != "UNESCO-MAB Biosphere Reserve")

shp.PA<-st_buffer(st_combine(shp.PA), 0) # Combine and fix some polygon problems

PAs.inter<-st_intersection(st_as_sf(pts), shp.PA)
chlist$PA<-revalue(as.factor(chlist$Liste %in% PAs.inter$`chlist$Liste`), c("TRUE"="PA", "FALSE"="unPA"))


# Calculate protection in a 1 km buffer (test in Supplementary Information)
bf.sf<-st_buffer(st_as_sf(pts),1000)
bf.sf$area<-st_area(bf.sf)
chlist$Buff.Area<-as.numeric(bf.sf$area[match(chlist$Liste, bf.sf$`chlist$Liste`)])

PAs.inter1000<-st_intersection(bf.sf, shp.PA)
PAs.inter1000$area<-st_area(PAs.inter1000) 
chlist$area.PA1000<-as.numeric(PAs.inter1000$area[match(chlist$Liste, PAs.inter1000$`chlist$Liste`)])
chlist$area.PA1000<-replace(chlist$area.PA1000, is.na(chlist$area.PA1000)==TRUE, 0)
chlist$PA1000<-chlist$area.PA1000/chlist$Buff.Area
chlist$Buff.Area<-chlist$area.PA1000<-NULL
hist(chlist$PA1000)


### Add raster values

# Create buffer
bf <- spTransform(gBuffer(pts, width = 1000, byid=TRUE, quadsegs=24, id=pts$`chlist$Liste`), CRSobj = CRS("+init=epsg:4238"))


### Altitude
alt<-raster("D:/Victor/PA GIS/globe_v1_all/w001001.adf")
alt<-crop(alt, extent(bf), snap="out") # reduce range
bf$alt<-extract(alt, bf, fun=median, df=T, na.rm=T)$w001001


### Land cover 2015
glc<-raster("D:/Victor/PA GIS/Land Cover/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif")
glcFor<-crop(glc, extent(bf), snap="out")
glcFor<-replace(glcFor, glcFor==210, NA) # Oceans are NA (not included in proportion calculations)
glcFor2<-replace(glcFor, glcFor %in% c(50:90, 160,170), 1) # Forest are one, the rest is zero so that the mean corresponds to the proportion
glcFor2<-replace(glcFor2, glcFor2 !=1, 0)

bf$For.prop<-extract(glcFor2, bf, fun=median, na.rm=T) 
bf$LUbin<-revalue(as.character(cut(bf$For.prop, breaks=c(-0.1,0.1,0.6,1.1))), c("(0.6,1.1]"="Forest", "(-0.1,0.1]"="Non-Forest", "(0.1,0.6]"="Intermed"))


### Human footprint
hfp<-raster("D:/Victor/PA GIS/Human Footprint/Footprint2.tif")
hfp<-crop(hfp, extent(bf), snap=  "out") 
bf$hfp<-extract(hfp, bf, fun=median, df=T, na.rm=T)

### Agricultural suitability
agri<-raster("D:/Victor/PA GIS/Agricultural suitability Zabel/overall_cropsuit_i_1981-2010.tif")
agri<-crop(agri, extent(bf), snap=  "out") 
bf$agri<-extract(agri, bf, fun=median, df=T, na.rm=T)$overall_cropsuit_i_1981.2010

### Remoteness
remote<-raster("D:/Victor/PA GIS/Accessibility.to.cities.Weiss2018/accessibility_to_cities_2015_v1.0.tif")
remote<-crop(remote, extent(bf), snap=  "out") 
remote<-replace(remote, remote==(-9999), NA)
bf$remoteness<-extract(remote, bf, fun=median, df=T, na.rm=T)


### Canopy height
Canopy<-raster("D:/Victor/PA GIS/PA South America/CanopyHeight/Simard_Pinto_3DGlobalVeg_L3C.tif")
Canopy<-crop(Canopy, extent(bf), snap=  "out") 
bf$canopy<-extract(Canopy, bf, fun=median, df=T, na.rm=T)


### Add in chlist
chlist$alt<-bf$alt[match(chlist$Liste, bf$`chlist$Liste`)]
chlist$Forbin<-bf$LUbin[match(chlist$Liste, bf$`chlist$Liste`)]
chlist$Forprop<-bf$For.prop[match(chlist$Liste, bf$`chlist$Liste`)]
chlist$hfp<-bf$hfp[,2][match(chlist$Liste, bf$`chlist$Liste`)]
chlist$remoteness<-bf$remoteness[,2][match(chlist$Liste, bf$`chlist$Liste`)]
chlist$canopy<-bf$canopy[,2][match(chlist$Liste, bf$`chlist$Liste`)]
chlist$agri<-bf$agri[,2][match(chlist$Liste, bf$`chlist$Liste`)]




### SAVE CHLIST VALUES
chlist<-subset(chlist, chlist$Forbin %in% c("Forest","Non-Forest")) # Remove intermediate forests
saveRDS(chlist, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script11.", version, ".", zone, ".rds"))

# Remove checklists outside the zone
obs<-obs[rownames(obs) %in% chlist$Liste,]
# Remove species that don't belong to the new subset (ie not observed in the checklists we've kept)
obs.zone<-replace(obs, is.na(obs)==F, 1)
obs.zone<-replace(obs.zone, is.na(obs.zone)==T, 0)
for(i in 1:ncol(obs.zone)){obs.zone[,i]<-as.numeric(obs.zone[,i])}
ObsSums<-apply(obs.zone, 2, sum)
obs<-obs[, ObsSums>0]
# Save
saveRDS(obs, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script11.", version, ".", zone, ".rds"))











