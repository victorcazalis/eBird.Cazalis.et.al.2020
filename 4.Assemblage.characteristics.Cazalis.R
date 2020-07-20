library(dplyr) ; library(mgcv) ; library(gtools) ; library(sf)

chlist<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script11.", version, ".", zone, ".rds"))
obs<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script20.", version, ".", zone, ".rds"))



# Remove China (no data on protected areas)
chlist<-subset(chlist, chlist$Country != "CN")
obs<-subset(obs, rownames(obs) %in% chlist$Liste) 

cat("Should be 100% TRUE")
table(rownames(obs)==chlist$Liste)



# Charge BirdLife International distributions
setwd("D:/Victor/eBird/BirdLife Distributions V7.0")
distributions <- st_read(dsn = "a00000009.gdbtable")


### Species list
species<-data.frame(Species=colnames(obs))

cat("Should be 100% '1'")
table(table(species$Species))


### Calculate the number of observations in the hotspot per species
for(i in 1:ncol(obs)){
  obs[,i]<-as.numeric(replace(obs[,i], is.na(obs[,i])==F, 1))
  obs[,i]<-as.numeric(replace(obs[,i], is.na(obs[,i])==T, 0))
}

species$N<-colSums(obs, na.rm=T)



### Endemism and range
endem_limits<-st_transform(read_sf("D:/Victor/eBird/Data/20.Endemism.limits/Endemism.limits.Hotspots.shp"), st_crs(distributions)) # Corresponds to the limits of hotspots without restrictions on biomes
endem_limits<-endem_limits[endem_limits$NAME==zone,]

species$range_IUCN<-NA
species$range_zone<-NA
for(i in 1:nrow(species)){
  sp<-subset(distributions, distributions$SCINAME==sub("[.]", " ", as.character(species$Species[i])))
  species$range_IUCN[i]<-sum(st_area(sp))
  species$range_zone[i]<-sum(st_area(st_intersection(st_buffer(sp,0), endem_limits)))
}
species$Endem90<-species$range_zone/species$range_IUCN>=0.9


### Red List status
iucn<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/global/RL status birds.csv", sep=";") 
species$IUCN<-iucn[match(species$Species, sub(" ", ".", iucn$Scientific.name)), "RL.2017"]
species$IUCN<-factor(droplevels(species$IUCN), c("CR", "EN", "VU", "NT", "LC"))


### Add Forest dependence 
for.dep<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/global/Forest dependency 2017.csv", sep=";")
species$for.dep<-for.dep[match(species$Species, gsub(" ",".", for.dep$Scientific.name)), "Forest.dependency"]



### Save 
cat("Should be 100% TRUE")
table(species$Species==colnames(obs))

write.csv(species, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/species.characteristics.", version, ".", zone, ".csv"), row.names = FALSE)

saveRDS(chlist, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script21.", version, ".", zone, ".rds"))
obs<-replace(obs, is.na(obs)==T, 0)
saveRDS(obs, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script21.", version, ".", zone, ".rds"))








####################################
### CALCULATE ASSEMBLAGE INDICES ###
####################################
library(dplyr) ; library(mgcv) ; library(gtools)
source_lines <- function(file, lines){source(textConnection(readLines(file)[lines]))}


##############################
### CHARGE DATA AND FORMAT ###
##############################

# Charge chlist, obs and list of species
species<-read.csv(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/species.characteristics.", version, ".", zone, ".csv"))
chlist<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script21.", version, ".", zone, ".rds"))
obs<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script21.", version, ".", zone, ".rds"))


cat("Should be 100% TRUE")
table(rownames(obs)==chlist$Liste)
cat("Should be 100% TRUE")
table(colnames(obs)==species$Species)


#################################
### CREATE ASSEMBLAGE INDEXES ###
#################################

# Richness in endemic species
obs.endem<-obs
for(i in 1:ncol(obs.endem)){
  obs.endem[,i]<-as.character(obs.endem[,i])
  obs.endem[,i]<-replace(obs.endem[,i], obs.endem[,i]==1, as.character(species$Endem90[i]))
}

chlist$endem.index<-NA
for(i in 1:nrow(obs.endem)){
  chlist$endem.index[i]<-table(factor(obs.endem[i,] %in% c("TRUE"), c("TRUE","FALSE")))["TRUE"]
}


### Richness in threatened and Near Threatened species
obs.iucn<-obs
for(i in 1:ncol(obs.iucn)){
  obs.iucn[,i]<-as.character(obs.iucn[,i])
  obs.iucn[,i]<-replace(obs.iucn[,i], obs.iucn[,i]==1, as.character(species$IUCN[i]))
}

chlist$iucn.index<-NA
for(i in 1:nrow(obs.iucn)){
  chlist$iucn.index[i]<-table(factor(obs.iucn[i,] %in% c("NT","VU","EN","CR"), c("TRUE","FALSE")))["TRUE"]
}


### Richness in forest dependent species
obs.for<-obs

for(i in 1:ncol(obs.for)){
  obs.for[,i]<-as.character(obs.for[,i])
  obs.for[,i]<-replace(obs.for[,i], obs.for[,i]==1, as.character(species$for.dep[i]))
}

chlist$Ass.specialisation2<-chlist$Ass.specialisation1<-NA
for(i in 1:nrow(obs.for)){
  chlist$Ass.specialisation2[i]<-table(factor(obs.for[i,] %in% c("High", "Medium"), c("TRUE","FALSE")))["TRUE"]
}


saveRDS(chlist, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script22.", version, ".", zone, ".rds"))








