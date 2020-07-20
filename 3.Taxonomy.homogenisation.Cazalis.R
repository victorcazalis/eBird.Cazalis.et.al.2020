library(dplyr) ; library(mgcv) ; library(gtools)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_)) # Function I'll use in subset to exclude data (opposite to %in%)



### Charge data
chlist<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script11.", version, ".", zone, ".rds"))
obs<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script11.", version, ".", zone, ".rds"))

# Remove China
chlist<-subset(chlist, chlist$Country != "CN")
obs<-subset(obs, rownames(obs) %in% chlist$Liste) 

# Remove columns when the only obs were in China
x<-c()
for(i in 1:ncol(obs)){x[i]<-table(factor(is.na(obs[,i]), c("TRUE", "FALSE")))["FALSE"]}
obs<-obs[, x!=0]

cat("Should be 100% TRUE", "\n")
print(table(rownames(obs)==chlist$Liste))





###########################################################
### CHANGE SPECIES NAMES WHEN SPECIES HAVE BEEN RENAMED ###
###########################################################
taxo_rename<-read.csv("Taxonomy_changes_Rename.csv", sep=";") # This is the second sheet of the excel file "Species taxonomy changes supplementary table 02.12.19"


for(i in 1:ncol(obs)){
  if(names(obs)[i] %in% taxo_rename$Name.eBird)   {names(obs)[i] <- as.character(taxo_rename$Name.Birdlife[taxo_rename$Name.eBird==names(obs)[i]])}
}



###########################################################################
### MERGE SPECIES WHEN 2 eBird SPECIES CORRESPOND TO 1 BirdLife species ###
###########################################################################
taxo_merge<-read.csv("Taxonomy_changes_BLLump.csv", sep=";") # This is the third sheet of the excel file "Species taxonomy changes supplementary table 02.12.19"
taxo_merge<-subset(taxo_merge, taxo_merge$Name.eBird1 %in% names(obs) | taxo_merge$Name.eBird2 %in% names(obs))

for(i in 1:nrow(taxo_merge)){
  names(obs)[which(names(obs)==taxo_merge$Name.eBird1[i])]<-"Lump1"    # First eBird name
    if("Lump1" %not in% names(obs)){obs$Lump1<-0} # Add 0 if the species has not been recorded (helps avoiding further bugs)
  names(obs)[which(names(obs)==taxo_merge$Name.eBird2[i])]<-"Lump2"   # Second eBird name
    if("Lump2" %not in% names(obs)){obs$Lump2<-0} # Add 0 if the species has not been recorded (helps avoiding further bugs)
  
  obs$Lump1<-as.numeric(replace(obs$Lump1, obs$Lump1=="X", 1000000)) # Replace X by a million to sum
  obs$Lump2<-as.numeric(replace(obs$Lump2, obs$Lump2=="X", 1000000))
  
  obs$Lumped<-rowSums(data.frame(Lump1=obs$Lump1,Lump2=obs$Lump2), na.rm=T) # Sum both eBird species
  obs$Lumped<-replace(as.character(obs$Lumped), obs$Lumped>=1000000, "X")   # Bring back Xs
  obs$Lumped<-replace(obs$Lumped, obs$Lumped==0, NA)
  
  # Change names
  obs$Lump1<-obs$Lump2<-NULL # Remove both species
  names(obs)[which(names(obs)=="Lumped")]<-as.character(taxo_merge$Name.Birdlife[i])  # Rename with Birdlife name (will be linked directly to birdlife characteristics)

}





###################################################################################
### SPLIT SPECIES WHEN 2 SPECIES IN BirdLife ARE MERGED IN ONE SPECIES IN eBIRD ###
###################################################################################

library(DescTools) ; library(sf) ; library(sp) ; library(rgdal) ; library(rgeos) ; library(raster)

taxo_split<-read.csv("Taxonomy_changes_eBirdlump.csv", sep=";") # This is the fourth sheet of the excel file "Species taxonomy changes supplementary table 02.12.19"
taxo_split<-subset(taxo_split, taxo_split$Hotspot == zone)

# Charge distributions if not done yet
setwd("D:/Victor/eBird/BirdLife Distributions V7.0")
if("distributions" %not in% ls()){distributions <- st_read(dsn = "a00000009.gdbtable")} # Don't charge if already charged


# Make points out of chlist
pts<-st_as_sf(spTransform(SpatialPointsDataFrame(coords=data.frame(chlist$lon, chlist$lat), data=as.data.frame(chlist$Liste), proj4string=CRS("+init=epsg:4238")), CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs")))

for(i in 1:nrow(taxo_split)){ # For each eBird species to split
  
  Name_eBird<-as.character(taxo_split$Name.eBird[i])
  Names_Birdlife<-as.character(unlist(c(taxo_split[i, c(3:7)])))
  Names_Birdlife<-gsub("[.]", " ", subset(Names_Birdlife, Names_Birdlife != ""))
  
  Distrib<-distributions[distributions$SCINAME %in% Names_Birdlife,]
  
  # Make a plot of species distribution around the hotspot, add points where the species was detected
  G_distr<-ggplot()+
    geom_sf(data=Distrib, aes(fill=SCINAME), col=NA)+
    ggtitle("Distribution and hotspot")+
    geom_sf(data=pts, alpha=0.3)+
    geom_sf(data=pts[is.na(obs[,which(names(obs)==Name_eBird)])==FALSE,], col="red") # Points where the species was detected
  
  
  ### Choose the closest distribution
  dist.SF<-st_distance(pts[is.na(obs[,which(names(obs)==Name_eBird)])==FALSE,], Distrib)
  SCINAME<-as.character(NA)
  for(n.pt in 1:nrow(dist.SF)){
    SCINAME[n.pt]<-as.character(Distrib$SCINAME[which(dist.SF[n.pt,]==min(dist.SF[n.pt,]))])
  }
  
  
  ### Apply changes in obs
  nb_spc<-nlevels(as.factor(SCINAME))
  
  df<-as.data.frame(matrix(rep(obs[,names(obs)==Name_eBird], nb_spc), byrow=F, ncol=nb_spc)) # Create a matrix with nb_spc col and all abundances in lines (fitting checklists order)
  names(df)<-levels(as.factor(SCINAME))
    
  for (SP in 1:nb_spc){
      df[is.na(df[,SP])==FALSE,SP]<-replace(df[is.na(df[,SP])==FALSE,SP], SCINAME != names(df)[SP], NA)
    }
    
    names(df)<-sub(" ", ".", names(df))
    
    obs[, which(names(obs)==Name_eBird)]<-NULL
    obs[, (ncol(obs)+1):(ncol(obs)+nb_spc)]<-df
    
    
    
    
    ## Make a plot to check 
    obs_melt<-reshape2::melt(obs[, names(obs)%in% sub(" ", ".",Names_Birdlife)], measure.vars=names(obs)[names(obs)%in% sub(" ", ".",Names_Birdlife)])
    obs_melt$lon<-rep(chlist$lon, nb_spc)
    obs_melt$lat<-rep(chlist$lat, nb_spc)
    obs_melt<-obs_melt[is.na(obs_melt$value)==F,]
    if(nb_spc==1){obs_melt$variable<-names(df[1])} # Not written in melt if no variation in the names variable
    
    obs_melt$variable<-as.factor(as.character(obs_melt$variable))
    
    N_obs2<-as.numeric(table(is.na(obs_melt$value))["FALSE"]) 
    
    G_obs<-ggplot()+
      geom_sf(data=Distrib, aes(col=SCINAME), alpha=0.7)+
      xlim(extent(pts)[1:2])+ylim(extent(pts)[3:4])+
      geom_point(data=obs_melt, aes(lon, lat, fill=variable), col="black", shape=21)+
      ggtitle(paste0("Nb of observations removed (out of distributions): ", (N_obs1-N_obs2), " over ", N_obs1, " observations"))

  

}




########################
### SAVE THE NEW OBS ###
########################

### Few last manual changes due to remaining mismatches
names(obs)[names(obs)=="Trogon.violaceous"]<-"Trogon.violaceus"
if(zone=="Atlantic.Forest"){names(obs)[names(obs)=="Ramphastos.vitellinus"]<-"Ramphastos.ariel"}
if(zone=="Tumbes"){names(obs)[names(obs)=="Cacicus.cela"]<-"Cacicus.flavicrissus"}



# Remove species not observed
obs.zone<-obs
for(i in 1:ncol(obs.zone)){obs.zone[,i]<-as.character(obs.zone[,i])}
obs.zone<-replace(obs.zone, is.na(obs.zone)==F, 1)
obs.zone<-replace(obs.zone, is.na(obs.zone)==T, 0)
for(i in 1:ncol(obs.zone)){obs.zone[,i]<-as.numeric(obs.zone[,i])}
ObsSums<-apply(obs.zone, 2, sum)
ObsSums<-replace(ObsSums, is.na(ObsSums), 0)
obs<-obs[, ObsSums>0]




saveRDS(obs, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script20.", version, ".", zone, ".rds"))



