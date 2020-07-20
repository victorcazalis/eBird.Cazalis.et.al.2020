library(dplyr) ; library(mgcv) ; library(gtools)
source_lines <- function(file, lines){source(textConnection(readLines(file)[lines]))}


##############################
### CHARGE DATA AND FORMAT ###
##############################

# Charger chlist, obs and list of species
species<-read.csv(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/species.characteristics.", version, ".", zone, ".csv"))
chlist<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/chlist.Script22.", version, ".", zone, ".rds"))
obs<-readRDS(paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obs.Script21.", version, ".", zone, ".rds"))


cat("Should be 100% TRUE")
table(rownames(obs)==chlist$Liste)
cat("Should be 100% TRUE")
table(colnames(obs)==species$Species)




#######################################
### ANALYSIS IIIa (forest presence) ###
#######################################

# Change factor order
chlist$Forbin<-factor(chlist$Forbin, c("Non-Forest", "Forest"))



for(IND in 1:4){
  Nom_indice<-c("Species richness", "Forest specialisation2", "Endemism", "IUCN")[IND]


if(Nom_indice=="Species richness"){chlist$INDICE<-chlist$rich ; fam="gaussian"}
if(Nom_indice=="Forest specialisation1"){chlist$INDICE<-chlist$Ass.specialisation1 ; fam="nb"}
if(Nom_indice=="Forest specialisation2"){chlist$INDICE<-chlist$Ass.specialisation2 ; fam="nb"}
if(Nom_indice=="IUCN"){chlist$INDICE<-chlist$iucn.index ; fam="nb"}
if(Nom_indice=="Endemism"){chlist$INDICE<-chlist$endem.index ; fam="nb"}


if(Nom_indice %in% c("Species richness")){
  M.INDICE  <-gam(INDICE ~ Forbin   + te(day, lat, lon) + s(duration, k=4) + s(alt, k=6) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
}

if(Nom_indice %in% c("Forest specialisation1", "Forest specialisation2", "IUCN", "Endemism")){
  M.INDICE  <-gam(INDICE ~ Forbin + log(rich)  + te(day, lat, lon) + s(duration, k=4) + s(alt, k=6) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
}
  
  
## Save results
  res.glob<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/III.Forest.results.csv", sep=",")
  
  INDI<-c("Richness", "ForDep", "Endemism", "IUCN")[IND]
  res.glob$Coef[res.glob$zone==zone & res.glob$Indice==INDI]    <-coefficients(M.INDICE)["ForbinForest"]
  res.glob$SE[res.glob$zone==zone & res.glob$Indice==INDI]      <-summary(M.INDICE)$se["ForbinForest"]
  res.glob$P.value[res.glob$zone==zone & res.glob$Indice==INDI] <-summary(M.INDICE)$p.pv["ForbinForest"]
  
  
  write.csv(res.glob, "D:/Victor/eBird/Data/Expand.Analyses.tables/III.Forest.results.csv", row.names=F)
  
} ### END OF THE LOOP



######################################
### ANALYSIS IIIb (forest quality) ###
######################################

# Keep forest checklists only
chlist<-subset(chlist, chlist$Forbin=="Forest")
obs<-obs[rownames(obs) %in% chlist$Liste,]
obs<-obs[,colSums(obs)>0]
species<-subset(species, species$Species %in% colnames(obs))
species$Species<-droplevels(species$Species)

# Inverse HFP (to talk about habitat quality rather than degradation)
chlist$hfp.inv<- (-1)*chlist$hfp


cat("Should be 100% TRUE")
table(rownames(obs)==chlist$Liste)
cat("Should be 100% TRUE")
table(colnames(obs)==species$Species)



### MODELS

# Stock plots
GT<-list()

for(IND in 1:4){
  Nom_indice<-c("Species richness", "Forest specialisation2", "Endemism", "IUCN")[IND]
  
  
  if(Nom_indice=="Species richness"){chlist$INDICE<-chlist$rich ; fam="gaussian"}
  if(Nom_indice=="Forest specialisation1"){chlist$INDICE<-chlist$Ass.specialisation1 ; fam="nb"}
  if(Nom_indice=="Forest specialisation2"){chlist$INDICE<-chlist$Ass.specialisation2 ; fam="nb"}
  if(Nom_indice=="IUCN"){chlist$INDICE<-chlist$iucn.index ; fam="nb"}
  if(Nom_indice=="Endemism"){chlist$INDICE<-chlist$endem.index ; fam="nb"}
  
  
  
  ### MODELE
  if(Nom_indice %in% c("Species richness")){
    M.INDICE  <-gam(INDICE ~ PA + scale(hfp.inv) + scale(Forprop) + scale(canopy) + te(day, lat, lon) + s(duration, k=4) + s(alt, k=6) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
  }
  
  if(Nom_indice %in% c("Forest specialisation1", "Forest specialisation2", "IUCN", "Endemism")){
    M.INDICE  <-gam(INDICE ~ PA  + scale(hfp.inv) + scale(Forprop) + scale(canopy) + log(rich)  + te(day, lat, lon) + s(duration, k=4) + s(alt, k=6) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
  }
  
  
  ### PREDICT
  res<-data.frame(Variable=c("PA","hfp.inv", "Forprop", "canopy"),
                  Coef=c(coefficients(M.INDICE)["PAPA"], coefficients(M.INDICE)["scale(hfp.inv)"], coefficients(M.INDICE)["scale(Forprop)"], coefficients(M.INDICE)["scale(canopy)"]),
                  SD=c(summary(M.INDICE)$se["PAPA"], summary(M.INDICE)$se["scale(hfp.inv)"], summary(M.INDICE)$se["scale(Forprop)"],summary(M.INDICE)$se["scale(canopy)"]),
                  Pv=c(summary(M.INDICE)$p.pv["PAPA"], summary(M.INDICE)$p.pv["scale(hfp.inv)"], summary(M.INDICE)$p.pv["scale(Forprop)"],summary(M.INDICE)$p.pv["scale(canopy)"]))
  
  res$Variable<-factor(res$Variable, c("PA", "hfp.inv", "Forprop", "canopy")) 
  
  res$Min<-res$Coef-1.96*res$SD 
  res$Max<-res$Coef+1.96*res$SD
  
  res$coul<-cut(res$Pv, breaks=c(-0.1,0.001,0.01,0.05,1))
  res$coul<-revalue(res$coul, c("(-0.1,0.001]"="firebrick4", "(0.001,0.01]"="firebrick3", "(0.01,0.05]"="tomato", "(0.05,1]"="darkgray"))
  res$coul<-replace(as.character(res$coul), is.na(res$coul)==TRUE, "black")
  
  
  titre<-Nom_indice ; if(IND>1){titre<-paste0(Nom_indice, "in: ", zone)}
  
  GT[[IND]]<-ggplot(res)+
    geom_pointrange(aes(x=Variable, y=Coef, ymin=Min, ymax=Max), col=as.character(res$coul))+
    geom_abline(intercept=0, slope=0, size=1.5)+
    xlab("")+
    ggtitle(paste(Nom_indice, "in: ", zone))
  
  
  
  ### EXPORT RESULTS TO MAKE A GLOBAL PLOT
  
  res.glob<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/AA.Degradation.Results.csv", sep=",")
  res.glob<-subset(res.glob, res.glob$Variable %in% res$Variable)
  res.glob$Variable<-as.character(res.glob$Variable)
  res$Variable<-as.character(res$Variable)
  
  # Write results
  res.glob$coul<-as.character(res.glob$coul)
  
  res$SD<-NULL
  
  res.glob[res.glob$zone==zone & res.glob$Indice==c("Richness", "ForDep", "Endemism", "IUCN")[IND],][,4:9]<-res[,1:6]
  
  
  write.csv(res.glob, "D:/Victor/eBird/Data/Expand.Analyses.tables/AA.Degradation.Results.csv", row.names=F)

} ### END OF THE LOOP


cowplot::save_plot(cowplot::plot_grid(GT[[1]], GT[[2]], GT[[3]], GT[[4]]), filename=paste0("D:/Victor/eBird/Figures/Expand.plots/Hab.degradation/Assemblage.degradation/Degradation.Assemblage.", version, ".", zone, ".png"), base_height=10, base_aspect_ratio=1.4) 
