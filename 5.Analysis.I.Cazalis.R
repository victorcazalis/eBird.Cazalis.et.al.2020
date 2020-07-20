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





#################################
### MODELS ASSEMBLAGE INDEXES ###
#################################
for(IND in 1:4){
  Nom_indice<-c("Species richness", "Forest specialisation2", "Endemism", "IUCN")[IND]


if(Nom_indice=="Species richness"){chlist$INDICE<-chlist$rich ; fam="gaussian"}
if(Nom_indice=="Forest specialisation1"){chlist$INDICE<-chlist$Ass.specialisation1 ; fam="nb"}
if(Nom_indice=="Forest specialisation2"){chlist$INDICE<-chlist$Ass.specialisation2 ; fam="nb"}
if(Nom_indice=="IUCN"){chlist$INDICE<-chlist$iucn.index ; fam="nb"}
if(Nom_indice=="Endemism"){chlist$INDICE<-chlist$endem.index ; fam="nb"}


if(Nom_indice %in% c("Species richness")){
  M.INDICE1  <-gam(INDICE ~ PA                                                   + te(day, lat, lon) + s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
  M.INDICE2  <-gam(INDICE ~ PA + s(remoteness, k=6) + s(agri, k=6) + s(alt, k=6) + te(day, lat, lon) + s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
}

if(Nom_indice %in% c("Forest specialisation1", "Forest specialisation2", "IUCN", "Endemism")){
  M.INDICE1  <-gam(INDICE ~ PA                                                   + log(rich)  + te(day, lat, lon) + s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
  M.INDICE2  <-gam(INDICE ~ PA + s(remoteness, k=6) + s(agri, k=6) + s(alt, k=6) + log(rich)  + te(day, lat, lon) + s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlist, family=fam)
}
  
  
## Save results
  res.glob<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/I.PAs.effects.birds.results.csv", sep=",")
  
  # Write the time of last change
  res.glob$Date<-as.character(res.glob$Date)
  res.glob$Date[res.glob$zone==zone]<-as.character(Sys.time())
  
  # Write results

  INDI<-c("Richness", "ForDep", "Endemism", "IUCN")[IND]
  res.glob$Coef   [res.glob$zone==zone & res.glob$Indice==INDI & res.glob$Test=="Loc.Impl"]      <-coefficients(M.INDICE1)["PAPA"]
  res.glob$SE     [res.glob$zone==zone & res.glob$Indice==INDI & res.glob$Test=="Loc.Impl"]      <-summary(M.INDICE1)$se["PAPA"]
  res.glob$P.value[res.glob$zone==zone & res.glob$Indice==INDI & res.glob$Test=="Loc.Impl"]      <-summary(M.INDICE1)$p.pv["PAPA"]
  
  res.glob$Coef   [res.glob$zone==zone & res.glob$Indice==INDI & res.glob$Test=="Impl.only"]      <-coefficients(M.INDICE2)["PAPA"]
  res.glob$SE     [res.glob$zone==zone & res.glob$Indice==INDI & res.glob$Test=="Impl.only"]      <-summary(M.INDICE2)$se["PAPA"]
  res.glob$P.value[res.glob$zone==zone & res.glob$Indice==INDI & res.glob$Test=="Impl.only"]      <-summary(M.INDICE2)$p.pv["PAPA"]
  
  
  write.csv(res.glob, "D:/Victor/eBird/Data/Expand.Analyses.tables/I.PAs.effects.birds.results.csv", row.names=F)
  
  

} ### END OF THE LOOP






#######################################################################
### PLOT THE RESULTS (once it has been calculated for each hotspot) ###
#######################################################################

res<-read.csv("D:/Victor/eBird/Data/Expand.Analyses.tables/I.PAs.effects.birds.results.csv", sep=",")

res$Indice<-revalue(res$Indice, c("Richness"="Richness", "ForDep"="Forest-dependent species", "Endemism"="Endemic species", "IUCN"="Threatened species"))
res$Indice<-factor(res$Indice, c("Richness", "Forest-dependent species", "Endemic species", "Threatened species"))


res$Test<-factor(revalue(res$Test, c("Impl.only"="I.Impl", "Loc.Impl"="I.Full")), c("I.Impl", "I.Full"))

res$Initial<-revalue(res$zone,c(
  "Eastern.Afromontane"="EAS", "Ghats.Lanka"="GHA", "Tumbes"="TUM", "Tropical.Andes"="AND", "Atlantic.Forest"="ATL", "Indo-Burma"="IND", "Mesoamerica"="MES", "Sundaland"="SUN"))
res$Initial<-factor(res$Initial, c("ATL", "AND", "TUM", "MES", "EAS", "GHA", "IND", "SUN"))

theme_min_VCa <- theme_minimal() %+replace%  theme(strip.text=element_text(hjust=0.5, size=13, face="bold", colour="black", margin=margin(c(8,0,4,0))),
                                                   axis.title.y=element_text(angle=90, size=13),
                                                   panel.spacing.y=unit(-0.4, "cm"))





G<-ggplot(res[res$Test=="I.Impl",], aes(x=Initial))+
  geom_bar(aes(y=Coef), fill="#cb831fff", stat="identity", show.legend = FALSE)+
  geom_hline(yintercept=0, size=0.7, col="gray40")+
  geom_pointrange(aes(y=Coef, ymin=(Coef-1.96*SE), ymax=(Coef+1.96*SE)), fatten=0.5)+
  geom_text(aes(y=(Coef*1.2), label=gtools ::stars.pval(P.value)))+
  guides(fill = guide_legend(reverse=TRUE))+
  ylab("Coefficient")+
  xlab("")+
  facet_wrap(~Indice, ncol=1, scale="free_y")+
  theme_min_VCa




# With labels

cowplot::save_plot("D:/Victor/eBird/Figures/I.PAs.effects.birds.svg", 
                   G, base_width=5.5, base_height = 6)






