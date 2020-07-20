setwd("D:/Victor/eBird/Data/To export")
version="dec18"

library(auk) ; library(plyr) ; library(readr) ; library(nlme)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_)) 

#################################
### Setup region of expertise ###
#################################

### Can be Africa, Asia, America
continent="Asia"

if(continent=="America"){countries<-c("AG","AI","AR","BO","BR","BS","BZ","CL","CO","CR","CU","DM","DO","EC","FK","GD","GF","GP","GT","GY","HN","HT","JM","KN","KY","LC","MQ","MS","MX","NI","PA","PE","PR","PY","SR","SV","TC","UY","VC","VE","VG","VI")} # 42 countries
if(continent=="Asia")   {countries<-c("BD","BN","BT","CN","HK","ID","IN","KH","LA","LK","MM","MY","NP","PG","PH","PK","TH","TW","VN")} # 19 countries
if(continent=="Africa") {countries<-c("AO","BF","BI","BJ","BW","CD","CF","CG","CI","CM","DJ","DZ","EG","EH","ER","ET","GA","GH","GM","GN","GQ","GW","KE","LR","LS","LY","MA","MG","ML","MR","MW","MZ","NA","NE","NG","RW","SD","SL","SN","SO","SS","ST","SZ","TD","TG","TN","TZ","UG","ZA","ZM","ZW")} #51 countries

# All raw files must be put in a common folder called "To export"
files.countries<-list.files("D:/Victor/eBird/Data/To export", recursive=TRUE)
files.countries<-subset(files.countries, substr(files.countries, nchar(files.countries)-17, nchar(files.countries)-16) %in% countries)

## Check that I have the right number of files
length(files.countries) # Should be 42 (America), 19 (Asia) and 51 (Africa)
table(table(substr(files.countries, nchar(files.countries)-17, nchar(files.countries)-16)))         # Should be 1 (unique countries)
table(substr(files.countries, nchar(files.countries)-17, nchar(files.countries)-16) %in% countries) # Should be TRUE (countries of interest)




##################################
### EXPORT WITHOUT RESTRICTION ###
##################################

### Boucle pour tous les pays 
for (i in 1:length(files.countries)){
  
  eb<-read_ebd(files.countries[i], rollup=T, unique=F)

  # Checklist id (not the same than in script 2 because here we did not use the auk_unique function as we do not want to merge checklists of different observers)
  eb$checklist<-eb$sampling_event_identifier
  eb$checklist_id<-eb$sampling_event_identifier<-NULL
  
  ### Reduce the dataset
  # Remove useless columns
  eb[,c("state_province", "subnational1_code","subspecies_scientific_name","state","state_code","global_unique_identifier", "last_edited_date", "subspecies_common_name", "breeding_bird_atlas_category", "effort_area_ha", "county", "subnational2_code", "iba_code", "bcr_code", "usfws_code", "atlas_block", "locality", "locality_id", "locality_type", "first_name", "last_name", "has_media", "reviewed", "x")]<-NULL
  # Remove disapproved observations
  eb<-subset(eb, eb$approved==T)
  # Remove checklists that did not report all species detected
  eb<-subset(eb, eb$all_species_reported==T)
  # Remove old observations
  eb$year<-as.numeric(format(eb$observation_date, "%Y"))
  eb<-subset(eb, eb$year>=2005)
  
  
  saveRDS(eb, file=paste("D:/Victor/eBird/Data/Expertise.calculation/Exported.OBSERVER/", continent,"/ebird.",version,".Export.", substr(files.countries[i], nchar(files.countries[i])-17, nchar(files.countries[i])-16), ".rds", sep=""))
  
  cat(i)
  
}




#####################
### MERGING FILES ###
#####################

setwd(paste0("D:/Victor/eBird/Data/Expertise.calculation/Exported.OBSERVER/", continent))


df<-readRDS(list.files()[1]) # Stock the first table

if(length(list.files())>1){
for(i in 2:length(list.files())){
  df<-rbind(df, readRDS(list.files()[i])) # Add countries starting at 2
}
}

saveRDS(df, paste0("D:/Victor/eBird/Data/Expertise.calculation/Exported.OBSERVER/0.MERGED/Merged.data.Observer.", version, ".", continent, ".rds"))





###########################
### TABLE FOR EXPERTISE ###
###########################

### Checklists

chlist<-ddply(df, .(df$checklist), function(x){data.frame(
  date=x$observation_date[1],
  time=x$time_observations_started[1],
  distance=x$effort_distance_km[1],
  n.observers=x$number_observers[1],
  duration=x$duration_minutes[1],
  observer=x$observer_id[1],
  rich=nrow(x),
  lon=x$longitude[1],
  lat=x$latitude[1],
  protocol=x$protocol_type[1]
)}, .progress = "win")


# Format data
names(chlist)<-replace(names(chlist), names(chlist)=="df$checklist", "checklist")

chlist$day<-as.numeric(format(chlist$date, "%j"))

chlist$time.min<-sapply(strsplit(as.character(chlist$time),":"),
       function(x) {
         x <- as.numeric(x)
         x[1]*60+x[2] })


### Observer
observer.var<-ddply(df, .(df$observer_id), function(x){data.frame(Nb_obs=nrow(x), Nb_spc=nlevels(droplevels(as.factor(x$scientific_name))), Nb_checklist=nlevels(droplevels(as.factor(x$checklist))))})
names(observer.var)[1]<-"observer_id"

# Remove observers with less than 100 species or less than 10 checklists and 30 species per checklist
obs.to.exclude<-subset(observer.var, observer.var$Nb_spc<100 | observer.var$Nb_checklist<10 | (observer.var$Nb_obs/observer.var$Nb_checklist)<30)
chlist2<-subset(chlist, chlist$observer %not in% obs.to.exclude$observer_id)
chlist2$observer<-droplevels(chlist2$observer)


# Remove NAs
cat(paste(100*round(table(is.na(chlist2$duration))["TRUE"]/nrow(chlist2),2), " % of duration values are NA"))
chlist2<-chlist2[,c("checklist", "rich", "protocol", "duration", "n.observers", "time.min", "lon", "lat", "day", "observer")]
chlist3<-chlist2[complete.cases(chlist2),]
chlist3$observer<-droplevels(chlist3$observer)




#################################
### MAKE MODELS FOR EXPERTISE ###
#################################


# Modele
mod.obs<-mgcv::gamm(
  rich ~ protocol + n.observers + s(duration) + s(time.min) + te(lon, lat, day),
  random=list(observer= ~1),
  data=chlist3,
  family="poisson"
)




### Predict
df<-data.frame(
  observer=levels(droplevels(chlist3$observer)),
  protocol="Stationary",
  n.observers=median(chlist3$n.observers, na.rm=T),
  duration=median(chlist3$duration, na.rm=T),
  time.min=median(chlist3$time.min, na.rm=T),
  lon=median(chlist3$lon, na.rm=T),
  lat=median(chlist3$lat, na.rm=T),
  day=median(chlist3$day, na.rm=T)
)



# Extract fixed effect predictions
df$Pred2<-predict(mod.obs$gam, newdata=df)

# Add random effects (cf help of gamm)
rand2<-ranef(mod.obs$lme)

rand.observer2<-as.data.frame(rand2$observer)
rand.observer2$observer<-substr(rownames(rand.observer2), 7, 50)

df$rand.intercept2<-rand.observer2[match(df$observer, rand.observer2$observer), "(Intercept)"]

df$expertise<-df$Pred2+df$rand.intercept2




### SAVE THE EXPERTISE SCORE
saveRDS(df, paste0("D:/Victor/eBird/Data/Expand.Analyses.tables/obsKelling.",continent, ".", version, ".rds"))




