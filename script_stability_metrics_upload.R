# Application of invariability-area relationship (Wang et al 2017) and stability-synchrony framework across 
# ecological hierarchies (Wang et al 2019) to spatio-temporal data sets (Iowa lakes, MZB North Sea)
# Doro Hodapp
# May 2019


rm(list = ls())

# load packages
#library(codyn)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(abind)
library(pracma)
library(mgcv)


#### Invariability-area relationship (IAR) ----


# I(A) = 1/(CV^2(A)) = I_1 * (A/(A-1)*rho_A + 1) --> IAR is governed by patterns of spatial synchrony 
# across the landscape, which determines rho_A.

# Calculate I(A) for areas of increasing size and plot log-log IAR curve 


#### PHYTOPLANKTON SPECIES ----
# Load data
phyto <- read.csv(paste(in_path, "/Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
dim(phyto)

lakes <- unique(phyto$lake_name)
sampleID <- unique(phyto$sample_ID)

# function to calculate invariability
invar <- function(x){
  1/(sd(x, na.rm=T)/mean(x, na.rm=T))^2
}

#### a) Calculate invariability for areas of increasing extent by adding the next closer lake biomass ----
# Load spatial distance information
dist <- read.csv(paste(in_path, "/Iowa lakes/Analyses/lake_dist.csv", sep = ""))
dist_sub <- dist %>% 
  filter(lake1 %in% lakes & lake2 %in% lakes)
  
# Load geographic coordinate information
coord <- read.csv(paste(in_path, "/Iowa lakes/Analyses/lake_dist.csv", sep = ""))
coord_sub <- coord %>% 
  filter(lake1 %in% lakes & lake2 %in% lakes)

invar_dist_ls <- list() # create list for analysis output

# Calculate invariability(I) for each lake, then add the next closest lake, pool their species and 
# calculate I for the pooled community, for i.e. increasing extent of samples area 
for ( i in 1:length(lakes)){
  
  print(i)
  
  temp <- phyto[phyto$lake_name== lakes[i],]
  # calculate I for one lake
  I <- temp %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass, na.rm=T)) %>%
    dplyr::summarise(I = invar(totbiom))
  
  # select closest lake
  lake_i <- dist_sub[dist_sub$lake1 %in% lakes[i] | dist_sub$lake2 %in% lakes[i],]
  lakes_in <- as.factor(unlist(lake_i[which.min(lake_i$distance),c("lake1","lake2")]))
  
  inv <- phyto[phyto$lake_name %in% lakes_in,] %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
    dplyr::summarise(I = invar(totbiom))
  
  I <- append(I, inv)
  
  # select additional lake with the shortest average distance from all other lakes in the subset
  lakes_remain <- lakes[!(lakes %in% lakes_in)]
  mean_dist_remain <- data.frame(lake_name = lakes_remain)
  
  while(length(lakes_remain) > 1){
    length(lakes_remain)
    lakes_remain <- lakes[!(lakes %in% lakes_in)] # update what lakes are still out there
    mean_dist_remain <- data.frame(lake_name = lakes_remain) # create data frame
    
    for (j in 1:length(lakes_remain)){
      mean_dist_remain$mean_dist[j] <- mean(dist_sub$distance[(dist_sub$lake1 %in% lakes_remain[j] & dist_sub$lake2 %in% lakes_in)|(dist_sub$lake1 %in% lakes_in & dist_sub$lake2 %in% lakes_remain[j])],na.rm=T) 
    }    
    lakes_in <- as.factor(c(as.character(lakes_in),as.character(mean_dist_remain$lake_name[which.min(mean_dist_remain$mean_dist)]))) # add next lake to lake subset
    
    # calculate invariability (I) for increased lake sample
    
    inv <- phyto[phyto$lake_name %in% lakes_in,] %>% 
      gather(., key = species, value = biomass, ends_with(".bio")) %>%
      group_by(sample_ID) %>% 
      dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
      dplyr::summarise(I = invar(totbiom))
    
    I <- append(I, inv)
  }
  
  invar_dist_ls[[i]] <- I
  
}

saveRDS(invar_dist_ls, paste(out_path, "/Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_phyto.rds", sep = ""))

  
#### b) Calculate invariability for areas of increasing extent by adding lakes randomly ----

lakes <- unique(phyto$lake_name)

invar_rand_ls <- list() # create list for analysis output

# Calculate invariability(I) for each lake, then add a randomly selected lake, pool their species and 
# calculate I for the pooled community, for i.e. increasing extent of sampled area 
for ( i in 1:length(lakes)){
  
  print(i)
  
  temp <- phyto[phyto$lake_name == lakes[i],]
  # calculate I for one lake
  I <- temp %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
    dplyr::summarise(I = invar(totbiom))
  I <- as.data.frame(I)
  
  lakes_in <- lakes[i] # first lake
  lakes_remain <- lakes[!(lakes %in% lakes_in)] # all other lakes
  
  while(length(lakes_remain) > 0){
    # randomly select next lake to be added
    ind <- sample(length(lakes_remain),1)
    lakes_in <- as.factor(c(as.character(lakes_in),as.character(lakes_remain[ind])))
    
    inv <- phyto[phyto$lake_name %in% lakes_in,] %>% 
      gather(., key = species, value = biomass, ends_with(".bio")) %>%
      group_by(sample_ID) %>% 
      dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
      dplyr::summarise(I = invar(totbiom))
    
    I <- append(I, as.data.frame(inv))
    lakes_remain <- lakes[!(lakes %in% lakes_in)]
  }
  
  invar_rand_ls[[i]] <- I
  
}

saveRDS(invar_rand_ls,paste(out_path, "/Iowa lakes/Analyses/stability metrics/invar_rand_Iowa_phyto.rds", sep = ""))


#### ZOOPLANKTON SPECIES ----
# load data
zoo <- read.csv(paste(in_path,"/Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
dim(zoo)


# Load spatial distance information
dist <- read.csv(paste(in_path,"/Iowa lakes/Analyses/lake_dist.csv", sep = ""))
# Load geographic coordinate information
coord <- read.csv(paste(in_path,"/Iowa lakes/Analyses/lake_coords.csv", sep = ""))


# function to calculate invariability
invar <- function(x){
  1/(sd(x, na.rm=T)/mean(x, na.rm=T))^2
}


#### a) Calculate invariability for areas of increasing extent by adding the next closer lake biomass ----

lakes <- unique(zoo$lake_name)

invar_dist_ls <- list() # create list for analysis output

# Calculate invariability(I) for each lake, then add the next closest lake, pool their species and 
# calculate I for the pooled community, for i.e. increasing extent of samples area 
for ( i in 1:length(lakes)){
  
  print(i)
  
  temp <- zoo[zoo$lake_name== lakes[i],]
  # calculate I for one lake
  I <- temp %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
    dplyr::summarise(I = invar(totbiom))
  
  # select closest lake
  lake_i <- dist_sub[dist_sub$lake1 %in% lakes[i] | dist_sub$lake2 %in% lakes[i],]
  lakes_in <- as.factor(unlist(lake_i[which.min(lake_i$distance),c("lake1","lake2")]))
  
  inv <- zoo[zoo$lake_name %in% lakes_in,] %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
    dplyr::summarise(I = invar(totbiom))
  
  I <- append(I, inv)
  
  # select additional lake with the shortest average distance from all other lakes in the subset
  lakes_remain <- lakes[!(lakes %in% lakes_in)]
  mean_dist_remain <- data.frame(lake_name = lakes_remain)
  
  while(length(lakes_remain) > 1){
    lakes_remain <- lakes[!(lakes %in% lakes_in)] # update what lakes are still out there
    mean_dist_remain <- data.frame(lake_name = lakes_remain) # create data frame
    
    for (j in 1:length(lakes_remain)){
      mean_dist_remain$mean_dist[j] <- mean(dist_sub$distance[(dist_sub$lake1 %in% lakes_remain[j] & dist_sub$lake2 %in% lakes_in)|(dist_sub$lake1 %in% lakes_in & dist_sub$lake2 %in% lakes_remain[j])],na.rm=T) 
    }    
    lakes_in <- as.factor(c(as.character(lakes_in),as.character(mean_dist_remain$lake_name[which.min(mean_dist_remain$mean_dist)]))) # add next lake to lake subset
    
    
    # calculate invariability (I) for additional lake sample
    
    inv <- zoo[zoo$lake_name %in% lakes_in,] %>% 
      gather(., key = species, value = biomass, ends_with(".bio")) %>%
      group_by(sample_ID) %>% 
      dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
      dplyr::summarise(I = invar(totbiom))
    
    I <- append(I, inv)
  }
  
  invar_dist_ls[[i]] <- I
  
}

saveRDS(invar_dist_ls,paste(out_path,"/Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_zoo.rds", sep = ""))

 
#### b) Calculate invariability for areas of increasing extent by adding lakes randomly ----

lakes <- unique(zoo$lake_name)

invar_rand_ls <- list() # create list for analysis output

# Calculate invariability(I) for each lake, then add a randomly selected lake, pool their species and 
# calculate I for the pooled community, for i.e. increasing extent of samples area 
for ( i in 1:length(lakes)){
  
  print(i)
  
  temp <- zoo[zoo$lake_name == lakes[i],]
  # calculate I for one lake
  I <- temp %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
    dplyr::summarise(I = invar(totbiom))
  I <- as.data.frame(I)
  
  lakes_in <- lakes[i] # first lake
  lakes_remain <- lakes[!(lakes %in% lakes_in)] # all other lakes
  
  while(length(lakes_remain) > 0){
    # randomly select next lake to be added
    ind <- sample(length(lakes_remain),1)
    lakes_in <- as.factor(c(as.character(lakes_in),as.character(lakes_remain[ind])))
    
    inv <- zoo[zoo$lake_name %in% lakes_in,] %>% 
      gather(., key = species, value = biomass, ends_with(".bio")) %>%
      group_by(sample_ID) %>% 
      dplyr::summarise(totbiom = sum(biomass,na.rm=T)) %>%
      dplyr::summarise(I = invar(totbiom))
    
    I <- append(I, as.data.frame(inv))
    lakes_remain <- lakes[!(lakes %in% lakes_in)]
  }
  
  invar_rand_ls[[i]] <- I
  
}

saveRDS(invar_rand_ls,paste(out_path,"/Iowa lakes/Analyses/stability metrics/invar_rand_Iowa_zoo.rds", sep  = ""))


#### MACROZOBENTHOS SPECIES ----
## Load data
dat <- read.csv(paste(in_path,"/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""))
dim(dat)
head(dat)

stations <- unique(dat$Posno)
dat$Posno <- as.character(dat$Posno)

# Aggregate biomass across species

mzb <- dat %>% 
  group_by(Posno, Year) %>% 
  dplyr:: summarise(totbiom = sum(biom, na.rm = T))

#### a) Calculate invariability for areas of increasing extent by adding the next closer station biomass ----
# Load spatial distance information
dist <- read.csv(paste(in_path, "/dist_geo.csv", sep = ""))
dist_sub <- dist %>% 
  filter(Posno1 %in% stations & Posno2 %in% stations)

invar_dist_ls <- list() # create list for analysis output

# Calculate invariability(I) for each lake, then add the next closest lake, pool their species and 
# calculate I for the pooled community, for i.e. increasing extent of samples area 
for (i in 1:length(stations)){
  
  print(i)
  
  temp <- mzb[mzb$Posno == stations[i],]
  # calculate I for one station
  I <- temp %>% 
    summarise(I = invar(totbiom))
  
  # select closest lake
  station_i <- dist_sub[dist_sub$Posno1 %in% stations[i] | dist_sub$Posno2 %in% stations[i],]
  stations_in <- as.factor(unlist(station_i[which.min(station_i$distance),c("Posno1","Posno2")]))
  
  inv <- mzb[mzb$Posno %in% stations_in,] %>% 
    group_by(Year) %>% 
    dplyr::summarise(totbiom = sum(totbiom)) %>% 
    dplyr::summarise(I = invar(totbiom))
  
  I <- append(I, inv)
  
  # select additional lake with the shortest average distance from all other lakes in the subset
  stations_remain <- stations[!(stations %in% stations_in)]
  mean_dist_remain <- data.frame(Posno = stations_remain)
  
  while(length(stations_remain) > 1){
    stations_remain <- stations[!(stations %in% stations_in)] # update what lakes are still out there
    mean_dist_remain <- data.frame(Posno = stations_remain) # create data frame
    
    for (j in 1:length(stations_remain)){
      mean_dist_remain$mean_dist[j] <- mean(dist_sub$distance[(dist_sub$Posno1 %in% stations_remain[j] & dist_sub$Posno2 %in% stations_in)|(dist_sub$Posno1 %in% stations_in & dist_sub$Posno2 %in% stations_remain[j])],na.rm=T) 
    }    
    stations_in <- as.factor(c(as.character(stations_in),as.character(mean_dist_remain$Posno[which.min(mean_dist_remain$mean_dist)]))) # add next lake to lake subset
    
    
    # calculate invariability (I) for additional lake sample
    
    inv <- mzb[mzb$Posno %in% stations_in,] %>% 
      group_by(Year) %>% 
      dplyr::summarise(totbiom = sum(totbiom)) %>% 
      dplyr::summarise(I = invar(totbiom))
    
    I <- append(I, inv)
  }
  
  invar_dist_ls[[i]] <- I
  
}

saveRDS(invar_dist_ls,paste(out_path,"/stability metrics/invar_dist_NS_mzb.rds", sep =""))

#### b) Calculate invariability for areas of increasing extent by adding stations randomly ----
stations <- unique(mzb$Posno)

invar_rand_ls <- list() # create list for analysis output

# Calculate invariability(I) for each lake, then add a randomly selected lake, pool their species and 
# calculate I for the pooled community, for i.e. increasing extent of samples area 
for ( i in 1:length(stations)){
  
  print(i)
  
  temp <- mzb[mzb$Posno == stations[i],]
  # calculate I for one lake
  I <- temp %>% 
    summarise(I = invar(totbiom))
  
  stations_in <- stations[i] # first station
  stations_remain <- stations[!(stations %in% stations_in)] # all other stations
  
  while(length(stations_remain) > 0){
    
    # randomly select next station to be added
    ind <- sample(length(stations_remain),1)
    stations_in <- as.factor(c(as.character(stations_in),as.character(stations_remain[ind])))
    
    inv <- mzb[mzb$Posno %in% stations_in,] %>% 
      group_by(Year) %>% 
      dplyr::summarise(totbiom = sum(totbiom)) %>% 
      dplyr::summarise(I = invar(totbiom))
    
    I <- append(I, inv)
    stations_remain <- stations[!(stations %in% stations_in)]
  }
  
  invar_rand_ls[[i]] <- I
  
}

saveRDS(invar_rand_ls, paste(out_path, "/stability metrics/invar_rand_NS_mzb.rds", sep =""))










#---------------------------------------------------------------------------------------------------------------
###### Stability-synchrony analysis (Wang et al 2019) ######----------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# The var.partition command (Wang et al. 2019, Appendix 2) calculates synchrony and stability/variability metrics at 4 levels of 
# ecological hierarchy (population, community, metapopulation, metacommunity)
# As input it needs a 3-dimensional array (N*T*M, where N is the number of species, T is the number of time points, and M is the 
# number of local communities)

## PHYTOPLANKTON SPECIES---------------------------------------------------------------------------------------------------
# load data
phyto <- read.csv(paste(in_path, "/Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
dim(phyto)

lakes <- unique(phyto$lake_name)

# calculate data set characteristics
n.of.spec <- length(which(str_detect(names(phyto),".bio"))) # 132 species
n.of.sampl <- max(phyto$sample_ID)                          #  21 sampling occasions
n.of.lakes <- length(unique(phyto$lake_name))               #  49 lakes


# transform data frame to format that can then be convertet to 3d array
phyto_spec <- phyto %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

# convert data frame to 3d array
spec_array <- abind(split(phyto_spec[,3:ncol(phyto_spec)],phyto_spec$lake_name), along = 3)  


# calculate variability and synchrony measures across levels of organization

source(paste(in_path, "/Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/function_var_partition.R", sep =""))

sync_var_spec <- var.partition(spec_array)
sync_var_spec

## ZOOPLANKTON SPECIES--------------------------------------------------------------------------------------------------

zoo <- read.csv(paste(in_path, "/Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
dim(zoo)


# calculate biomass values for the different functional groups
# transform data frame to format that can then be convertet to 3d array
zoo_spec <- zoo %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

# convert data frame to 3d array
zoo_array <- abind(split(zoo_spec[,3:ncol(zoo_spec)],zoo_spec$lake_name), along = 3) 

# calculate variability and synchrony measures across levels of organization

source(paste(in_path, "/Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/function_var_partition.R", sep =""))

sync_var_zoo <- var.partition(zoo_array)
sync_var_zoo


## MACROZOOBENTHOS SPECIES --------------------------------------------------------------------------------------------------
## SPECIES
# select relevant information from data file
data <- mzb_sub %>% 
  select(Posno, Year, species, biom) %>% 
  spread(., key = Year, value = biom)



data_spec <- as.data.frame(data)
data_spec[is.na(data_spec) == T] <- 0
data_spec$Posno <- data_spec$Posno[,drop=T]

# convert data frame to 3d array
MZB_array <- abind(split(data_spec[,3:ncol(data_spec)],data_spec$Posno), along = 3)

# calculate variability and synchrony measures across levels of organization

source(paste(in_path,"/Analyses/stability metrics/function_var_partition.R", sep = ""))

sync_var_spec <- var.partition(MZB_array)
sync_var_spec

#### Prepare data files for upload to repository ####
spec_names_phyto <- paste0(c(rep("Species")),1:87)   
names(phyto) <- c("lake_ID","sample_ID", spec_names_phyto)

# Convert to long format, because Dryad always gave out an error message when trying to upload the data files in wide format.. 
## PHYTOPLANKTON
phyto_long <- phyto %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Species", values_to = "Biomass")

write.csv(phyto_long, paste(in_path, "/Iowa lakes/Data/Iowa_lake_phyto_upload.csv", sep = ""), row.names = F)

## ZOOPLANKTON
spec_names_zoo <- paste0(c(rep("Species")),1:21)
names(zoo) <- c("lake_ID","sample_ID", spec_names_zoo)

zoo_long <- zoo %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Species", values_to = "Biomass")

write.csv(zoo_long, paste(in_path, "/Iowa lakes/Data/Iowa_lake_zoo_upload.csv", sep = ""), row.names = F)

## MACROZOOBENTHOS
mzb <- dat %>% 
  select(Posno, Year, species, biom) %>% 
  pivot_wider(names_from = species, values_from = biom)

spec_names_mzb <- paste0(c(rep("Species")),1:195)
names(mzb) <- c("Station_ID","Year", spec_names_mzb)

mzb_long <- mzb %>% 
  pivot_longer(cols = 3:ncol(.), names_to = "Species", values_to = "Biomass")

write.csv(mzb_long,paste(in_path, "/Iowa lakes/Data/North_Sea_mzb_upload.csv", sep = ""), row.names = F)




