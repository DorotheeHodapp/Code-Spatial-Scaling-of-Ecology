#### Analysis on contribution of single sites to overall synchrony-stability relationships, IAR, etc.. 
#### i.e. How much does each single site contribute to overall stability of regional ecosystem functioning?
#### Doro Hodapp 

rm(list = ls())

# load packages

library(plyr)
#library(codyn)
library(tidyr)
library(stringr)
library(ggplot2)
library(abind)
library(pracma)
library(vegan)
library(gridExtra)
library(RColorBrewer)
library(dplyr)

##### Choose what species group to run...
type <- "phyto" # choose "phyto", "zoo" or "mzb" here

##### Load data #####
if(type == "phyto"){
  data <- read.csv(paste(in_path,"Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
  }else if(type == "zoo"){
  data <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
  }else{
  data <- read.csv(paste(in_path, "Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""))
}


##### Data preparation for analyses ----
if (type == "phyto" | type == "zoo"){
  
    lakes <- unique(data$lake_name)
    
    data_long <- data %>% 
      gather(., key = species, value = biomass, ends_with(".bio"))
    
    data_long <- as.data.frame(data_long)
    
    species_all <- unique(data_long$species)
  
  }else if(type == "mzb"){
  
  dim(data)
  head(data)
  
  stations <- unique(data$Posno)
  
  species_all <- unique(data$species)
  
} 



#------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### (1a) CONTRIBUTION OF EACH LAKE/ STATION TO THE SHAPE OF THE IAR #####

# function to calculate invariability
invar <- function(x){
  1/(sd(x, na.rm=T)/mean(x, na.rm=T))^2
}

# FOR EACH LAKE/ STATION: Calculate invariability for areas of increasing extent by adding the next closer lake/station biomass


if(type == "phyto" | type == "zoo"){
  
  # Load spatial distance information
  dist <- read.csv(paste(in_path,"Iowa lakes/Analyses/lake_dist.csv", sep = ""))
  dist <- dist %>% 
    filter(lake1 %in% lakes &lake2 %in% lakes)
  
  # Load geographic coordinate information
  coord <- read.csv(paste(in_path, "Iowa lakes/Analyses/lake_coords.csv", sep = ""))
  
  invar_dist_ls <- list() # create list for analysis output
  
  data_long <- data %>%
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    dplyr::group_by(lake_name, sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass, na.rm = T))
  
  
  data_long <- as.data.frame(data_long)
  lakes <- unique(data_long$lake_name)

  for (k in 1:length(lakes)){
    
    print(k)
    data_sub <- data_long[data_long$lake_name != lakes[k],]
    #data_sub$lake_name <- data_sub$lake_name[drop=T,]
    lakes_sub <- unique(data_sub$lake_name)
    dist_sub <- dist[dist$lake1 %in% lakes_sub & dist$lake2 %in% lakes_sub,]
    
    # Calculate invariability(I) for each lake, then add the next closest lake, pool their species and 
    # calculate I for the pooled community, for i.e. increasing extent of samples area 
    for (i in 1:length(lakes_sub)){
      
      print(i)
      temp <- data_sub[data_sub$lake_name == lakes_sub[i],]
      
      # calculate I for one lake
      I <- temp %>% 
        summarise(I = invar(totbiom))
      
      # select closest lake
      lake_i <- dist_sub[dist_sub$lake1 %in% lakes_sub[i] | dist_sub$lake2 %in% lakes_sub[i],]
      lakes_in <- as.factor(unlist(lake_i[which.min(lake_i$distance),c("lake1","lake2")]))
      
      inv <- data_sub[data_sub$lake_name %in% lakes_in,] %>% 
        group_by(sample_ID) %>% 
        dplyr::summarise(tbiom = sum(totbiom,na.rm=T)) %>%
        dplyr::summarise(I = invar(tbiom))
      
      I <- append(I, inv)
      
      # select additional lake with the shortest average distance from all other lakes in the subset
      lakes_remain <- lakes_sub[!(lakes_sub %in% lakes_in)]
      
      while(length(lakes_remain) > 1){
        lakes_remain <- lakes_sub[!(lakes_sub %in% lakes_in)] # update what lakes are still out there
        mean_dist_remain <- data.frame(lake_name = lakes_remain) # create data frame
        
        for (j in 1:length(lakes_remain)){
          mean_dist_remain$mean_dist[j] <- mean(dist_sub$distance[(dist_sub$lake1 %in% lakes_remain[j] & dist_sub$lake2 %in% lakes_in)|(dist_sub$lake1 %in% lakes_in & dist_sub$lake2 %in% lakes_remain[j])],na.rm=T) 
        }    
        lakes_in <- as.factor(c(as.character(lakes_in),as.character(mean_dist_remain$lake_name[which.min(mean_dist_remain$mean_dist)]))) # add next lake to lake subset
        
        # calculate invariability (I) for increased lake sample
        
        inv <- data_sub[data_sub$lake_name %in% lakes_in,] %>% 
          group_by(sample_ID) %>% 
          dplyr::summarise(tbiom = sum(totbiom,na.rm=T)) %>%
          dplyr::summarise(I = invar(tbiom))
        
        I <- append(I, inv)
      }
      
      invar_dist_ls[[i]] <- I
    }
    
    saveRDS(invar_dist_ls,paste(out_path, "Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_", type, "_",lakes[k],".rds", sep=""))
  }
  
#### MZB ----
}else if(type == "mzb"){
  
  # Load spatial distance information
  
  dist <- read.csv(paste(in_path, "Armonies data/dist_geo.csv", sep = ""))
  dist <- dist %>% 
    filter(Posno1 %in% stations & Posno2 %in% stations) 
  # Load geographic coordinate information
  #coord <- read.csv("C:/Users/dhodapp/Documents/HIFMB/Projects/Spatio-temporal Turnover/Iowa lakes/Analyses/lake_coords.csv")
  
  dat_mzb <- data %>% 
    group_by(Posno, Year) %>% 
    dplyr::summarise(totbiom = sum(biom,na.rm=T))
  
  
  dat_mzb <- as.data.frame(dat_mzb)
  dat_mzb$Posno <- as.factor(dat_mzb$Posno)
  
  stations <- unique(dat_mzb$Posno)
  
  invar_dist_ls <- list() # create list for analysis output
  
  
  for (k in 1:length(stations)){   
    
    print(k)
    mzb <- dat_mzb[dat_mzb$Posno != stations[k],]
    mzb$Posno <- mzb$Posno[drop=T,]
    stations_sub <- unique(mzb$Posno)
    dist_sub <- dist[dist$Posno1 %in% stations_sub & dist$Posno2 %in% stations_sub,]
    
    # Calculate invariability(I) for each lake, then add the next closest lake, pool their species and 
    # calculate I for the pooled community, for i.e. increasing extent of samples area 
    for (i in 1:length(stations_sub)){
      
      print(i)

      temp <- mzb[mzb$Posno == stations_sub[i],]
      
      # calculate I for one station
      I <- temp %>%
        summarise(I = invar(totbiom))
      
      # select closest lake
      station_i <- dist_sub[dist_sub$Posno1 %in% stations_sub[i] | dist_sub$Posno2 %in% stations_sub[i],]
      stations_in <- as.factor(unlist(station_i[which.min(station_i$distance),c("Posno1","Posno2")]))
      
      inv <- mzb[mzb$Posno %in% stations_in,] %>% 
        group_by(Year) %>% 
        dplyr::summarise(totbiom = sum(totbiom)) %>% 
        dplyr::summarise(I = invar(totbiom))
      
      I <- append(I, inv)
      
      # select additional station with the shortest average distance from all other stations in the subset
      stations_remain <- stations_sub[!(stations_sub %in% stations_in)]
      mean_dist_remain <- data.frame(station = stations_remain)
      
      while(length(stations_remain) > 1){
        stations_remain <- stations_sub[!(stations_sub %in% stations_in)] # update names of stations that have not been included yet
        mean_dist_remain <- data.frame(station = stations_remain) # create data frame
        
        for (j in 1:length(stations_remain)){
          mean_dist_remain$mean_dist[j] <- mean(dist_sub$distance[(dist_sub$Posno1 %in% stations_remain[j] & dist_sub$Posno2 %in% stations_in)|(dist_sub$Posno1 %in% stations_in & dist_sub$Posno2 %in% stations_remain[j])],na.rm=T) 
        }
        stations_in <- as.factor(c(as.character(stations_in),as.character(mean_dist_remain$station[which.min(mean_dist_remain$mean_dist)]))) # add next lake to lake subset
        
        # calculate invariability (I) for increased lake sample
        
        inv <- mzb[mzb$Posno %in% stations_in,] %>% 
          group_by(Year) %>% 
          dplyr::summarise(totbiom = sum(totbiom)) %>% 
          dplyr::summarise(I = invar(totbiom))
        
        I <- append(I, inv)
      }
      
      invar_dist_ls[[i]] <- I
      
    }
    
    saveRDS(invar_dist_ls, paste(out_path,"Armonies data/Analyses/stability metrics/invar_dist_NS_", type, "_",stations[k],".rds", sep=""))
  }
  
}


#### compare IAR curves between all lakes and lake subsets ----

library('readr')
library('dplyr')
library('ggplot2')
theme_set(theme_bw())
library('mgcv')

if(type == "phyto" | type =="zoo"){

# Load IAR all lakes results
all_contr <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_", type,".rds", sep = ""))
all_m <- do.call("rbind", all_contr) # convert list entries to matrix by rbinding all list elements
all_m <- matrix((as.numeric(all_m)), nrow = length(lakes), ncol = length(lakes))

all_means <- colMeans(all_m) # calculate mean invariability values across all lakes


# load IAR lake subset results (always one lake missing)

single_means <- list()

for (k in 1:length(lakes)){
single_contr <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_", type, "_",lakes[k],".rds", sep=""))

single_m <- do.call("rbind", single_contr) # convert list entries to matrix by rbinding all list elements
single_m <- matrix((as.numeric(single_m)),nrow = length(lakes)-1, ncol = length(lakes)-1)

single_means[[k]] <- colMeans(single_m)
}


# Calculate correlation between IAR based on all lakes and IAR based on lake subsets (always on lake missing) 
# in order to estimate influence of each lake on overall IAR
IAR_corr <- data.frame(lake_name = lakes, corr = NA)

for (i in 1:length(lakes)){
  IAR_corr$corr[i] <- cor(all_means[1:(length(lakes)-1)],single_means[[i]],method = "spearman")
  }

write.csv(IAR_corr, paste(out_path, "Iowa lakes/Analyses/stability metrics/IAR_corr_patch_", type, ".csv", sep = ""))

lake_corr_lower_1 <- IAR_corr$lake_name[IAR_corr$corr < 0.99]

# create data frame for plotting
IAR_means <- data.frame(lake = c(rep("all",length(lakes)), rep(as.character(lakes), each = length(lakes)-1)), area = c(1:length(lakes), rep(1:(length(lakes)-1), length(lakes))), I = c(all_means, as.numeric(do.call("cbind", single_means))))

IAR_means_sub <- IAR_means %>% 
  dplyr::filter(lake %in% lake_corr_lower_1)

# plot data 
png(paste(out_path,"Iowa lakes/Analyses/Graphs/IAR_patch_contr_", type, ".png", sep =""), width = 800, height = 800, res = 150)
ggplot(IAR_means, aes(x = area, y = I, color = lake)) +
#  geom_point(aes(color = lake)) +
  geom_point(colour = "black", shape=21, size = 2, alpha = 0.8,
             aes(fill = lake)) +
  geom_smooth(method = 'loess', se = FALSE) +
  theme(legend.position = "none") +
  #theme(legend.text = element_text(size=10)) + 
  #theme(legend.title = element_text(size = 12)) + 
  ylab("Invariability") + 
  xlab("Area (number of samples)") +
  theme(plot.title = element_text(size=16), axis.title = element_text(size=14)) +
  ggtitle(paste("Single lake contribution to IAR \n(", type,"plankton subset)", sep = ""))
dev.off()

}else if(type == "mzb"){
# Load IAR all stations results
all_contr <- readRDS(paste(in_path, "Armonies data/Analyses/stability metrics/invar_dist_NS_mzb.rds", sep = ""))

all_m <- do.call("rbind", all_contr) # convert list entries to matrix by rbinding all list elements
all_m <- matrix((as.numeric(all_m)), nrow = length(stations), ncol = length(stations))

all_means <- colMeans(all_m) # calculate mean invariability values across all lakes


# Load IAR station subset results (always one station missing)
single_means <- list()

for (k in 1:length(stations)){  
  single_contr <- readRDS(paste(out_path, "Armonies data/Analyses/stability metrics/invar_dist_NS_", type, "_", stations[k],".rds", sep=""))
  
  single_m <- do.call("rbind", single_contr) # convert list entries to matrix by rbinding all list elements
  single_m <- matrix((as.numeric(single_m)), nrow = length(stations)-1, ncol = length(stations)-1)
  
  single_means[[k]] <- colMeans(single_m, na.rm=T)
}

station_name <- rep(NA,length(stations))
for (i in 1:length(stations)){
  station_name[i] <- paste("Station",stations[i],sep="")
}



# Calculate correlation between IAR based on all lakes and IAR based on lake subsets (always on station missing) 
# in order to estimate influence of each lake on overall IAR
IAR_corr <- data.frame(station = station_name, corr = NA)

for (i in 1:length(stations)){
  IAR_corr$corr[i] <- cor(all_means[1:(length(stations)-1)],single_means[[i]],method = "spearman")
}

write.csv(IAR_corr, paste(out_path, "Armonies data/Analyses/stability metrics/IAR_corr_patch_", type,".csv", sep = ""))

station_corr_lower_1 <- IAR_corr$station[IAR_corr$corr < 0.4]

# create data frame for plotting

IAR_means <- data.frame(station = c(rep("all",length(stations)), rep(as.character(station_name), each = length(stations)-1)), area = c(1:length(stations), rep(1:(length(stations)-1), length(stations))), I = c(all_means, as.numeric(do.call("cbind", single_means))))

IAR_means_sub <- IAR_means %>% 
  filter(station %in% station_corr_lower_1)

# plot data 
png(paste(out_path,"Armonies data/Analyses/Graphs/IAR_patch_contr_", type, ".png", sep =""), width = 800, height = 800, res = 150)
ggplot(IAR_means, aes(x = area, y = I, colour = station)) +
  geom_point(colour = "black", shape=21, size = 2, alpha = 0.8,
             aes(fill = station)) +
  geom_smooth(method = 'loess', se = FALSE) +
  theme(legend.position = "none") +
  #theme(legend.text = element_text(size=10)) + 
  #theme(legend.title = element_text(size = 12)) + 
  ylab("Invariability") + 
  xlab("Area (number of samples)") +
  theme(plot.title = element_text(size=16), axis.title = element_text(size=14)) + 
  ggtitle(paste("Single station contribution to IAR \n(", type,")", sep = ""))
dev.off()


}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### (1b) CONTRIBUTION OF EACH SPECIES TO THE SHAPE OF THE IAR #####

# function to calculate invariability
invar <- function(x){
  1/(sd(x, na.rm=T)/mean(x, na.rm=T))^2
}


if(type == "phyto" | type == "zoo"){
  
  # Load spatial distance information
  dist <- read.csv(paste(in_path, "Iowa lakes/Analyses/lake_dist.csv", sep = ""))
  dist_sub <- dist %>% 
    filter(lake1 %in% lakes & lake2 %in% lakes)
  
  # Load geographic coordinate information
  #coord <- read.csv(paste(in_path, "Iowa lakes/Analyses/lake_coords.csv", sep = ""))
  
  data_long <- data %>% 
    gather(., key = species, value = biomass, ends_with(".bio"))
  
  data_long <- as.data.frame(data_long)
  lakes <- unique(data_long$lake_name)
  
  species_all <- unique(data_long$species)

  invar_dist_ls <- list() # create list for analysis output
  
  for (k in 1:length(species_all)){
    
    print(k)
    
    data_sub <- data_long %>% 
      dplyr::filter(species != species_all[k]) %>% 
      group_by(lake_name, sample_ID) %>% 
      dplyr::summarise(totbiom = sum(biomass, na.rm=T)) %>% 
      ungroup()
    

    # Calculate invariability(I) for each lake, then add the next closest lake, pool their species and 
    # calculate I for the pooled community, for i.e. increasing extent of samples area 
    for (i in 1:(length(lakes))){
      
      print(i)
      temp <- data_sub[data_sub$lake_name == lakes[i],]
      
      # calculate I for one lake
      I <- temp %>% 
      #  dplyr::select(sample_ID, totbiom) %>% 
        dplyr::summarise(I = invar(totbiom))
      
      # select closest lake
      lake_i <- dist_sub[dist_sub$lake1 %in% lakes[i] | dist_sub$lake2 %in% lakes[i],]
      lakes_in <- as.factor(unlist(lake_i[which.min(lake_i$distance),c("lake1","lake2")]))
      
      inv <- data_sub[data_sub$lake_name %in% lakes_in,] %>% 
        group_by(sample_ID) %>% 
        dplyr::summarise(tbiom = sum(totbiom,na.rm=T)) %>%
        dplyr::summarise(I = invar(tbiom))
      
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
        
        # calculate invariability (I) for increased lake sample
        
        inv <- data_sub[data_sub$lake_name %in% lakes_in,] %>% 
          group_by(sample_ID) %>% 
          dplyr::summarise(tbiom = sum(totbiom,na.rm=T)) %>%
          dplyr::summarise(I = invar(tbiom))
        
        I <- append(I, inv)
      }
      
      invar_dist_ls[[i]] <- I
    }
    
    saveRDS(invar_dist_ls,paste(out_path, "Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_", type, "_",species_all[k],".rds", sep=""))
  }
  
}else if(type == "mzb"){
  
  # Load spatial distance information
  
  dist <- read.csv(paste(in_path, "Armonies data/dist_geo.csv", sep = ""))
  dist_sub <- dist %>% 
  filter(Posno1 %in% stations & Posno2 %in% stations)
  # Load geographic coordinate information
  #coord <- read.csv("C:/Users/dhodapp/Documents/HIFMB/Projects/Spatio-temporal Turnover/Iowa lakes/Analyses/lake_coords.csv")
  
  
  invar_dist_ls <- list() # create list for analysis output
  
  species_all <- unique(data$species)
  
  for (k in 1:length(species_all)){
    
    print(k)
    
    dat_mzb <- data %>% 
      filter(species != species_all[k]) %>%   # exclude one species from species assemblage
      group_by(Posno, Year) %>% 
      dplyr::summarise(totbiom = sum(biom,na.rm=T)) # calculate total biomass per station per sampling occasion

    dat_mzb <- as.data.frame(dat_mzb)
    stations <- unique(dat_mzb$Posno)
    
    # Calculate invariability(I) for one lake, then add the next closest lake, pool their species and 
    # calculate I for the pooled community, for i.e. increasing extent of sampling area 
    for (i in 1:length(stations)){
      
      print(i)
      
      temp <- dat_mzb[dat_mzb$Posno == stations[i],]
      
      # calculate I for one station
      I <- temp %>%
        summarise(I = invar(totbiom))
      
      # select closest lake
      station_i <- dist_sub[dist_sub$Posno1 %in% stations[i] | dist_sub$Posno2 %in% stations[i],]
      stations_in <- as.factor(unlist(station_i[which.min(station_i$distance),c("Posno1","Posno2")]))
      
      inv <- dat_mzb[dat_mzb$Posno %in% stations_in,] %>% 
        group_by(Year) %>% 
        dplyr::summarise(totbiom = sum(totbiom,na.rm = T)) %>% 
        dplyr::summarise(I = invar(totbiom))
      
      I <- append(I, inv)
      
      # select additional station with the shortest average distance from all other stations in the subset
      stations_remain <- stations[!(stations %in% stations_in)]
      
      while(length(stations_remain) > 1){
        stations_remain <- stations[!(stations %in% stations_in)] # update names of stations that have not been included yet
        mean_dist_remain <- data.frame(station = stations_remain) # create data frame
        
        for (j in 1:length(stations_remain)){
          mean_dist_remain$mean_dist[j] <- mean(dist_sub$distance[(dist_sub$Posno1 %in% stations_remain[j] & dist_sub$Posno2 %in% stations_in)|(dist_sub$Posno1 %in% stations_in & dist_sub$Posno2 %in% stations_remain[j])],na.rm=T) 
        }
        stations_in <- as.factor(c(as.character(stations_in),as.character(mean_dist_remain$station[which.min(mean_dist_remain$mean_dist)]))) # add next lake to lake subset
        
        # calculate invariability (I) for increased lake sample
        
        inv <- dat_mzb[dat_mzb$Posno %in% stations_in,] %>% 
          group_by(Year) %>% 
          dplyr::summarise(totbiom = sum(totbiom)) %>% 
          dplyr::summarise(I = invar(totbiom))
        
        I <- append(I, inv)
      }
      
      invar_dist_ls[[i]] <- I
      
    }
    
    saveRDS(invar_dist_ls, paste(out_path, "Armonies data/Analyses/stability metrics/invar_dist_NS_", type, "_",species_all[k],".rds", sep=""))
  }
  
}

##-----------------------------------------------------------------------------------------
#### compare IAR curves between all species and species subsets ----

library('readr')
library('dplyr')
library('ggplot2')
theme_set(theme_bw())
library('mgcv')

if(type == "phyto" | type =="zoo"){
  # Load IAR all lakes results
  all_contr <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_", type,".rds", sep = ""))
  all_m <- do.call("rbind", all_contr) # convert list entries to matrix by rbinding all list elements
  all_m <- matrix((as.numeric(all_m)), nrow = length(lakes), ncol = length(lakes))
  
  all_means <- colMeans(all_m) # calculate mean invariability values across all lakes
  
  
  
  # Load IAR lake subset results (always one species missing)
  
  single_means <- list()
  
  for (k in 1: length(species_all)){   
    single_contr <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/invar_dist_Iowa_", type, "_",species_all[k],".rds", sep=""))

    single_m <- do.call("rbind", single_contr) # convert list entries to matrix by rbinding all list elements
    single_m <- matrix((as.numeric(single_m)), nrow = length(lakes), ncol = length(lakes))
    single_means[[k]] <- colMeans(single_m)
  }
  
  lake_name <- rep(NA,length(lakes))
  for (i in 1:length(lakes)){
    lake_name[i] <- paste(lakes[i],"_",sep="")
  }
  
  # Calculate correlation between IAR based on all species and IAR based on species subsets (always one species missing) 
  # in order to estimate influence of each lake on overall IAR
  IAR_corr <- data.frame(species = species_all, corr = NA)
  
  for (i in 1:length(species_all)){
    IAR_corr$corr[i] <- cor(all_means,single_means[[i]],method = "spearman")
  }
  
  write.csv(IAR_corr, paste(out_path, "Iowa lakes/Analyses/stability metrics/IAR_corr_species_", type, ".csv", sep = ""))
  
  spec_corr_lower_1 <- IAR_corr$species[IAR_corr$corr < 0.9999999]
  
  # create data frame for plotting
  # for readability reasons: only plot species with correlation between full data set and species subset < 1 
  IAR_means <- data.frame(species = c(rep("all",length(lakes)), rep(as.character(species_all), each = length(lakes))), area = c(1:length(lakes), rep(1:length(lakes), length(species_all))), I = c(all_means, as.numeric(do.call("cbind", single_means))))
  
  IAR_means_sub <- IAR_means %>% 
    filter(IAR_means$species %in% spec_corr_lower_1)
  
  # plot data 
  png(paste(out_path, "Iowa lakes/Analyses/Graphs/IAR_species_contr_", type, ".png", sep =""), width = 800, height = 800, res = 150)
  ggplot(IAR_means, aes(x = area, y = I, colour = species)) +
    geom_point(colour = "black", shape=21, size = 2, alpha = 0.8,
               aes(fill = species)) +
    geom_smooth(method = 'loess', se = FALSE) +
    theme(legend.position = "none") +
    #theme(legend.text = element_text(size=10)) + 
    #theme(legend.title = element_text(size = 12)) + 
    ylab("Invariability") + 
    xlab("Area (number of samples)") +
    theme(plot.title = element_text(size=16), axis.title = element_text(size=14)) + 
    ggtitle(paste("Single species contribution to IAR \n(", type,"plankton subset)", sep = ""))
  dev.off()
  
  
  
  
}else if(type == "mzb"){
  
  # load IAR all stations results
  all_contr <- readRDS(paste(in_path, "Armonies data/Analyses/stability metrics/invar_dist_NS_mzb.rds", sep = ""))
  
  all_m <- do.call("rbind", all_contr) # convert list entries to matrix by rbinding all list elements
  all_m <- matrix((as.numeric(all_m)), nrow = length(stations), ncol = length(stations))
  
  all_means <- colMeans(all_m) # calculate mean invariability values across all lakes
  
  
  # load IAR station subset results (always one station missing)
  single_means <- list()
  
  for (k in 1:length(species_all)){  
    single_contr <- readRDS(paste(in_path, "Armonies data/Analyses/stability metrics/invar_dist_NS_", type, "_", species_all[k],".rds", sep=""))
    
    single_m <- do.call("rbind", single_contr) # convert list entries to matrix by rbinding all list elements
    single_m <- matrix((as.numeric(single_m)), nrow = length(stations), ncol = length(stations))
    
    single_means[[k]] <- colMeans(single_m, na.rm=T)
  }
  
  station_name <- rep(NA,length(stations))
  for (i in 1:length(stations)){
    station_name[i] <- paste("Station", stations[i], sep="")
  }
  
  
  # Calculate correlation between IAR based on all species and IAR based on species subsets (always one species missing) 
  # in order to estimate influence of each lake on overall IAR
  IAR_corr <- data.frame(species = species_all, corr = NA)
  
  for (i in 1:length(species_all)){
    IAR_corr$corr[i] <- cor(all_means, single_means[[i]], method = "spearman")
  }
  
  
  write.csv(IAR_corr, paste(out_path, "Armonies data/Analyses/stability metrics/IAR_corr_species_", type, ".csv", sep = ""))
  
  spec_corr_lower_1 <- IAR_corr$species[IAR_corr$corr < 0.9999999]
  
  # create data frame for plotting
  
  IAR_means <- data.frame(species = c(rep("all",length(stations)), rep(as.character(species_all), each = length(stations))), area = c(1:length(stations), rep(1:(length(stations)), length(species_all))), I = c(all_means, as.numeric(do.call("cbind", single_means))))
  
  IAR_means_sub <- IAR_means %>% 
    filter(IAR_means$species %in% spec_corr_lower_1)
  
  # plot data 
  png(paste(out_path, "Armonies data/Analyses/Graphs/IAR_species_contr_", type, ".png", sep =""), width = 800, height = 800, res = 150)
  ggplot(IAR_means, aes(x = area, y = I, colour = species)) +
    geom_point(colour = "black", shape=21, size = 2, alpha = 0.8,
               aes(fill = species)) +
    geom_smooth(method = 'loess', se = FALSE) +
    #theme(legend.position = "none") +
    #theme(legend.text = element_text(size=10)) + 
    #theme(legend.title = element_text(size = 12)) + 
    ylab("Invariability") + 
    xlab("Area (number of samples)") +
    theme(plot.title = element_text(size=16), axis.title = element_text(size=14)) + 
    ggtitle(paste("Single species contribution to IAR \n(", type,")", sep = ""))
  dev.off()
    
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### (3) CONTRIBUTION OF EACH LAKE/ STATION TO STAB-SYNC FRAMEWORK (Wang et al 2019) #####
## --> calculate stab-syn framework metrics for subsets excluding one species, respectively site from analysis and calculate 
## log(phi'/phi) --> positive values indicate stabilizing effects (Lamy et al. 2019)

## PHYTOPLANKTON SPECIES ----------------
# load data
phyto <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
dim(phyto)

source("function_var_partition.R")

## all species and all sitees
# transform data frame to format that can then be converted to 3d array
phyto_spec <- phyto %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

# convert data frame to 3d array
spec_array <- abind(split(phyto_spec[,3:ncol(phyto_spec)],phyto_spec$lake_name), along = 3)  


# calculate variability and synchrony measures across levels of organization
sync_var_all <- var.partition(spec_array)

#### A) SINGLE SPECIES CONTRIBUTION
## all species minus 1
contr_each_spec_phyto <- list()

# transform data frame to format that can then be convertet to 3d array
phyto_spec <- phyto %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

spec_phyto <- unique(phyto_spec$species)

for (i in 1:length(spec_phyto)){
  
  print(i)
  phyto_spec_sub <- phyto_spec %>% 
    filter(species != spec_phyto[i])
  
  
  # convert data frame to 3d array
  spec_array <- abind(split(phyto_spec_sub[,3:ncol(phyto_spec_sub)],phyto_spec_sub$lake_name), along = 3)  
  
  # calculate variability and synchrony measures across levels of organization
  
  sync_var_spec_sub <- var.partition(spec_array)

  contr_each_spec_phyto[[i]] <- log(sync_var_spec_sub[8]/sync_var_all[8])     # 8: phi_S2C_R 
  
}


## B) SINGLE SITE CONTRIBUTION

## all lakes minus 1
contr_each_lake_phyto <- list()
lakes <- unique(phyto$lake_name)

# transform data frame to format that can then be convertet to 3d array
phyto_lake <- phyto %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

for (i in 1:length(lakes)){
  
  print(i)
  phyto_lake_sub <- phyto_lake %>% 
    filter(lake_name != lakes[i])
  phyto_lake_sub$lake_name <- phyto_lake_sub$lake_name[,drop = T]
  
  # convert data frame to 3d array
  lake_array <- abind(split(phyto_lake_sub[,3:ncol(phyto_lake_sub)],phyto_lake_sub$lake_name), along = 3)  
 
  
  # calculate variability and synchrony measures across levels of organization
  
  sync_var_lake_sub <- var.partition(lake_array)
  
  contr_each_lake_phyto[[i]] <- log(sync_var_lake_sub[6]/sync_var_all[6])     # 6: phi_C_L2R 
  
}

contr_each_lake_phyto

## ZOOPLANKTON SPECIES -----------------

source("function_var_partition.R")

# load data
zoo <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep  = ""))
dim(zoo)

## all species
# transform data frame to format that can then be converted to 3d array
zoo_spec <- zoo %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

# convert data frame to 3d array
spec_array <- abind(split(zoo_spec[,3:ncol(zoo_spec)],zoo_spec$lake_name), along = 3)  


# calculate variability and synchrony measures across levels of organization
sync_var_all <- var.partition(spec_array)

## A) SINGLE SPECIES CONTRIBUTION
## all species minus 1
contr_each_spec_zoo <- list()

# transform data frame to format that can then be converted to 3d array
zoo_spec <- zoo %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)

spec_zoo <- unique(zoo_spec$species)

for (i in 1:length(spec_zoo)){
  
  print(i)
  zoo_spec_sub <- zoo_spec %>% 
    filter(species != spec_zoo[i])
  
  
  # convert data frame to 3d array
  spec_array <- abind(split(zoo_spec_sub[,3:ncol(zoo_spec_sub)],zoo_spec_sub$lake_name), along = 3)  
  
  # calculate variability and synchrony measures across levels of organization
  
  sync_var_spec_sub <- var.partition(spec_array)
  
  contr_each_spec_zoo[[i]] <- log(sync_var_spec_sub[8]/sync_var_all[8])     # 8: phi_S2C_R 
  
}

contr_each_spec_zoo

## B) SINGLE SITE CONTRIBUTION

## all lakes minus 1
contr_each_lake_zoo <- list()
lakes <- unique(zoo$lake_name)

# transform data frame to format that can then be convertet to 3d array
zoo_lake <- zoo %>% 
  gather(., key = species, value = biomass, ends_with(".bio")) %>% 
  spread(.,key = sample_ID, value = biomass)


for (i in 1:length(lakes)){
  
  print(i)
  zoo_lake_sub <- zoo_lake %>% 
    filter(lake_name != lakes[i])
  zoo_lake_sub$lake_name <- zoo_lake_sub$lake_name[,drop = T]
  
  # convert data frame to 3d array
  lake_array <- abind(split(zoo_lake_sub[,3:ncol(zoo_lake_sub)],zoo_lake_sub$lake_name), along = 3)  
  
  # calculate variability and synchrony measures across levels of organization
  
  sync_var_lake_sub <- var.partition(lake_array)
  
  contr_each_lake_zoo[[i]] <- log(sync_var_lake_sub[6]/sync_var_all[6])     # 6: phi_C_L2R 
  
}

contr_each_lake_zoo

## MZB SPECIES -----------------

source(paste(in_path, "Armonies data/Analyses/stability metrics/function_var_partition.R", sep = ""))

mzb <- read.csv(paste(in_path, "Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""))

# all species
# select relevant information from data file and transform data frame to format that can then be converted to 3d array
data <- mzb %>% 
  select(Posno, Year, species, biom) %>% 
  spread(., key = Year, value = biom)

Year <- unique(data$Year)

data_spec <- as.data.frame(data)
data_spec[is.na(data_spec) == T] <- 0
data_spec$Posno <- data_spec$Posno[,drop=T]

# convert data frame to 3d array
MZB_array <- abind(split(data_spec[,3:ncol(data_spec)],data_spec$Posno), along = 3)

# calculate variability and synchrony measures across levels of organization
sync_var_all <- var.partition(MZB_array)


## A) SINGLE SPECIES CONTRIBUTION
# all species minus 1

contr_each_spec_mzb <- list() # create list for output
spec_mzb <- unique(data$species)

# select relevant information from data file and transform data frame to format that can then be converted to 3d array

for (i in 1:length(spec_mzb)){ 
  
  print(i)
  
  data <- mzb %>% 
    select(Posno, Year, species, biom) %>% 
    filter(species != spec_mzb[i]) %>% 
    spread(., key = Year, value = biom)
  
  
  data_spec <- as.data.frame(data)
  data_spec[is.na(data_spec) == T] <- 0
  data_spec$Posno <- data_spec$Posno[,drop=T]
  
  # convert data frame to 3d array
  MZB_array <- abind(split(data_spec[,3:ncol(data_spec)],data_spec$Posno), along = 3)
  
  # calculate variability and synchrony measures across levels of organization
  sync_var_spec_sub <- var.partition(MZB_array)
  
  contr_each_spec_mzb[[i]] <- log(sync_var_spec_sub[8]/sync_var_all[8])     # 8: phi_S2C_R 
  
}

contr_each_spec_mzb
#spec[contr_each_spec_mzb == 0]


## B) SINGLE SITE CONTRIBUTION
# all stations minus 1

contr_each_station_mzb <- list() # create list for output
stations <- unique(mzb$Posno)

# select relevant information from data file and transform data frame to format that can then be converted to 3d array

for (i in 1:length(stations)){ 
  
  print(i)
  
  data <- mzb %>% 
    select(Posno, Year, species, biom) %>% 
    filter(Posno != stations[i]) %>% 
    spread(., key = Year, value = biom)
  
  
  data_stat <- as.data.frame(data)
  data_stat[is.na(data_stat) == T] <- 0
  data_stat$Posno <- data_stat$Posno[,drop=T]
  
  # convert data frame to 3d array
  MZB_array <- abind(split(data_stat[,3:ncol(data_stat)],data_stat$Posno), along = 3)
  
  # calculate variability and synchrony measures across levels of organization
  sync_var_station_sub <- var.partition(MZB_array)
  
  contr_each_station_mzb[[i]] <- log(sync_var_station_sub[6]/sync_var_all[6])     # 6: phi_C_L2R 
  
}

contr_each_station_mzb



#### Combine and save results  for all species groups ----

contr_spec_phyto <- unlist(contr_each_spec_phyto)      #  87
contr_lake_phyto <- unlist(contr_each_lake_phyto)      #  49
contr_spec_zoo <- unlist(contr_each_spec_zoo)          #  21
contr_lake_zoo <- unlist(contr_each_lake_zoo)          #  49
contr_spec_mzb <- unlist(contr_each_spec_mzb)          #  195
contr_station_mzb <- unlist(contr_each_station_mzb)    #  72


contr_spec <- data.frame(spec_group = c(rep("phyto", length(contr_spec_phyto)), rep("zoo",length(contr_spec_zoo)), rep("mzb", length(contr_spec_mzb))), contr = c(contr_spec_phyto, contr_spec_zoo, contr_spec_mzb), species = c(as.character(spec_phyto), as.character(spec_zoo), as.character(spec_mzb)))
contr_space <- data.frame(spec_group = c(rep("phyto", length(contr_lake_phyto)), rep("zoo", length(contr_lake_zoo)), rep("mzb", length(contr_station_mzb))), contr = c(contr_lake_phyto, contr_lake_zoo, contr_station_mzb), patch = c(as.character(lakes),as.character(lakes),as.character(stations)))

write.csv(contr_spec, paste(out_path, "Iowa lakes/Analyses/stability metrics/contr_spec.csv", sep = ""))
write.csv(contr_space, paste(out_path, "Iowa lakes/Analyses/stability metrics/contr_patch.csv", sep = ""))



