## Calculate site and species characteristics of the three data subsets Iowa lakes phytoplankton, zooplankton, and North Sea macrozoobenthos
## in order to relate those to the contribution size of single sites and species to overall stabilizing effects with increasing area aggregation
## Doro Hodapp, May 2020


rm(list = ls())

## Load packages
library(dplyr)
library(tidyr)
library(betapart)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(fuzzyjoin)
library(vegan)

# diversity function from vegan package adjusted for mZB analysis
diversity_dh <- function (x, index = "shannon", MARGIN = 1, base = exp(1)) 
{
  x <- drop(as.matrix(x))
  if (!is.numeric(x)) 
    stop("input data must be numeric")
  if (any(x < 0, na.rm = TRUE)) 
    stop("input data must be non-negative")
  INDICES <- c("shannon", "simpson", "invsimpson")
  index <- match.arg(index, INDICES)
  if (length(dim(x)) > 1) {
    total <- apply(x, MARGIN, sum, na.rm = T)
    x <- sweep(x, MARGIN, total, "/")
  }
  else {
    x <- x/(total <- sum(x, na.rm = T))
  }
  if (index == "shannon") 
    x <- -x * log(x, base)
  else x <- x * x
  if (length(dim(x)) > 1) 
    H <- apply(x, MARGIN, sum, na.rm = TRUE)
  else H <- sum(x, na.rm = TRUE)
  if (index == "simpson") 
    H <- 1 - H
  else if (index == "invsimpson") 
    H <- 1/H
  if (any(NAS <- is.na(total))) 
    H[NAS] <- NA
  H
}


# specnumber function from vegan package adjusted for MZB analysis
specnumber_dh <- function (x, groups, MARGIN = 1) 
{
  if (!missing(groups)) {
    if (length(groups) == 1) 
      groups <- rep(groups, nrow(x))
    x <- aggregate(x, list(groups), max, na.rm = T)
    rownames(x) <- x[, 1]
    x <- x[, -1]
  }
  if (length(dim(x)) > 1) 
    apply(x > 0, MARGIN, sum, na.rm = T)
  else sum(x > 0, na.rm = T)
}

## Choose species group to run analysis for..
type = "phyto"  #("phyto", "zoo", "mzb")


#### Calculate species and patch characteristics ----
if(type == "phyto"){
  data <- read.csv(paste(in_path,"Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
  # Convert to long format
  data_long <- data %>%
    gather(., key = species, value = biom, ends_with(".bio"))
  
  patches <- unique(data_long$lake_name)
  names(data_long) <- c("patch", "sample", "species", "biom")
  names(data)[1:2] <- c("patch","sample")
  
  # Load information on single species contribution (IAR)
  contr_spec_IAR_phyto <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/IAR_corr_species_phyto.csv", sep = ""))
  # Load information on single patch contribution (IAR)
  contr_patch_IAR_phyto <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/IAR_corr_patch_phyto.csv", sep = ""))
  names(contr_patch_IAR_phyto)[2] <- "patch" 
  
  #### Patch characteristics ----
  
  ## Biomass over time --
  stats_biom <- data_long %>%
    group_by(patch) %>% 
    dplyr::summarise(tot.biom = sum(biom, na.rm = T), mean.biom = mean(biom, na.rm = T), var.biom = var(biom, na.rm = T))
  
  
  ## Community uniqueness -- calculated as average compositional dissimilarity from all other patches
  
  patch_spec_biom <- data_long %>% 
    group_by(patch, species) %>% 
    dplyr::summarise(mean.biom = sum(biom, na.rm = T)) %>%       # calculate mean biomass per species and patch
    spread(., key = species, value = mean.biom)           # convert to wide format for further analyses
  
  dissim_pair <- beta.pair.abund(patch_spec_biom[,c(2:ncol(patch_spec_biom))], index.family = "bray")   # calculate pairwise dissimilarity indices 
  
  dissim.m.bray <- as.matrix(dissim_pair$beta.bray)  # overall dissimilarity in abundances
  dissim.m.bal <- as.matrix(dissim_pair$beta.bray.bal)    # balanced variation in abundances
  dissim.m.gra <- as.matrix(dissim_pair$beta.bray.gra)    # abundance gradients
  
  # Calculate average dissimilarity to all other patches
  
  dist.bray <- colSums(dissim.m.bray)/(nrow(dissim.m.bray)-1)
  dist.bal <- colSums(dissim.m.bal)/(nrow(dissim.m.bal)-1)
  dist.gra <- colSums(dissim.m.gra)/(nrow(dissim.m.gra)-1)
  
  dissim <- data.frame(patch = patch_spec_biom$patch, dissim.bal = dist.bal, dissim.gra = dist.gra, dissim.bray = dist.bray)
  
  ## Diversity (richness, effective species number (ESN) derived from Hulbert's probability of interspecifc encounter (PIE)) --
  
  data_m <- as.matrix(data[,c(3:ncol(data))])
  
  p_i <- data_m/rowSums(data_m, na.rm = T)
  
  p_i_2 <- p_i^2 
  
  ESN <- apply(p_i_2, 1, function(x) 1/sum(x, na.rm = T))
  sample.esn <- cbind(data[,c(1,2)],ESN)
  
  patch_esn <- sample.esn %>% 
    group_by(patch) %>% 
    dplyr::summarise(esn.mean = mean(ESN, na.rm = T), esn.var = var(ESN, na.rm = T))
  
  
  patch_info_phyto <- cbind(stats_biom, dissim[,c(2:ncol(dissim))], patch_esn[,c(2:ncol(patch_esn))])
  
  #### Species characteristics ----
  
  ## Species overall biomass --
  
  spec_tot_biom <- data_long %>%
    group_by(species) %>% 
    dplyr::summarise(spec.tot.biom = sum(biom, na.rm = T))
  
  hist(log(spec_tot_biom$spec.tot.biom))
  
  
  ## Species mean temporal variability across patches, mean temporal biomass across patches, and spatial variability of temporal mean biomass --
  
  spec_biom_stats_phyto <- data_long %>% 
    group_by(patch, species) %>% 
    dplyr::summarise(spec.biom.var = var(biom, na.rm = T), spec.biom.mean = mean(biom, na.rm = T)) %>% 
    group_by(species) %>% 
    dplyr::summarise(mean.temp.var = mean(spec.biom.var, ne.rm = T), mean.biom.patch = mean(spec.biom.mean,na.rm = T), biom.spatial.var = var(spec.biom.mean, na.rm = T))
  
  
  #### Overall sampling area characteristics ----
  # Average and variability of evenness of local communities
  
  # function to calculate evenness per sample
  func.even <- function(x){
    H <- diversity(x)
    even <- H/log(specnumber(x))
  }
  
  even_sample <- apply(data[,3:ncol(data)], 1, func.even)
  even_sample_df <- cbind(data[,1:2], even_sample)
  
  mean_var_even <- even_sample_df %>% 
    summarise(even_mean = mean(even_sample, na.rm = T), even_var = var(even_sample, na.rm = T))
  
  # Average and variability of effective species number across sample
  
  mean_var_esn <- patch_esn %>% 
    dplyr::summarise(esn_mean = mean(esn.mean, na.rm = T), esn_var = var(esn.mean, na.rm = T))
  
  # Average and variability of pairwise dissimilarity across samples
  
  mean_var_dissim <- data.frame(dissim_mean = NA, dissim_var = NA)
  mean_var_dissim$dissim_mean <- mean(dissim.m.bray, na.rm = T)
  mean_var_dissim$dissim_var <- var(c(dissim.m.bray), na.rm = T)
  
  
}else if(type == "zoo"){
  data <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
  # Convert to long format
  data_long <- data %>%
    gather(., key = species, value = biom, ends_with(".bio"))
  
  patches <- unique(data_long$lake_name)
  names(data_long) <- c("patch", "sample", "species", "biom")
  names(data)[1:2] <- c("patch","sample")
  
  # Load information on single species contribution (IAR)
  contr_spec_IAR_zoo <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/IAR_corr_species_zoo.csv", sep = ""))
  # Load information on single patch contribution (IAR)
  contr_patch_IAR_zoo <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/IAR_corr_patch_zoo.csv", sep = ""))
  names(contr_patch_IAR_zoo)[2] <- "patch" 
  
  #### Patch characteristics ----
  
  ## Total patch biomass over time --
  stats_biom <- data_long %>%
    group_by(patch) %>% 
    dplyr::summarise(tot.biom = sum(biom, na.rm = T), mean.biom = mean(biom, na.rm = T), var.biom = var(biom, na.rm = T))
  
  
  ## Community uniqueness -- calculated as average compositional dissimilarity from all other patches
  
  patch_spec_biom <- data_long %>% 
    group_by(patch, species) %>% 
    dplyr::summarise(mean.biom = sum(biom, na.rm = T)) %>%       # calculate mean biomass per species and patch
    spread(., key = species, value = mean.biom)           # convert to wide format for further analyses
  
  dissim_pair <- beta.pair.abund(patch_spec_biom[,c(2:ncol(patch_spec_biom))], index.family = "bray")   # calculate pairwise dissimilarity indices 
  
  dissim.m.bray <- as.matrix(dissim_pair$beta.bray)  # overall dissimilarity in abundances
  dissim.m.bal <- as.matrix(dissim_pair$beta.bray.bal)    # balanced variation in abundances
  dissim.m.gra <- as.matrix(dissim_pair$beta.bray.gra)    # abundance gradients
  
  # Calculate average dissimilarity to all other patches
  
  dist.bray <- colSums(dissim.m.bray)/(nrow(dissim.m.bray)-1)
  dist.bal <- colSums(dissim.m.bal)/(nrow(dissim.m.bal)-1)
  dist.gra <- colSums(dissim.m.gra)/(nrow(dissim.m.gra)-1)
  
  dissim <- data.frame(patch = patch_spec_biom$patch, dissim.bal = dist.bal, dissim.gra = dist.gra, dissim.bray = dist.bray)
  
  ## Diversity (richness, effective species number (ESN) derived from Hulbert's probability of interspecfic encounter (PIE)) --
  
  data_m <- as.matrix(data[,c(3:ncol(data))])
  
  p_i <- data_m/rowSums(data_m, na.rm = T)
  
  p_i_2 <- p_i^2 
  
  ESN <- apply(p_i_2, 1, function(x) 1/sum(x, na.rm = T))
  sample.esn <- cbind(data[,c(1,2)],ESN)
  
  patch_esn <- sample.esn %>% 
    group_by(patch) %>% 
    dplyr::summarise(esn.mean = mean(ESN, na.rm = T), esn.var = var(ESN, na.rm = T))
  
  mean_var_esn <- patch_esn %>% 
    dplyr::summarise(esn_mean = mean(esn.mean, na.rm = T), esn_var = var(esn.mean, na.rm = T))
  
  patch_info_zoo <- cbind(stats_biom, dissim[,c(2:ncol(dissim))], patch_esn[,c(2:ncol(patch_esn))])
  
  
  #### Species characteristics ----
  
  ## Species overall biomass --
  
  spec_tot_biom <- data_long %>%
    group_by(species) %>% 
    dplyr::summarise(spec.tot.biom = sum(biom, na.rm = T))
  
  hist(log(spec_tot_biom$spec.tot.biom))
  
  
  ## Species mean temporal variability across patches, mean temporal biomass across patches, and spatial variability of temporal mean biomass --
  
  spec_biom_stats_zoo <- data_long %>% 
    group_by(patch, species) %>% 
    dplyr::summarise(spec.biom.var = var(biom, na.rm = T), spec.biom.mean = mean(biom, na.rm = T)) %>% 
    group_by(species) %>% 
    dplyr::summarise(mean.temp.var = mean(spec.biom.var, ne.rm = T), mean.biom.patch = mean(spec.biom.mean,na.rm = T), biom.spatial.var = var(spec.biom.mean, na.rm = T))
  
  #### Overall sampling area characteristics ----
  # Average and variability of evenness of local communities
  
  # function to calculate evenness per sample
  func.even <- function(x){
    H <- diversity(x)
    even <- H/log(specnumber(x))
  }
  
  even_sample <- apply(data[,3:ncol(data)], 1, func.even)
  even_sample_df <- cbind(data[,1:2],even_sample)
  
  mean_var_even <- even_sample_df %>% 
    summarise(even_mean = mean(even_sample, na.rm = T), even_var = var(even_sample, na.rm = T))
  # even_mean = 0.5647366, even_var = 0.02384553
  
  # Average and variability of effective species number across sample
  
  mean_var_esn <- patch_esn %>% 
    dplyr::summarise(esn_mean = mean(esn.mean, na.rm = T), esn_var = var(esn.mean, na.rm = T))
  
  # Average and variability of pairwise dissimilarity across samples
  
  mean_var_dissim <- data.frame(dissim_mean = NA, dissim_var = NA)
  mean_var_dissim$dissim_mean <- mean(dissim.m.bray, na.rm = T)
  mean_var_dissim$dissim_var <- var(c(dissim.m.bray), na.rm = T)
  
}else{
  data_long <- read.csv(paste(in_path, "Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""))
  
  patches <- unique(data$Posno)
  
  spec_tot_biom <- data_long %>%
    group_by(species) %>% 
    dplyr::summarise(spec.tot.biom = sum(biom, na.rm = T)) %>% 
    filter(spec.tot.biom > 0)
  
  spec.mzb <- unique(spec_tot_biom$species) 
  
  data_long <- data_long %>% 
    select(Posno, Year, species, biom) %>% 
    filter(species %in% spec.mzb)
  
  
  names(data_long) <- c("patch", "sample", "species", "biom")
  
  data <- data_long %>% 
    spread(., key = species, value = biom) 
  
  names(data)[1:2] <- c("patch","sample")
  
  # Load information on single species contribution (IAR)
  contr_spec_IAR_mzb <- read.csv(paste(in_path, "Armonies data/Analyses/stability metrics/IAR_corr_species_mzb.csv", sep = ""))
  # Load information on single patch contribution (IAR)
  contr_patch_IAR_mzb <- read.csv(paste(in_path, "Armonies data/Analyses/stability metrics/IAR_corr_patch_mzb.csv", sep = ""))
  contr_patch_IAR_mzb$patch <- substr(contr_patch_IAR_mzb$station, 8, 10)
  contr_patch_IAR_mzb <- contr_patch_IAR_mzb[, c(1,4,3)]
  
  #### Patch characteristics ----
  
  ## Total patch biomass over time --
  stats_biom <- data_long %>%
    group_by(patch) %>% 
    dplyr::summarise(tot.biom = sum(biom, na.rm = T), mean.biom = mean(biom, na.rm = T), var.biom = var(biom, na.rm = T))
  
  
  ## Community uniqueness -- calculated as average compositional dissimilarity from all other patches
  
  patch_spec_biom <- data_long %>% 
    group_by(patch, species) %>% 
    dplyr::summarise(mean.biom = sum(biom, na.rm = T)) %>%       # calculate mean biomass per species and patch
    spread(., key = species, value = mean.biom)           # convert to wide format for further analyses
  
  dissim_pair <- beta.pair.abund(patch_spec_biom[,c(2:ncol(patch_spec_biom))], index.family = "bray")   # calculate pairwise dissimilarity indices 
  
  dissim.m.bray <- as.matrix(dissim_pair$beta.bray)  # overall dissimilarity in abundances
  dissim.m.bal <- as.matrix(dissim_pair$beta.bray.bal)    # balanced variation in abundances
  dissim.m.gra <- as.matrix(dissim_pair$beta.bray.gra)    # abundance gradients
  
  # Calculate average dissimilarity to all other patches
  
  dist.bray <- colSums(dissim.m.bray)/(nrow(dissim.m.bray)-1)
  dist.bal <- colSums(dissim.m.bal)/(nrow(dissim.m.bal)-1)
  dist.gra <- colSums(dissim.m.gra)/(nrow(dissim.m.gra)-1)
  
  dissim <- data.frame(patch = patch_spec_biom$patch, dissim.bal = dist.bal, dissim.gra = dist.gra, dissim.bray = dist.bray)
  
  ## Diversity (richness, effective species number (ESN) derived from Hulbert's probability of interspecfic encounter (PIE)) --
  
  data_m <- as.matrix(data[,c(3:ncol(data))])
  
  p_i <- data_m/rowSums(data_m, na.rm = T)
  
  p_i_2 <- p_i^2 
  
  ESN <- apply(p_i_2, 1, function(x) 1/sum(x, na.rm = T))
  sample.esn <- cbind(data[,c(1,2)],ESN)
  
  patch_esn <- sample.esn %>% 
    group_by(patch) %>% 
    dplyr::summarise(esn.mean = mean(ESN, na.rm = T), esn.var = var(ESN, na.rm = T))
  
  mean_var_esn <- patch_esn %>% 
    dplyr::summarise(esn_mean = mean(esn.mean, na.rm = T), esn_var = var(esn.mean, na.rm = T))
  
  patch_info_mzb <- cbind(stats_biom, dissim[,c(2:ncol(dissim))], patch_esn[,c(2:ncol(patch_esn))])
  
  #### Species characteristics ----
  
  ## Species overall biomass --
  
  spec_tot_biom <- data_long %>%
    group_by(species) %>% 
    dplyr::summarise(spec.tot.biom = sum(biom, na.rm = T)) %>% 
    filter(spec.tot.biom > 0)
  
  
  hist(log(spec_tot_biom$spec.tot.biom))
  
  
  ## Species mean temporal variability across patches, mean temporal biomass across patches, and spatial variability of temporal mean biomass --
  
  spec_biom_stats_mzb <- data_long %>% 
    group_by(patch, species) %>% 
    dplyr::summarise(spec.biom.var = var(biom, na.rm = T), spec.biom.mean = mean(biom, na.rm = T)) %>% 
    group_by(species) %>% 
    dplyr::summarise(mean.temp.var = mean(spec.biom.var, ne.rm = T), mean.biom.patch = mean(spec.biom.mean,na.rm = T), biom.spatial.var = var(spec.biom.mean, na.rm = T))
  
  #### Overall sampling area characteristics ----
  # Average and variability of evenness of local communities
  
  # function to calculate evenness per sample
  func.even <- function(x){
    H <- diversity_dh(x)
    even <- H/log(specnumber_dh(x))
  }
  
  even_sample <- apply(data[,3:ncol(data)], 1, func.even)
  even_sample_df <- cbind(data[,1:2],even_sample)
  
  mean_var_even <- even_sample_df %>% 
    summarise(even_mean = mean(even_sample, na.rm = T), even_var = var(even_sample, na.rm = T))
  
  # Average and variability of effective species number across sample
  
  mean_var_esn <- patch_esn %>% 
    dplyr::summarise(esn_mean = mean(esn.mean, na.rm = T), esn_var = var(esn.mean, na.rm = T))
  
  # Average and variability of pairwise dissimilarity across samples
  
  mean_var_dissim <- data.frame(dissim_mean = NA, dissim_var = NA)
  mean_var_dissim$dissim_mean <- mean(dissim.m.bray, na.rm = T)
  mean_var_dissim$dissim_var <- var(c(dissim.m.bray), na.rm = T)
  
}

#### Combine data sets and save ----
# Combine species info from all species groups
spec_biom_stats_phyto$spec_group <- "phyto"
spec_biom_stats_zoo$spec_group <- "zoo"
spec_biom_stats_mzb$spec_group <- "mzb"

spec_biom_stats_all <- rbind(spec_biom_stats_phyto, spec_biom_stats_zoo, spec_biom_stats_mzb)

# Save species info
write.csv(spec_biom_stats_all, paste(out_path, "Iowa lakes/Analyses/stability metrics/spec_info_all.csv", sep =""), row.names =F)


# Combine patch info from all species groups
patch_info_phyto$spec_group <- "phyto"
patch_info_zoo$spec_group <- "zoo"
patch_info_mzb$spec_group <- "mzb"

patch_info_all <- rbind(patch_info_phyto, patch_info_zoo, patch_info_mzb)

# Save patch info
write.csv(patch_info_all, paste(out_path, "Iowa lakes/Analyses/stability metrics/patch_info_all.csv", sep =""), row.names =F)



#### Distribution of biomass across all sites ####

phyto <- read.csv(paste(in_path,"Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
zoo <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
mzb_long <- read.csv(paste(in_path, "Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""))

# PHYTO
phyto_long <- phyto %>%
  gather(., key = species, value = biom, ends_with(".bio"))

patches <- unique(phyto_long$lake_name)
names(phyto_long) <- c("patch", "sample", "species", "biom")

phyto_totbiom_spatial <- phyto_long %>% 
  group_by(patch) %>% 
  dplyr::summarise(totbiom = sum(biom))


# ZOO
zoo_long <- zoo %>%
  gather(., key = species, value = biom, ends_with(".bio"))

patches <- unique(zoo_long$lake_name)
names(zoo_long) <- c("patch", "sample", "species", "biom")

zoo_totbiom_spatial <- zoo_long %>% 
  group_by(patch) %>% 
  dplyr::summarise(totbiom = sum(biom))


# MZ
mzb_totbiom_spatial <- mzb_long %>% 
  group_by(Posno) %>% 
  dplyr::summarise(totbiom = sum(biom))


par(mfrow = c(1,3))
hist(phyto_totbiom_spatial$totbiom, breaks = 15)
hist(zoo_totbiom_spatial$totbiom, breaks = 15)
hist(mzb_totbiom_spatial$totbiom, breaks = 15)



#### Calculate synchrony between single species or sites and the combined biomass of remaining community ----
#### Additional analysis suggested by reviewer, synchrony measure taken from Gross et al. 2014 - AmNat (https://doi.org/10.1086/673915) ----

## PHYTO species
data <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
# Convert to long format
data_long <- data %>%
  gather(., key = species, value = biom, ends_with(".bio"))

nspec <- unique(data_long$species)

corr_phyto_df <- data.frame(species = nspec, corr = c(rep(NA,length(nspec))))

# Calculate correlation of each single species with the summed biomass of the rest of the community
for (i in 1:length(nspec)){
  group_biom <- data_long %>% 
    filter(species != nspec[i]) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(group_biom = sum(biom, na.rm = T))
  
  spec_biom <- data_long %>% 
    filter(species == nspec[i]) %>% 
    group_by(sample_ID) %>% 
    dplyr::summarise(spec_biom = sum(biom, na.rm = T))
  
  corr_biom <- rcorr(group_biom$group_biom, spec_biom$spec_biom)
  corr_phyto_df$corr [i] <- corr_biom$r[1,2]
}

## PHYTO sites

nsites <- unique(data_long$lake_name)

corr_phyto_sites_df <- data.frame(sites = nsites, corr = c(rep(NA,length(nsites))))

# Calculate correlation of each single species with the summed biomass of the rest of the community
for (i in 1:length(nsites)){
  group_biom <- data_long %>% 
    filter(lake_name != nsites[i]) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(group_biom = sum(biom, na.rm = T))
  
  site_biom <- data_long %>% 
    filter(lake_name == nsites[i]) %>% 
    group_by(sample_ID) %>% 
    dplyr::summarise(site_biom = sum(biom, na.rm = T))
  
  corr_biom <- rcorr(group_biom$group_biom, site_biom$site_biom)
  corr_phyto_sites_df$corr [i] <- corr_biom$r[1,2]
}


## ZOO species 
data <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
# Convert to long format
data_long <- data %>%
  gather(., key = species, value = biom, ends_with(".bio"))

nspec <- unique(data_long$species)

corr_zoo_df <- data.frame(species = nspec, corr = c(rep(NA,length(nspec))))

# Calculate correlation of each single species with the summed biomass of the rest of the community
for (i in 1:length(nspec)){
  group_biom <- data_long %>% 
    filter(species != nspec[i]) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(group_biom = sum(biom, na.rm = T))
  
  spec_biom <- data_long %>% 
    filter(species == nspec[i]) %>% 
    group_by(sample_ID) %>% 
    dplyr::summarise(spec_biom = sum(biom, na.rm = T))
  
  corr_biom <- rcorr(group_biom$group_biom, spec_biom$spec_biom)
  corr_zoo_df$corr [i] <- corr_biom$r[1,2]
} 

## ZOO sites

nsites <- unique(data_long$lake_name)

corr_zoo_sites_df <- data.frame(sites = nsites, corr = c(rep(NA,length(nsites))))

# Calculate correlation of each single species with the summed biomass of the rest of the community
for (i in 1:length(nsites)){
  group_biom <- data_long %>% 
    filter(lake_name != nsites[i]) %>%
    group_by(sample_ID) %>% 
    dplyr::summarise(group_biom = sum(biom, na.rm = T))
  
  site_biom <- data_long %>% 
    filter(lake_name == nsites[i]) %>% 
    group_by(sample_ID) %>% 
    dplyr::summarise(site_biom = sum(biom, na.rm = T))
  
  corr_biom <- rcorr(group_biom$group_biom, site_biom$site_biom)
  corr_zoo_sites_df$corr [i] <- corr_biom$r[1,2]
}


## MZB species
data_long <- read.csv(paste(in_path, "Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""))

nspec <- unique(data_long$species)

corr_mzb_df <- data.frame(species = nspec, corr = c(rep(NA,length(nspec))))

# Calculate correlation of each single species with the summed biomass of the rest of the community
for (i in 1:length(nspec)){
  group_biom <- data_long %>% 
    filter(species != nspec[i]) %>%
    group_by(Year) %>% 
    dplyr::summarise(group_biom = sum(biom, na.rm = T))
  
  spec_biom <- data_long %>% 
    filter(species == nspec[i]) %>% 
    group_by(Year) %>% 
    dplyr::summarise(spec_biom = sum(biom, na.rm = T))
  
  corr_biom <- rcorr(group_biom$group_biom, spec_biom$spec_biom)
  corr_mzb_df$corr [i] <- corr_biom$r[1,2]
}

## MZB sites

nsites <- unique(data_long$Posno)

corr_mzb_sites_df <- data.frame(sites = nsites, corr = c(rep(NA,length(nsites))))

# Calculate correlation of each single species with the summed biomass of the rest of the community
for (i in 1:length(nsites)){
  group_biom <- data_long %>% 
    filter(Posno != nsites[i]) %>%
    group_by(Year) %>% 
    dplyr::summarise(group_biom = sum(biom, na.rm = T))
  
  site_biom <- data_long %>% 
    filter(Posno == nsites[i]) %>% 
    group_by(Year) %>% 
    dplyr::summarise(site_biom = sum(biom, na.rm = T))
  
  corr_biom <- rcorr(group_biom$group_biom, site_biom$site_biom)
  corr_mzb_sites_df$corr [i] <- corr_biom$r[1,2]
}

## Combine species groups and save

corr_info_spec <- rbind(corr_phyto_df, corr_zoo_df, corr_mzb_df)
write.csv(corr_info_spec, paste(out_path, "Iowa lakes/Analyses/stability metrics/spec_info_corr.csv", sep =""), row.names = F)

corr_info_sites <- rbind(corr_phyto_sites_df, corr_zoo_sites_df, corr_mzb_sites_df)
write.csv(corr_info_sites, paste(out_path, "Iowa lakes/Analyses/stability metrics/patch_info_corr.csv", sep =""), row.names = F)



