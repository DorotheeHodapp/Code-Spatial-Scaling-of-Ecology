#### Create figures
#### Doro Hodapp  

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())
library(sf)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(devtools)
# install.packages("rnaturalearthhires",
#                  repos = "http://packages.ropensci.org",
#                  type = "source")
library(rnaturalearthhires)
library(maps)
library(gridExtra)


#### Figure 1: Sampling sites ----

#### Iowa lakes ----

data_geo <- read.csv(paste(in_path,"Iowa lakes/Analyses/lake_coords.csv", sep = ""), header = T)
iowa_dat <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""), header = T)
iowa_dat_zoo <- read.csv(paste(in_path, "Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""), header = T)
lakes <- unique(iowa_dat$lake_name)

data_geo_sub <- data_geo %>% 
  filter(lake_name %in% lakes)

# Convert data frame to sf object
lake_coords <- st_as_sf(data_geo_sub, coords = c("longitude","latitude"))
lake_coords <- st_set_crs(lake_coords, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


# Download data for world map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# Plot world map 
png(paste(out_path, "Iowa lakes/Analyses/Graphs/worldmap_countries.png", sep = ""), width = 2000, height = 1400, res = 300)
ggplot(data = world) +
  geom_sf(color = "dark grey", fill = "light grey") +
  coord_sf(xlim = c(-120, 30), ylim = c(0, 70), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")
dev.off()

us_states <- ne_states(country = 'united states of america', geounit = NULL, iso_a2 = NULL, spdf = NULL,
          returnclass = "sf")

states <- ne_states(country = NULL, geounit = NULL, iso_a2 = NULL, spdf = NULL,
                       returnclass = "sf")

# Plot Iowa map with sampling locations
png(paste(out_path, "Iowa lakes/Analyses/Graphs/sites_iowa.png", sep = ""), width = 2000, height = 1400, res = 300)
ggplot(data = us_states) +
  geom_sf(fill = "white") +
  geom_sf(data = lake_coords, size = 0.9) +
  coord_sf(xlim = c(-97, -90), ylim = c(40, 44), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("Iowa (US)", subtitle = "Phyto- and Zooplankton")
dev.off()


# Plot Iowa lakes with total biomass values (Phytoplankton)
# First calculate total biomass per lake, then merge with geo information
phyto_long <- iowa_dat %>%
  gather(., key = species, value = biom, ends_with(".bio"))

patches <- unique(phyto_long$lake_name)
names(phyto_long) <- c("patch", "sample", "species", "biom")

phyto_totbiom_spatial <- phyto_long %>% 
  filter(biom > 0) %>% 
  group_by(patch) %>% 
  dplyr::summarise(totbiom = sum(biom), meanbiom = mean(biom, na.rm = T), medbiom = median(biom, na.rm = T), varbiom = var(biom,na.rm = T))

lake_coords_biom <- merge(data_geo_sub, phyto_totbiom_spatial, by.x = c("lake_name"), by.y = c("patch"))

# Plot
png(paste(out_path, "Iowa lakes/Analyses/Graphs/sites_iowa_phyto_biom.png", sep = ""), width = 2000, height = 1400, res = 300)
ggplot(data = us_states) +
  geom_sf(fill = "white") +
  geom_sf(data = lake_coords, size = 0.9) +
  geom_point(data=lake_coords_biom,
             aes(x=longitude, y= latitude, color = varbiom)) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  coord_sf(xlim = c(-97, -90), ylim = c(40, 44), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("Iowa (US)", subtitle = "Total phytoplankton biomass per site")
dev.off()



# Plot Iowa lakes with total biomass values (Zooplankton)
# First calculate total biomass per lake, then merge with geo information
zoo_long <- iowa_dat_zoo %>%
  gather(., key = species, value = biom, ends_with(".bio"))

patches <- unique(zoo_long$lake_name)
names(zoo_long) <- c("patch", "sample", "species", "biom")

zoo_totbiom_spatial <- zoo_long %>% 
  filter(biom > 0) %>% 
  group_by(patch) %>% 
  dplyr::summarise(totbiom = sum(biom), meanbiom = mean(biom, na.rm = T), medbiom = median(biom, na.rm = T), varbiom = var(biom,na.rm = T))

lake_coords_biom <- merge(data_geo_sub, zoo_totbiom_spatial, by.x = c("lake_name"), by.y = c("patch"))

# Plot
png(paste(out_path, "Iowa lakes/Analyses/Graphs/sites_iowa_zoo_biom.png", sep = ""), width = 2000, height = 1400, res = 300)
ggplot(data = us_states) +
  geom_sf(fill = "white") +
  geom_sf(data = lake_coords, size = 0.9) +
  geom_point(data=lake_coords_biom,
             aes(x=longitude, y= latitude, color = varbiom)) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  coord_sf(xlim = c(-97, -90), ylim = c(40, 44), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("Iowa (US)", subtitle = "Total zooplankton biomass per site")
dev.off()



# Download coastline data
coastline <- ne_coastline(returnclass = "sf")

png(paste(out_path, "Iowa lakes/Analyses/Graphs/worldmap.png", sep = ""), width = 1200, height = 800, res = 300)
ggplot(data = coastline) +
  geom_sf(color = "dark grey", fill = "grey") +
  coord_sf(xlim = c(-120, 30), ylim = c(0, 70), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")
dev.off()



#### North Sea stations ----
# load data
dat_sub <- read.csv(paste(in_path, "Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv", sep = ""), header = T)
dat <- read.csv("/Users/dhodapp/Documents/Hodapp, Dorothee/HIFMB/Projects/Spatio-temporal Turnover/Armonies data/Data/Armonies Original Data.csv", header = T)
stations <- unique(dat_sub$Posno)

## GEOGRAPHIC DISTANCE
dat_geo <- dat %>% 
  mutate(Date_new = as.POSIXct(as.character(Date), format = "%m/%d/%Y %H:%M")) %>% 
  mutate(Year = format(Date_new,format = "%Y"), Month = format(Date_new,format = "%m"), Day = format(Date_new,format = "%d")) %>% 
  mutate(Season = ifelse(Month %in% c("01","02","12"), "W",ifelse(Month %in% c("03","04","05"), "Sp", ifelse(Month %in% c("06","07","08"),"Su",ifelse(Month %in% c("09","10","11"),"F","NA"))))) %>% 
  filter(Season=="F") %>% 
  select(Posno, StNBr, StELae) %>% 
  group_by(Posno) %>% 
  summarise_all(mean) %>%  # calculate mean longitude and latitude for each station
  filter(Posno %in% stations)



station_coords <- st_as_sf(dat_geo, coords = c("StELae","StNBr"))
station_coords <- st_set_crs(station_coords, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


# Plot North Sea site map with sampling locations
png(paste(out_path, "Iowa lakes/Analyses/Graphs/sites_north_sea.png", sep = ""), width = 2000, height = 1400, res = 300)
ggplot(data = states) +
  geom_sf(fill = "white") +
  geom_sf(data = station_coords, size = 0.9) +
  coord_sf(xlim = c(5, 10), ylim = c(53, 56), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("North Sea", subtitle = "Macrozoobenthos")
dev.off()

# Plot North Sea site map with total biomass per site
mzb_totbiom_spatial <- dat_sub %>% 
  filter(Posno %in% stations) %>% 
  filter(biom > 0) %>% 
  group_by(Posno) %>% 
  dplyr::summarise(totbiom = sum(biom), meanbiom = mean(biom, na.rm = T), medbiom = median(biom, na.rm = T), varbiom = var(biom,na.rm = T))

station_coords_biom <- merge(dat_geo, mzb_totbiom_spatial, by = c("Posno"))

# Plot
png(paste(out_path, "Iowa lakes/Analyses/Graphs/sites_NS_mzb_biom.png", sep = ""), width = 2000, height = 1400, res = 300)
ggplot(data = states) +
  geom_sf(fill = "white") +
  geom_sf(data = station_coords, size = 0.9) +
  geom_point(data=station_coords_biom,
             aes(x=StELae, y= StNBr, color = varbiom)) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  coord_sf(xlim = c(5, 10), ylim = c(53, 56), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("North Sea", subtitle = "Total macrozoobenthos biomass per site")
dev.off()

# Plot spatial biomass distribution per year 
station_coords_biom_annual <- merge(dat_geo, dat_sub, by = c("Posno"))


ggplot(data = states) +
  geom_sf(fill = "white") +
  geom_sf(data = station_coords, size = 0.9) +
  geom_point(data=station_coords_biom_annual[station_coords_biom_annual$Year == 2011,],
             aes(x=StELae, y= StNBr, color = totbiom)) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  coord_sf(xlim = c(5, 10), ylim = c(53, 56), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle("North Sea", subtitle = "Macrozoobenthos biomass")



#### Figure 2: IAR for the three species groups ----
iar_phyto_dist <- read.csv(".../invar_dist_average_phyto.csv")
iar_phyto_dist$spec_group <-"phyto"
iar_phyto_rand <- read.csv(".../invar_rand_average_phyto.csv")
iar_phyto_rand$spec_group <-"phyto"
iar_zoo_dist <- read.csv(".../invar_dist_average_zoo.csv")
iar_zoo_dist$spec_group <-"zoo"
iar_zoo_rand <- read.csv(".../invar_rand_average_zoo.csv")
iar_zoo_rand$spec_group <-"zoo"
iar_mzb_dist <- read.csv(".../invar_dist_average_mzb.csv")
iar_mzb_dist$spec_group <-"mzb"
iar_mzb_rand <- read.csv(".../invar_rand_average_mzb.csv")
iar_mzb_rand$spec_group <-"mzb"

iar_dist <- rbind(iar_phyto_dist, iar_zoo_dist, iar_mzb_dist)
iar_rand <- rbind(iar_phyto_rand, iar_zoo_rand, iar_mzb_rand)


p1 <- ggplot(iar_dist, aes(log(area), log(I_median))) +
  geom_line(aes(col = spec_group)) + 
  geom_ribbon(aes(ymin = log(I_lower), ymax = log(I_upper), fill = spec_group), alpha = 0.2) + 
#  ylim(-0.7, 1.8) +
  scale_colour_manual(values = c("#969696","#008080","#99CCFF"), name = "Species group") + 
  scale_fill_manual(values = c("#969696","#008080","#99CCFF"), name = "Species group") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  annotate(geom="text", x=0.2, y=1.7, label="(A)") +
  theme(legend.position = "bottom") + 
  xlab("log('Area' = number of sites)") + 
  ylab("log(Invariability)")

p2 <- ggplot(iar_rand, aes(log(area), log(I_median))) +
  geom_line(aes(col = spec_group)) + 
  geom_ribbon(aes(ymin = log(I_lower), ymax = log(I_upper), fill = spec_group), alpha = 0.2) + 
#  ylim(-0.7, 1.8) +
  scale_colour_manual(values = c("#969696","#008080","#99CCFF"), name = "Species group") + 
  scale_fill_manual(values = c("#969696","#008080","#99CCFF"), name = "Species group") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  annotate(geom="text", x=0.2, y=1.8, label="(B)") +
  xlab("log('Area' = number of sites)") + 
  ylab("log(Invariability)")

mylegend <- g_legend(p1)

png(paste(out_path,"Iowa lakes/Manuscript/Figures/Figure2.png", sep = ""), width = 2200, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"), 
                         p2 + theme(legend.position = "none"), nrow = 1), 
             mylegend, nrow=2, heights=c(10, 1))
dev.off()

tiff(paste(out_path,"Iowa lakes/Manuscript/Figures/Figure2.tiff", sep = ""), width = 2200, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"), 
                         p2 + theme(legend.position = "none"), nrow = 1), 
             mylegend, nrow=2, heights=c(10, 1))
dev.off()




#### Figure 3: created with powerpoint... ----
#### Figure 4: Plotting observed and null model distribution of synchrony measures between species and patches ----

#stab_sync_phyto <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/stab_sync_phyto.RDS", sep = ""))
stab_sync_test_phyto <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/stab_sync_test_phyto.RDS", sep = ""))

#stab_sync_zoo <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/stab_sync_zoo.RDS", sep = ""))
stab_sync_test_zoo <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/stab_sync_test_zoo.RDS", sep = ""))

#stab_sync_mzb <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/stab_sync_mzb.RDS", sep =""))
stab_sync_test_mzb <- readRDS(paste(in_path, "Iowa lakes/Analyses/stability metrics/Synchrony-stability-analysis/Synchrony-stability-analysis/stab_sync_test_mzb.RDS", sep =""))


null_distr <- stab_sync_test_phyto$null_model

null_distr_long_phyto <- null_distr %>% 
  filter(row_number()!=501) %>% 
  gather(., key = sync_measure, value = phi, 1:4) %>% 
  mutate(spec_group = "phyto")

obs_phyto <- null_distr %>%
  filter(row_number()== 501) %>% 
  gather(., key = sync_measure, value = phi, 1:4) %>% 
  mutate(spec_group = "phyto")

null_distr <- stab_sync_test_zoo$null_model

null_distr_long_zoo <- null_distr %>% 
  filter(row_number()!=501) %>% 
  gather(., key = sync_measure, value = phi, 1:4) %>% 
  mutate(spec_group = "zoo")

obs_zoo <- null_distr %>%
  filter(row_number()== 501) %>% 
  gather(., key = sync_measure, value = phi, 1:4) %>% 
  mutate(spec_group = "zoo")

null_distr <- stab_sync_test_mzb$null_model

null_distr_long_mzb <- null_distr %>% 
  filter(row_number()!=501) %>% 
  gather(., key = sync_measure, value = phi, 1:4) %>% 
  mutate(spec_group = "mzb")

obs_mzb <- null_distr %>%
  filter(row_number()== 501) %>% 
  gather(., key = sync_measure, value = phi, 1:4) %>% 
  mutate(spec_group = "mzb")

null_distr_long <- rbind(null_distr_long_phyto, null_distr_long_zoo, null_distr_long_mzb)
obs <- rbind(obs_phyto, obs_zoo, obs_mzb)


tiff(paste(out_path, "Iowa lakes/Manuscript/Figures/Figure4.tiff", sep = ""), width = 2800, height = 2000, res = 300)
ggplot(null_distr_long, aes(sync_measure, phi)) +
  geom_violin(aes(fill = spec_group)) + 
  geom_point(data = obs, aes(sync_measure, phi, colour = spec_group, fill = spec_group), position = position_dodge(width = 0.9), shape = 21, size = 4, colour = "black") +
  theme(axis.text.x = element_text(size = 14)) + 
  scale_fill_manual(name = "Species group", values=c("#440154FF","#21908CFF","#FDE725FF")) + 
  scale_color_manual(name = "Species group", values=c("#440154FF","#21908CFF","#FDE725FF")) +
  xlab("Synchrony measure") + 
  ylab("Synchrony") + 
  theme(axis.title = element_text(size = 18)) + 
  theme(legend.title = element_text(size=14))
dev.off()



ggplot(obs, aes(sync_measure, phi)) + 
  geom_point(aes(colour = spec_group), position = position_dodge(width = 1)) +
  



#### Figure 5: Correlation between single species or patch characteristics and contribution strength (neg/pos) to stability increase/decrease ----
## Load single species and patch contribution (stab-sync framework) --
# stab-sync
contr_spec_stab_sync_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/contr_spec.csv", sep  = ""))
contr_patch_stab_sync_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/contr_patch.csv", sep  = ""))

# species and patch info
spec_info_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/spec_info_all.csv", sep = ""))
patch_info_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/patch_info_all.csv", sep = "")) 
spec_info_corr <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/spec_info_corr.csv", sep = ""))
patch_info_corr <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/patch_info_corr.csv", sep = ""))


# Combine species, patch and stability contribution data

cc_spec <- merge(spec_info_all, contr_spec_stab_sync_all, by = "species")

cc_patch <- merge(patch_info_all, contr_patch_stab_sync_all, by = "patch")

cc_spec_corr <-  merge(spec_info_corr, contr_spec_stab_sync_all, by = "species")

cc_patch_corr <-  merge(patch_info_corr, contr_patch_stab_sync_all, by.x = "sites", by.y = "patch")

# helper function to extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


cc_spec <- cc_spec %>% 
  mutate(contr_abs = abs(contr))

cc_patch <- cc_patch %>% 
  mutate(contr_abs = abs(contr))

## Phytoplankton
cc_spec_phyto <- cc_spec %>% 
  filter(spec_group.x == "phyto")

cc_patch_phyto <- cc_patch %>% 
  filter(spec_group.x == "phyto")

# species characteristics 
cor.test(log(cc_spec_phyto$mean.biom.patch), log(cc_spec_phyto$contr_abs), method = "kendall")  # tau = 0.80; p-value < 2.2e-16  
cor.test(log(cc_spec_phyto$mean.temp.var), log(cc_spec_phyto$contr_abs), method = "kendall")  # tau = 0.81; p-value < 2.2e-16
cor.test(log(cc_spec_phyto$biom.spatial.var), log(cc_spec_phyto$contr_abs), method = "kendall")  # tau = 0.82; p-value < 2.2e-16

# patch characteristics
cor.test(log(cc_patch_phyto$mean.biom), log(cc_patch_phyto$contr_abs), method = "kendall")  # tau = 0.48; p-value = 2.061e-12
cor.test(log(cc_patch_phyto$var.biom), log(cc_patch_phyto$contr_abs), method = "kendall")  # tau = 0.50; p-value = 4.945e-13

# community characteristics
cor.test(log(cc_patch_phyto$dissim.bray), log(cc_patch_phyto$contr_abs), method = "kendall")  # tau = 0.04; p-value = 0.5271  
cor.test(cc_patch_phyto$esn.mean, cc_patch_phyto$contr_abs, method = "kendall")  # tau = -0.19; p-value = 0.0046
cor.test(cc_patch_phyto$esn.var, cc_patch_phyto$contr_abs, method = "kendall")  # tau = -0.09; p-value = 0.1671

## Zooplankton
cc_spec_zoo <- cc_spec %>% 
  filter(spec_group.x == "zoo")

cc_patch_zoo <- cc_patch %>% 
  filter(spec_group.x == "zoo")


cor.test(log(cc_spec_zoo$mean.biom.patch), log(cc_spec_zoo$contr_abs), method = "kendall")  # tau = 0.76; p-value = 5.824e-08  
cor.test(log(cc_spec_zoo$mean.temp.var), log(cc_spec_zoo$contr_abs), method = "kendall")  # tau = 0.8; p-value = 5.99e-09
cor.test(log(cc_spec_zoo$biom.spatial.var), log(cc_spec_zoo$contr_abs), method = "kendall")  # tau = 0.81; p-value = 3.209e-09

# patch characteristics
cor.test(log(cc_patch_zoo$mean.biom), log(cc_patch_zoo$contr_abs), method = "kendall")  # tau = 0.33; p-value = 1.527e-06
cor.test(log(cc_patch_zoo$var.biom), log(cc_patch_zoo$contr_abs), method = "kendall")  # tau = 0.32; p-value = 2.975e-06

# community characteristics
cor.test(log(cc_patch_zoo$dissim.bray), log(cc_patch_zoo$contr_abs), method = "kendall")  # tau = 0.27; p-value = 7.88e-05  
cor.test(cc_patch_zoo$esn.mean, cc_patch_zoo$contr_abs, method = "kendall")  # tau = 0.13; p-value = 0.0562
cor.test(cc_patch_zoo$esn.var, cc_patch_zoo$contr_abs, method = "kendall")  # tau = 0.04; p-value = 0.5932


## Macrozoobenthos
cc_spec_mzb <- cc_spec %>% 
  filter(spec_group.x == "mzb")

cc_patch_mzb <- cc_patch %>% 
  filter(spec_group.x == "mzb")


cor.test(log(cc_spec_mzb$mean.biom.patch), log(cc_spec_mzb$contr_abs), method = "kendall")  # tau = 0.82; p-value < 2.2e-16  
cor.test(log(cc_spec_mzb$mean.temp.var), log(cc_spec_mzb$contr_abs), method = "kendall")  # tau = 0.84; p-value < 2.2e-16
cor.test(log(cc_spec_mzb$biom.spatial.var), log(cc_spec_mzb$contr_abs), method = "kendall")  # tau = 0.83; p-value < 2.2e-16

# patch characteristics
cor.test(log(cc_patch_mzb$mean.biom), log(cc_patch_mzb$contr_abs), method = "kendall")  # tau = 0.37; p-value = 4.453e-06
cor.test(log(cc_patch_mzb$var.biom), log(cc_patch_mzb$contr_abs), method = "kendall")  # tau = 0.42; p-value = 1.442e-07

# community characteristics
cor.test(log(cc_patch_mzb$dissim.bray), log(cc_patch_mzb$contr_abs), method = "kendall")  # tau = 0.22; p-value = 0.0063  
cor.test(cc_patch_mzb$esn.mean, cc_patch_mzb$contr_abs, method = "kendall")  # tau = -0.11; p-value = 0.1395
cor.test(cc_patch_mzb$esn.var, cc_patch_mzb$contr_abs, method = "kendall")  # tau = -0.01; p-value = 0.8611

# helper function to extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


## Contribution of species and patches to stab-sync framework --  
# Species
p1 <- ggplot(cc_spec, aes(log(mean.biom.patch), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  geom_smooth(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  theme_bw() +
  ylab("Single species contribution\n to meta-community stability") + 
  xlab("log(Mean total species biomass)") + 
  annotate(geom="text", x=-13, y=0.6, label="(A)") +
  theme(legend.position = "bottom", legend.direction = "horizontal")

p2 <- ggplot(cc_spec, aes(log(mean.temp.var), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  geom_smooth(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  ylab("") +
  xlab("log(Temporal species biomass variability)") +
  annotate(geom="text", x=-22, y=0.19, label="(B)") +
  theme_bw()

p3 <- ggplot(cc_spec, aes(log(biom.spatial.var), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  geom_smooth(aes(color = spec_group.x)) +  
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  ylab("") + 
  xlab("log(Spatial species biomass variability)") +
  annotate(geom="text", x=-23.5, y=0.19, label="(C)") +
  theme_bw()

mylegend <- g_legend(p1)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/spec_biom_corr_contr_abs.png", sep = ""), width = 3300, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"), 
                         p2 + theme(legend.position = "none"), 
                         p3 + theme(legend.position = "none"), nrow = 1))#, 
           #  mylegend, nrow=2, heights=c(10, 1))
dev.off()


# Patch
p4 <-ggplot(cc_patch, aes(log(tot.biom), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  geom_smooth(aes(color = spec_group.x)) + 
  #  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() + 
  ylab("Single site contribution\n to meta-community stability") +
  xlab("log(Total biomass)") + 
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom="text", x=2.5, y=0.145, label="(D)") +
  theme(legend.position = "bottom", legend.direction = "horizontal")

p5 <- ggplot(cc_patch, aes(log(mean.biom), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  geom_smooth(aes(color = spec_group.x)) + 
  xlab("log(Mean biomass)") + 
  ylab("Single site contribution\n to meta-community stability") +
  #  geom_hline(yintercept = 0, linetype="dashed") +
  annotate(geom="text", x=-4.7, y=0.145, label="(D)") +
  theme_bw() #+ 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p6 <- ggplot(cc_patch, aes(log(var.biom), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"),name = "Species group") +
  geom_smooth(aes(color = spec_group.x)) +  
  xlab("log(Temporal biomass variability)") + 
  ylab("") +
  #  geom_hline(yintercept = 0, linetype="dashed") +
  annotate(geom="text", x=-3.9, y=0.145, label="(E)") +
  theme_bw()# + 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



mylegend<-g_legend(p5)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/patch_biom_corr_contr_abs.png", sep = ""), width = 2200, height = 1000, res = 300)
grid.arrange(arrangeGrob(p5 + theme(legend.position="none"),
                         p6 + theme(legend.position="none"),
                         nrow=1))#, 
            # mylegend, nrow=2, heights=c(10, 1))
dev.off()



p9 <- ggplot(cc_patch, aes(log(dissim.bray), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  geom_smooth(aes(color = spec_group.x)) +  
  #  ylim(-0.15,0.18) +
  #  xlim(0.22,0.9) +
  xlab("Compositional Uniqueness") +
  ylab("Single site contribution\n to meta-community stability") + 
  annotate(geom="text", x = -1, y = 0.17, label="(F)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = "bottom")


p10 <- ggplot(cc_patch, aes(log(esn.mean), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  geom_smooth(aes(color = spec_group.x)) + 
  #  ylim(-0.15,0.18) +
  #  xlim(0,6) +
  xlab("Mean ESN") +
  ylab("") +
  theme_bw() + 
  annotate(geom="text", x=0.25, y=0.17, label="(G)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(legend.position = "bottom")


p11 <- ggplot(cc_patch, aes(log(esn.var), log(contr_abs))) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  geom_smooth(aes(color = spec_group.x)) + 
  #  ylim(-0.15,0.18) +
  xlab("Temporal variability of ESN") +
  ylab("") +
  annotate(geom="text", x= -2.3, y=0.17, label="(H)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw()

mylegend<-g_legend(p9)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/patch_biodiv_corr_contr_abs.png", sep = ""), width = 3300, height = 1000, res = 300)
grid.arrange(arrangeGrob(p9 + theme(legend.position="none"),
                         p10 + theme(legend.position="none"),
                         p11 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))
dev.off()




#### Figure S3: Geographic versus environmental distance ----
## Iowa
# Load spatial distance information
dist <- read.csv(".../lake_dist.csv")

dist_sub <- dist %>% 
  filter(lake1 %in% lakes & lake2 %in% lakes)

# North Sea
# Load spatial distance information
dist <- read.csv(".../dist_geo.csv")

dist_sub <- dist %>% 
  filter(Posno1 %in% station & Posno2 %in% station)

## Calculate average temporal correlation between any two sampling locations

dist <- read.csv(paste(in_path, "/Iowa lakes/Analyses/lake_dist.csv", sep = ""))
dist_sub <- dist %>% 
  dplyr::filter(lake1 %in% lakes & lake2 %in% lakes)


dist_env <- read.csv(paste(in_path, "/Iowa lakes/Analyses/lake_dist_env.csv", sep = ""))
dist_env_sub <- dist_env %>% 
  dplyr::filter(lake1 %in% lakes & lake2 %in% lakes)


for (i in 1:nrow(dist_sub)){   
  
  print(i)
  
  lake1 <- dist_sub$lake1[i]
  lake2 <- dist_sub$lake2[i]
  
  lake_pair <- phyto %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    filter(as.character(lake_name) %in% c(as.character(lake1),as.character(lake2))) %>% 
    group_by(lake_name, sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass, na.rm = T))
  
  lake_pair <- as.data.frame(lake_pair)
  
  lake_pair_wide <- merge(lake_pair[as.character(lake_pair$lake_name) == as.character(lake1),], lake_pair[as.character(lake_pair$lake_name) == as.character(lake2),], by = "sample_ID")
  
  dist_sub$corr_phyto[i] <- cor(lake_pair_wide$totbiom.x, lake_pair_wide$totbiom.y)
  
  lake_pair <- zoo %>% 
    gather(., key = species, value = biomass, ends_with(".bio")) %>%
    filter(as.character(lake_name) %in% c(as.character(lake1),as.character(lake2))) %>% 
    group_by(lake_name, sample_ID) %>% 
    dplyr::summarise(totbiom = sum(biomass, na.rm = T))
  
  lake_pair <- as.data.frame(lake_pair)
  
  lake_pair_wide <- merge(lake_pair[as.character(lake_pair$lake_name) == as.character(lake1),], lake_pair[as.character(lake_pair$lake_name) == as.character(lake2),], by = "sample_ID")
  
  dist_sub$corr_zoo[i] <- cor(lake_pair_wide$totbiom.x, lake_pair_wide$totbiom.y)
  
  
}

dist_geo_env_sync <- merge(dist_sub, dist_env_sub, by = c("lake1", "lake2"))

dist_long <- dist_geo_env_sync %>% 
  gather(., key = species_group, value = cor, corr_phyto, corr_zoo)


## Plot: synchrony vs distance

png(paste(out_path, "/Iowa lakes/Manuscript/Figures/FigureS3_phyto_zoo.png", sep = ""), width = 1400, height = 1300, res = 200)
ggplot(dist_long, aes(distance/1000,cor)) +
  geom_point(aes(color = species_group), size = 2, alpha = 0.6) + 
  geom_smooth(aes(color = species_group), size = 1) + 
  ylim(-1,1) + 
  xlab("Distance [km]") +
  ylab("Synchrony (rho)") + 
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18)) + 
  scale_color_manual(name="Species group",
                     breaks=c("corr_phyto", "corr_zoo"),
                     values =  c("#21908CFF","#FDE725FF"),
                     labels=c("Phytoplankton", "Zooplankton")) + 
  annotate(geom="text", x=0, y=1, label="(A)", size = 8) +
  theme(legend.position = c(0.8, 0.1), legend.title = element_text(size = 16), legend.text = element_text(size = 16)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()


#### Figure S4: Median biomass vs temporal variability across sites ----
## Data preparation see Figure 1

p1 <-ggplot(phyto_totbiom_spatial, aes(medbiom,varbiom)) + 
  geom_point(aes(size = totbiom), col = "#008080") + 
  xlab("Median biomass") + 
  ylab("Temporal biomass variability") +
  annotate(geom="text", x=0.03, y=350000, label="(A)", size = 4) +
  ggtitle("Phytoplankton")

p2 <-ggplot(zoo_totbiom_spatial, aes(medbiom,varbiom)) + 
  geom_point(aes(size = totbiom), col = "#99CCFF") + 
  xlab("Median biomass") + 
  ylab("Temporal biomass variability") +
  annotate(geom="text", x=1, y=520000, label="(B)", size = 4) +
  ggtitle("Zooplankton")

p3 <-ggplot(mzb_totbiom_spatial, aes(medbiom,varbiom)) + 
  geom_point(aes(size = totbiom), col = "#969696") + 
  xlab("Median biomass") + 
  ylab("Temporal biomass variability") +
  annotate(geom="text", x=0.004, y=340, label="(C)", size = 4) +
  ggtitle("Macrozoobenthos")


# helper function to extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


mylegend <- g_legend(p1)

png(paste(out_path,"Iowa lakes/Manuscript/Figures/FigureS4.png", sep = ""), width = 3600, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 , 
                         p2 , 
                         p3 , nrow = 1))
dev.off()

tiff(paste(out_path,"Iowa lakes/Manuscript/Figures/FigureS4.tiff", sep = ""), width = 3600, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 , 
                         p2 , 
                         p3 , nrow = 1))
dev.off()
#### Figure S5: Biomass dynamics of phyto- and zooplankton ----

## Phyto
phyto <- read.csv(paste(in_path, "/Iowa lakes/Data/Iowa_lakes_phyto_biom_sample_rounds.csv", sep = ""))
dim(phyto)

lakes <- unique(phyto$lake_name)
sampleID <- unique(phyto$sample_ID)

phyto_long <- phyto %>% 
  gather(., key = species, value = biom, ends_with(".bio")) %>% 
  group_by(lake_name, sample_ID) %>%   
  dplyr::summarise(tot_biom = sum(biom, na.rm = T))

p1 <- ggplot(phyto_long, aes(sample_ID, log(tot_biom))) + 
  geom_line(aes(col = lake_name, group = lake_name)) + 
#  scale_colour_gradient(low = "#00C9C4", high = "#005654") +
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = c(2,5,8,11,14,17,20),
                                    labels = c(2002,2002,2003,2004,2005,2006,2007)) +
  annotate(geom="text", x=1, y=9, label="(A)", size = 4) +
  xlab("Time (Years)") + 
  ylab("log(Total Biomass)") +
  ggtitle("Phytoplankton")
  

## Zoo
zoo <- read.csv(paste(in_path,"/Iowa lakes/Data/Iowa_lakes_zoo_biom_sample_rounds.csv", sep = ""))
dim(zoo)

zoo_long <- zoo %>% 
  gather(., key = species, value = biom, ends_with(".bio")) %>% 
  group_by(lake_name, sample_ID) %>%   
  dplyr::summarise(tot_biom = sum(biom, na.rm = T))

p2 <- ggplot(zoo_long, aes(sample_ID, log(tot_biom))) + 
  geom_line(aes(col = lake_name)) + 
  theme(legend.position = "none") + 
#  scale_colour_gradient(low = "#E5F2FF", high = "#379BFF") +
  scale_x_continuous(breaks = c(2,5,8,11,14,17,20),
                     labels = c(2002,2002,2003,2004,2005,2006,2007)) +
  annotate(geom="text", x=1, y=9.3, label="(B)", size = 4) +
  xlab("Time (Years)") + 
  ylab("log(Total Biomass)") +
  ggtitle("Zooplankton")


## MZB
dat <- read.csv("/Users/dhodapp/Documents/Hodapp, Dorothee/HIFMB/Projects/Spatio-temporal Turnover/Armonies data/Data/NorthSea_mzb_biom_annual_complete.csv")
dim(dat)
head(dat)

stations <- unique(dat$Posno)
dat$Posno <- as.character(dat$Posno)

mzb_long <- dat %>% 
  group_by(Posno, Year) %>% 
  dplyr::summarise(tot_biom = sum(biom, na.rm = T))

p3 <- ggplot(mzb_long, aes(Year, log(tot_biom))) + 
  geom_line(aes(col = Posno)) + 
  theme(legend.position = "none") + 
  #  scale_colour_gradient(low = "#E5F2FF", high = "#379BFF") +
  scale_x_continuous(breaks = c(2005,2006,2007,2008,2009,2010,2011),
                     labels = c(2005,2006,2007,2008,2009,2010,2011)) +
  annotate(geom="text", x=2005, y=5.2, label="(C)", size = 4) +
  xlab("Time (Years)") + 
  ylab("log(Total Biomass)") + 
  ggtitle("Macrozoobenthos")

png(paste(out_path,"Iowa lakes/Manuscript/Figures/FigureS5.png", sep = ""), width = 3300, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"), 
                         p2 + theme(legend.position = "none"), 
                         p3 + theme(legend.position = "none"),nrow = 1))
dev.off()

tiff(paste(out_path,"Iowa lakes/Manuscript/Figures/FigureS5.tiff", sep = ""), width = 3300, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"), 
                         p2 + theme(legend.position = "none"), 
                         p3 + theme(legend.position = "none"),nrow = 1))
dev.off()










#### Figure S6: Contribution of single species and patches to synchrony ----

## Load single species and patch contribution (stab-sync framework) --
# stab-sync
contr_spec_stab_sync_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/contr_spec.csv", sep  = ""))
contr_patch_stab_sync_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/contr_patch.csv", sep  = ""))

# species and patch info
spec_info_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/spec_info_all.csv", sep = ""))
patch_info_all <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/patch_info_all.csv", sep = "")) 
spec_info_corr <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/spec_info_corr.csv", sep = ""))
patch_info_corr <- read.csv(paste(in_path, "Iowa lakes/Analyses/stability metrics/patch_info_corr.csv", sep = ""))


# Combine species, patch and stability contribution data

cc_spec <- merge(spec_info_all, contr_spec_stab_sync_all, by = "species")

cc_patch <- merge(patch_info_all, contr_patch_stab_sync_all, by = "patch")

cc_spec_corr <-  merge(spec_info_corr, contr_spec_stab_sync_all, by = "species")

cc_patch_corr <-  merge(patch_info_corr, contr_patch_stab_sync_all, by.x = "sites", by.y = "patch")

# helper function to extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


## Contribution of species and patches to stab-sync framework --  
# Species
p1 <- ggplot(cc_spec, aes(log(mean.biom.patch), contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  geom_smooth(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  theme_bw() +
  ylab("Single species contribution\n to meta-community stability") + 
  xlab("log(Mean total species biomass)") + 
  annotate(geom="text", x=-13, y=0.19, label="(A)") +
  theme(legend.position = "bottom", legend.direction = "horizontal")

p2 <- ggplot(cc_spec, aes(log(mean.temp.var), contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  geom_smooth(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  ylab("") +
  xlab("log(Temporal species biomass variability)") +
  annotate(geom="text", x=-22, y=0.19, label="(B)") +
  theme_bw()

p3 <- ggplot(cc_spec, aes(log(biom.spatial.var), contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  geom_smooth(aes(color = spec_group.x)) +  
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  ylab("") + 
  xlab("log(Spatial species biomass variability)") +
  annotate(geom="text", x=-23.5, y=0.19, label="(C)") +
  theme_bw()

# Additional plot: single species synchrony with rest of metacommunity vs stability contribution
p3.1 <- ggplot(cc_spec_corr, aes(corr, contr)) + 
  geom_point(aes(color = spec_group)) + 
  geom_smooth(aes(color = spec_group)) +  
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  ylab("") + 
  xlab("Synchrony with all other species") +
  #  annotate(geom="text", x=-23.5, y=0.19, label="(C)") +
  theme_bw() +
  theme(legend.position="bottom")

# Additional plot: single species synchrony with rest of metacommunity vs stability contribution
p3.2 <- ggplot(cc_patch_corr, aes(corr, contr)) + 
  geom_point(aes(color = spec_group)) + 
  geom_smooth(aes(color = spec_group)) +  
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  ylab("") + 
  xlab("Synchrony with all other sites") +
  #  annotate(geom="text", x=-23.5, y=0.19, label="(C)") +
  theme_bw()


mylegend <- g_legend(p1)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/spec_biom_sync_contr.png", sep = ""), width = 3300, height = 1000, res = 300)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"), 
                         p2 + theme(legend.position = "none"), 
                         p3 + theme(legend.position = "none"), nrow = 1))#, 
#  mylegend, nrow=2, heights=c(10, 1))
dev.off()

mylegend <- g_legend(p3.1)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/spec_patch_biom_sync_contr.png", sep = ""), width = 2200, height = 1000, res = 300)
grid.arrange(arrangeGrob(p3.1 + theme(legend.position = "none"), 
                         p3.2 + theme(legend.position = "none"), nrow = 1), 
             mylegend, nrow=2, heights=c(10, 1))
dev.off()


# Patch
p4 <-ggplot(cc_patch, aes(log(tot.biom), contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  # geom_smooth(aes(color = spec_group.x)) + 
  #  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() + 
  ylab("Single site contribution\n to meta-community stability") +
  xlab("log(Total biomass)") + 
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom="text", x=2.5, y=0.145, label="(D)") +
  theme(legend.position = "bottom", legend.direction = "horizontal")

p5 <- ggplot(cc_patch, aes(log(mean.biom), contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  # geom_smooth(aes(color = spec_group.x)) + 
  xlab("log(Mean biomass)") + 
  ylab("Single site contribution\n to meta-community stability") +
  #  geom_hline(yintercept = 0, linetype="dashed") +
  annotate(geom="text", x=-4.7, y=0.145, label="(D)") +
  theme_bw() #+ 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p6 <- ggplot(cc_patch, aes(log(var.biom), contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"),name = "Species group") +
  # geom_smooth(aes(color = spec_group.x)) +  
  xlab("log(Temporal biomass variability)") + 
  ylab("") +
  #  geom_hline(yintercept = 0, linetype="dashed") +
  annotate(geom="text", x=-3.9, y=0.145, label="(E)") +
  theme_bw()# + 
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



mylegend<-g_legend(p5)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/patch_biom_sync_contr.png", sep = ""), width = 2200, height = 1000, res = 300)
grid.arrange(arrangeGrob(p5 + theme(legend.position="none"),
                         p6 + theme(legend.position="none"),
                         nrow=1))#, 
#   mylegend, nrow=2, heights=c(10, 1))
dev.off()


p9 <- ggplot(cc_patch, aes(dissim.bray, contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  #  geom_smooth(aes(color = spec_group.x)) +  
  #  ylim(-0.15,0.18) +
  #  xlim(0.22,0.9) +
  xlab("Compositional Uniqueness") +
  ylab("Single site contribution\n to meta-community stability") + 
  annotate(geom="text", x=0.24, y=0.17, label="(F)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = "bottom")


p10 <- ggplot(cc_patch, aes(esn.mean, contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  #  geom_smooth(aes(color = spec_group.x)) + 
  #  ylim(-0.15,0.18) +
  #  xlim(0,6) +
  xlab("Mean ESN") +
  ylab("") +
  theme_bw() + 
  annotate(geom="text", x=0.3, y=0.17, label="(G)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(legend.position = "bottom")


p11 <- ggplot(cc_patch, aes(esn.var, contr)) + 
  geom_point(aes(color = spec_group.x)) + 
  scale_colour_manual(values = c("#440154FF","#21908CFF","#FDE725FF"), name = "Species group") +
  #  geom_smooth(aes(color = spec_group.x)) + 
  #  ylim(-0.15,0.18) +
  xlab("Temporal variability of ESN") +
  ylab("") +
  annotate(geom="text", x=0.22, y=0.17, label="(H)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw()

mylegend<-g_legend(p9)

png(paste(out_path,"Iowa lakes/Analyses/Graphs/patch_biodiv_sync_contr.png", sep = ""), width = 3300, height = 1000, res = 300)
grid.arrange(arrangeGrob(p9 + theme(legend.position="none"),
                         p10 + theme(legend.position="none"),
                         p11 + theme(legend.position="none"),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))
dev.off()