# Exploring influence of filtering and thinning data (fig s1, table s1)

library(devtools)
install_github('mcooper/moranfast')
library(moranfast)  
library(ff)
library(sf)
library(dismo)
library(rgeos)
library(rgdal)
library(maps)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)

setwd("/vast/palmer/home.mccleary/jc3893/")

# Species names
list <- read_csv("30x30/species_list.csv")

# white-breasted nuthatch
sci_name <- as.character(list[468,2])
sci.name <- gsub("_"," ",sci_name)

# season (do both)
# season="breeding" 
season="nonbreeding" 

# get data
all_preds <- c(
  "year","day_of_year","time_observations_started","duration_minutes",
  "effort_distance_km","number_observers","bio1","bio12","bio15","cloudsd",
  "evisum","tri","elev","pland_10_cropland_rainfed", 
  "pland_100_mosaic_tree_shrub","pland_11_cropland_rainfed", 
  "pland_110_mosaic_herbacious","pland_12_cropland_rainfed", 
  "pland_120_shrubland","pland_121_shrubland","pland_122_shrubland",
  "pland_130_grassland","pland_140_lichens_mosses","pland_150_sparse",
  "pland_152_sparse","pland_153_sparse","pland_160_flooded_freshwater",
  "pland_170_flooded_saltwater","pland_180_flooded_shrub", 
  "pland_190_urban","pland_20_cropland_irrigated","pland_200_barren",     
  "pland_201_barren","pland_202_barren","pland_210_water",
  "pland_220_ice","pland_30_mosaic_cropland","pland_40_mosaic_natural_veg",
  "pland_50_evergreen_broadleaf","pland_60_deciduous_broadleaf", 
  "pland_61_deciduous_broadleaf","pland_62_deciduous_broadleaf", 
  "pland_70_evergreen_needleleaf","pland_71_evergreen_needleleaf",
  "pland_72_evergreen_needleleaf","pland_80_deciduous_needleleaf",
  "pland_81_deciduous_needleleaf","pland_82_deciduous_needleleaf",
  "pland_90_mixed_forest")
ebird_full <- read.csv.ffdf(file=paste0("30x30/ebird_",season,"_new.csv"),
                            colClasses=c(rep("factor", 2), rep("numeric",2),
                                         rep("integer",2), "numeric", "integer",
                                         "numeric","integer",rep("integer",677),
                                         rep("numeric",44)))
ebird <- data.frame(ebird_full[,c(gsub(" ",".",sci.name),all_preds,
                                 "latitude","longitude","checklist_id")])

# check limits
summary(ebird$duration_minutes)
summary(ebird$effort_distance_km)
summary(ebird$number_observers)

### unfiltered, unthinned data
# n
n.full = nrow(ebird)
pres = ebird[ebird$Sitta.carolinensis==1,]
abs = ebird[ebird$Sitta.carolinensis==0,]
np.full = nrow(pres)
na.full = nrow(abs)
# centroids
cent.full = c(mean(ebird$latitude), mean(ebird$longitude))
centp.full = c(mean(pres$latitude), mean(pres$longitude))
centa.full = c(mean(abs$latitude), mean(abs$longitude))
# autocorrelation
out.full = c()
for (i in 1:50){
ebirdsub = sample_n(ebird, 50000)
auto.full = moranfast(ebirdsub$Sitta.carolinensis, 
                      ebirdsub$longitude, ebirdsub$latitude)
out.full = c(out.full, auto.full$observed)
print(i)
}
mean.full = mean(out.full)
sd.full = sd(out.full)

{ # visualization
  world <- ne_countries(scale='medium', returnclass = 'sf') 
  pres_st = st_as_sf(pres, coords = c("longitude", "latitude"), 
                     crs = crs(world))
  lim <- st_bbox(pres_st)            
  # plot (all panels form figure s1)
  map <- ggplot() +
    geom_sf(data=world, bg="gray95") + 
    coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)], expand=F) +
    theme_bw() +
    theme(text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=13),
          axis.text.y = element_text(color = "black", size=13),
          plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
    geom_point(abs, mapping=aes(x=longitude, y=latitude), col="gray50", size=.2) +
    geom_point(pres, mapping=aes(x=longitude, y=latitude), col="black", bg="white", size=.5) +
    xlab("") +
    ylab("")
  ggsave(paste0("30x30/",sci_name,"_point_map_with_abs_full_",season,".jpeg"), map, "jpeg",
         height=5, width=7, units="in")
  }


### filtered data
ebird.filt = ebird[ebird$duration_minutes <= 300,]
ebird.filt = ebird.filt[ebird.filt$effort_distance_km <= 1,]
# sample size
n.filt = nrow(ebird.filt)
pres = ebird.filt[ebird.filt$Sitta.carolinensis==1,]
abs = ebird.filt[ebird.filt$Sitta.carolinensis==0,]
np.filt = nrow(pres)
na.filt = nrow(abs)
# centroids
cent.filt = c(mean(ebird.filt$latitude), mean(ebird.filt$longitude))
centp.filt = c(mean(pres$latitude), mean(pres$longitude))
centa.filt = c(mean(abs$latitude), mean(abs$longitude))
# autocorrelation
out.filt = c()
for (i in 1:50){
  ebirdsub = sample_n(ebird.filt, 50000)
  auto.filt = moranfast(ebirdsub$Sitta.carolinensis, 
                        ebirdsub$longitude, ebirdsub$latitude)
  out.filt = c(out.filt, auto.filt$observed)
  print(i)
}
mean.filt = mean(out.filt)
sd.filt = sd(out.filt)

{ # visualization
  world <- ne_countries(scale='medium', returnclass = 'sf') 
  pres_st = st_as_sf(pres, coords = c("longitude", "latitude"), 
                     crs = crs(world))
  lim <- st_bbox(pres_st)            
  # plot
  map <- ggplot() +
    geom_sf(data=world, bg="gray95") + 
    coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)], expand=F) +
    theme_bw() +
    theme(text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=13),
          axis.text.y = element_text(color = "black", size=13),
          plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
    geom_point(abs, mapping=aes(x=longitude, y=latitude), col="gray50", size=.2) +
    geom_point(pres, mapping=aes(x=longitude, y=latitude), col="black", bg="white", size=.5) +
    xlab("") +
    ylab("")
  ggsave(paste0("30x30/",sci_name,"_point_map_with_abs_filt_",season,".jpeg"), map, "jpeg",
         height=5, width=7, units="in")
}



### Spatiotemporally thinned data
ebird$week = ceiling(ebird$day_of_year/7)
ebird$lat = round(ebird$latitude, 2)
ebird$lon = round(ebird$longitude, 2)
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird.thin <- ebird %>% 
  group_by(Sitta.carolinensis, year, week, lat, lon) %>% 
  sample_n(size = 1) %>% 
  ungroup()
# n
n.thin = nrow(ebird.thin)
pres = ebird.thin[ebird.thin$Sitta.carolinensis==1,]
abs = ebird.thin[ebird.thin$Sitta.carolinensis==0,]
np.thin = nrow(pres)
na.thin = nrow(abs)
# centroids
cent.thin = c(mean(ebird.thin$latitude), mean(ebird.thin$longitude))
centp.thin = c(mean(pres$latitude), mean(pres$longitude))
centa.thin = c(mean(abs$latitude), mean(abs$longitude))
# autocorrelation
out.thin = c()
for (i in 1:50){
  ebirdsub = sample_n(ebird.thin, 50000)
  auto.thin = moranfast(ebirdsub$Sitta.carolinensis, 
                        ebirdsub$longitude, ebirdsub$latitude)
  out.thin = c(out.thin, auto.thin$observed)
  print(i)
}
mean.thin = mean(out.thin)
sd.thin = sd(out.thin)

{ # visualization
  world <- ne_countries(scale='medium', returnclass = 'sf') 
  pres_st = st_as_sf(pres, coords = c("longitude", "latitude"), 
                     crs = crs(world))
  lim <- st_bbox(pres_st)            
  # plot
  map <- ggplot() +
    geom_sf(data=world, bg="gray95") + 
    coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)], expand=F) +
    theme_bw() +
    theme(text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=13),
          axis.text.y = element_text(color = "black", size=13),
          plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
    geom_point(abs, mapping=aes(x=longitude, y=latitude), col="gray50", size=.2) +
    geom_point(pres, mapping=aes(x=longitude, y=latitude), col="black", bg="white", size=.5) +
    xlab("") +
    ylab("")
  ggsave(paste0("30x30/",sci_name,"_point_map_with_abs_thin_",season,".jpeg"), map, "jpeg",
         height=5, width=7, units="in")
}




### Filtered and thinned data
ebird.final = ebird.thin[ebird.thin$duration_minutes <= 300,]
ebird.final = ebird.final[ebird.final$effort_distance_km <= 1,]
# n
n.final = nrow(ebird.final)
pres = ebird.final[ebird.final$Sitta.carolinensis==1,]
abs = ebird.final[ebird.final$Sitta.carolinensis==0,]
np.final = nrow(pres)
na.final = nrow(abs)
# centroids
cent.final = c(mean(ebird.final$latitude), mean(ebird.final$longitude))
centp.final = c(mean(pres$latitude), mean(pres$longitude))
centa.final = c(mean(abs$latitude), mean(abs$longitude))
# autocorrelation
out.final = c()
for (i in 1:50){
  ebirdsub = sample_n(ebird.final, 50000)
  auto.final = moranfast(ebirdsub$Sitta.carolinensis, 
                        ebirdsub$longitude, ebirdsub$latitude)
  out.final = c(out.final, auto.final$observed)
  print(i)
}
mean.final = mean(out.final)
sd.final = sd(out.final)

{ # visualization
  world <- ne_countries(scale='medium', returnclass = 'sf') 
  pres_st = st_as_sf(pres, coords = c("longitude", "latitude"), 
                     crs = crs(world))
  lim <- st_bbox(pres_st)            
  # plot
  map <- ggplot() +
    geom_sf(data=world, bg="gray95") + 
    coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)], expand=F) +
    theme_bw() +
    theme(text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=13),
          axis.text.y = element_text(color = "black", size=13),
          plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
    geom_point(abs, mapping=aes(x=longitude, y=latitude), col="gray50", size=.2) +
    geom_point(pres, mapping=aes(x=longitude, y=latitude), col="black", bg="white", size=.5) +
    xlab("") +
    ylab("")
  ggsave(paste0("30x30/",sci_name,"_point_map_with_abs_final_",season,".jpeg"), map, "jpeg",
         height=5, width=7, units="in")
}


# collect all - table S1
table = rbind(c(n.full, np.full, na.full,
  paste0(round(cent.full[1],3),", ", round(cent.full[2],3)),
                paste0(round(centp.full[1],3),", ", round(centp.full[2],3)),
                paste0(round(centa.full[1],3),", ", round(centa.full[2],3)), 
                paste0(round(mean.full,3), " (", round(sd.full,3),")")),
              c(n.filt, np.filt, na.filt,
                paste0(round(cent.filt[1],3),", ", round(cent.filt[2],3)),
                paste0(round(centp.filt[1],3),", ", round(centp.filt[2],3)),
                paste0(round(centa.filt[1],3),", ", round(centa.filt[2],3)), 
                paste0(round(mean.filt,3), " (", round(sd.filt,3),")")),
              c(n.thin, np.thin, na.thin,
                paste0(round(cent.thin[1],3),", ", round(cent.thin[2],3)),
                paste0(round(centp.thin[1],3),", ", round(centp.thin[2],3)),
                paste0(round(centa.thin[1],3),", ", round(centa.thin[2],3)), 
                paste0(round(mean.thin,3), " (", round(sd.thin,3),")")),
              c(n.final, np.final, na.final,
                paste0(round(cent.final[1],3),", ", round(cent.final[2],3)),
                paste0(round(centp.final[1],3),", ", round(centp.final[2],3)),
                paste0(round(centa.final[1],3),", ", round(centa.final[2],3)), 
                paste0(round(mean.final,3), " (", round(sd.final,3),")")))
table = data.frame(table)
names(table) = c("Total Observations", "Presence Observations", "Absence Observations", 
                 "Centroid (all data)", "Centroid (presence)",
                 "Centroid (absence)", "Moran's I (Mean, SD)") 
write_csv(table, paste0("30x30/ebird_exploration_",season,".csv"))


