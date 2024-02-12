# Richness and range size rarity based on SDMs

library(tidyverse)
library(raster)
library(rgdal)
library(RColorBrewer)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)

# Tables with species, model, area
cols = c("sciname","landcover","season","range_area")
path = "/gpfs/loomis/pi/jetz/data/results_30x30/sdm_birds"
splist = read_csv(paste0(path,"/sdm_birds_V2.1.csv"))[,cols]

# Plotting
countries = readOGR("30x30/shapefiles", "countries")
blues = c("white", brewer.pal(9, "Blues"))

# wd
setwd(paste0(path,"/V2.1_30x30/"))

# # parallel computing
# n.cores <- 3
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "PSOCK"
# )
# doParallel::registerDoParallel(cl = my.cluster)
# 
# # Season
# foreach(season = c("nonbreeding","breeding","resident")) %dopar% {
#   #for (season in c("nonbreeding","breeding","resident")){
#   library(tidyverse)
#   library(raster)
#   library(rgdal)
#   library(RColorBrewer)
#   library(stringr)
#   library(parallel)
#   library(doParallel)
#   library(foreach)
#   # With/out landcover
#   for (lc in c("landcover","no_landcover")){
#     # Subset
#     splist_sub = splist[splist$landcover==lc &
#                           splist$season==season,]
#     # Files
#     files = list.files(pattern=paste0("_",season,"_",lc,"_00.tif"),
#                        recursive=T, full.names=T)
#     files = str_subset(files, "PA")
#     # Richness
#     if(file.exists(paste0("agg/",season,"_",lc,"_rich.tif"))==F){
#     for (file in 1:length(files)){
#       try({
#       tif = files[file]
#       sciname = gsub("./","",dirname(dirname(dirname(tif))))
#       print(sciname)
#         r_tif = raster(tif)
#         if(file==1){rich = r_tif
#         }else{
#           rich = sum(rich, r_tif, na.rm=T)
#         }
#       })
#     }
#     writeRaster(rich, paste0("agg/",season,"_",lc,"_rich.tif"),
#                 overwrite=T)
#     jpeg(paste0("agg/",season,"_",lc,"_rich.jpeg"),
#          width=1000, height=1000, units="px")
#     plot(rich, col=blues)
#     plot(countries, col="transparent", border="black",
#          usePolypath=F, add=T)
#     dev.off()
#     }
#     # Range size rarity
#     if(file.exists(paste0("agg/",season,"_",lc,"_trsr.tif"))==F){
#     for (file in 1:length(files)){
#       try({
#       tif = files[file]
#       sciname = gsub("./","",dirname(dirname(dirname(tif))))
#       print(sciname)
#       area = splist_sub[splist_sub$sciname==sciname,]$range_area
#       if(length(area)>0){
#         w_tif = raster(tif)*(1/area)
#         if(file==1){rsr = w_tif
#         }else{
#           rsr = sum(rsr, w_tif, na.rm=T)
#         }
#       }
#       })
#     }
#     writeRaster(rsr, paste0("agg/",season,"_",lc,"_trsr.tif"),
#                 overwrite=T)
#     jpeg(paste0("agg/",season,"_",lc,"_trsr.jpeg"),
#          width=1000, height=1000, units="px")
#     plot(rsr, col=blues)
#     plot(countries, col="transparent", border="black",
#          usePolypath=F, add=T)
#     dev.off()
#     }
#   }
# }

# Nonbreeding/breeding + residents for actual richness/rarity
for (lc in c("landcover", "no_landcover")){
  for(season in c("breeding","nonbreeding")){
    for (type in c("rich","trsr")){
      res <- raster(paste0("agg/resident_",lc,"_",type,".tif"))
      r <- raster(paste0("agg/",season,"_",lc,"_",type,".tif"))
      r <- r + res
      writeRaster(r, paste0("agg/",season,"_",lc,"_",type,"_total.tif"),
                  overwrite=T)
      jpeg(paste0("agg/",season,"_",lc,"_",type,"_total.jpeg"),
           width=1000, height=1000, units="px")
      plot(r, col=blues)
      plot(countries, col="transparent", border="black",
           usePolypath=F, add=T)
      dev.off()
      if(type=="rich"){rich = r}else{rsr = r}
    }
    # Average range-size rarity
    arsr = rsr/rich
    writeRaster(arsr, paste0("agg/",season,"_",lc,"_arsr_total.tif"),
                overwrite=T)
    jpeg(paste0("agg/",season,"_",lc,"_arsr_total.jpeg"),
         width=1000, height=1000, units="px")
    plot(arsr, col=blues)
    plot(countries, col="transparent", border="black",
         usePolypath=F, add=T)
    dev.off()
  }
}

# Delta nolc-lc
for (type in c("rich","trsr")){
    for(season in c("breeding","nonbreeding")){
      lc <- raster(paste0("agg/",season,"_landcover_",type,"_total.tif"))
      nolc <- raster(paste0("agg/",season,"_no_landcover_",type,"_total.tif"))
      dif <- nolc - lc
      writeRaster(dif, paste0("agg/",season,"_delta_",type,".tif"),
                  overwrite=T)
    }
  }
