# Biodiversity estimations: Richness and range size rarity

library(tidyverse)
library(raster)
library(rgdal)
library(RColorBrewer)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)
library(rasterVis)
library(sf)

# Tables with species, model, area
cols = c("sciname","scale","season","range_area")
splist = read_csv("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/sdm_bird_outputs_VS3.1.csv")[,cols]

# Plotting
# colors
blues = c("white", brewer.pal(9, "Blues"))
# geographic borders
countries = st_read("/vast/palmer/home.mccleary/jc3893/30x30/shapefiles/", "countries") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
canada = st_read("/vast/palmer/home.mccleary/jc3893/30x30/shapefiles/", "Canada") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
# coordinate reference system for plotting
plotcrs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83"


# version with just countries
countries_outline <- countries[c(9,38),] %>%
  st_transform(crs=plotcrs, allow_ballpark=T) %>%
  st_crop(extent(-5771203, 2991315, -1690735, 4518085))

countries = countries[c(9,44:54,56:80,82:87,89:97),] 

# version with all states/provs
countries_proj <- countries %>%
  st_union(canada) %>%
     st_transform(crs=plotcrs)

# working directory
setwd("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/VS3.1_30x30/")
dir.create("agg/")

# parallel computing
n.cores <- 3
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)


### Additive biodiversity
# Season
#foreach(season = c("nonbreeding","breeding","resident")) %dopar% {
  for (season in c("nonbreeding","breeding","resident")){
  library(tidyverse)
  library(raster)
  library(rgdal)
  library(RColorBrewer)
  library(stringr)
  library(parallel)
  library(doParallel)
  library(foreach)
  # Across scales
  for (scale in c(1,3,5,10,50)){
    # Subset species list
    splist_sub = splist[splist$scale==scale &
                          splist$season==season,]
    # Files
    files = list.files(pattern=paste0("_",season,"_",scale,"_00.tif"),
                       recursive=T, full.names=T)
    files = str_subset(files, "PA")
    # remove species without predictions
    files = files[grepl(paste0(splist_sub[,1][[1]],collapse = "|"), files)]
    # Richness - add all species predictions
    if(file.exists(paste0("agg/",season,"_",scale,"_rich.tif"))==F){
    for (file in 1:length(files)){
      try({
      tif = files[file]
      sciname = gsub("./","",dirname(dirname(dirname(tif))))
      print(sciname)
        r_tif = raster(tif)
        if(file==1){rich = r_tif
        }else{
          rich = sum(rich, r_tif, na.rm=T)
        }
      })
    }
      # write raster and plot
    writeRaster(rich, paste0("agg/",season,"_",scale,"_rich.tif"),
                overwrite=T)
    jpeg(paste0("agg/",season,"_",scale,"_rich.jpeg"),
         width=1000, height=1000, units="px")
    plot(rich, col=blues)
    plot(countries, col="transparent", border="black",
         usePolypath=F, add=T)
    dev.off()
    }
    # Range size rarity
    if(file.exists(paste0("agg/",season,"_",scale,"_rsr.tif"))==F){
    for (file in 1:length(files)){
      try({
      tif = files[file]
      sciname = gsub("./","",dirname(dirname(dirname(tif))))
      print(sciname)
      area = splist_sub[splist_sub$sciname==sciname,]$range_area
      # add weighted range for each species
      if(length(area)>0){
        w_tif = raster(tif)*(1/log(area))
        if(file==1){rsr = w_tif
        }else{
          rsr = sum(rsr, w_tif, na.rm=T)
        }
      }
      })
    }
    writeRaster(rsr, paste0("agg/",season,"_",scale,"_rsr.tif"),
                overwrite=T)
    jpeg(paste0("agg/",season,"_",scale,"_rsr.jpeg"),
         width=1000, height=1000, units="px")
    plot(rsr, col=blues)
    plot(countries, col="transparent", border="black",
         usePolypath=F, add=T)
    dev.off()
    }
  }
}

# Nonbreeding/breeding + resident predictions for seasonal richness/rarity
for (type in c("rich","rsr")){
  for (scale in c(1,3,5,10,50)){
    res <- raster(paste0("agg/resident_",scale,"_",type,".tif"))
    for(season in c("breeding","nonbreeding")){
      r <- raster(paste0("agg/",season,"_",scale,"_",type,".tif"))
      r <- r + res
      writeRaster(r, paste0("agg/",season,"_",scale,"_",type,"_total.tif"),
                  overwrite=T)
      jpeg(paste0("agg/",season,"_",scale,"_",type,"_total.jpeg"),
           width=1000, height=1000, units="px")
      plot(r, col=blues)
      plot(countries, col="transparent", border="black",
           usePolypath=F, add=T)
      dev.off()
    }
  }
}


# Prep rasters for nicer plots - remove predictions outside US/Canada, project
for (type in c("rich", "rsr")){
  for(season in c("breeding","nonbreeding")){
    for (scale in c(1,3,5,10,50)){
      r <- raster(paste0("agg/",season,"_",scale,"_",type,"_total.tif")) %>%
        crop(countries) %>%
        mask(mask = countries) %>% 
        projectRaster(crs=plotcrs) %>%
        trim()
      writeRaster(r, paste0("agg/",season,"_",scale,"_",type,"_v3.tif"),
                  overwrite=T)
    }
  }
}

# Nicer plots
# colors
bp <- rev(brewer.pal(9, "Spectral"))
colors = c("#0F70DE" ,bp[1:2], bp[3], bp[5:9], "#B30707")
# plot with new colors
for (type in c("rich","rsr")){
  if(type=="rich"){zrange=c(0,184)}else{zrange=c(0,13)}
  jpeg(paste0("agg/agg_",type,"_v3.jpeg"),
       width=1200, height=400, units="px")
  par(mfrow=c(2,5), mar=c(1.1,.6,1.1,0))
  for(season in c("breeding","nonbreeding")){
    for (scale in c(1,3,5,10,50)){
  r <- raster(paste0("agg/",season,"_",scale,"_",type,"_v3.tif"))
  plot(r, col=colors, zlim=zrange, legend=F, xaxt='n', yaxt='n', box=F)
  plot(countries_outline, col="transparent", border="black",
       usePolypath=F, add=T)
}
}
  dev.off()
}


# Example species (fig. 2a-c)
# indigo bunting = steelblue1 "Passerina_cyanea"
# sedge wren = burlywood4 "Cistothoru_stellaris"
# lesser goldfinch = gold1 "Spinus_psaltria"
sp_name = c("Passerina_cyanea", "Cistothoru_stellaris", "Spinus_psaltria")
col = c("steelblue1", "coral4", "gold1")
season = c("breeding","breeding","breeding")
for (sp in 1:3){
jpeg(paste0("agg/",sp_name[sp],"_scale_comp3.jpeg"),
     width=1000, height=250, units="px")
par(mfrow=c(1,5), mar=c(1.1,1.1,1.1,0))
# loop scales
for (scale in c(1,3,5,10,50)){
  # get raster, crop, mask, project
r <- raster(paste0(sp_name[sp],"/predictions/PA/Thresholded_predictions_",
                   sp_name[sp],"_",season[sp],"_",scale,"_00.tif")) %>%
  crop(countries) %>%
  mask(mask = countries) %>%
  projectRaster(crs=plotcrs)
# if scale is 1, trim, otherwise disaggregate and process (to match extent/grain for stacking)
if(scale==1){r1 = r %>% trim()
r1test <- reclassify(r1, cbind(-Inf, 1.5, 2))
r1test[is.na(r1test[])] <- 0 
plot(r1, col=col[sp], legend=F, xaxt='n', yaxt='n', box=F, axes=F)
}else{
rda <- r %>% 
  disaggregate(scale) %>%
  projectRaster(r1)
rda[is.na(rda[])] <- 0 
rda <- reclassify(rda, cbind(.5, 1.5, 1))
r <- r1test-rda %>%
  trim()
# plotting
cpal <- c('gray80','transparent', col[sp], 'black')
r2 <- ratify(r)
# # for legend
# levelplot(r2,col.regions=cpal,att='ID', 
#           scales=list(x=list(at=NULL), y=list(at=NULL)),
#           par.settings = list(axis.line = list(col = "transparent")), xaxt='n', yaxt='n', box=F, axes=F)
plot(r, col=cpal, legend=F, xaxt='n', yaxt='n', box=F, axes=F)}
plot(countries_proj, col="transparent", border="black",
     usePolypath=F, add=T)
assign(paste0("r.",scale),r)
}
dev.off()
}


# Delta biodiversity maps
bp <- brewer.pal(9, "BrBG")
colors = rev(c(rep(bp[1],1), bp[2:4], rep(bp[5],2), bp[6:8], rep(bp[9],1)))

for (type in c("rich", "rsr")){ 
  if(type=="rich"){sequ <- seq(-50,50,10)}else{sequ <- seq(-5,5,1)}
  for(season in c("breeding","nonbreeding")){ 
    r1 <- raster(paste0("agg/",season,"_1_",type,"_v3.tif"))
    for (scale in c(3,5,10,50)){ 
      r <- raster(paste0("agg/",season,"_",scale,"_",type,"_v3.tif"))
      rda <- r %>% 
        disaggregate(scale)
      r1 <- projectRaster(r1, rda)
      rd <- rda-r1
      writeRaster(rd, paste0("agg/",season,"_",scale,"_",type,"_delta_v3.tif"),
                  overwrite=T)
      rd <- aggregate(rd, 4)
      if(type=="rich"){
      rd[rd>50] <- NA; rd[rd<(-50)] <- NA
      }else{
        rd[rd>5] <- NA; rd[rd<(-5)] <- NA }
      jpeg(paste0("agg/",season,"_",scale,"_",type,"_delta_v3.jpeg"),
           width=500, height=400, units="px")
      plot(rd, breaks=sequ, col=colors)
      plot(countries_proj, col="transparent", border="black",
           usePolypath=F, add=T)
      dev.off()
    }
  }
}


# Creating figure 3a-b
# colors
bp <- rev(brewer.pal(9, "Spectral"))
colors = c("#0F70DE" ,bp[1:2], bp[3], bp[5:9], "#B30707")
bp2 <- brewer.pal(9, "BrBG")
colors2 = rev(c(rep(bp2[1],1), bp2[2:4], rep(bp2[5],2), bp2[6:8], rep(bp2[9],1)))

jpeg("agg/fig1.jpeg", width=1200, height=400, units="px")
par(mfrow=c(2,5), mar=c(1.1,.6,1.1,0))
# loop seasons, combine plots across seasons and grains
  for(season in c("breeding","nonbreeding")){
      r <- raster(paste0("agg/",season,"_1_rich_v3.tif"))
      plot(r, col=colors, zlim=c(0,184), legend=F, 
           xaxt='n', yaxt='n', box=F, axes=F)
      plot(countries_outline, col="transparent", border="black",
           usePolypath=F, add=T, box=F)
      for (scale in c(3,5,10,50)){
       d <- raster(paste0("agg/",season,"_",scale,"_rich_delta_v3.tif")) 
       plot(d, col=colors2, breaks=seq(-50,50,10), legend=F, 
            xaxt='n', yaxt='n', box=F, axes=F)
       plot(countries_outline, col="transparent", border="black",
            usePolypath=F, add=T, box=F)
    }
  }
  dev.off()

  # Creating figure s5
  
  jpeg("agg/figs5.jpeg",
       width=1200, height=400, units="px")
  par(mfrow=c(2,5), mar=c(1.1,.6,1.1,0))
  
  for(season in c("breeding","nonbreeding")){
    for (scale in c(1,3,5,10,50)){
      r <- raster(paste0("agg/",season,"_",scale,"_rich_v3.tif")) 
      plot(r, col=colors, zlim=c(0,184), legend=F, 
           xaxt='n', yaxt='n', box=F, axes=F)
      plot(countries_outline, col="transparent", border="black",
           usePolypath=F, add=T, box=F)
    }
  }
  dev.off()

  
  # Creating figure s6
  
  jpeg("agg/figs6.jpeg",
       width=1200, height=400, units="px")
  par(mfrow=c(2,5), mar=c(1.1,.6,1.1,0))
  
  for(season in c("breeding","nonbreeding")){
    r <- raster(paste0("agg/",season,"_1_rsr_v3.tif"))
    plot(r, col=colors, zlim=c(0,13), legend=F, 
         xaxt='n', yaxt='n', box=F, axes=F)
    plot(countries_outline, col="transparent", border="black",
         usePolypath=F, add=T, box=F)
    for (scale in c(3,5,10,50)){
      d <- raster(paste0("agg/",season,"_",scale,"_rsr_delta_v3.tif")) 
      plot(d, col=colors2, breaks=seq(-5,5,1), legend=F, 
           xaxt='n', yaxt='n', box=F, axes=F)
      plot(countries_outline, col="transparent", border="black",
           usePolypath=F, add=T, box=F)
    }
  }
  dev.off() 
  
  
  # Aggregated biodiversity
  
  # Season
  foreach(season = c("nonbreeding","breeding","resident")) %dopar% {
    #for (season in c("nonbreeding","breeding","resident")){
    library(tidyverse)
    library(raster)
    library(rgdal)
    library(RColorBrewer)
    library(stringr)
    library(parallel)
    library(doParallel)
    library(foreach)
    # Across scales
    files = list.files(pattern=paste0("_",season,"_1_00.tif"),
                       recursive=T, full.names=T)
    files = str_subset(files, "PA")
    for (scale in c(3,5,10,50)){
      # Subset
      splist_sub = splist[splist$scale==scale &
                            splist$season==season,]
      # remove sp not included
      files_sub = files[grepl(paste0(splist_sub[,1][[1]],collapse = "|"), files)]
      # Richness - upscale predictions from 1km
      if(file.exists(paste0("agg/",season,"_",scale,"_rich_agg.tif"))==F){
        for (file in 1:length(files_sub)){
          try({
            tif = files_sub[file]
            sciname = gsub("./","",dirname(dirname(dirname(tif))))
            print(sciname)
            r_tif = raster(tif) %>%
              aggregate(scale, modal)
            if(file==1){rich = r_tif
            }else{
              rich = sum(rich, r_tif, na.rm=T)
            }
          })
        }
        writeRaster(rich, paste0("agg/",season,"_",scale,"_rich_agg.tif"),
                    overwrite=T)
        jpeg(paste0("agg/",season,"_",scale,"_rich_agg.jpeg"),
             width=1000, height=1000, units="px")
        plot(rich, col=blues)
        plot(countries, col="transparent", border="black",
             usePolypath=F, add=T)
        dev.off()
      }
      # # Range size rarity
      # if(file.exists(paste0("agg/",season,"_",scale,"_rsr_agg.tif"))==F){
      #   for (file in 1:length(files_sub)){
      #     try({
      #       tif = files_sub[file]
      #       sciname = gsub("./","",dirname(dirname(dirname(tif))))
      #       print(sciname)
      #       area = splist_sub[splist_sub$sciname==sciname,]$range_area
      #       if(length(area)>0){
      #         w_tif = raster(tif)*(1/log(area))
      #         w_tif = aggregate(w_tif, scale, modal)
      #         if(file==1){rsr = w_tif
      #         }else{
      #           rsr = sum(rsr, w_tif, na.rm=T)
      #         }
      #       }
      #     })
      #   }
      #   writeRaster(rsr, paste0("agg/",season,"_",scale,"_rsr_agg.tif"),
      #               overwrite=T)
      #   jpeg(paste0("agg/",season,"_",scale,"_rsr_agg.jpeg"),
      #        width=1000, height=1000, units="px")
      #   plot(rsr, col=blues)
      #   plot(countries, col="transparent", border="black",
      #        usePolypath=F, add=T)
      #   dev.off()
      # }
    }
  }
  
  # Nonbreeding/breeding + residents for seasonal richness/rarity
  for (type in c("rich")){
    for (scale in c(3,5,10,50)){
      res <- raster(paste0("agg/resident_",scale,"_",type,"_agg.tif"))
      for(season in c("breeding","nonbreeding")){
        if(file.exists(paste0("agg/",season,"_",scale,"_",type,"_total_agg.tif"))==F){
          r <- raster(paste0("agg/",season,"_",scale,"_",type,"_agg.tif"))
          r <- r + res
          writeRaster(r, paste0("agg/",season,"_",scale,"_",type,"_total_agg.tif"),
                      overwrite=T)
        }
      }
    }
  }
  
  # Nicer plots
  for (type in c("rich")){
    for (scale in c(3,5,10,50)){
      for(season in c("breeding","nonbreeding")){
        if(file.exists(paste0("agg/",season,"_",scale,"_",type,"_agg_v3.tif"))==F){
          r <- raster(paste0("agg/",season,"_",scale,"_",type,"_total_agg.tif")) %>%
            crop(countries) %>%
            mask(mask = countries) %>%
            projectRaster(crs=plotcrs) %>%
            trim()
          writeRaster(r, paste0("agg/",season,"_",scale,"_",type,"_agg_v3.tif"),
                      overwrite=T)
        }
      }
    }
  }
  
  # Additive vs Aggregated figure (fig s7)
colors <- brewer.pal(9, "YlOrRd")
  
  jpeg("agg/fig_iva.jpeg",
       width=1000, height=400, units="px")
  par(mfrow=c(2,4), mar=c(1.1,.6,1.1,0))
  for (type in c("rich")){ 
    if(type=="rich"){sequ <- seq(0,180,20)}else{sequ <- seq(0,18,2)}
    for(season in c("breeding","nonbreeding")){
      for (scale in c(3,5,10,50)){ 
        # get deltas between additive and aggregated methods
        r_agg <- raster(paste0("agg/",season,"_",scale,"_",type,"_agg_v3.tif"))
        r <- raster(paste0("agg/",season,"_",scale,"_",type,"_v3.tif"))
        r <- projectRaster(r, r_agg)
        rd <- r_agg-r
        writeRaster(rd, paste0("agg/",season,"_",scale,"_",type,"_comparison_delta.tif"),
                    overwrite=T)
        if(type=="rich"){
          rd[rd>180] <- NA; rd[rd<(-200)] <- NA
        }else{
          rd[rd>18] <- NA; rd[rd<(-20)] <- NA }
        rd <- rd %>% 
          projectRaster(crs=plotcrs) %>%
          trim()
        plot(rd, breaks=sequ, col=colors, legend=F, 
             xaxt='n', yaxt='n', box=F, axes=F)
        plot(countries_outline, col="transparent", border="black",
             usePolypath=F, add=T)
      }
    }
  }
  dev.off()
  
  