# Validating models against site lists

library(lubridate)
library(sf)
library(rgdal)
library(gridExtra)
library(tidyverse)
library(raster)
library(RColorBrewer)
library(MetBrewer)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(Metrics)
library(data.table)
setwd("~/ranges/")
memory.limit(300000)

# grid100 <- raster("grids/Global_reference_raster_100km_CEA.tif")
# grid200 <- aggregate(grid100, 2)
# writeRaster(grid200, "grids/Global_reference_raster_200km_forval_CEA.tif")

countries <- readOGR("shapefiles", "countries")
bbox <- readOGR("shapefiles", "30x30BBox_NAmerica") %>%
  spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
  st_as_sf()
grid200 <- raster("grids/Global_reference_raster_200km_CEA.tif") %>%
  crop(bbox)
grid200 <- setValues(grid200, 1:ncell(grid200))

# richness and species at the top birded site per cell
for (season in c("breeding", "nonbreeding")){
  # season = "nonbreeding"
  ebird <- read_csv(paste0("ebird_",season,"_CEA.csv"))
  for (cell in 1:ncell(grid200)){
    if(cell==1 & season=="breeding"){rm(out)}
    print(cell)
    rcell <- rasterFromCells(grid200, cell)
    ecell <- ebird[ebird$x > xmin(rcell) &
                     ebird$x < xmax(rcell) &
                     ebird$y > ymin(rcell) &
                     ebird$y < ymax(rcell),]
    if (nrow(ecell) >= 100){
      top <- modal(as.factor(ecell$locality_id))
      esite <- ecell[ecell$locality_id==top,]
      if (nrow(esite) >= 100){
        #   for (thresh in c(0, .02, .05, .1)){
        thresh = 0
        sp_pres <- ifelse(colMeans(esite[,c(11:687)]) > thresh, 1, 0)
        siterich <- sum(as.numeric(sp_pres))
        #   }
        total <- unlist(c(cell, season, esite[1,"x"], esite[1,"y"], siterich, esite$locality_id[[1]],
                          #, rich.0.02, rich.0.05, rich.0.1 #reinstate assigning names
                          sp_pres))
        if(exists("out")){out <- rbind(out,total)}else{out <- total}
      }
    }
  }
}
out <- data.frame(out)
colnames(out) <- c("cell","season","x","y","rich.0", "locality_id",
                   #"rich.02","rich.05","rich.1",
                   gsub(" ","_",names(sp_pres)))
rownames(out) <- NULL
write_csv(out, "site_richness_validation_withsp0.csv")

# Add estimated richness across scales 
out <- read_csv("site_richness_validation_withsp0.csv")
coordinates(out) <- c("x","y")
for (season in c("breeding", "nonbreeding")){
  stack <- stack()
  for (scale in c(1,3,5,10,50)){
    r <- raster(paste0("agg/",season,"_",scale,"_rich_total_cropped.tif"))
    if(scale>1){r <- projectRaster(r, stack)}
    stack <- stack(stack, r)
  }
  ests <- raster::extract(stack, out)
  colnames(ests) <- paste0("est.",c(1,3,5,10,50))
  assign(paste0("ests.",season), ests)
}
val <- cbind(data.frame(out), rbind(ests.breeding, ests.nonbreeding))
val$optional <- NULL
val <- val[complete.cases(val),]
val[,684:688] <- round(val[,684:688],0)
val$delta.1 <- val$est.1 - val$rich.0
val$delta.3 <- val$est.3 - val$rich.0
val$delta.5 <- val$est.5 - val$rich.0
val$delta.10 <- val$est.10 - val$rich.0
val$delta.50 <- val$est.50 - val$rich.0
write_csv(val, "site_richness_validation_withsp_comp0.csv")


#### At this point, get richness estimates for all scales (hpc)

# # Evaluate species
# 
# # hotspots
# season="breeding"
# list <- read_csv("species_list.csv")
# list$Scientific <- gsub(" ","_",list$Scientific)
# val <- read_csv("site_richness_validation_withsp_comp05.csv")
# val <- cbind(val[,c(1:5)],val[,list$Scientific],val[,683:692])
# names(val)[6:677] <-  gsub(" ","_",list$Common)
# View(val)
# vals <- val[val$season==season,]
# coordinates(vals) <- c("x","y")
# crs(vals) <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# vals <- spTransform(vals, "+proj=longlat +datum=WGS84 +nodefs")
# vals <- data.frame(vals)
# sppt <- sort(names(vals)[(which(vals[120,4:675]==1)) +3])
# vals[120,c("x","y")]
# 
# # grid cells
# scale=5
# sp_est <- read_csv(paste0("patterns/sp_in_cells_",season,"_",scale,".csv"))
# est_points <- sp_est
# coordinates(est_points) <- c("x","y")
# crs(est_points) <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
# est_points <- spTransform(est_points, "+proj=longlat +datum=WGS84 +nodefs")
# est_points <- data.frame(est_points)
# wh <- which(list$Scientific %in% names(est_points))
# est_points <- cbind(est_points[,list$Scientific[wh]],est_points[,603:604])
# names(est_points)[1:597] <-  list$Common[wh]
# spr <- sort(names(est_points)[(which(est_points[120,1:597]==1))])
# est_points[120,c("x","y")]
# 
# tpres <- length(which(spr %in% sppt))
# fpres <- length(setdiff(spr, sppt))
# tabs <- length(which(spr_abs %in% sppt_abs))
# fabs <- length(setdiff(spr_abs, sppt_abs))
# pres_rate <- tpres/(tpres+fpres)
# abs_rate <- tabs/(tabs+fabs)
# 


# Compare hotspot richness to estimated richness
# and sp to sp
val <- read_csv("site_richness_validation_withsp_comp0.csv")
# cut any added points for large hotspots
val <- val[val$type=="auto",]
spnames <- names(val[7:683])
for (season in c("breeding", "nonbreeding")){
  vals <- val[val$season==season,]
  for (scale in c(1,3,5,10,50)){
    # get richness and sp list per cell
    vals$est <- vals[,paste0("est.",scale)][[1]] 
    # lists of sp observed 
    sp_est <- read_csv(paste0("patterns/sp_in_cells_",season,"_",scale,".csv"))
    est_points <- sp_est
    coordinates(est_points) <- c("x","y")
    for (i in 1:nrow(vals)){
      cell <- vals[i,"cell"]
      #find closest cell to point, get sp
      val_coords <- vals[i, c("x","y")]
      coordinates(val_coords) <- c("x","y")
      crs(val_coords) <- crs(est_points) <- 
        "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
      wm <- which.min(pointDistance(val_coords, est_points))
      spr <- sort(names(sp_est)[(which(sp_est[wm,3:ncol(sp_est)]==1)) +2])
      # hotspot species list
      sppt <- sort(names(vals)[(which(vals[i,7:683]==1)) +6])
      # same for absent species
      spr_abs <- setdiff(names(vals[7:683]), spr)
      sppt_abs <- setdiff(names(vals[7:683]), sppt)
      # get true pres, false pres, true abs, false abs
      tpres <- length(which(spr %in% sppt))
      fpres <- length(setdiff(spr, sppt))
      tabs <- length(which(spr_abs %in% sppt_abs))
      fabs <- length(setdiff(spr_abs, sppt_abs))
      pres_rate <- tpres/(tpres+fpres)
      abs_rate <- tabs/(tabs+fabs)
      fpres_rate <- fpres/(tpres+fpres)
      fabs_rate <- fabs/(tabs+fabs)
      ### species-level accuracy by location
      # first find which species are correctly or incorrectly pres/abs
      sp_corr_pres <- spr[which(spr %in% sppt)]
      sp_corr_abs <- spr_abs[which(spr_abs %in% sppt_abs)]
      sp_false_pres <- spr[which(spr %in% sppt_abs)]
      sp_false_abs <- spr_abs[which(spr_abs %in% sppt)]
      # then assign sp-level values
      loc_tpres <- which(spnames %in% sp_corr_pres)
      loc_tabs <- which(spnames %in% sp_corr_abs)
      loc_fpres <- which(spnames %in% sp_false_pres)
      loc_fabs <- which(spnames %in% sp_false_abs)
      pred <- rep(NA, length(spnames))
      pred[loc_tpres] <- "tpres"
      pred[loc_tabs] <- "tabs"
      pred[loc_fpres] <- "fpres"
      pred[loc_fabs] <- "fabs"
      # collect loc-level info
      pa <- c(cell[[1]], season, scale, tpres, fpres, tabs, fabs, 
              pres_rate, abs_rate, fpres_rate, fabs_rate)
      if(i==1 & scale==1){coll <- pa}else{coll <- rbind(coll, pa)}
      # collect sp-level info
      spl <- c(cell[[1]], season, scale, pred)
      if(i==1 & scale==1){splev <- spl}else{splev <- rbind(splev, spl)}
    }
    # species-level metrics
    splev <- data.frame(splev)
    colnames(splev) <- c("cell","season","scale",spnames)
    sub_splev <- splev[,4:ncol(splev)]
    foo = function(x){
      len = length(x)
      t_pres = sum(x == "tpres")/len
      t_abs = sum(x == "tabs")/len
      f_pres = sum(x == "fpres")/len
      f_abs = sum(x == "fabs")/len
      return(c(t_pres, t_abs, f_pres, f_abs))
    }
    sp_sum <- t(sapply(sub_splev, foo))
    sp_sum <- cbind(rownames(sp_sum),season, scale, sp_sum)
    sp_sum <- data.frame(sp_sum)
    rownames(sp_sum) <- NULL
    colnames(sp_sum) <- c("spname","season","scale","tpres","tabs","fpres","fabs")
    # collect info
    if(scale==1){sptab <- sp_sum}else{sptab <- rbind(sptab, sp_sum)}
    
    # richness scatter plot
    if(scale==1){ylab="Estimated richness\nfrom stacked SDMs"}else{ylab=NULL}
    if(season=="breeding"){xlab=NULL}else{xlab="Hotspot richness"}
    gg <-  ggplot(vals, aes(rich.0, est)) +
      geom_point(size=2) +
      geom_smooth(color="black", method="lm") +
      theme(text = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size=12),
            axis.text.x = element_text(color = "black", size=12)) +
      ylim(0,160) +
      xlim(0,250) +
      xlab(xlab) +
      ylab(ylab)
    assign(paste0("gg.",season,".",scale), gg)
    
    # positives and negatives
    if(scale==1){ylab="True negatives rate"}else{ylab=NULL}
    if(season=="breeding"){xlab=NULL}else{xlab="True positives rate"}
    coll.df <- data.frame(coll)
    colnames(coll.df) <- c("cell","season","scale","true_pres","false_pres",
                           "true_abs","false_abs","pres_rate","abs_rate",
                           "fpres_rate","fabs_rate")
    # fix columns
    coll.df <- mutate_at(coll.df, c('pres_rate', 'abs_rate',
                                    'fpres_rate', 'fabs_rate'), as.numeric)
    coll.df <- mutate_at(coll.df, 'scale', as.factor)
    coll.df$scale <- factor(coll.df$scale, levels=c('1', '3', '5', '10', '50'))
    coll1 <- coll.df[coll.df$scale==1,]
    collc <- coll.df[coll.df$scale!=1,]
    collc$delta_pres <- collc$pres_rate - coll1$pres_rate
    collc$delta_abs <- collc$abs_rate - coll1$abs_rate
    collc$delta_fpres <- collc$fpres_rate - coll1$fpres_rate
    collc$delta_fabs <- collc$fabs_rate - coll1$fabs_rate
    
    coll_sub <- coll.df[coll.df$scale==scale,] 
    pos <-  ggplot(coll_sub, aes(pres_rate, abs_rate)) +
      geom_point(size=2) +
      geom_smooth(color="black", method="lm") +
      theme(text = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size=12),
            axis.text.x = element_text(color = "black", size=12)) +
      ylim(0.9,1) +
      xlim(0,1) +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(paste0(scale,"km, ",season))
    assign(paste0("pos.",season,".",scale), pos)
    
    # plot trends during the last scale run
    if(scale==50){
      
      # delta comparisons
      rows <- nrow(vals)
      valp <- cbind(rep("points",rows),vals$rich.0)
      colnames(valp) <- c("scale","value")
      valp <- data.frame(valp)
      valp$value <- as.numeric(valp$value)
      val2 <- cbind(c(rep(1,rows),rep(3,rows),
                      rep(5,rows),rep(10,rows),rep(50,rows)),
                    c(vals$delta.1,vals$delta.3,vals$delta.5,
                      vals$delta.10,vals$delta.50))
      colnames(val2) <- c("scale","value")
      val2 <- data.frame(val2)
      val2$value <- as.numeric(val2$value)
      val2$scale <- as.factor(val2$scale)
      
      # fix columns - sp level data
      sptab <- data.frame(sptab)         
      sptab <- mutate_at(sptab, c('tpres', 'tabs', 'fpres', 'fabs'), as.numeric)   
      # rates (denominator set at minimum 1 to avoid div by 0)
      sptab$pres_rate <- sptab$tpres/ifelse((sptab$tpres+sptab$fpres)==0,
                                            1,(sptab$tpres+sptab$fpres))
      sptab$abs_rate <- sptab$tabs/ifelse((sptab$tabs+sptab$fabs)==0,
                                          1,(sptab$tabs+sptab$fabs))
      sptab$fpres_rate <- sptab$fpres/ifelse((sptab$tpres+sptab$fpres)==0,
                                             1,(sptab$tpres+sptab$fpres))
      sptab$fabs_rate <- sptab$fabs/ifelse((sptab$tabs+sptab$fabs)==0,
                                           1,(sptab$tabs+sptab$fabs))
      sptab <- mutate_at(sptab, c('pres_rate', 'abs_rate',
                                  'fpres_rate', 'fabs_rate'), as.numeric)
      sptab <- mutate_at(sptab, 'scale', as.factor)
      sptab$scale <- factor(sptab$scale, levels=c('1', '3', '5', '10', '50'))
      # divide to 1km and other scale dfs
      sp1 <- sptab[sptab$scale==1,]
      spc <- sptab[sptab$scale!=1,]
      spc$delta_pres <- spc$pres_rate - sp1$pres_rate
      spc$delta_abs <- spc$abs_rate - sp1$abs_rate
      spc$delta_fpres <- spc$fpres_rate - sp1$fpres_rate
      spc$delta_fabs <- spc$fabs_rate - sp1$fabs_rate
      
      # violin plot, delta richness between points and 5 scales
      if(season=="breeding"){col="#F8766D"}else{col="#00BFC4"}
      # if(season=="nonbreeding"){xlab="grain size (km)"}else{xlab=NULL}
      vv <- ggplot(data=valp, aes(x=scale, y=value)) + 
        geom_hline(yintercept = mean(valp$value), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Java", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('Richness') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      vv2 <- ggplot(data=val2, aes(x=scale, y=value, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " Richness"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("vv.",season), vv)
      assign(paste0("vv2.",season), vv2)
      
      # violin plot - change in percent positive or negative
      
      pv <- ggplot(data=coll1, aes(x=scale, y=pres_rate)) + 
        geom_hline(yintercept = mean(coll1$pres_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Cross", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('True positives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      pv2 <- ggplot(data=collc, aes(x=scale, y=delta_pres, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col, 4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " True positives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("pv.",season), pv)
      assign(paste0("pv2.",season), pv2)
      
      
      nv <- ggplot(data=coll1, aes(x=scale, y=abs_rate)) + 
        geom_hline(yintercept = mean(coll1$abs_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Cross", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('True negatives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      nv2 <- ggplot(data=collc, aes(x=scale, y=delta_abs, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " True negatives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("nv.",season), nv)
      assign(paste0("nv2.",season), nv2)
      
      # violin plot - change in percent positive or negative of false pres/abs
      
      fpv <- ggplot(data=coll1, aes(x=scale, y=fpres_rate)) + 
        geom_hline(yintercept = mean(coll1$fpres_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Demuth", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('False positives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      fpv2 <- ggplot(data=collc, aes(x=scale, y=delta_fpres, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " False positives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("fpv.",season), fpv)
      assign(paste0("fpv2.",season), fpv2)
      
      fnv <- ggplot(data=coll1, aes(x=scale, y=fabs_rate)) + 
        geom_hline(yintercept = mean(coll1$fabs_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Demuth", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('False negatives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      fnv2 <- ggplot(data=collc, aes(x=scale, y=delta_fabs, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " False negatives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("fnv.",season), fnv)
      assign(paste0("fnv2.",season), fnv2)
      
      
      #### sp-level validation plots
      
      # violin plot - change in percent positive or negative
      
      sppv <- ggplot(data=sp1, aes(x=scale, y=pres_rate)) + 
        geom_hline(yintercept = mean(sp1$pres_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Cross", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('True positives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      sppv2 <- ggplot(data=spc, aes(x=scale, y=delta_pres, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col, 4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " True positives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("sppv.",season), sppv)
      assign(paste0("sppv2.",season), sppv2)
      
      
      spnv <- ggplot(data=sp1, aes(x=scale, y=abs_rate)) + 
        geom_hline(yintercept = mean(sp1$abs_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Cross", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('True negatives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      spnv2 <- ggplot(data=spc, aes(x=scale, y=delta_abs, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " True negatives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("spnv.",season), spnv)
      assign(paste0("spnv2.",season), spnv2)
      
      # violin plot - change in percent positive or negative of false pres/abs
      
      spfpv <- ggplot(data=sp1, aes(x=scale, y=fpres_rate)) + 
        geom_hline(yintercept = mean(sp1$fpres_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Demuth", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('False positives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      spfpv2 <- ggplot(data=spc, aes(x=scale, y=delta_fpres, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " False positives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("spfpv.",season), spfpv)
      assign(paste0("spfpv2.",season), spfpv2)
      
      spfnv <- ggplot(data=sp1, aes(x=scale, y=fabs_rate)) + 
        geom_hline(yintercept = mean(sp1$fabs_rate), size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = met.brewer("Demuth", n = 5)[1])  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('False negatives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      spfnv2 <- ggplot(data=spc, aes(x=scale, y=delta_fabs, fill=scale)) + 
        geom_hline(yintercept = 0, size = 1.5, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = rep(col,4)) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " False negatives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = 15))
      assign(paste0("spfnv.",season), spfnv)
      assign(paste0("spfnv2.",season), spfnv2)
      
      
      
    }
    
    
    
    # stats for validation
    beta <- coef(summary(lm(est~rich.0, data=vals)))[2,1]
    r2 <- summary(lm(est~rich.0, data=vals))$r.squared
    rmse <- sqrt(mean((vals$rich.0 - vals$est)^2))
    smape <- smape(vals$rich.0, vals$est)
    row <- c(season, scale, beta, r2, rmse, smape)
    if(scale==1 & season=="breeding"){table <- row}else{table <- rbind(table, row)}
  }
  colnames(coll) <- c("cell","season","scale","true_pres","false_pres",
                      "true_abs","false_abs","pres_rate","abs_rate",
                      "fpres_rate","fabs_rate")
  write_csv(data.frame(coll), paste0("patterns/point_performance_",season,".csv"))
}
gga <- ggarrange(gg.breeding.1, gg.breeding.3, gg.breeding.5,
                 gg.breeding.10, gg.breeding.50,
                 gg.nonbreeding.1, gg.nonbreeding.3, gg.nonbreeding.5,
                 gg.nonbreeding.10, gg.nonbreeding.50,
                 nrow=2, ncol=5, widths = c(1.15,1,1,1,1))
ggsave(paste0("patterns/validation_scatter.jpeg"), gga, "jpeg",
       height=6, width=12, units="in")

posa <- ggarrange(pos.breeding.1, pos.breeding.3, pos.breeding.5,
                  pos.breeding.10, pos.breeding.50,
                  pos.nonbreeding.1, pos.nonbreeding.3, pos.nonbreeding.5,
                  pos.nonbreeding.10, pos.nonbreeding.50,
                  nrow=2, ncol=5, widths = c(1.15,1,1,1,1))
ggsave("patterns/validation_pos_neg_rate_scatter.jpeg", posa, "jpeg",
       height=6, width=12, units="in")

# vva <- ggarrange(vv.breeding, vv2.breeding, vv.nonbreeding, vv2.nonbreeding, 
#                  nrow=2, ncol=2, widths = c(1.5,3.5))
# ggsave("patterns/validation_deltas.jpeg", vva, "jpeg",
#        height=6, width=6, units="in")
# 
# colnames(table) <- c("season","scale","coef","r2","rmse","smape")
# write_csv(data.frame(table), "patterns/validation_table.csv")

pngg <- ggarrange(pv.breeding, pv2.breeding, nv.breeding, nv2.breeding, 
                  pv.nonbreeding, pv2.nonbreeding, nv.nonbreeding, nv2.nonbreeding, 
                  nrow=2, ncol=4, widths = c(1.5,3.5), 
                  labels=c("a)","","b)","","c)","","d)",""))
ggsave("patterns/validation_posneg_rate_violin2.jpeg", pngg, "jpeg",
       height=6, width=12, units="in")

fgg <- ggarrange(fpv.breeding, fpv2.breeding, fnv.breeding, fnv2.breeding, 
                 fpv.nonbreeding, fpv2.nonbreeding, fnv.nonbreeding, fnv2.nonbreeding, 
                 nrow=2, ncol=4, widths = c(1.5,3.5), 
                 labels=c("a)","","b)","","c)","","d)",""))
ggsave("patterns/validation_fposneg_rate_violin2.jpeg", fgg, "jpeg",
       height=6, width=12, units="in")

sppngg <- ggarrange(sppv.breeding, sppv2.breeding, spnv.breeding, spnv2.breeding, 
                    sppv.nonbreeding, sppv2.nonbreeding, spnv.nonbreeding, spnv2.nonbreeding, 
                  nrow=2, ncol=4, widths = c(1.5,3.5), 
                  labels=c("a)","","b)","","c)","","d)",""))
ggsave("patterns/validation_posneg_rate_violin_splevel.jpeg", sppngg, "jpeg",
       height=6, width=12, units="in")

spfgg <- ggarrange(spfpv.breeding, spfpv2.breeding, spfnv.breeding, spfnv2.breeding, 
                   spfpv.nonbreeding, spfpv2.nonbreeding, spfnv.nonbreeding, spfnv2.nonbreeding, 
                 nrow=2, ncol=4, widths = c(1.5,3.5), 
                 labels=c("a)","","b)","","c)","","d)",""))
ggsave("patterns/validation_fposneg_rate_violin_splevel.jpeg", spfgg, "jpeg",
       height=6, width=12, units="in")

colnames(table) <- c("season","scale","coef","r2","rmse","smape")
write_csv(data.frame(table), "patterns/validation_table.csv")


# map of validation sites 
vals_sp <- st_as_sf(vals, coords=c("x","y")) %>%
  st_set_crs("+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  st_transform(crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
  sf:::as_Spatial() %>%
  data.frame()
vals_sub <- vals_sp[vals_sp$cell %in% c(808, 727),]
world <- ne_countries(scale='medium', returnclass = 'sf') %>%
  st_transform(crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
map <- ggplot() +
  geom_sf(data=world) + 
  coord_sf(xlim=c(min(vals_sp$coords.x1)-50000, max(vals_sp$coords.x1)+50000),
           ylim=c(min(vals_sp$coords.x2)-50000, max(vals_sp$coords.x2)+50000), expand=F) +
  theme_bw() +
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8),
        plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
  geom_point(vals_sp, mapping=aes(x=coords.x1, y=coords.x2)) +
  geom_point(vals_sub, mapping=aes(x=coords.x1, y=coords.x2), col="red") +
  xlab("") +
  ylab("")
ggsave(paste0("validation_map100.jpeg"), map, "jpeg",
       height=6, width=6, units="in")




