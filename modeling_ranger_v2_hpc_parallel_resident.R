# Species Distribution Models

{library(maps)
  library(rgdal)
  library(raster)
  library(dismo)
  library(data.table)	
  library(scales)
  library(ranger)
  library(scam)
  library(fields)  
  library(PresenceAbsence) 
  library(pdp)
  library(plyr)
  library(tidyverse)
  library(gbm)
  library(sf)
  library(rgeos)
  library(SDMtune)
  library(ff)
  library(ffbase)
  library(foreach)
  library(parallel)
  library(doParallel)
  
  sessionInfo()
  capabilities() 
}

###############################
#MODELS (RF)
###############################
{ version <- "V2.1"
set.seed(1)
setwd("/gpfs/ysm/home/jc3893/")

# Species names
list <- read_csv("30x30/species_list.csv")

# shapefile
countries <- readOGR("30x30/shapefiles", "countries")
countries_proj <- spTransform(countries, CRS="+init=epsg:26978")
aba <- countries[c(9,38),] %>%
  st_as_sf()

{# grid for thinning
  bb <- bbox(countries)
  cs <- c(.05,.05)
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  grd
  sp_grd <- SpatialGridDataFrame(grd,
                                 data=data.frame(id=1:prod(cd)),
                                 proj4string=CRS(proj4string(countries)))
}

# define season - change as needed
#season="breeding"; months=6:8; mindate=153; maxdate=243; scode=2
#season <- "nonbreeding"; months <- c(12,1,2); mindate=335; maxdate=59; scode=3
#season <- "prebreeding_migration"; months <- 3:5; mindate=60; maxdate=151
#season <- "postbreeding_migration"; months <- 9:11; mindate=244; maxdate=334
season="resident"; months=1:12; mindate=1; maxdate=365; scode=1

# ebird data
for (s in c("breeding","nonbreeding",
           "prebreeding_migration","postbreeding_migration")){
ebird_ff <- read.csv.ffdf(file=paste0("30x30/ebird_",s,"_new.csv"),
                            colClasses=c(rep("factor", 2), rep("numeric",2),
                                         rep("integer",2), "numeric", "integer",
                                         "numeric","integer",rep("integer",677),
                                         rep("numeric",44)))
assign(paste0("ebird_ff_",s), ebird_ff)
}

# prediction surface
pred_surface_full <- read.csv.ffdf(file="30x30/habitat_prediction-surface-lc.csv",
                                   colClasses=rep("numeric",46))

# standard grid
bbox <- read_sf("30x30/shapefiles", "30x30BBox_NAmerica")
grid <- raster("30x30/grids/Global_reference_raster_1km_WGS84.tif") %>%
  crop(bbox)

# Predictors (with weather data)
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
# For removing landcover
nolc_preds <- all_preds[1:13]

# mol, birdlife ranges
mol.range <- readOGR("/gpfs/ysm/home/jc3893/30x30/range_polygons", "jetz_ranges_sub")
bl.range <- readOGR("/gpfs/ysm/home/jc3893/30x30/range_polygons", "birdlife")

# directories
dir <- paste0("/gpfs/loomis/pi/jetz/data/results_30x30/sdm_birds/",version,"_30x30/")
dir.create(dir)
dir.create(paste0(dir, "plots/"))
setwd(dir)

# # parallel computing
# #n.cores <- parallel::detectCores()
# n.cores <- 4
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "PSOCK"
# )
# doParallel::registerDoParallel(cl = my.cluster)

}

# ----------------------------------------------------------------------
# RANGER
# ----------------------------------------------------------------------		
# find a species
# which(list[[1]]=="Botteri's Sparrow")
{# loop species
   #series <- 1:nrow(list)
  series <- c(190,195,199,223,224,225,245,246,421,522,523,524,525,526,548)
  #foreach (i=1:nrow(list)) %do% { 
   # foreach (i=series) %dopar% { 
   for(i in series){
    try({
      # {library(maps)
      #   library(rgdal)
      #   library(raster)
      #   library(dismo)
      #   library(data.table)	
      #   library(scales)
      #   library(ranger)
      #   library(scam)
      #   library(fields)  
      #   library(PresenceAbsence) 
      #   library(pdp)
      #   library(tidyverse)
      #   library(gbm)
      #   library(sf)
      #   library(rgeos)
      #   library(plyr)
      #   library(SDMtune)
      #   library(ff)
      #   library(ffbase)
      #   library(foreach)
      #   library(parallel)
      #   library(doParallel)
      # }
      # cleanup
      rm()
      #removeTmpFiles(0)
      gc()
      # Establish names
      common <- as.character(list[i,1])
      sci_name <- as.character(list[i,2])
      sci.name <- gsub("_"," ",sci_name)
      code <- as.character(list[i,5])  
      print(common)
      # check if species has been done
      if(!file.exists(paste0(sci_name,"_complete_",season,".csv"))){ 
         if(!file.exists(paste0(sci_name,"_complete_breeding.csv"))){
         if(!file.exists(paste0(sci_name,"_complete_nonbreeding.csv"))){
        rm(sprange)
        # Expert range = modeling domain
        try({
          total.range <- readOGR(unzip(
            paste0("/gpfs/ysm/home/jc3893/30x30/range_polygons/",
                   code,"-range-2020.gpkg.zip"),
            paste0(code,"-range-mr-2020.gpkg")), "range")
        # Get seasonal range of species
        sprange <- total.range[total.range$season_name==season,]
        })
        if(!exists("sprange")){
          try({
            total.range <- readOGR(unzip(
              paste0("/gpfs/ysm/home/jc3893/30x30/range_polygons/",
                     code,"-range-2021.gpkg.zip"),
              paste0(code,"-range-mr-2021.gpkg")), "range")
            # Get seasonal range of species
            sprange <- total.range[total.range$season==season,]
          })
        }
           if(!exists("sprange")){
             sprange <- mol.range[mol.range$sciname==sci.name & 
                                    mol.range$season==scode,]
             sprange$season <- season
           }
           if(!exists("sprange")){ 
             sprange <- bl.range[bl.range$sci_name==sci.name,]
             sprange$season <- season
             
           }
        if(nrow(sprange)==0){print("no range exists"); STOP}
        # Unproj, box range
        sprange <- spTransform(sprange, CRS="+proj=longlat +datum=WGS84 +no_defs") 
        if(bbox(sprange)[4]<78){
          # Proj range, buffer, back to unproj
        sprange <- spTransform(sprange, CRS="+init=epsg:26978")
          # if(area(sprange)>200000){
          #   sprange <- buffer(sprange, -5000) # remove small bits
          # }
        sprange <- buffer(sprange, 200000) %>%
          spTransform(CRS="+proj=longlat +datum=WGS84 +no_defs")
        }
        box <- bbox(sprange)
          #plot(sprange, col="bisque")
          # Back to SP dataframe
          range.df <- data.frame(ID=1, row.names="buffer")
          sprange <- SpatialPolygonsDataFrame(sprange, range.df, match.ID=F)
          # ST range
          rangest <- st_as_sf(sprange)
          # check that range is in our bbox
          sf::sf_use_s2(FALSE)
          intersect <- st_intersects(aba, rangest, sparse=T)
          if(sum(lengths(intersect))>0){
          print(paste("START", i, common, season))
          # directories
          dir.create(paste0(sci_name,"/"))
          dir.create(paste0(sci_name,"/predictions/"))
          dir.create(paste0(sci_name,"/predictions/ROR/"))
          dir.create(paste0(sci_name,"/predictions/PA/"))
          dir.create(paste0(sci_name,"/points/"))
          dir.create(paste0(sci_name,"/accuracy/"))
          dir.create(paste0(sci_name,"/rangemap/"))

          # Save range
          writeOGR(sprange, paste0(sci_name,"/rangemap"),
                   paste0(sci_name,"_",season), driver="ESRI Shapefile",
                   overwrite_layer=T)
          
          # ebird - get cropped rows, columns needed
          for (s in c("breeding","nonbreeding",
                          "prebreeding_migration","postbreeding_migration")){
          ebird_full <- get(paste0("ebird_ff_",s))
          ebird_index <- ffwhich(ebird_full, 
                                   latitude<box[2,2] &
                                     latitude>box[2,1] &
                                     longitude<box[1,2] &
                                     longitude>box[1,1])
          ebird_full <- data.frame(ebird_full[ebird_index,
                                         c(gsub(" ",".",sci.name),all_preds,
                                           "latitude","longitude","checklist_id")])
          # NA removal (only a few rows)
          ebird_full <- ebird_full[complete.cases(ebird_full),]
          
          # Cut ebird to range
          print("processing ebird data")
          ebird_full <- st_as_sf(ebird_full, coords = c("longitude", "latitude"), 
                            crs = crs(rangest))
          ebird.int <- st_intersects(ebird_full, rangest, sparse=F)
          ebird_full <- ebird_full[as.vector(ebird.int),]
          #plot(ebird[sample.int(nrow(ebird), nrow(ebird)/20),], max.plot=1)
          ebird_full <- data.frame(sf:::as_Spatial(ebird_full))
          ebird_full$optional <- NULL
          names(ebird_full)[(ncol(ebird_full)-1):ncol(ebird_full)] <- 
            c("longitude", "latitude")
          rm(ebird.int)
          
          # Impose data size limit
          if(nrow(ebird_full) > 1000000){
          ebird_full <- ebird_full[sample(nrow(ebird_full),1000000),]
          }
          
          # add to other seasons
          if(s=="breeding"){ebird <- ebird_full}else{
            ebird <- rbind(ebird, ebird_full)
          }
          }
          
          # randomize rows
          ebird <- ebird[sample(nrow(ebird)),]
          
          # Column for presence/absence
          ebird$species_observed <- ebird[,gsub(" ",".",sci.name)]
    
          # Cut pred surface to range
          print("processing prediction surface")
          surface_index <- ffwhich(pred_surface_full, 
                                   latitude<box[2,2] &
                                     latitude>box[2,1] &
                                     longitude<box[1,2] &
                                     longitude>box[1,1])
          pred_surface <- data.frame(pred_surface_full[surface_index,])
          pred_surface <- st_as_sf(pred_surface, coords = c("longitude", "latitude"), 
                                   crs = crs(sprange))
          surface.int <- st_intersects(pred_surface, rangest, sparse=F)
          pred_surface <- pred_surface[as.vector(surface.int),]
          pred_surface <- data.frame(sf:::as_Spatial(pred_surface))
          pred_surface$optional <- NULL
          names(pred_surface)[(ncol(pred_surface)-1):ncol(pred_surface)] <- 
            c("longitude", "latitude")
          pred_surface$twi <- NULL
          pred_surface <- pred_surface[complete.cases(pred_surface),]
          rm(surface.int)
          gc()
          # Populate with values
          pred_surface$duration_minutes <- 60
          pred_surface$effort_distance_km <- 1
          pred_surface$number_observers <- 1
          pred_surface$year <- 2020
          # Find max hr of day for spatial predictions
          hrs <- ebird[,c("time_observations_started","species_observed")] %>%
            group_by(round(time_observations_started)) %>%
            summarize_at(vars(species_observed), mean)
          pred_surface$time_observations_started <- as.numeric(
            hrs[which.max(hrs$species_observed),1])
          # set day for spatial predictions
          if(season=="nonbreeding"){
            pred_surface$day_of_year <- sample(
              c(1:maxdate,mindate:365), nrow(pred_surface), replace=T)
          }else{
            pred_surface$day_of_year <- sample(
              mindate:maxdate, nrow(pred_surface), replace=T)
          }
          
          
          # Splits
          print("splitting and filtering data")
          {#training and testing samples
            index <- runif(nrow(ebird))
            train <- ebird[which(index>.4),]
            test <- ebird[which(index<=.4 & index>.2),]
            train.oob <- ebird[which(index<=.2),]
            rm(index)
          }
          
          {# Spatiotemporal thinning
            train$id <- NULL
            coordinates(train) <- c("longitude", "latitude")
            crs(train) <- crs(sp_grd)
            over <- over(train, sp_grd)
            train <- cbind(data.frame(train), over)
            # get cell id and week number for each checklist
            train$optional <- NULL
            checklist_cell <- train %>% 
              mutate(cell = id,
                     week = ceiling(day_of_year/7))
            # sample one checklist per grid cell per week
            # sample detection/non-detection independently 
            train <- checklist_cell %>% 
              group_by(species_observed, year, week, cell) %>% 
              sample_n(size = 1) %>% 
              ungroup()
          }
          {# Assemble TrainOOB
            train.oob$id <- NULL
            coordinates(train.oob) <- c("longitude","latitude")
            crs(train.oob) <- crs(sp_grd)
            over <- over(train.oob, sp_grd)
            train.oob <- cbind(data.frame(train.oob), over)
            train$optional <- NULL
            checklist_cell_oob <- train.oob %>% 
              mutate(cell = id,
                     week = ceiling(day_of_year/7))
            train.oob <- checklist_cell_oob %>% 
              group_by(species_observed, year, week, cell) %>% 
              sample_n(size = 1) %>% 
              ungroup()
          }
          
          # check there is enough data
          if(nrow(train[train[,"species_observed"]>0,])<10){print("too few rows"); STOP}
          
          # remove oversampling for occurrence model, limit size
          train <- train[!duplicated(train$checklist_id), ]
          # Binary Occurrence Response, coded as numeric
          train$pres_abs <- as.numeric(train[[which(names(train) == "species_observed")]] > 0)
          # balanced sampling
          pos_fraction <- mean(as.numeric(train$pres_abs))
          # Check that positives ARE the minority class
          if (pos_fraction > 0.5) pos_fraction <- 1 - pos_fraction
          # Ranger binary response model requires code response as factor
          train$pres_abs <- as.factor(as.numeric(train$pres_abs))
          
          # Save points alone
          write_csv(train[,c("latitude","longitude","pres_abs")], 
                    paste0(sci_name,"/points/pts_binary_occurrences_",
                           sci_name,"_",season,"_00.csv"))
          
          # Loop with or without landcover
          for (lc in c("landcover","no_landcover")){
            print(paste("model", lc))
            if(lc=="landcover"){pred_names=all_preds}else{
              pred_names=nolc_preds}
            
            {#MODEL
              # Model Formula, threads
              m_formula <- paste("pres_abs ~", paste(pred_names, collapse = "+"))
              n_threads <- 7
              # Balanced Ranger Occurrence Model
              rf_occ <- ranger::ranger(
                formula =  m_formula, 
                num.trees = 100, 
                #max.depth = 0, 
                importance = "impurity",
                num.threads = n_threads,
                respect.unordered.factors = "order",
                always.split.variables = NULL,
                probability = TRUE,
                replace = TRUE, 
                sample.fraction = c(pos_fraction, pos_fraction),
                data = train)
            }
            
            # Test Set PPMS
            test$occ_pred <- predict(
              rf_occ, data = test, 
              type = "response",
              num.threads = n_threads)$predictions[,2]
            
            # # Calibration Plots
            # calibration_model <- scam(obs ~ s(occ_pred, k=5, bs="mpi"),
            #                           data=test, family=binomial)
            # seq <- seq(0, 1, by=0.01)
            # pred_calibration <- predict(calibration_model, data.frame(occ_pred = seq),
            #                             type="response")
            # jpeg(paste0("~/ranges/calibration/calibration_",
            #             common,"_",season,'.jpeg'))
            # plot(test$occ_pred[1:10000], jitter(test$obs[1:10000]), ylim = c(-.25, 1.25),
            #      pch = 20, cex = 0.2)
            # lines(seq, pred_calibration, col=2)
            # abline(0,1, col="blue");abline(h=0, col="blue");abline(h=1, col="blue")
            # dev.off()
            
            {#Select Optimal Threshold on train OOB data
              trainOOB.pred <- predict(rf_occ, data=train.oob, type="response")
              pa.df <- data.frame("space.holder.tag",
                                  obs = train.oob$species_observed, 
                                  ppp = trainOOB.pred$predictions[,2] )
              pa.metrics <- presence.absence.accuracy(
                pa.df,
                threshold = quantile(
                  pa.df$ppp,
                  probs = seq(from=0, to=1, length=1000), na.rm =T),
                na.rm = T, st.dev = F)
              pa.metrics$PCC <- pa.metrics$Kappa <- NULL
              pa.metrics$TSS <- (pa.metrics$sensitivity + pa.metrics$specificity) - 1
              
              # Save confusion matrix
              write_csv(pa.metrics, 
                        paste0(sci_name,"/accuracy/accuracy_matrix_",
                               sci_name,"_",season,"_",lc,"_00.csv"))
              
              # optimal_thresh_position <- which.max(pa.metrics$Kappa) 
              optimal_thresh_position <- which.max(pa.metrics$sensitivity +
                                                     pa.metrics$specificity)
            }
            {# Compute Test set PPMs
              test.pred <- predict(rf_occ, data=test, type="response")
              pa.df <- data.frame("space.holder.tag",
                                  obs = test$species_observed, 
                                  ppp = test.pred$predictions[,2] )
              pa.metrics <- presence.absence.accuracy(
                pa.df,
                threshold = pa.metrics$threshold[ optimal_thresh_position ] ,
                na.rm = T, st.dev = F)
              thresh <- pa.metrics$threshold
              pa.metrics["PCC"] <- pa.metrics["Kappa"] <- NULL
              pa.metrics["TSS"] <- (pa.metrics$sensitivity + pa.metrics$specificity) - 1
            }
            
            # Spatial predictions
            pred_surface$occ_pred <- predict(
              rf_occ, data = pred_surface, type = "response",
              num.threads = n_threads)$predictions[,2]
            pred_surface$occ_thresh <- ifelse(
              pred_surface$occ_pred>pa.metrics$threshold,1,0)
            
            # Save rasterized predictions
            coordinates(pred_surface) <- c("longitude","latitude")
            ras <- rasterize(pred_surface, grid, pred_surface$occ_pred, fun='first')
            writeRaster(ras, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                    sci_name,"_",season,"_",lc,"_00.tif"),
                        overwrite=T)
            pred_ones <- pred_surface[pred_surface$occ_thresh==1,]
            ras <- rasterize(pred_ones, grid, pred_ones$occ_thresh, fun='first')
            writeRaster(ras, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                    sci_name,"_",season,"_",lc,"_00.tif"),
                        overwrite=T)
            pred_surface <- data.frame(pred_surface)
            pred_ones <- data.frame(pred_ones)
            
            occ.probs <- pred_surface[,c("longitude","latitude","occ_pred")]
            write_csv(occ.probs, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                        sci_name,"_",season,"_",lc,"_00.csv"))
            occ.probs <- pred_ones[,c("longitude","latitude","occ_thresh")]
            write_csv(occ.probs, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                        sci_name,"_",season,"_",lc,"_00.csv"))
            rm(pred_ones)
            
            # Predictions projected
            if(lc=="landcover"){
              pred_surface_proj <- pred_surface
              coordinates(pred_surface_proj) <- c("longitude","latitude")
              crs(pred_surface_proj) <- crs(countries)
              pred_surface_proj <- spTransform(pred_surface_proj, CRS="+init=epsg:26978")
              pred_surface_proj <- data.frame(pred_surface_proj)
            }
            pred_surface_proj$occ_pred <- predict(
              rf_occ, data = pred_surface_proj, type = "response",
              num.threads = n_threads)$predictions[,2]
            pred_surface_proj$occ_thresh <- ifelse(
              pred_surface_proj$occ_pred>pa.metrics$threshold,1,0)

            {# Plot maps
              jpeg(paste0("plots/PO_spatial_predictions_",
                          sci_name,"_",season,"_",lc,"_00.jpeg"),
                   width=500, height=500, units="px")
              # Plot prob of occurrence
              plotres <- 500
              quilt.plot(
                x = pred_surface_proj$longitude,
                y = pred_surface_proj$latitude,
                z = pred_surface_proj$occ_pred,
                nrow = plotres,
                ncol = plotres,
                na.rm = T,
                main = paste(common,"-",sci.name,"\n",season,"-",lc))
              plot(countries_proj, col="transparent", border="black",
                   usePolypath=F, add=T)
              dev.off()
              # Plot thresholded occurrence
              jpeg(paste0("plots/Thresholded_spatial_predictions_",
                          sci_name,"_",season,"_",lc,"_00.jpeg"),
                   width=500, height=500, units="px")
              quilt.plot(
                x = pred_surface_proj$longitude,
                y = pred_surface_proj$latitude,
                z = pred_surface_proj$occ_thresh,
                nrow = plotres,
                ncol = plotres,
                na.rm = T,
                add.legend = F,
                col = c("transparent", "darkgreen"),
                main = paste(common,"-",sci.name,"\n",season,"-",lc))
              plot(countries_proj, col="transparent", border="black",
                   usePolypath=F, add=T)
              dev.off()
              }
            
            # Species paths
            path <- paste0("/mnt/Work/Yale/MOL/results_30x30/sdm_birds/",
                           version,"_30x30/",sci_name)
            modName <- paste0(sci_name,"_",season,"_",lc)
            modURL <- modPDF <- bad_ok_good <- comments <- 
              rangeOffset <- elevOffset <- deltaAUC <- NA
            modPathROR <- paste0(path,"/predictions/ROR/PO_spatial_predictions_",
                                 sci_name,"_",season,"_",lc,"_00.jpeg")
            modPathPA <- paste0(path,"/predictions/PA/Thresholded_spatial_predictions_",
                                sci_name,"_",season,"_",lc,"_00.jpeg")
            rangePath <- paste0(path,"/rangemap/",sci_name,"_",season,".shp")
            ptsPath <- paste0(path,"/points/pts_occurrences_",season,".csv")
            ptsBgPath <- NA
            confMatPath <- paste0(path,"/accuracy/accuracy_matrix_",
                                  sci_name,"_",season,"_",lc,"_00.csv")
            POpredsPath <- paste0(path,"/predictions/ROR/PO_predictions_",
                                  sci_name,"_",season,"_",lc,"_00.csv")
            ThreshpredsPath <- paste0(path,"/predictions/PA/Thresholded_predictions_",
                                      sci_name,"_",season,"_",lc,"_00.csv")
            POpredsrasPath <- paste0(path,"/predictions/ROR/PO_predictions_",
                                     sci_name,"_",season,"_",lc,"_00.tif")
            ThreshpredsrasPath <- paste0(path,"/predictions/PA/Thresholded_predictions_",
                                         sci_name,"_",season,"_",lc,"_00.tif")
            envVars <- RDataPath <- rorOrigPath <- paOrigPath <- rangeOrigPath <-
              ptsBgOrigPath <- statsOrigPath <- NA
            
            # Species values
            if(lc=="landcover"){imp.scores <- rf_occ$variable.importance}else{
              imp.scores <- c(rf_occ$variable.importance,rep(NA,36))}
            model.output <- unlist(c(common, sci_name, season, lc, version, 
                                     modName, modURL, modPDF, bad_ok_good, 
                                     comments, rangeOffset, elevOffset,
                                     nrow(train), nrow(occ.probs), pa.metrics[5],
                                     deltaAUC,pa.metrics[c(2,6,3,4)],
                                     modPathROR,modPathPA,rangePath,ptsPath,
                                     ptsBgPath,confMatPath,POpredsPath,ThreshpredsPath,
                                     POpredsrasPath,ThreshpredsrasPath,
                                     envVars,imp.scores,RDataPath,rorOrigPath,paOrigPath,
                                     rangeOrigPath,ptsBgOrigPath,statsOrigPath))
            names(model.output)<- c("common_name","sciname", 
                                    "season","landcover","version","modName",
                                    "modURL","modPDF","bad_ok_good","comments",
                                    "rangeOffset","elevOffset","noPts","range_area",
                                    "AUC","deltaAUC","SPSthreshold","TSS",
                                    "Sensitivity","Specificity",
                                    "modPathROR","modPathPA",
                                    "rangePath","ptsPath","ptsBgPath","confMatPath",
                                    "POPredsPath","threshPredsPath","POPredsRasPath",
                                    "threshPredsRasPath","envVars",
                                    paste0(all_preds,"_importance"),
                                    "RDataPath","rorOrigPath","paOrigPath",	
                                    "rangeOrigPath","ptsBgOrigPath","statsOrigPath")
            if(lc=="landcover"){sp.output <- model.output
            }else{sp.output <- data.frame(rbind(sp.output, model.output))}
          } #end landcover inclusion loop
          # get delta AUC
          deltaAUC <- round(max(as.numeric(as.character(sp.output$AUC))) - 
                              as.numeric(as.character(sp.output$AUC)), 4)
          sp.output$deltaAUC <- deltaAUC
          write_csv(data.frame(sp.output), 
                    paste0(sci_name,"_complete_",season,".csv"))
        } # end range exists in bbox check
      }}} #end species completion check
    }) # end try
  }#end sp
} #end models



