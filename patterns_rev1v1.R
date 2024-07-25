# cross-scale patterns
# for season, species attributes, and ecoregions

{# fusing output CVSs of distribution models
  library(tidyverse)
  library(RColorBrewer)
  library(MetBrewer)
  library(gridExtra)
  library(sf)
 # library(rgdal)
  library(ggpubr)
  library(car)
  library(raster)
  # version
  ver <- "VS3.1"}


# load data
table <- read_csv(paste0("~/ranges/sdm_bird_outputs_",ver,".csv"))

# cut species that don't have all scales modeled
rcounts <- table[table$season=="resident",] %>%
  count(common_name)
rsp <- rcounts[rcounts$n==5,1][[1]]
bcounts <- table[table$season=="breeding",] %>%
  count(common_name)
bsp <- bcounts[bcounts$n==5,1][[1]]
nbcounts <- table[table$season=="nonbreeding",] %>%
  count(common_name)
nbsp <- nbcounts[nbcounts$n==5,1][[1]]

# new file with only completed species
b <- table[table$common_name %in% bsp & table$season=="breeding",]
nb <- table[table$common_name %in% nbsp & table$season=="nonbreeding",]
r <- table[table$common_name %in% rsp & table$season=="resident",]
full <- rbind(b, nb, r)

write_csv(full, paste0("~/ranges/sdm_bird_outputs_",ver,"_complete.csv"))


######### MODEL METRICS & AREA OF OCCURRENCE & PREDICTOR IMPORTANCES

data <- read_csv(paste0("~/ranges/sdm_bird_outputs_",ver,"_complete.csv"))

metrics <- c("AUC","SPSthreshold",
             "TSS","Sensitivity","Specificity")
# Get species-level delta values of model performance, range size between 1 and 50km scales
rm(out)
for (season in c("breeding", "nonbreeding", "resident")){
  for (i in 1:length(unique(data$sciname))){
    sci_name <- unique(data$sciname)[i]
    common <- unique(data$common_name)[i]
    sub <- data[data$sciname==sci_name & data$season==season,
                c("scale","range_area",metrics)]
    if(nrow(sub>0)){
      for (scale in c(3,5,10,50)){
        difs <- round(sub[sub$scale==scale,]-sub[1,], 4)[3:7]
        aocdif <- round(((sub[sub$scale==scale,"range_area"]/round(sub[1,"range_area"], 4) - 1) * 100),4)
        collect <- c(common, sci_name, season, scale, as.numeric(aocdif), as.numeric(difs))
        if(scale>3){spcollect <- rbind(spcollect, collect)}else{spcollect <- collect}
      }
      if(!exists("out")){out <- spcollect}else{out <- rbind(out, spcollect)}
    }
  }
}
colnames(out) <- c("common_name", "sciname", "season", "scale",
                   "delta_range_area", paste0("delta_",metrics))
rownames(out) <- NULL
out <- data.frame(out)
write_csv(out, "~/ranges/birds_vs3.1_scale_comp.csv")


# Range size estimates
out <- read_csv("~/ranges/sdm_bird_outputs_VS1_complete.csv")
  out$season <- plyr::revalue(out$season, c(
    "breeding"="summer","nonbreeding"="winter"))
out2 <- out[out$scale==1 & out$season!="resident",]
outr <- outr2 <- out[out$scale==1 & out$season=="resident",]
outr$season = "summer"
outr2$season = "winter"
out2 <- rbind(out2, outr, outr2)
out2$range_area <- log(out2$range_area)
  aoc.1 <- ggplot(out2, aes(x=range_area, fill=season)) +
    geom_density(alpha=.5) +
    geom_vline(xintercept=0, color="black", linetype="dashed", cex=1) +
    theme(legend.position = c(.23,.8),
          legend.title = element_blank(),
          legend.text=element_text(size=16),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=15),
          plot.title = element_text(size = 15, face = "bold", hjust=.5)) + 
    xlim(5,18) +
    theme(text = element_text(size=16)) +
    xlab(expression(Range ~ size ~ (log(km^2)))) +
    ylab("Density")
  
out <- read_csv("~/ranges/birds_vs3.1_scale_comp.csv")
out$season <- plyr::revalue(out$season, c(
  "breeding"="summer","nonbreeding"="winter"))
for (scale in c(3,5,10,50)){
    # maneuvering to get two sets of resident rows, one under each season
    out2 <- out[out$scale==scale & out$season!="resident",]
    outr <- outr2 <- out[out$scale==scale & out$season=="resident",]
    outr$season = "summer"
    outr2$season = "winter"
    out2 <- rbind(out2, outr, outr2)
    
    # get under/over percentages
    nrow(out2[out2$season=="summer" & out2$delta_range_area<0,])/nrow(
      out2[out2$season=="summer",])
    nrow(out2[out2$season=="winter" & out2$delta_range_area<0,])/nrow(
      out2[out2$season=="winter",])
    
    # plot (fig. 2d)
    aocgg <- ggplot(out2, aes(x=delta_range_area, fill=season)) +
      geom_density(alpha=.5) +
      geom_vline(xintercept=0, color="black", linetype="dashed", cex=1) +
      xlim(-100,300) +
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(color = "black"),
            axis.text.x = element_text(color = "black", size=15),
            plot.title = element_text(size = 15, face = "bold", hjust=.5)) + 
      theme(text = element_text(size=15)) +
      xlab(expression("% Difference in Range Size")) +
      ylab("")
    assign(paste0("aoc.",scale), aocgg)
}
gg.aoc <- ggarrange(aoc.1, aoc.3, aoc.5, aoc.10, aoc.50, nrow=1, ncol=5)
ggsave("~/ranges/scale_aoc_density_rev1.jpeg", gg.aoc, "jpeg", height=3, width=17, units="in")

# Violin plot for AUC, TSS metrics (Fig. s4)
data <- read_csv("~/ranges/sdm_bird_outputs_VS3.1_complete.csv")
 metrics <- c("AUC","SPSthreshold",
              "TSS","Sensitivity","Specificity")
for(season in c("breeding","nonbreeding")){
  if(season=="breeding"){col="#F8766D"}else{col="#00BFC4"}
  for(metric in metrics[c(1,3)]){
#    if(season=="nonbreeding" & metric=="TSS"){xlab="grain size (km)"}else{xlab=""}
    data1km <- data[data$scale==1 & (data$season==season | data$season=="resident"),]
    data1km$scale <- as.factor(data1km$scale)
    data1km$metric <- data1km[,metric][[1]]
    gg <- ggplot(data=data1km, aes(x=scale, y=metric)) + 
      geom_hline(yintercept = mean(data1km$metric), size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .75, fill=col) +
      scale_fill_manual(values = col)  +
      stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
      xlab(NULL) +
      ylab(metric) +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 15))
    
    out2 <- out[out$season==season | out$season=="resident",]
    out2$scale <- as.factor(out2$scale)
    out2$metric <- out2[,paste0("delta_",metric)][[1]]
    if(metric=="range_area"){out2$metric <- ifelse(out2$metric>0,out2$metric^(1/3),(-(-out2$metric)^(1/3)))}
    if(metric=="range_area"){out2$metric <- log(out2$metric)}
    gg2 <- ggplot(data=out2, aes(x=scale, y=metric, fill=scale)) + 
      geom_hline(yintercept = 0, size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .75) +
      scale_fill_manual(values = rep(col,4)) +
      stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
      ylab(expr(paste(Delta, " ", !!metric))) +
      xlab(NULL) +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 15))
    gg2
    assign(paste0("gg.",metric,".",season),gg)
    assign(paste0("gg2.",metric,".",season),gg2)
  }
}
gg.ppm <- ggarrange(gg.AUC.breeding, gg2.AUC.breeding, gg.TSS.breeding, gg2.TSS.breeding,
             gg.AUC.nonbreeding, gg2.AUC.nonbreeding, gg.TSS.nonbreeding, gg2.TSS.nonbreeding,
             nrow=2, ncol=4, widths = c(1.5,3.5), labels=c("a)","","b)","","c)","","d)",""))
ggsave("~/ranges/scale_ppms_violin_rev1.jpeg", gg.ppm, "jpeg", height=6, width=12, units="in")



{##### SPECIES-LEVEL VARIATION
  setwd("~/ranges/")
  
  # Species level data -combine deltas with model outputs
  output <- read_csv("sdm_bird_outputs_VS3.1_complete.csv")
  output2 <- output %>%
    dplyr::select(common_name,season,scale,noPts,range_area)
  data <- read_csv("birds_vs3.1_scale_comp.csv") %>%
    left_join(output2, by=c("common_name","season","scale"))
  write_csv(data, "birds_vs3.1_sp_level.csv")
  
  # Create habitat diversity index
  output <- output[output$scale==1,]
  ldis <- matrix(nrow=0, ncol=3)
  for (row in 1:nrow(output)){
  hab <- output[row, 46:81]
  hab <- hab/sum(hab, na.rm = T)
  hab[ hab == 0 ] <- NA
  ldi <- -1* sum(hab * log(hab), na.rm=T) / log(length(hab))
  ldis <- rbind(ldis, c(as.character(output[row, c(1,3)]), ldi))
  }
  ldis <- data.frame(ldis)
  names(ldis) <- c("common_name","season","ldi")
  data <- left_join(data, ldis, by=c("common_name","season"))
  data$ldi <- as.numeric(data$ldi)
  data$area <- sqrt(data$range_area)/1000
  data <- data[complete.cases(data),]
  write_csv(data, "birds_vs3.1_sp_level_corrected_with_ldi.csv")
  
  # restrict comparison to 1 vs 5 km models
  data <- read_csv("birds_vs3.1_sp_level_corrected_with_ldi.csv")
  comp <- data[data$scale==5,]
  
    # Model
    comp.st <- comp
    comp.st$ldi <- scale(comp.st$ldi)[,1]
    comp.st$area <- scale(sqrt(comp.st$area))[,1]
    glm <- glm(delta_AUC ~ ldi + area + season, data=comp.st)
    table <- data.frame(summary(glm)$coef)
    write.csv(table, file=paste0("glm_table_v2.csv")) # table S4
}

# plotting (fig. 6)
  gg1 = ggplot(comp.st, aes(x=area, y=delta_AUC)) +
    geom_point() +
    geom_smooth(method="lm") +
    theme(text = element_text(size=15),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color="gray70"),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=12),
          axis.text.y = element_text(color = "black", size=12)) +
    xlab("Estimated Range Size (sqrt(km^2)) \n(standardized)") +
    ylab(expr(paste(Delta, " AUC"))) +
    ggtitle(NULL)
  gg2 = ggplot(comp.st, aes(x=ldi, y=delta_AUC)) +
    geom_point() +
    geom_smooth(method="lm") +
    theme(text = element_text(size=15),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color="gray70"),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=12),
          axis.text.y = element_text(color = "black", size=12)) +
    xlab("Landcover Diversity Index \n(standardized)") +
    ylab(expr(paste(Delta, " AUC"))) +
    ggtitle(NULL)
  gga <- ggarrange(gg1, gg2, nrow=1, ncol=2, labels=c("a)","b)"))
  ggsave(paste0("sp_attributes.jpeg"), gga, "jpeg",height=5, width=9, units="in")


##### Site-level ecoregion vs validation performance

# connect coords, ecoregion to sites
setwd("~/ranges/")
b = read_csv(paste0("point_performance_breeding_rev1.csv"))
nb = read_csv(paste0("point_performance_nonbreeding_rev1.csv"))
data = rbind(b, nb)
out <- read_csv("site_richness_validation_withsp0.csv")
data = base::merge(data, out[,1:4], by=c("cell", "season"))

eco <- st_read("shapefiles", "Ecoregions2017")
eco <- eco[eco$REALM=="Nearctic" | eco$REALM=="Neotropic",]
eco <- st_transform(eco, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

data <- st_as_sf(data, coords = c("x", "y"), crs = "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
data <- st_transform(data, crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
plot(eco, max.plot=1)
plot(data, add=T)

sf_use_s2(FALSE)
data.int <- data.frame(st_intersects(data, eco, sparse=F))
ers = c()
for (i in 1:nrow(data)){
  ww = which(data.int[i,]==T) ## NEXT STEP GET ECOREGION OF THIS POLYGON NUMBER
  er = eco[ww,]$BIOME_NAME
  print(er)
  er = ifelse(identical(er, character(0)) == T, NA, er)
  ers = c(ers, er)
}
data <- data.frame(sf:::as_Spatial(data))
data <- rename(data, "longitude"="coords.x1", "latitude"="coords.x2")
data$optional = NULL
data$biome = ers

b = data[data$season=="breeding",]
write_csv(b, "patterns/point_performance_breeding_coords.csv")
nb = data[data$season=="nonbreeding",]
write_csv(nb, "patterns/point_performance_nonbreeding_coords.csv")

# models/plots (fig. 5)
for(season in c("breeding","nonbreeding")){
  if(season=="breeding"){col="#F8766D"}else{col="#00BFC4"}
  comp = read_csv(paste0("point_performance_",season,"_coords.csv"))
  comp1 = comp[comp$scale==1,]
  comp5 = comp[comp$scale==5,]
  comp1$delta_pres_rate = comp1$pres_rate - comp5$pres_rate
  comp1 = comp1[comp1$biome!="Tropical & Subtropical Grasslands, Savannas & Shrublands" &
                   comp1$biome!="Mediterranean Forests, Woodlands & Scrub" &
                   comp1$biome!="Tropical & Subtropical Coniferous Forests" &
                  comp1$biome!="Tundra",]
  comp1 = comp1[complete.cases(comp1$biome),]
  comp1$biome <- plyr::revalue(comp1$biome, c("Temperate Conifer Forests"="Temperate Broadleaf & Mixed Forests",
                                              "Flooded Grasslands & Savannas"="Temperate Grasslands, Savannas & Shrublands"))
  write_csv(comp1, paste0("point_performance_",season,"_deltas.csv"))
bio <- ggplot(data=comp1, aes(x=biome, y=delta_pres_rate)) + 
  geom_hline(yintercept = 0, linewidth = 1.5, color = "gray") +
  geom_violin(trim=F, alpha = .75, fill=col) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
  ylab(expr(paste(Delta, " True positives rate"))) +
  xlab("") +
  theme_light() +
  theme(legend.position = "none",
        text = element_text(size = 16),
        axis.text.x=element_blank())
assign(paste0("bio.",season), bio)
}

gga <- ggarrange(bio.breeding, bio.nonbreeding,
                 nrow=2, ncol=1)
ggsave(paste0("sp_level_biomes_bysite_rev1.jpeg"), gga, "jpeg",
       height=7, width=6, units="in")


b = read_csv(paste0("point_performance_breeding_deltas.csv"))
nb = read_csv(paste0("point_performance_nonbreeding_deltas.csv"))
data = rbind(b, nb)

data$biome = factor(data$biome)
data$season = factor(data$season)
# glm
glm <- glm(delta_pres_rate ~ biome + season, data=data)
table <- data.frame(summary(glm)$coef)
write.csv(table, file=paste0("glm_table_biome_rev1.csv"))
# anova (table s3)
anova <- aov(delta_pres_rate ~ biome + season, data=data)
table <- summary(anova)$
write.csv(table, file=paste0("anova_table_biome_rev1.csv"))
# tukeys (table s3)
tukey = TukeyHSD(anova, conf.level=.95)
tukey$biome
write.csv(tukey$biome, file=paste0("tukey_table_biome_rev1.csv"))


###### Histogram of sample sizes by season (fig. s2)
setwd("~/ranges/")
data <- read_csv(paste0("sdm_bird_outputs_",ver,"_complete.csv"))
data$season <- plyr::revalue(data$season, c(
  "breeding"="summer","nonbreeding"="winter"))
data2 <- data[data$scale==1 & data$season!="resident",]
datar <- datar2 <- data[data$scale==1 & data$season=="resident",]
datar$season = "summer"
datar2$season = "winter"
data2 <- rbind(data2, datar, datar2)
gghist <- ggplot(data2, aes(x=noPts, fill=season)) +
  geom_density(alpha=.5) +
  geom_vline(xintercept=0, color="black", linetype="dashed", cex=1) +
  theme(legend.position = c(.7,.8),
        legend.title = element_blank(),
        legend.text=element_text(size=16),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=15),
        plot.title = element_text(size = 15, face = "bold", hjust=.5)) + 
  theme(text = element_text(size=16)) +
  xlab("Number of observation points") +
  ylab("Density")
ggsave(paste0("patterns/sp_level_sample_size.jpeg"), gghist, "jpeg",
       height=5, width=5, units="in")
# min observations by season
min(data2[data2$season=="summer",]$noPts)
min(data2[data2$season=="winter",]$noPts)
