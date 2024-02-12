# cross-scale patterns across species

{# fusing output csvs of SDMs
  library(tidyverse)
  library(RColorBrewer)
  library(MetBrewer)
  library(gridExtra)
  library(sf)
  library(rgdal)
  library(ggpubr)
  library(car)
  library(raster)
  
  ver <- "VS1"}

# files <- list.files(paste0("~/ranges/",
#                            ver,"_outputs/"), full.names=T,  recursive=F)
# #files <- files[- grep(".gpkg", files)]
# for (i in 1:length(files)){
#   file <- read_csv(files[i])
#   if(i==1){table <- file}else{table <- rbind(table, file)}
# }
# # cut redundant rows
# table <- table[!duplicated(table[c("common_name","scale","season","water_cover")]),]
# # save
# write_csv(table, paste0("~/ranges/sdm_bird_outputs_",ver,".csv"))
# 
# # load data
# table <- read_csv(paste0("~/ranges/sdm_bird_outputs_",ver,".csv"))
# 
# # cut species that don't have all scales
# rcounts <- table[table$water_cover=="cont" &
#                    table$season=="resident",] %>%
#   count(common_name)
# rsp <- rcounts[rcounts$n==5,1][[1]]
# bcounts <- table[table$water_cover=="cont" &
#                    table$season=="breeding",] %>%
#   count(common_name)
# bsp <- bcounts[bcounts$n==5,1][[1]]
# nbcounts <- table[table$water_cover=="cont" &
#                     table$season=="nonbreeding",] %>%
#   count(common_name)
# nbsp <- nbcounts[nbcounts$n==5,1][[1]]
# 
# # new file with only completed species
# b <- table[table$common_name %in% bsp & table$season=="breeding",]
# nb <- table[table$common_name %in% nbsp & table$season=="nonbreeding",]
# r <- table[table$common_name %in% rsp & table$season=="resident",]
# full <- rbind(b, nb, r)
# 
# # data for lukas
# lukas <- full[full$status=="waterbird",]
# write_csv(lukas, paste0("~/ranges/sdm_bird_outputs_",ver,"_lukas.csv"))
# 
# # functional traits for lukas
# # water affinity added in excel
# lukas <- read_csv(paste0("~/ranges/sdm_bird_outputs_",ver,"_lukas.csv"))
# traits <- read_csv("~/stoat/functional/sp_matrix_traits.csv") %>% distinct() %>%
#   dplyr::select(common, BodyMass.Value, diet, distance)
# lukas <- left_join(lukas, traits, by=c("common_name"="common"))
# # also remove duplicates
# lukas <- lukas[!((lukas$common_name=="Black Guillemot" |
#                     lukas$common_name=="Northern Fulmar" |
#                     lukas$common_name=="Black-legged Kittiwake") &
#                    lukas$season=="resident"),]
# write_csv(lukas, paste0("~/ranges/sdm_bird_outputs_",ver,"_lukas.csv"))
# 
# # remove binary water cover models
# full <- full[full$water_cover=="cont",]
# #remove duplicates for 3 species with all 3 seasons
# full <- full[!((full$common_name=="Black Guillemot" |
#                   full$common_name=="Northern Fulmar" |
#                   full$common_name=="Black-legged Kittiwake") &
#                  full$season=="resident"),]
# write_csv(full, paste0("~/ranges/sdm_bird_outputs_",ver,"_complete.csv"))


######### MODEL METRICS & AREA OF OCCURRENCE & PREDICTOR IMPORTANCES

data <- read_csv(paste0("~/ranges/sdm_bird_outputs_",ver,"_complete.csv"))
# % importance of clim, top, lc suites
data$clim <- rowSums(data[,40:44])/rowSums(data[,40:83])
data$top <- rowSums(data[,45:47])/rowSums(data[,40:83])
data$lc <- rowSums(data[,48:83])/rowSums(data[,40:83])

metrics <- c("AUC","SPSthreshold",
             "TSS","Sensitivity","Specificity")
# Get change between 1 and 50km scales
rm(out)
for (season in c("breeding", "nonbreeding", "resident")){
  for (i in 1:length(unique(data$sciname))){
    sci_name <- unique(data$sciname)[i]
    common <- unique(data$common_name)[i]
    sub <- data[data$sciname==sci_name & data$season==season,
                c("scale","range_area",metrics,"clim","top","lc")]
    if(nrow(sub>0)){
      for (scale in c(3,5,10,50)){
        difs <- round(sub[sub$scale==scale,]-sub[1,], 4)[3:10]
        aocdif <- round(((sub[sub$scale==scale,"range_area"]/round(sub[1,"range_area"], 4) - 1) * 100),4)
        collect <- c(common, sci_name, season, scale, as.numeric(aocdif), as.numeric(difs))
        if(scale>3){spcollect <- rbind(spcollect, collect)}else{spcollect <- collect}
      }
      if(!exists("out")){out <- spcollect}else{out <- rbind(out, spcollect)}
    }
  }
}
colnames(out) <- c("common_name", "sciname", "season", "scale",
                   "delta_range_area", paste0("delta_",metrics),
                   "delta_clim", "delta_top","delta_lc")
rownames(out) <- NULL
out <- data.frame(out)
write_csv(out, "~/ranges/birds_vs1_scale_comp.csv")

# Explore PPM change between 1 and 50km scales (old plot)
out <- read_csv("~/ranges/birds_vs1_scale_comp.csv")
for(season in c("breeding","nonbreeding")){
  if(season=="breeding"){col="darkorange1"}else{col="mediumorchid"}
  jpeg(paste0("~/ranges/scale_ppms_",season,".jpeg"), width=1200, height=800, units="px")
  par(mfrow=c(3,4), cex=1, mar=c(4.1, 4.1, 2.1, 2.1))
  for (scale in c(5,10,50)){
    out2 <- out[(out$season==season | out$season=="resident") &
                  out$scale==scale,]
    nbins <- seq(-.42, .14, by = .02)
    if(scale == 5){title="Change in AUC"}else{title=NULL}
    hist(out2$delta_AUC, breaks=nbins, col=col, border=col, main=title, 
         xlim=c(-.42, .14), xlab=expr(paste(Delta, " AUC")), ylab="Density")
    abline(v=0, col="gray40", lty="dashed")
    abline(v=median(out2$delta_AUC), col="black")
    nbins <- seq(-.76,.6, by = .04)
    if(scale == 5){title="Change in TSS"}else{title=NULL}
    hist(out2$delta_TSS, breaks=nbins, col=col, border=col, main=title, 
         xlim=c(-.76, .6), xlab=expr(paste(Delta, " TSS")), ylab=NULL)
    abline(v=0, col="gray40", lty="dashed")
    abline(v=median(out2$delta_TSS), col="black")
    nbins <- seq(-.72,.52, by = .04)
    if(scale == 5){title="Change in Sensitivity"}else{title=NULL}
    hist(out2$delta_Sensitivity, breaks=nbins, col=col, border=col, main=title, 
         xlim=c(-.72, .52), xlab=expr(paste(Delta, " Sensitivity")), ylab=NULL)
    abline(v=0, col="gray40", lty="dashed")
    abline(v=median(out2$delta_Sensitivity), col="black")
    nbins <- seq(-.48,.28, by = .04)
    if(scale == 5){title="Change in Specificity"}else{title=NULL}
    hist(out2$delta_Specificity, breaks=nbins, col=col, border=col, main=title, 
         xlim=c(-.48, .28), xlab=expr(paste(Delta, " Specificity")), ylab=NULL)
    abline(v=0, col="gray40", lty="dashed")
    abline(v=median(out2$delta_Specificity), col="black")
  }
  dev.off()
}

# AOC
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
  
out <- read_csv("~/ranges/birds_vs1_scale_comp.csv")
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
    
    # plot
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
      xlab(expression("% Change in Range Size")) +
      ylab("")
    assign(paste0("aoc.",scale), aocgg)
}
gg.aoc <- ggarrange(aoc.1, aoc.3, aoc.5, aoc.10, aoc.50, nrow=1, ncol=5)
ggsave("~/ranges/scale_aoc_density2.jpeg", gg.aoc, "jpeg", height=3, width=17, units="in")

# Violin plot for PPMs
data <- read_csv("~/ranges/sdm_bird_outputs_VS1_complete.csv")
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
ggsave("~/ranges/scale_ppms_violin2.jpeg", gg.ppm, "jpeg", height=6, width=12, units="in")

# Predictor importances
# Violin plot
for(season in c("breeding","nonbreeding")){
  for(metric in c("clim","top","lc")){
    if(metric=="clim"){metricname="Climate vars"}
    if(metric=="top"){metricname="Topographic vars"}
    if(metric=="lc"){metricname="Landcover vars"}
    data1km <- data[data$scale==1 & (data$season==season | data$season=="resident"),]
    data1km$scale <- as.factor(data1km$scale)
    data1km$metric <- data1km[,metric][[1]]
    gg <- ggplot(data=data1km, aes(x=scale, y=metric)) + 
      geom_hline(yintercept = mean(data1km$metric), size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .75, fill=met.brewer("Java", n = 5)[1]) +
      scale_fill_manual(values = met.brewer("Johnson", n = 5)[5])  +
      stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
      xlab('grain size (km)') +
      ylab(metricname) +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 16))
    
    out2 <- out[out$season==season | out$season=="resident",]
    out2$scale <- as.factor(out2$scale)
    out2$metric <- out2[,paste0("delta_",metric)][[1]]
    if(metric=="range_area"){out2$metric <- ifelse(out2$metric>0,out2$metric^(1/3),(-(-out2$metric)^(1/3)))}
    if(metric=="range_area"){out2$metric <- log(out2$metric)}
    gg2 <- ggplot(data=out2, aes(x=scale, y=metric, fill=scale)) + 
      geom_hline(yintercept = 0, size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .75) +
      scale_fill_manual(values = met.brewer("Johnson", n = 5, direction=-1)[2:5]) +
      stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
      ylab(expr(paste(Delta, " ", !!metricname))) +
      xlab('grain size (km)') +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 16))
    gg2
    assign(paste0("gg.",metric),gg)
    assign(paste0("gg2.",metric),gg2)
    
  }
  jpeg(paste0("~/ranges/scale_predimps_violin_",season,".jpeg"), width=1200, height=800)
  grid.arrange(gg.clim, gg2.clim, gg.top, gg2.top, gg.lc, gg2.lc,
               nrow=2, widths = c(1.5,3.5,1.5,3.5))
  dev.off()
}


{##### SPECIES-LEVEL VARIATION
  # runs AFTER ecoregions sheet
  
  ##### manually check biome assignments
  setwd("~/ranges/")
  
  # Species level data
  data <- read_csv("birds_vs1_sp_level_corrected.csv")
  
  # Create habitat diversity index
  output <- read_csv("sdm_bird_outputs_VS1_complete.csv")
  output <- output[output$scale==1,]
  ldis <- matrix(nrow=0, ncol=3)
  for (row in 1:nrow(output)){
  hab <- output[row, 48:83]
  hab <- hab/sum(hab, na.rm = T)
  hab[ hab == 0 ] <- NA
  ldi <- -1* sum(hab * log(hab), na.rm=T) / log(length(hab))
  ldis <- rbind(ldis, c(as.character(output[row, c(1,4)]), ldi))
  }
  ldis <- data.frame(ldis)
  names(ldis) <- c("common_name","season","ldi")
  data <- left_join(data, ldis, by=c("common_name","season"))
  data$ldi <- as.numeric(data$ldi)
  write_csv(data, "birds_vs1_sp_level_corrected_with_ldi.csv")
  
  # Process data
  data <- read_csv("birds_vs1_sp_level_corrected_with_ldi.csv")
  data$area <- sqrt(data$range_area)/1000
  data <- data[complete.cases(data),]
  data$biome <- recode_factor(data$biome, multiple = "Multiple") 
  comp <- data[data$scale==5,]
  
  
  # plotting
  for(season in c("breeding","nonbreeding")){
    comp.s <- comp[comp$season==season | comp$season=="resident",]
    for (var in c("ldi","area")){
      comp.s$pred <- comp.s[,var][[1]]
      if(var=="ldi"){xlab="Landcover Diversity Index"}else{
        xlab="grain size (km)"}
      if(season=="breeding"){xlab=""}
      # scatter plot
      gg <- ggplot(comp.s, aes(pred, delta_AUC)) +
        geom_point(size=2) +
        geom_smooth(color="black", method="lm") +
        theme(text = element_text(size=12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.y = element_text(color = "black", size=13),
              axis.text.x = element_text(color = "black", size=13)) +
        geom_hline(yintercept=0, lty=2) +
        xlab(xlab) +
        ylab(expr(paste(Delta, " AUC"))) 
      assign(paste0("gg.",var,".",season), gg)
    }
    
    ### figure coding region by delta AUC
    
    if(season=="breeding"){col="#F8766D"}else{col="#00BFC4"}
    comp2 = comp.s[comp.s$biome!="Multiple" &
                   comp.s$biome!="Temperate Conifer Forests"   ,]
    comp2$biome <- as.factor(comp2$biome)
    
    comp2$biome <- droplevels(comp2$biome)
    
    
    bio <- ggplot(data=comp2, aes(x=biome, y=delta_AUC)) + 
      geom_hline(yintercept = 0, size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .75, fill=col) +
      stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
      ylab(expr(paste(Delta, " AUC"))) +
      xlab("") +
      theme_light() +
      theme(legend.position = "none",
            text = element_text(size = 12),
            axis.text.x=element_blank())
    assign(paste0("bio.",season), bio)
  }
    # Model
    comp.st <- comp
    comp.st = comp.st[comp.st$biome!="Temperate Conifer Forests",]
    comp.st$biome <- droplevels(comp.st$biome)
    
    comp.st$ldi <- scale(comp.st$ldi)[,1]
    comp.st$area <- scale(comp.st$area)[,1]
    glm <- glm(delta_AUC ~ ldi + area + biome + season, data=comp.st)
    table <- data.frame(summary(glm)$coef)
    write.csv(table, file=paste0("patterns/glm_table.csv"))
    anova <- Anova(glm, test="F")
    table <- data.frame(anova)
    write.csv(table, file=paste0("patterns/anova_table.csv"))
  
  
  gga <- ggarrange(gg.ldi.breeding, gg.area.breeding,
                   gg.ldi.nonbreeding, gg.area.nonbreeding,
                   nrow=2, ncol=2)
  ggsave(paste0("patterns/sp_level_cont.jpeg"), gga, "jpeg",
         height=6, width=6, units="in")

  gga <- ggarrange(bio.breeding, bio.nonbreeding,
                   nrow=2, ncol=1)
  ggsave(paste0("patterns/sp_level_biomes.jpeg"), gga, "jpeg",
         height=7, width=7, units="in")
  
  levels(comp2$biome)
}














