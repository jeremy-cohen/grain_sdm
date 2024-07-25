# Simulate prediction accuracy in heterogenous and homogenous landscapes
# with fine and coarse grain layers

library(devtools)
install_github("syanco/checkyourself")
library(checkyourself)
library(raster)

# 100, 10 for homogenous
# 10, 100 for heterogenous
# parameter matrix
pm = rbind(c(100,10),c(10,100))
type = c("hom","het")

# 1000 iterations
table = matrix(ncol=5, nrow=0)
for (sim in 1:1000){
  if(sim>1){dev.off()}
  try({
  # loop to create two replicate "environmental layers" 
  # that are either homogenous or heterogenous
  for (g in 1:2){
    # loop to create homogenous and heterogenous rasters
    for (f in 1:2){
      mat.side <- 60 #set starting matrix size
      x <- matrix(1, mat.side, mat.side) #create the matrix
      #Declarations
      size.clusters <- pm[f,1] #define target size of each cluster to be grown (in # of cells)
      n.clusters <- pm[f,2] #define approximate number of clusters of Habitat A to "grow"
      count.max <- 100 #set maximum number of iterations throught he while-loop
      #Required Initial Objects
      n <- mat.side * mat.side #total number of cells in matrix
      cells.left <- 1:n #create 'cells.left' object
      cells.left[x!=1] <- -1 # Indicates occupancy of cells
      i <- 0 #i counts clusters created and should start at 0 always
      indices <- c() #create empty vector for indices
      ids <- c() #create empty vector for ids
      
      while(i < n.clusters && length(cells.left) >= size.clusters && count.max > 0) {
        count.max <- count.max-1 #countdown against max number of loops
        xy <- sample(cells.left[cells.left > 0], 1) #randomly draw an unoccupied cell
        cluster <- expand(x, size.clusters, xy) #run expand function to grow that cluster
        if (!is.na(cluster[1]) && length(cluster)==size.clusters) {
          i <- i+1 #add to cluster count
          ids <- c(ids, rep(i, size.clusters)) #add cluster to id list
          indices <- c(indices, cluster) #add cluster to indices list
          cells.left[indices] <- -1 #remove all cells in the cluster grown from the available list
        }
      }
      
      y <- matrix(NA, mat.side, mat.side) #create blank matrix of the same size as `x`.
      
      #Add the cluster ids to the matrix at locations in indices - this adds each cluster id to the cells indicated by the
      #vector 'indices' and leaves the rest of the cells as 'NA'
      y[indices] <- ids
      
      # constrain two landcover values per raster
      y[y<=(n.clusters/2)] = 0 # set to one of two landcover values
      y[y>(n.clusters/2)] = 1 # set to one of two landcover values
      y[is.na(y)] = 1 # resolve NAs
      ras = raster(y)
      # plot(ras)
      assign(paste0(type[f],".ras.",g), ras)
    }
  }
  
  # stack each set of environmental covariates
  hom.stack = stack(hom.ras.1, hom.ras.2)
  het.stack = stack(het.ras.1, het.ras.2)
  # vector versions of fine grain landscapes
  hom.dat.fine = getValues(hom.stack)
  het.dat.fine = getValues(het.stack)
  # coarse grain version of each landscape
  hom.stack.coarse = raster::aggregate(hom.stack, 3, modal)
  hom.dat.coarse = getValues(hom.stack.coarse)
  het.stack.coarse = raster::aggregate(het.stack, 3, modal)
  het.dat.coarse = getValues(het.stack.coarse)
  
  # simulate 100 points with probabilities tied to landscape
  # same pts for hom/het landscapes
  hom.pts = het.pts = data.frame(cbind(runif(100),runif(100)))
  # get pres/abs based on presence of both landscape features
  ex = raster::extract(hom.stack, hom.pts)
  hom.pts$ex = ex[,1]*ex[,2]
  # flip 10% of features to simulate randomness of occurrence
  samp = sample(100,10)
  hom.pts$ex[samp] = 1 - hom.pts$ex [samp]
  # same for heterogenous landscapes
  ex = raster::extract(het.stack, het.pts)
  het.pts$ex = ex[,1]*ex[,2]
  het.pts$ex[samp] = 1 - het.pts$ex [samp]
  
  # SDMs (4 combos of landscape x grain)
  # homogenous landscape
  land = raster::extract(hom.stack, hom.pts[,1:2])
  hom.pts$land1 = land[,1]
  hom.pts$land2 = land[,2]
  sdm.hom.fine = lm(ex ~ land1 + land2, data = hom.pts)
  land = raster::extract(hom.stack.coarse, hom.pts[,1:2])
  hom.pts$land1c = land[,1]
  hom.pts$land2c = land[,2]
  sdm.hom.coarse = lm(ex ~ land1c + land2c, data = hom.pts)
  # heterogenous landscape
  land = raster::extract(het.stack, het.pts[,1:2])
  het.pts$land1 = land[,1]
  het.pts$land2 = land[,2]
  sdm.het.fine = lm(ex ~ land1 + land2, data = het.pts)
  land = raster::extract(het.stack.coarse, het.pts[,1:2])
  het.pts$land1c = land[,1]
  het.pts$land2c = land[,2]
  sdm.het.coarse = lm(ex ~ land1c + land2c, data = het.pts)
  
  {# predictions
    
    # blank rasters
    fine = raster(nrows=60, ncols=60, xmn=0, xmx=1, ymn=0, ymx=1)
    coarse = raster(nrows=20, ncols=20, xmn=0, xmx=1, ymn=0, ymx=1)
    
    # plot (fig. 1b)
    jpeg(paste0("~/ranges/simulation/plots/simulation_",sim,".jpeg"), 
         width=1000, height=250, units="px", quality=100)
    par(mfrow=c(1,4), mar=c(2,2,2,1))
    coefs = c()
    
    het.dat2 = data.frame(het.dat.fine)
    names(het.dat2) = c("land1","land2")
    pred = predict(sdm.het.fine, het.dat2)
    pred = round(pred, 0)
    pred[pred==-1] = 0
    pred.het.fine = setValues(fine, pred)
    plot(pred.het.fine, legend=F, axes=F)
    points(het.pts[het.pts$ex==1,], pch=16, cex=2)
    points(het.pts[het.pts$ex==0,], cex=2)
    het.pts$pre = raster::extract(pred.het.fine, het.pts[,1:2])
    coefs = c(coefs, coef(lm(ex ~ pre, data=het.pts))[2])
    
    het.dat2 = data.frame(het.dat.coarse)
    names(het.dat2) = c("land1c","land2c")
    pred = predict(sdm.het.coarse, het.dat2)
    pred = round(pred, 0)
    pred[pred==-1] = 0
    pred.het.coarse = setValues(coarse, pred)
    plot(pred.het.coarse, legend=F, axes=F)
    points(het.pts[het.pts$ex==1,], pch=16, cex=2)
    points(het.pts[het.pts$ex==0,], cex=2)
    het.pts$pre.coarse = raster::extract(pred.het.coarse, het.pts[,1:2])
    coefs = c(coefs, coef(lm(ex ~ pre.coarse, data=het.pts))[2])
    
    hom.dat2 = data.frame(hom.dat.fine)
    names(hom.dat2) = c("land1","land2")
    pred = predict(sdm.hom.fine, hom.dat2)
    pred = round(pred, 0)
    pred[pred==-1] = 0
    pred.hom.fine = setValues(fine, pred)
    plot(pred.hom.fine, legend=F, axes=F)
    points(hom.pts[hom.pts$ex==1,], pch=16, cex=2)
    points(hom.pts[hom.pts$ex==0,], cex=2)
    hom.pts$pre = raster::extract(pred.hom.fine, hom.pts[,1:2])
    coefs = c(coefs, coef(lm(ex ~ pre, data=hom.pts))[2])
    
    hom.dat2 = data.frame(hom.dat.coarse)
    names(hom.dat2) = c("land1c","land2c")
    pred = predict(sdm.hom.coarse, hom.dat2)
    pred = round(pred, 0)
    pred[pred==-1] = 0
    pred.hom.coarse = setValues(coarse, pred)
    plot(pred.hom.coarse, legend=F, axes=F)
    points(hom.pts[hom.pts$ex==1,], pch=16, cex=2)
    points(hom.pts[hom.pts$ex==0,], cex=2)
    hom.pts$pre.coarse = raster::extract(pred.hom.coarse, hom.pts[,1:2])
    coefs = c(coefs, coef(lm(ex ~ pre.coarse, data=hom.pts))[2])
    
  } # end plot
  
  # compile coefs
  coefs = c(sim, coefs)
  table = rbind(table, coefs)
  
  print(paste("end simulation",sim))
  })
} # end iteration
dev.off()
# name and save table
table = data.frame(table[complete.cases(table),])
names(table) = c("sim", "het.fine", "het.coarse", "hom.fine", "hom.coarse")
means = colMeans(table[,2:5])
means
write_csv(table, "~/ranges/simulation/simulation_table.csv")
write.csv(means, "~/ranges/simulation/simulation_means.csv")
