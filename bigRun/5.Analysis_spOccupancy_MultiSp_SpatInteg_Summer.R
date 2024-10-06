##############################################################################
##############################################################################
########                  Run spOccupancy - A-BioTrack           #############
########                           June 2024                     #############
########       Single-species spatial integrated  occupancy      #############
########                Multiple species aggregations            #############
########                           Warm period                   #############
##############################################################################

## high performance compute settings
N_THREADS = 32

## Upload package and organize A-BioTrack data 
library(spOccupancy)
library(coda)
library(stars)
library(ggplot2)
set.seed(102)

## Upload Grid data for mapping
FullGrid <- st_read("Shapes/FullDataGrid.shp")

SharkSpp <- c("Alopias vulpinus", "Carcharhinus acronotus", "Carcharhinus brevipinna", "Carcharhinus leucas",
              "Carcharhinus limbatus", "Carcharhinus obscurus", "Carcharhinus perezii", "Carcharhinus plumbeus",
              "Carcharias taurus", "Galeocerdo cuvier", "Ginglymostoma cirratum", "Isurus oxyrinchus",
              "Mustelus canis", "Negaprion brevirostris", "Prionace glauca", "Sphyrna mokarran")

TurtleSpp <- c("Caretta caretta", "Chelonia mydas", "Dermochelys coriacea", "Lepidochelys olivacea")

RaySpp <- c("Bathytoshia centroura", "Hypanus americanus", "Hypanus sabinus",  "Hypanus say",
            "Mobula birostris", "Rhinoptera bonasus", "Rostroraja eglanteria")

SealSpp <- c("Halichoerus grypus")

FishSpp <- c("Acipenser brevirostrum", "Acipenser oxyrinchus", "Albula vulpes",
             "Alosa aestivalis", "Alosa pseudoharengus", "Alosa sapidissima",
             "Archosargus probatocephalus", "Centropomus undecimalis", "Centropristis striata",
             "Cynoscion nebulosus", "Cynoscion regalis",
             "Epinephelus adscensionis", "Epinephelus guttatus", "Epinephelus itajara",
             "Epinephelus morio", "Epinephelus striatus", "Gadus morhua",
             "Lachnolaimus maximus", "Lutjanus analis", "Lutjanus apodus",
             "Lutjanus griseus", "Lutjanus jocu", "Megalops atlanticus",
             "Micropterus salmoides", "Morone saxatilis", "Mugil cephalus",
             "Mycteroperca bonaci", "Mycteroperca interstitialis", "Mycteroperca microlepis",
             "Mycteroperca phenax", "Mycteroperca venenosa", "Ocyurus chrysurus",
             "Paralichthys dentatus", "Paralichthys lethostigma", "Pseudopleuronectes americanus",
             "Pterois volitans", "Rachycentron canadum", "Sciaenops ocellatus",
             "Sphyraena barracuda", "Trachinotus falcatus")

AllSpp <- c(SharkSpp, TurtleSpp, RaySpp, SealSpp, FishSpp)
AllCladeNames <- c("Selachii", "Chelonioidea", "Batoidea", "Pinnipedia", "Actinopterygii")
# "Decapoda"
AllCladeList <- c(rep("Selachii", times = length(SharkSpp)),
                  rep("Chelonioidea", times = length(TurtleSpp)),
                  rep("Batoidea", times = length(RaySpp)),
                  rep("Pinnipedia", times = length(SealSpp)),
                  rep("Actinopterygii", times = length(FishSpp))
                  # ,
                  # rep("Decapoda", times = length(InvetSpp))
)

# Loop over spp indexes to import all matrices and add them to a list
spDetectList_Tel <- list()
spDetectList_Obis <- list()
GridDetEnvList_Obis <- list()

## Import detection histories for all species 
for (i in 1: length(AllSpp)){
  objName.tel <- paste ("sp",i ,".tel", sep = "")
  fileName.tel <- paste ("spOccupancy_MultiSpp_FullArea/MultiSpp_DetHist_Sp",i ,".txt", sep = "")
  file.tel <- read.table(file = fileName.tel , header = T)
  spDetectList_Tel[[i]] <- assign(objName.tel, value = file.tel)
  
  objName.obis <- paste ("sp",i ,".obis", sep = "")
  fileName.obis <- paste ("spOccupancy_MultiSpp_FullArea/MultiSpp_ObisHist_Sp",i ,".txt", sep = "")
  file.obis <- read.table(file = fileName.obis , header = T)
  spDetectList_Obis[[i]] <- assign(objName.obis, value = file.obis)
  
  objName.detcov.obis <- paste ("Grid_DetEnv_Obis_sp",i, sep = "")
  fileName.detcov.obis <- paste ("spOccupancy_MultiSpp_FullArea/MultiSpp_OBIS_DetectCovs_Sp",i ,".txt", sep = "")
  file.detcov.obis <- read.table(file = fileName.detcov.obis , header = T)
  GridDetEnvList_Obis[[i]] <- assign(objName.detcov.obis, value = file.detcov.obis)
  
}

## Import Detection Covs   *Add later OBIS det covs per species
Grid_DetEnv_Tel <- read.table(file = "spOccupancy_MultiSpp_FullArea/Grid_DetCovs.txt", header = T)
# Grid_DetEnv_Obis <- read.table(file = "spOccupancy_MultiSpp_FullArea/MultiSpp_OBIS_DetectCovs_Sp7.txt", header = T)
Grid_DetEnv <- read.table(file = "spOccupancy_MultiSpp_FullArea/Grid_DetCovs.txt", header = T)

## Import data on occupancy covariates and full detect covariates  *Figure out why depth is not in the matrix
# Grid_OccEnv <- read.table(file = "spOccupancy_MultiSpp_FullArea/Grid_OccEnv.txt", header = T)
Grid_OccEnv <- read.table(file = "Grid_OccEnv_Seasonal.txt", header = T)

## Code added by VZ: Save the objects we will need
save(AllSpp, spDetectList_Tel, spDetectList_Obis, GridDetEnvList_Obis, Grid_OccEnv, 
     Grid_DetEnv, Grid_DetEnv_Tel, FullGrid,
     file="DataSummer.RData")

## Code added by VZ: Remove everything from the environment
rm(list=ls())

## Code added by VZ: Load the saved objects back
load("DataSummer.RData")

## Create an empty matrix to compute Bayesian p-value and k-fold estimates
ModelValid <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
colnames(ModelValid) <- c("spID", "Species", "k-fold.integ.tel", "k-fold.integ.obis",
                          "k-fold.tel.alone", "k-fold.obis.alone")

## Loop over species
for(j in 1:length(AllSpp)){
  
  ## Fill summer months (May to October)
  SpDetectHistory_Tel <- spDetectList_Tel[[j]][,5:10]
  SpDetectHistory_Obis <- spDetectList_Obis[[j]][,5:10]
  Grid_DetEnv_Obis <- GridDetEnvList_Obis[[j]] # May need add seasonal detection covs in the future
  
  ## Filter data and keep just the grids where the species can be detected 
  GridCellsToRemove_Tel <- as.numeric(names(which((rowSums(is.na(SpDetectHistory_Tel)) == 6))))
  GridCellsToRemove_Obis <- as.numeric(names(which((rowSums(is.na(SpDetectHistory_Obis)) == 6))))
  
  # Telemetry
  sp.y.tel <- as.matrix(SpDetectHistory_Tel[!rowSums(is.na(SpDetectHistory_Tel)) == 6, ])
  coords.tel <- as.matrix(Grid_OccEnv[as.numeric(row.names(sp.y.tel)), c(1, 2)])
  occ.covs.tel <- as.matrix(Grid_OccEnv[as.numeric(row.names(sp.y.tel)), c(3:18)]) ##**Add depth later
  det.covs.tel <- as.matrix(Grid_DetEnv_Tel[as.numeric(row.names(sp.y.tel)), c(1,8)])
  # Depth.tel <- as.matrix((det.covs.tel[, 1]))
  Depth.tel_RAW <- as.matrix(Grid_OccEnv[as.numeric(row.names(sp.y.tel)), "Depth"]) ##Add depth into detection covs
  Depth.tel_RAW[Depth.tel_RAW >= 0] <- -0.1   # Remove positive values
  Depth.tel <- -Depth.tel_RAW  # Invert values to allow a log transformation
  NReceiv.tel <- as.matrix((det.covs.tel[, 2]))
  sites.tel <- as.numeric(rownames(sp.y.tel))
  
  # Obis
  sp.y.obis <- as.matrix(SpDetectHistory_Obis[!rowSums(is.na(SpDetectHistory_Obis)) == 6, ])
  # sp.y.obis[is.na(sp.y.obis)] <- 0
  coords.obis <- as.matrix(Grid_OccEnv[as.numeric(row.names(sp.y.obis)), c(1, 2)])
  occ.covs.obis <- as.matrix(Grid_OccEnv[as.numeric(row.names(sp.y.obis)), c(3:18)]) ##**Add depth later
  det.covs.obis <- as.matrix(Grid_DetEnv_Obis[!rowSums(is.na(SpDetectHistory_Obis)) == 6, c(1,2)])
  Depth.obis_RAW <- as.matrix((occ.covs.obis[, 1]))
  Depth.obis_RAW[Depth.obis_RAW >= 0] <- -0.1  # Remove positive values
  Depth.obis <- -Depth.obis_RAW  # Invert values to allow a log transformation
  NLists.obis <- as.matrix((det.covs.obis[, 1]))
  NSppDetec.obis <- as.matrix((det.covs.obis[, 2]))
  sites.obis <- as.numeric(rownames(sp.y.obis))
  
  ## Check the unique site IDs with data, detections and non-detections
  UniqueSites <- sort(unique(c(sites.obis, sites.tel)))
  length(UniqueSites)/dim(FullGrid)[1] # Proportion
  
  ## Rename unique and duplicated site IDs following a sequential order
  CellSeqIDs <- cbind("SeqID" = 1:length(UniqueSites), UniqueSites)
  OrigIDsTel <- as.numeric(rownames(sp.y.tel))
  OrigIDsObis <- as.numeric(rownames(sp.y.obis))
  TelIDs <- CellSeqIDs[CellSeqIDs[,2] %in% OrigIDsTel, 1]
  ObisIDs <- CellSeqIDs[CellSeqIDs[,2] %in% OrigIDsObis, 1]
  
  ## Add new names to original matrices
  row.names(sp.y.tel) <- TelIDs
  row.names(coords.tel) <- TelIDs
  row.names(occ.covs.tel) <- TelIDs
  row.names(det.covs.tel) <- TelIDs
  row.names(Depth.tel) <- TelIDs
  row.names(NReceiv.tel) <- TelIDs
  # row.names(Season.tel) <- TelIDs
  sites.tel <- TelIDs
  
  row.names(sp.y.obis) <- ObisIDs
  row.names(coords.obis) <- ObisIDs
  row.names(occ.covs.obis) <- ObisIDs
  row.names(det.covs.obis) <- ObisIDs
  row.names(Depth.obis) <- ObisIDs
  row.names(NLists.obis) <- ObisIDs
  row.names(NSppDetec.obis) <- ObisIDs
  # row.names(Season.obis) <- ObisIDs
  sites.obis <- ObisIDs
  
  
  ## Create an integrated occupancy covariates matrix
  # occ.covs.int <- Grid_OccEnv[UniqueSites, 3, drop = FALSE]
  occ.covs.int <- Grid_OccEnv[UniqueSites, 3:18]   # Add depth
  occ.covs.int$Depth[occ.covs.int$Depth >= 0] <- -0.1 # Remove positive values
  occ.covs.int$Depth <- -occ.covs.int$Depth # Convert depth to positve
  rownames(occ.covs.int) <- CellSeqIDs[,1]
  coords.int <- Grid_OccEnv[UniqueSites, c(1,2)]
  rownames(coords.int) <- CellSeqIDs[,1]
  # occ.covs.int <- Grid_OccEnv[ , 3, drop = FALSE]
  # coords.int <- Grid_OccEnv[ , c(1,2)]
  
  
  ## Merge datasets to run the integrated model
  y.int <- list(telemetry = sp.y.tel[,1:6], obis = sp.y.obis[,1:6])
  # y.int <- list(telemetry = sp.y.tel, obis = sp.y.obis)
  rownames(y.int$telemetry) <- 1:dim(y.int$telemetry)[1]
  rownames(y.int$obis) <- 1:dim(y.int$obis)[1]
  
  det.covs.int <- list(telemetry = list("NReceiv" = NReceiv.tel, "Depth" = Depth.tel),
                       obis = list("NLists" = NLists.obis, "Depth" = Depth.obis, "NSppDetec" = NSppDetec.obis))
  
  sites.int <- list(telemetry = sites.tel, obis = sites.obis)
  
  ## Add everything to one single list
  data.int <- list("y" = y.int, 
                   "occ.covs" = occ.covs.int, 
                   "det.covs" = det.covs.int, 
                   "sites" = sites.int,
                   "coords" = coords.int)
  
  str(data.int)
  
  
  ##### Run model
  occ.formula.int <- ~  scale(Depth) + scale(SST_s) +  I(scale(SST_s)^2) +
    scale(Chlor_s) + I(scale(Chlor_s)^2) + scale(TSM_s) + scale(SSH_s)
  
  det.formula.int <- list(telemetry = ~ scale(Depth) + scale(log(NReceiv)),
                          obis = ~ scale(log(NLists)))   #scale(Depth) +
  
  ## Initial values
  dist.int <- dist(data.int$coords)
  min.dist <- min(dist.int)
  max.dist <- max(dist.int)
  J <- nrow(data.int$occ.covs)
  
  # Exponential covariance model
  cov.model <- "exponential"
  
  ## Inits (abbreviated format) - Spatial
  inits.list <- list(alpha = list(0, 0),
                     beta = 0,
                     # sigma.sq.p = list(0.5, 0.5),
                     z = rep(1, J),
                     sigma.sq = 2,
                     phi = 3 / mean(dist.int),
                     w = rep(0, J))
  
  ## Priors - Spatial
  prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                     alpha.normal = list(mean = list(0, 0),
                                         var = list(2.72, 2.72)),
                     sigma.sq.ig = c(2, 1),
                     phi.unif = c(3 / max.dist, 3 / min.dist))
  
  ## Model settings - Spatial
  batch.length <- 80  # Samples per chain = batch length * n.batch
  n.batch <- 4000
  n.burn <- 150000
  n.thin <- 200
  tuning <- list(phi = 0.2)
  
  ## Run model - Spatial
  out.sp.int <- spIntPGOcc(occ.formula = occ.formula.int,
                           det.formula = det.formula.int,
                           data = data.int,
                           inits = inits.list,
                           priors = prior.list,
                           tuning = tuning,
                           cov.model = cov.model,
                           NNGP = TRUE, # NNGP Vs Full Gaussian Process
                           verbose = TRUE,
                           n.neighbors = 6,
                           n.batch = n.batch,
                           n.burn = n.burn,
                           n.thin = n.thin,
                           n.chains = 3,
                           n.omp.threads = N_THREADS, # Within chain parallel running
                           batch.length = batch.length,
                           n.report = 100)
  
  ## Export output
  sink(file = paste("Output_Summer_", AllSpp[j],".txt", sep = ""))
  print(AllSpp[j])
  summary(out.sp.int)
  sink(file = NULL)
  

  #####################  Model validation  #########
  
  ## Model validation - Goodness of Fit (GoF) assessment (Bayesian p-value)
  ppc.out.g1 <- ppcOcc(out.sp.int, fit.stat = 'freeman-tukey', group = 1)
  summary(ppc.out.g1)
  
  ## Export output
  sink(file = paste("Bayesian_p-valueFullsummer_", AllSpp[j],".txt", sep = ""))
  print(AllSpp[j])
  summary(ppc.out.g1)
  sink(file = NULL)
  
  pdf(file = paste("FullsummerGOF2_", AllSpp[j],".pdf", sep = ""))
  
  # create layout for figures
  par(mfrow=c(1,2))
  
  ppc.df <- data.frame(fit = ppc.out.g1$fit.y[[1]],
                       fit.rep = ppc.out.g1$fit.y.rep[[1]],
                       color = 'lightskyblue1')
  ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
  plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21,
       ylab = 'Fit', xlab = 'True', main = c(AllSpp[j], "Telemetry"))
  lines(ppc.df$fit, ppc.df$fit, col = 'black')
  
  ppc.df <- data.frame(fit = ppc.out.g1$fit.y[[2]],
                       fit.rep = ppc.out.g1$fit.y.rep[[2]],
                       color = 'lightskyblue1')
  ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
  plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21,
       ylab = 'Fit', xlab = 'True', main = c(AllSpp[j], "OBIS"))
  lines(ppc.df$fit, ppc.df$fit, col = 'black')
  
  dev.off()
  
  ModelValid[j,1] <- j
  ModelValid[j,2] <- AllSpp[j]
  #ModelValid[j,3] <- ppc.out.g1$fit.y.rep
  ppc.df$color[ppc.df$fit.rep > ppc.df$fit]
  
  diff.fit.tel <- ppc.out.g1$fit.y.rep.group.quants[[1]][3, ] - ppc.out.g1$fit.y.group.quants[[1]][3, ]
  plot(diff.fit.tel, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')
  
  diff.fit.obis <- ppc.out.g1$fit.y.rep.group.quants[[2]][3, ] - ppc.out.g1$fit.y.group.quants[[2]][3, ]
  plot(diff.fit.obis, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')
  
  
  boxplot(diff.fit.tel, diff.fit.obis)
  mean(diff.fit.tel)
  quantile(diff.fit.tel, probs = c(0.025, 0.975))
  
  mean(diff.fit.obis)
  quantile(diff.fit.obis, probs = c(0.025, 0.975))
  
  
  ######## k-fold cross-validation
  ## Integrated model (Spatial)
  kfold.sp.int <- spIntPGOcc(occ.formula = occ.formula.int,
                             det.formula = det.formula.int,
                             data = data.int,
                             inits = inits.list,
                             priors = prior.list,
                             tuning = tuning,
                             cov.model = cov.model,
                             NNGP = TRUE, # NNGP Vs Full Gaussian Process
                             verbose = TRUE,
                             n.neighbors = 6,
                             n.batch = n.batch,
                             n.burn = n.burn,
                             n.chains = 3,
                             n.omp.threads = 1, # Within chain parallel running
                             batch.length = batch.length,
                             n.report = 100,
                             k.fold = 4,
                             k.fold.only = T)
  
  ModelValid[j,3:4] <- kfold.sp.int$k.fold.deviance
  
  ## Telemetry alone (Spatial)
  kfold.sp.tel <- spIntPGOcc(occ.formula = occ.formula.int,
                             det.formula = det.formula.int,
                             data = data.int,
                             inits = inits.list,
                             priors = prior.list,
                             tuning = tuning,
                             cov.model = cov.model,
                             NNGP = TRUE, # NNGP Vs Full Gaussian Process
                             verbose = TRUE,
                             n.neighbors = 6,
                             n.batch = n.batch,
                             n.burn = n.burn,
                             n.chains = 3,
                             n.omp.threads = 1, # Within chain parallel running
                             batch.length = batch.length,
                             n.report = 100,
                             k.fold = 4,
                             k.fold.data = 1,
                             k.fold.only = T)
  
  ModelValid[j,5] <- kfold.sp.tel$k.fold.deviance
  
  ## OBIS alone (Spatial)
  kfold.sp.obis <- spIntPGOcc(occ.formula = occ.formula.int,
                              det.formula = det.formula.int,
                              data = data.int,
                              inits = inits.list,
                              priors = prior.list,
                              tuning = tuning,
                              cov.model = cov.model,
                              NNGP = TRUE, # NNGP Vs Full Gaussian Process
                              verbose = TRUE,
                              n.neighbors = 6,
                              n.batch = n.batch,
                              n.burn = n.burn,
                              n.chains = 3,
                              n.omp.threads = 1, # Within chain parallel running
                              batch.length = batch.length,
                              n.report = 100,
                              k.fold = 4,
                              k.fold.data = 2,
                              k.fold.only = T)
  
  ModelValid[j,6] <- kfold.sp.obis$k.fold.deviance
  
  write.table(ModelValid, file = "kFoldCrossResultsSummer.txt")
  
  
  ##################################################
  
  #### Predict  ##########################
  
  ## Predictions for the whole study area
  str(Grid_OccEnv)
  DepthNeg <- Grid_OccEnv$Depth
  DepthNeg[DepthNeg >= 0] <- -0.1
  DepthTrans <- -DepthNeg
  Depth.pred <- (DepthTrans - mean(data.int$occ.covs[, 1])) / sd(data.int$occ.covs[, 1])  ## Add depth
  SST.pred <- (Grid_OccEnv$SST_s - mean(data.int$occ.covs[, 3])) / sd(data.int$occ.covs[, 3])
  Chlor.pred <- (Grid_OccEnv$Chlor_s - mean(data.int$occ.covs[, 6])) / sd(data.int$occ.covs[, 6])
  TSM.pred <- (Grid_OccEnv$TSM_s - mean(data.int$occ.covs[, 9])) / sd(data.int$occ.covs[, 9])
  SSH.pred <- (Grid_OccEnv$SSH_s - mean(data.int$occ.covs[, 15])) / sd(data.int$occ.covs[, 15])
  # These are the new intercept and covariate data.
  # X.0 <- cbind(1, Depth.pred, SST.pred,  SST.pred^2, Chlor.pred, Chlor.pred^2, SSH.pred)  #Depth.pred
  X.0 <- cbind(1, Depth.pred, SST.pred,  SST.pred^2, Chlor.pred, Chlor.pred^2, TSM.pred, SSH.pred)  #Depth.pred
  coords.0 <- as.matrix(Grid_OccEnv[, c('X', 'Y')])
  out.sp.pred <- predict(out.sp.int, X.0, coords.0, verbose = FALSE) # Spatial
  # out.sp.pred <- predict(out.sp.int, X.0) # Non-spatial

  # Produce a species distribution map (posterior predictive means of occupancy)
  plot.dat <- data.frame(x = Grid_OccEnv$X,
                         y = Grid_OccEnv$Y,
                         mean.psi = apply(out.sp.pred$psi.0.samples, 2, mean),
                         sd.psi = apply(out.sp.pred$psi.0.samples, 2, sd))

  class(FullGrid)
  plot.grid <- cbind(st_as_sf(FullGrid), plot.dat$mean.psi, plot.dat$sd.psi)

  ## Export predictive map as shapefile
  spName_shape <- chartr(" ", "_", AllSpp[j])
  file_name_shape = paste("SDM_Shape_Summer_", spName_shape, ".shp", sep="")
  st_write(plot.grid, file_name_shape, append = FALSE)

  ## Plot mean occupancy
  file_name_mean = paste("Mean_Summer_", AllSpp[j], ".jpeg", sep="")
  mean.Psi.Plot <- ggplot() +
    ggtitle(AllSpp[j]) +
    geom_sf(data = plot.grid, mapping = aes(fill = plot.dat.mean.psi), color = NA) +
    scale_fill_viridis_c("mean psi", na.value = "transparent")
  mean.Psi.Plot
  # Export map
  ggsave(filename = paste(file_name_mean),
         plot =  mean.Psi.Plot, width = 10, height = 9, units = 'cm',
         scale = 2, dpi = 100)

  ## Plot SD
  file_name_sd = paste("SD_Summer_", AllSpp[j], ".jpeg", sep="")
  SD.Psi.Plot <- ggplot() +
    ggtitle(AllSpp[j]) +
    geom_sf(data = plot.grid, aes(fill = plot.dat.sd.psi), color = NA) +
    scale_fill_viridis_c("SD psi", na.value = "transparent")
  SD.Psi.Plot
  # Export map
  ggsave(filename = paste(file_name_sd),
         plot =  SD.Psi.Plot, width = 10, height = 9, units = 'cm',
         scale = 2, dpi = 100)


  ## Create an empty matrix to fill with z values per site (1000 posterior estimates)
  NSites <- dim(out.sp.pred$z.0.samples)[2]
  Z_sp <- matrix(NA, nrow = NSites, ncol = 1000)

  ## Fill out the matrix with posterior z estimates for every site
  for (i in 1:NSites){
    sampledZ <- sample(out.sp.pred$z.0.samples[,i], size = 100)
    Z_sp[i,] <- sampledZ
  }

  ## Export z posteriors
  write.table(Z_sp, file = paste("Zposterior_Summer_Sp", AllSpp[j],".txt", sep = ""))
  
  
  ## Clean up memory to run for other species
  rm(out.sp.int, out.sp.pred, plot.dat, plot.grid, Z_sp)

  print(j)
  
}






###############################################################################
##########   Generate multi-species maps based on posterior Z  ################
###############################################################################

# Loop over posterior Zs to estimate spp. richness
# ZList <- list()
# 
# ## Import data on posterior Zs (one matrix per sp)
# for (i in 1: length(AllSpp)){
#   objName.Z <- paste ("Zsp",i, sep = "")
#   fileName.Z <- paste ("spOccupancy_MultiSpp_FullArea/Zposterior_Sp",i ,".txt", sep = "")
#   file.Z <- read.table(file = fileName.Z , header = T)
#   ZList[[i]] <- assign(objName.Z, value = file.Z)
# }
# 
# 
# ## Sum matrices to generate species richness
# Zrich <- ZList[[1]] + ZList[[2]] + ZList[[3]] + ZList[[4]] + ZList[[5]] + ZList[[6]] +
#   ZList[[7]] + ZList[[8]] + ZList[[9]] + ZList[[10]] + ZList[[11]] + ZList[[12]] +
#   ZList[[13]] + ZList[[14]] + ZList[[15]] + ZList[[16]] + ZList[[17]] + ZList[[18]] +
#   ZList[[19]] + ZList[[20]] + ZList[[21]] + ZList[[22]] + ZList[[23]] + ZList[[24]] +
#   ZList[[25]] + ZList[[26]] + ZList[[27]] + ZList[[28]] + ZList[[29]] + ZList[[30]] +
#   ZList[[31]] + ZList[[32]]
# 
# # Produce a species distribution map (posterior predictive means of richness)
# plot.dat <- data.frame(x = Grid_OccEnv$X, 
#                        y = Grid_OccEnv$Y, 
#                        mean.rich = apply(Zrich, 1, mean), 
#                        sd.rich = apply(Zrich, 1, sd))
# 
# class(FullGrid)
# plot.grid <- cbind(st_as_sf(FullGrid), plot.dat$mean.rich, plot.dat$sd.rich)
# 
# ## Export predictive map as shapefile
# file_name_Richshape = paste("SDM_RichShape_", "Sharks", ".shp", sep="")
# st_write(plot.grid, file_name_Richshape, append = FALSE)
# 
# ## Plot Richness mean
# file_name_RichMean = "0MultiRichMean.jpeg"
# MeanRichMap <- ggplot() +
#   geom_sf(data = plot.grid, aes(fill = plot.dat.mean.rich), color = NA) +
#   scale_fill_viridis_c("N Spp", na.value = "transparent", option = "A", direction = -1)
# MeanRichMap
# # Export map
# ggsave(filename = paste(file_name_RichMean),
#        plot =  MeanRichMap, width = 10, height = 9, units = 'cm',
#        scale = 2, dpi = 100)
# 
# 
# # Plot Richness SD
# file_name_RichSD = "0MultiRichSD.jpeg"
# SDRichMap <- ggplot() +
#   geom_sf(data = plot.grid, aes(fill = plot.dat.sd.rich), color = NA) +
#   scale_fill_viridis_c("SD", na.value = "transparent", option = "A", direction = -1)
# SDRichMap
# # Export map
# ggsave(filename = paste(file_name_RichSD),
#        plot =  SDRichMap, width = 10, height = 9, units = 'cm',
#        scale = 2, dpi = 100)
# 
# 
# ###############################################################################
# ##########         Generate ecoregion-based indicators         ################
# ###############################################################################
# 
# # Loop over posterior Zs to estimate spp. richness
# ZList <- list()
# 
# ## Import data on posterior Zs (one matrix per sp)
# for (i in 1: 16){
#   objName.Z <- paste ("Zsp",i, sep = "")
#   fileName.Z <- paste ("spOccupancy_MultiSpp_FullArea/Zposterior_Sp",i ,".txt", sep = "")
#   file.Z <- read.table(file = fileName.Z , header = T)
#   ZList[[i]] <- assign(objName.Z, value = file.Z)
# }
# 
# 
# ## Sum matrices to generate species richness
# Zrich <- ZList[[1]] + ZList[[2]] + ZList[[3]] + ZList[[4]] + ZList[[5]] + ZList[[6]] +
#   ZList[[7]] + ZList[[8]] + ZList[[9]] + ZList[[10]] + ZList[[11]] + ZList[[12]] +
#   ZList[[13]] + ZList[[14]] + ZList[[15]] + ZList[[16]]
# 
# 
# # Produce a species distribution map (posterior predictive means of richness)
# plot.dat <- data.frame(x = Grid_OccEnv$X, 
#                        y = Grid_OccEnv$Y, 
#                        mean.rich = apply(Zrich, 1, mean), 
#                        sd.rich = apply(Zrich, 1, sd))
# 
# class(FullGrid)
# plot.grid <- cbind(st_as_sf(FullGrid), plot.dat$mean.rich, plot.dat$sd.rich)
# 
# ## Export predictive map as shapefile
# file_name_Richshape = paste("SDM_RichShape_", "Sharks", ".shp", sep="")
# st_write(plot.grid, file_name_Richshape, append = FALSE)
# 
# ## Plot Richness mean
# file_name_RichMean = "0SharkRichMean.jpeg"
# MeanRichMap <- ggplot() +
#   geom_sf(data = plot.grid, aes(fill = plot.dat.mean.rich), color = NA) +
#   scale_fill_viridis_c("N Spp", na.value = "transparent", option = "A", direction = -1)
# MeanRichMap
# # Export map
# ggsave(filename = paste(file_name_RichMean),
#        plot =  MeanRichMap, width = 10, height = 9, units = 'cm',
#        scale = 2, dpi = 100)
# 
# 
# # Plot Richness SD
# file_name_RichSD = "0SharkRichSD.jpeg"
# SDRichMap <- ggplot() +
#   geom_sf(data = plot.grid, aes(fill = plot.dat.sd.rich), color = NA) +
#   scale_fill_viridis_c("SD", na.value = "transparent", option = "A", direction = -1)
# SDRichMap
# # Export map
# ggsave(filename = paste(file_name_RichSD),
#        plot =  SDRichMap, width = 10, height = 9, units = 'cm',
#        scale = 2, dpi = 100)
# 
# 
# #### Wheighted by ecoregion
# ## Import grid file with ecoregion IDs
# Grid_Ecoregion_RAW <- st_read("Shapes/FullDataGrid_EcoregionID.shp")
# 
# ## Get the IDs of cells within two or more ecoregions
# SplitCellIds <- unique(Grid_Ecoregion_RAW$FID_1[which(duplicated(Grid_Ecoregion_RAW$FID_1))])
# 
# ## Create a matrix to be filled out with ecoregion IDs and max richness
# ecoreg.code <- rep(NA, times = dim(plot.dat)[1])
# plot.grid.ecoreg <- cbind(plot.dat, ecoreg.code)
# 
# 
# ## Loop over cell ids to fill out corresponding ecoregions
# for (j in 1: dim(plot.grid.ecoreg)[1]){
#   
#   ## Grid FID starts with 0
#   CellID <- j - 1
#   
#   if(CellID %in% SplitCellIds == F){
#     if(CellID == 27374 | CellID == 27441){
#       plot.grid.ecoreg[j, 5] <- NA
#     }
#     if(CellID != 27374 & CellID != 27441){
#       plot.grid.ecoreg[j, 5] <- Grid_Ecoregion_RAW$ECO_CODE[Grid_Ecoregion_RAW$FID_1 == CellID]
#     }
#   }
#   
#   if(CellID %in% SplitCellIds == T){
#     SplitCellsData <- Grid_Ecoregion_RAW[Grid_Ecoregion_RAW$FID_1 == CellID,]
#     plot.grid.ecoreg[j, 5] <- SplitCellData$ECO_CODE[SplitCellData$Shape_Area == max(SplitCellData$Shape_Area)]
#   }
#   
#   print(j)
# }
# 
# 
# ### Fill out with max mean richness per ecoregion
# 
# ## Check unique ecoregion codes
# unique(plot.grid.ecoreg$ecoreg.code)
# 
# ## Get max richness per ecorigion id
# max.mean.rich <- tapply(plot.grid.ecoreg$mean.rich, plot.grid.ecoreg$ecoreg.code, max)
# 
# ## Create a vector with ecoregion as factors and then replace levels by richness estimate 
# ecoregFact <- factor(plot.grid.ecoreg$ecoreg.code)
# levels(ecoregFact) <- as.character(max.mean.rich)
# max.mean.rich <- as.numeric(as.character(ecoregFact))
# ecoreg.rich <- plot.grid.ecoreg$mean.rich/max.mean.rich
# 
# # Produce a species distribution map with the weighted richness
# plot.grid <- cbind(st_as_sf(FullGrid), max.mean.rich, ecoreg.rich)
# 
# 
# ## Export predictive map as shapefile
# file_name_WRichshape = paste("SDM_WeightedRichShape_", "Sharks", ".shp", sep="")
# st_write(plot.grid, file_name_WRichshape, append = FALSE)
# 
# ## Plot Weighted richness mean
# file_name_RichMean = "0SharkWeightedRich.jpeg"
# MeanRichMap <- ggplot() +
#   geom_sf(data = plot.grid, aes(fill = ecoreg.rich), color = NA) +
#   scale_fill_viridis_c("", na.value = "transparent", option = "A", direction = -1)
# MeanRichMap
# # Export map
# ggsave(filename = paste(file_name_RichMean),
#        plot =  MeanRichMap, width = 10, height = 9, units = 'cm',
#        scale = 2, dpi = 100)
