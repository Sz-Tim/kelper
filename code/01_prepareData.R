# KELPER
# Empirical datasets and covariates
# Tim Szewczyk

# This script processes the datasets and covariates used to parameterize the regional model for *L. hyperborea* in the UK.


########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "tidyverse", "sf", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
options(mc.cores=4)

# switches
PAR_datasource <- c("MODIS", "POWER")[1]
lmType <- c("lm", "brms")[1]
refitRegressions <- T
gridRes <- 0.5 # currently in arc-seconds
remakeGrid <- T

# directories
gis.dir <- "..\\..\\00_gis\\"
data.dir <- "..\\data\\digitized\\"
lm.dir <- paste0("data\\", lmType)

# maximum extent for covariates
UK_bbox <- st_bbox(c(xmin=-11, xmax=2, ymin=49, ymax=61), crs=st_crs(4326))

# datasets
covars.ls <- loadCovariates(gis.dir, UK_bbox, saveFile="data\\covar_ls.rds")
data.ls <- compileDatasets(data.dir) %>%
  extractCovarsToDatasets(., covars.ls, PAR_datasource)

if(remakeGrid) {
  grid.domain <- UK_bbox %>%
    st_make_grid(cellsize=c(gridRes, gridRes)) %>%
    st_sf(id=1:length(.)) %>%
    extractCovarsToGrid(., covars.ls, PAR_datasource) %>%
    filter(!is.na(KD_mn) & !is.na(PAR_surface) & !is.na(fetch)) %>%
    mutate(id=row_number())
  st_write(grid.domain, paste0("data\\grid_", gridRes, ".shp"))
} else {
  grid.domain <- st_read(paste0("data\\grid_", gridRes, ".shp"))
}





########
##-- Regressions
# 1) StipeWeight (g) ~ StipeLength (mm)
# 2) FrondWeight (g) ~ StipeWeight (g) * kappa
# 3) FrondArea (cm2->m2) ~ FrondWeight (g) * PAR[depth]
# 4) FrondWeight (g) ~ FrondArea (cm2->m2) * PAR[depth]
# 5) CanopyHeight (mm) ~ SST
# 6) FrondAreaIndex ~ PAR[depth]
# 7) Density (N/m2) ~ PAR[depth]

if(lmType=="brms") {
  priors <- prior(normal(0,5))
  # 1)
  lenStipe_wtStipe <- data.ls$lengthStipe_weightStipe %>%
    mutate(logWtStipe=log(weightStipe), 
           logLenStipe=log(lengthStipe)) %>%
    brm(logWtStipe ~ logLenStipe, prior=priors,  data=.)
  
  # 2)
  wtStipe_wtFrond_kappa <- data.ls$weightStipe_weightFrond %>% 
    mutate(propClear=as.numeric(habitat=="clearing"),
           logWtFrond=log(weightFrond), 
           logWtStipe=log(weightStipe)) %>%
    brm(logWtFrond ~ logWtStipe*propClear, prior=priors, data=.)
  
  # 3)
  wtFrond_areaFrond_PARdepth <- data.ls$weightFrond_areaFrond %>%
    mutate(logWtFrond=log(weightFrond),
           logAreaFrond=log(areaFrond/1e4)) %>%
    brm(logAreaFrond ~ logWtFrond*PARdepth, prior=priors, data=.)
  
  # 4)
  areaFrond_wtFrond_PARdepth <- data.ls$weightFrond_areaFrond %>%
    mutate(logWtFrond=log(weightFrond),
           logAreaFrond=log(areaFrond/1e4)) %>%
    brm(logWtFrond ~ logAreaFrond*PARdepth, prior=priors, data=.)
  
  # 5)
  sst_canopyHeight <- data.ls$depth_maxStipeLen %>%
    brm(maxStipeLen ~ sstDay_mn, prior=priors, data=.)
  
  # 6)
  PARdepth_FAI <- data.ls$depth_FAI %>%
    mutate(logFAI=log(FAI),
           logPAR=log(PAR_atDepth)) %>%
    brm(logFAI ~ logPAR, prior=priors, data=.)
  
  # 7)
  PARdepth_N <- data.ls$lengthStipe_NperSqM %>%
    left_join(., data.ls$lengthStipe_NperSqM %>% 
                group_by(PAR_atDepth) %>% 
                summarise(stipeMin=min(lengthStipe), stipeMax=max(lengthStipe),
                          top33=(stipeMax-stipeMin)*2/3 + stipeMin), 
              by="PAR_atDepth") %>%
    mutate(stage=c("canopy", "subcanopy")[1 + (lengthStipe < top33)]) %>%
    group_by(PAR_atDepth, stage) %>%
    summarise(NperSqM=sum(NperSqM)) %>%
    mutate(logN=log(NperSqM)) %>%
    brm(logN ~ PAR_atDepth*stage, prior=priors, data=.)
  
} else if(lmType=="lm") {
  # 1)
  lenStipe_wtStipe <- data.ls$lengthStipe_weightStipe %>%
    mutate(logWtStipe=log(weightStipe), 
           logLenStipe=log(lengthStipe)) %>%
    lm(logWtStipe ~ logLenStipe,  data=.)
  
  # 2)
  wtStipe_wtFrond_kappa <- data.ls$weightStipe_weightFrond %>% 
    mutate(propClear=as.numeric(habitat=="clearing"),
           logWtFrond=log(weightFrond), 
           logWtStipe=log(weightStipe)) %>%
    lm(logWtFrond ~ logWtStipe*propClear, data=.)
  
  # 3)
  wtFrond_areaFrond_PARdepth <- data.ls$weightFrond_areaFrond %>%
    mutate(logWtFrond=log(weightFrond),
           logAreaFrond=log(areaFrond/1e4)) %>%
    lm(logAreaFrond ~ logWtFrond*PARdepth, data=.)
  
  # 4)
  areaFrond_wtFrond_PARdepth <- data.ls$weightFrond_areaFrond %>%
    mutate(logWtFrond=log(weightFrond),
           logAreaFrond=log(areaFrond/1e4)) %>%
    lm(logWtFrond ~ logAreaFrond*PARdepth, data=.)
  
  # 5)
  sst_canopyHeight <- data.ls$depth_maxStipeLen %>%
    lm(maxStipeLen ~ sstDay_mn, data=.)
  
  # 6)
  PARdepth_FAI <- data.ls$depth_FAI %>%
    mutate(logFAI=log(FAI),
           logPAR=log(PAR_atDepth)) %>%
    lm(logFAI ~ logPAR, data=.)
  
  # 7)
  PARdepth_N <- data.ls$lengthStipe_NperSqM %>%
    left_join(., data.ls$lengthStipe_NperSqM %>% 
                group_by(PAR_atDepth) %>% 
                summarise(stipeMin=min(lengthStipe), stipeMax=max(lengthStipe),
                          top33=(stipeMax-stipeMin)*2/3 + stipeMin), 
              by="PAR_atDepth") %>%
    mutate(stage=c("canopy", "subcanopy")[1 + (lengthStipe < top33)]) %>%
    group_by(PAR_atDepth, stage) %>%
    summarise(NperSqM=sum(NperSqM)) %>%
    mutate(logN=log(NperSqM)) %>%
    lm(logN ~ PAR_atDepth*stage, data=.)
}

saveRDS(lenStipe_wtStipe, paste0(lm.dir, "_lenStipe_wtStipe.rds"))
saveRDS(wtStipe_wtFrond_kappa, paste0(lm.dir, "_wtStipe_wtFrond_kappa.rds"))
saveRDS(wtFrond_areaFrond_PARdepth, paste0(lm.dir, "_wtFrond_areaFrond_PARdepth.rds"))
saveRDS(areaFrond_wtFrond_PARdepth, paste0(lm.dir, "_areaFrond_wtFrond_PARdepth.rds")) 
saveRDS(sst_canopyHeight, paste0(lm.dir, "_sst_canopyHeight.rds")) 
saveRDS(PARdepth_FAI, paste0(lm.dir, "_PARdepth_FAI_K.rds"))
saveRDS(PARdepth_N, paste0(lm.dir, "_PARdepth_N_K.rds"))

