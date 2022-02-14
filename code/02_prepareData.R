# KELPER
# Fit regressions
# Tim Szewczyk

# This script processes the datasets and covariates used to parameterize the
# regional model for *L. hyperborea* in the UK and fits the regressions used in
# the simulations.




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "tidyverse", "sf", "lme4", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
options(mc.cores=4)

# switches
PAR_datasource <- c("MODIS", "POWER")[1]
lmType <- c("lm", "brms")[1]

# directories
gis.dir <- "..\\..\\00_gis\\"
data.dir <- "..\\data\\digitized\\"
lm.dir <- paste0("data\\", lmType)
supp.f <- "..\\data\\collab\\collab_all.xlsx"

# maximum extent for covariates
UK_bbox <- st_bbox(c(xmin=-11, xmax=2, ymin=49, ymax=61), crs=st_crs(4326))

# datasets
covars.ls <- loadCovariates(gis.dir, UK_bbox, loadFile="data\\covar_ls.rds")
data.ls <- compileDatasets(data.dir, supp.f) %>%
  extractCovarsToDatasets(., covars.ls, PAR_datasource)








########
##-- Define equations

reg.forms <- list(
  lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * PAR_atDepth + (1|location)",
  lenSt_to_wtFr.lm="logWtFrond ~ logLenStipe + sstDay_mn + (1|location)",
  lenSt_to_wtTot.lm="logWtTotal ~ logLenStipe * PAR_atDepth + sstDay_mn + (1|location)",
  wtSt_to_wtFr.lm="logWtFrond ~ logWtStipe * propClear + (1|location)",
  wtFr_to_arFr.lm="logAreaFrond ~ logWtFrond * PARdepth + (1|location)",
  arFr_to_wtFr.lm="logWtFrond ~ logAreaFrond * PARdepth + (1|location)",
  canopyHeight.lm="maxStipeLen ~ sstDay_mn",
  FAI.lm="logFAI ~ logPAR + fetch",
  N_stage.lm="logN ~ stage * PAR_atDepth * sstDay_mn * fetch",
  N_canopy.lm="logN ~ PAR_atDepth * sstDay_mn * fetch"
)

priors <- c(prior(normal(0, 5), class=b),
            prior(normal(0, 5), class=Intercept),
            prior(cauchy(0, 1), class=sigma))





########
##-- Munge datasets

reg.dfs <- vector("list", length(reg.forms)) %>% setNames(names(reg.forms))

reg.dfs[[names(reg.forms)[1]]] <- data.ls$lengthStipe_weightStipe %>%
  mutate(logWtStipe=log(weightStipe), 
         logLenStipe=log(lengthStipe))

reg.dfs[[names(reg.forms)[2]]] <- data.ls$lengthStipe_weightFrond %>%
  mutate(logLenStipe=log(lengthStipe), 
         logWtFrond=log(weightFrond))

reg.dfs[[names(reg.forms)[3]]] <- data.ls$lengthStipe_weightTotal %>%
  mutate(logLenStipe=log(lengthStipe), 
         logWtTotal=log(weightTotal))

reg.dfs[[names(reg.forms)[4]]] <- data.ls$weightStipe_weightFrond %>% 
  mutate(propClear=as.numeric(habitat=="clearing"),
         logWtFrond=log(weightFrond), 
         logWtStipe=log(weightStipe))

reg.dfs[[names(reg.forms)[5]]] <- data.ls$weightFrond_areaFrond %>%
  mutate(logWtFrond=log(weightFrond),
         logAreaFrond=log(areaFrond/1e4)) %>%
  rename(PARdepth=PAR_atDepth)

reg.dfs[[names(reg.forms)[6]]] <- reg.dfs[[names(reg.forms)[5]]]

reg.dfs[[names(reg.forms)[7]]] <- data.ls$depth_maxStipeLen %>% mutate(location="a")

reg.dfs[[names(reg.forms)[8]]] <- data.ls$depth_FAI %>%
  mutate(logFAI=log(FAI),
         logPAR=log(PAR_atDepth))

reg.dfs[[names(reg.forms)[9]]] <- data.ls$lengthStipe_NperSqM %>%
  left_join(., data.ls$lengthStipe_NperSqM %>% 
              group_by(location, PAR_atDepth) %>% 
              summarise(stipeMin=min(lengthStipe), stipeMax=max(lengthStipe),
                        top33=(stipeMax-stipeMin)*2/3 + stipeMin), 
            by=c("PAR_atDepth", "location")) %>%
  mutate(stage=c("canopy", "subcanopy")[1 + (lengthStipe < top33)]) %>%
  group_by(location, sstDay_mn, PAR_atDepth, fetch, stage) %>%
  summarise(NperSqM=sum(NperSqM)) %>%
  bind_rows(data.ls$depth_NperSqM %>% 
              select(2:4, sstDay_mn, PAR_atDepth, fetch) %>%
              rename(canopy=NperSqM, subcanopy=N_subcanopy, recruits=N_recruits) %>%
              pivot_longer(1:3, names_to="stage", values_to="NperSqM")) %>%
  mutate(logN=log(NperSqM+1))

reg.dfs[[names(reg.forms)[10]]] <- reg.dfs[[names(reg.forms)[9]]] %>%
  filter(stage=="canopy")





########
##-- Fit regressions

if(lmType=="lm") {
  reg.fit <- vector("list", length(reg.forms)) %>% setNames(names(reg.forms))
  for(i in seq_along(reg.forms)) {
    if(grepl("location", reg.forms[[i]])) {
      reg.fit[[i]] <- lmer(reg.forms[[i]], reg.dfs[[i]])
    } else {
      reg.fit[[i]] <- lm(reg.forms[[i]], reg.dfs[[i]]) 
    }
  }
} else if(lmType=="brms") {
  reg.fit <- map2(reg.forms, reg.dfs, ~brm(.x, data=.y))
}




########
##-- Save output

saveRDS(reg.fit, paste0(lm.dir, "_fits.rds"))


