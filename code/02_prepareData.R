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
lmType <- c("lm", "brms")[2]

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
  lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * PAR_atDepth + ( 1 | location )",
  lenSt_to_wtFr.lm="logWtFrond ~ logLenStipe + sstDay_mn + ( 1 | location )",
  # lenSt_to_wtTot.lm="logWtTotal ~ logLenStipe * PAR_atDepth + sstDay_mn + ( 1 | location )",
  # wtSt_to_wtFr.lm="logWtFrond ~ logWtStipe * propClear + ( 1 | location )",
  wtFr_to_arFr.lm="logAreaFrond ~ logWtFrond * PARdepth + ( 1 | location )",
  arFr_to_wtFr.lm="logWtFrond ~ logAreaFrond * PARdepth + ( 1 | location )",
  canopyHeight.lm="maxStipeLen ~ sstDay_mn",
  FAI.lm="logFAI ~ logPAR * fetch",
  # N_stage.lm="logN ~ stage * PAR_atDepth * sstDay_mn + ( 1 | location )",
  N_canopy.lm="logN ~ PAR_atDepth * sstDay_mn + ( 1 | location )"
)

reg.full <- list(
  lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * PAR_atDepth + ( 1 | location )",
  lenSt_to_wtFr.lm="logWtFrond ~ logLenStipe * sstDay_mn + ( 1 | location )",
  # lenSt_to_wtTot.lm="logWtTotal ~ logLenStipe * PAR_atDepth * sstDay_mn + ( 1 | location )",
  # wtSt_to_wtFr.lm="logWtFrond ~ logWtStipe * propClear + ( 1 | location )",
  wtFr_to_arFr.lm="logAreaFrond ~ logWtFrond * PARdepth + ( 1 | location )",
  arFr_to_wtFr.lm="logWtFrond ~ logAreaFrond * PARdepth + ( 1 | location )",
  canopyHeight.lm="maxStipeLen ~ sstDay_mn",
  FAI.lm="logFAI ~ logPAR * fetch",
  # N_stage.lm="logN ~ stage * PAR_atDepth * sstDay_mn * fetch + ( 1 | location )",
  N_canopy.lm="logN ~ PAR_atDepth * sstDay_mn * fetch + ( 1 | location )"
)

priors <- c(prior(normal(0, 5), class=b),
            prior(normal(0, 5), class=Intercept),
            prior(cauchy(0, 1), class=sigma))





########
##-- Munge datasets

reg.dfs <- vector("list", length(reg.forms)) %>% setNames(names(reg.forms))

reg.dfs$lenSt_to_wtSt.lm <- data.ls$lengthStipe_weightStipe %>%
  mutate(logWtStipe=log(weightStipe), 
         logLenStipe=log(lengthStipe)) %>%
  select(any_of(str_split(reg.full$lenSt_to_wtSt.lm, " ")[[1]]), location)

reg.dfs$lenSt_to_wtFr.lm <- data.ls$lengthStipe_weightFrond %>%
  mutate(logLenStipe=log(lengthStipe), 
         logWtFrond=log(weightFrond)) %>%
  select(any_of(str_split(reg.full$lenSt_to_wtFr.lm, " ")[[1]]), location)

# reg.dfs$lenSt_to_wtTot.lm <- data.ls$lengthStipe_weightTotal %>%
#   mutate(logLenStipe=log(lengthStipe), 
#          logWtTotal=log(weightTotal)) %>%
#   select(any_of(str_split(reg.full$lenSt_to_wtTot.lm, " ")[[1]]), location)

# reg.dfs$wtSt_to_wtFr.lm <- data.ls$weightStipe_weightFrond %>% 
#   mutate(propClear=as.numeric(habitat=="clearing"),
#          logWtFrond=log(weightFrond), 
#          logWtStipe=log(weightStipe)) %>%
#   select(any_of(str_split(reg.full$wtSt_to_wtFr.lm, " ")[[1]]), location)

reg.dfs$wtFr_to_arFr.lm <- data.ls$weightFrond_areaFrond %>%
  mutate(logWtFrond=log(weightFrond),
         logAreaFrond=log(areaFrond/1e4)) %>%
  rename(PARdepth=PAR_atDepth) %>%
  select(any_of(str_split(reg.full$wtFr_to_arFr.lm, " ")[[1]]), location)

reg.dfs$arFr_to_wtFr.lm <- reg.dfs$wtFr_to_arFr.lm %>%
  select(any_of(str_split(reg.full$arFr_to_wtFr.lm, " ")[[1]]), location)

reg.dfs$canopyHeight.lm <- data.ls$depth_maxStipeLen %>% 
  mutate(location="a") %>%
  select(any_of(str_split(reg.full$canopyHeight.lm, " ")[[1]]), location)

reg.dfs$FAI.lm <- data.ls$depth_FAI %>%
  mutate(logFAI=log(FAI),
         logPAR=log(PAR_atDepth)) %>%
  select(any_of(str_split(reg.full$FAI.lm, " ")[[1]]), location)

# reg.dfs$N_stage.lm <- data.ls$lengthStipe_NperSqM %>%
#   left_join(., data.ls$lengthStipe_NperSqM %>% 
#               group_by(location, PAR_atDepth) %>% 
#               summarise(stipeMin=min(lengthStipe), stipeMax=max(lengthStipe),
#                         top33=(stipeMax-stipeMin)*2/3 + stipeMin), 
#             by=c("PAR_atDepth", "location")) %>%
#   mutate(stage=c("canopy", "subcanopy")[1 + (lengthStipe < top33)]) %>%
#   group_by(location, sstDay_mn, PAR_atDepth, fetch, stage) %>%
#   summarise(NperSqM=sum(NperSqM)) %>%
#   bind_rows(data.ls$depth_NperSqM %>% 
#               select(2:4, location, sstDay_mn, PAR_atDepth, fetch) %>%
#               rename(canopy=NperSqM, subcanopy=N_subcanopy, recruits=N_recruits) %>%
#               pivot_longer(1:3, names_to="stage", values_to="NperSqM")) %>%
#   mutate(logN=log(NperSqM+1)) %>%
#   ungroup %>%
#   select(any_of(str_split(reg.full$N_stage.lm, " ")[[1]]), location)

reg.dfs$N_canopy.lm <- data.ls$lengthStipe_NperSqM %>%
  left_join(., data.ls$lengthStipe_NperSqM %>% 
              group_by(location, PAR_atDepth) %>% 
              summarise(stipeMin=min(lengthStipe), stipeMax=max(lengthStipe),
                        top33=(stipeMax-stipeMin)*2/3 + stipeMin), 
            by=c("PAR_atDepth", "location")) %>%
  mutate(stage=c("canopy", "subcanopy")[1 + (lengthStipe < top33)]) %>%
  group_by(location, sstDay_mn, PAR_atDepth, fetch, stage) %>%
  summarise(NperSqM=sum(NperSqM)) %>%
  bind_rows(data.ls$depth_NperSqM %>% 
              select(2:4, location, sstDay_mn, PAR_atDepth, fetch) %>%
              rename(canopy=NperSqM, subcanopy=N_subcanopy, recruits=N_recruits) %>%
              pivot_longer(1:3, names_to="stage", values_to="NperSqM")) %>%
  mutate(logN=log(NperSqM+1)) %>%
  ungroup %>%
  filter(stage=="canopy") %>%
  select(any_of(str_split(reg.full$N_canopy.lm, " ")[[1]]), location)

reg.dfs <- map(reg.dfs, ~filter(.x, complete.cases(.x)))



########
##-- Fit regressions

compareModels <- F
if(lmType=="lm") {
  options(na.action="na.fail")
  reg.fit <- vector("list", length(reg.forms)) %>% setNames(names(reg.forms))
  for(i in seq_along(reg.forms)) {
    df_i <- reg.dfs[[i]]
    if(compareModels) {
      form_i <- reg.full[[i]]
    } else {
      form_i <- reg.forms[[i]]
    }
    if(grepl("location", form_i)) {
      full_i <- lmer(form_i, df_i)
    } else {
      full_i <- lm(form_i, df_i)
    }
    if(compareModels) {
      reg.fit[[i]] <- MuMIn::dredge(full_i)
    } else {
      reg.fit[[i]] <- full_i
    }
  }
  options(na.action="na.omit")
} else if(lmType=="brms") {
  reg.fit <- map2(reg.forms, reg.dfs, ~brm(.x, data=.y, 
                                           control=list(adapt_delta=0.9)))
}




########
##-- Save output

saveRDS(reg.fit, paste0(lm.dir, "_fits.rds"))


