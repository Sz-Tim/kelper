# KELPER
# Fit regressions
# Tim Szewczyk

# This script processes the datasets and covariates used to parameterize the
# regional model for *L. hyperborea* in the UK and fits the regressions used in
# the simulations.




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "glue", "lubridate", "tidyverse", "sf", "lme4", "glmmTMB", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
options(mc.cores=4)

# switches
PAR_datasource <- c("MODIS", "POWER")[1]
lmType <- c("lm", "brms")[2]
gridRes <- 0.25

# directories
gis.dir <- "..\\..\\00_gis\\"
data.dir <- "..\\data\\digitized\\"
lm.dir <- paste0("data\\", lmType)
supp.f <- "..\\data\\collab\\collab_all.xlsx"

# datasets
grid.sf <- st_read(glue("data\\grid_{gridRes}_{PAR_datasource}.gpkg")) %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn)
covars.ls <- loadCovariates(loadFile="data\\covar_ls.rds")
data.ls <- compileDatasets(data.dir, supp.f) %>%
  extractCovarsToDatasets(., grid.sf=grid.sf)







########
##-- Define equations

reg.full <- list(
  lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * PAR_atDepth * sstDay_mn * fetch + (1|location)",
  lenSt_to_wtFr.lm="logWtFrond ~ logLenStipe * PAR_atDepth * sstDay_mn * fetch + (1|location)",
  wtFr_to_arFr.lm="logAreaFrond ~ logWtFrond * PAR_atDepth * fetch + sstDay_mn",
  arFr_to_wtFr.lm="logWtFrond ~ logAreaFrond * PAR_atDepth * fetch + sstDay_mn",
  canopyHeight.lm="maxStipeLen ~ sstDay_mn * PAR_atDepth * fetch",
  FAI.lm="FAI ~ PAR_atDepth * sstDay_mn * fetch",
  N_canopy.lm="N ~ PAR_atDepth * sstDay_mn * fetch + (1|location)"
)

# Selected by AIC(c) using the unthinking MuMIn::dredge() 
# -- I think these overfit. They lead to wildly extreme predictions.
# -- Better to base each regression off of previous analyses & suggested
# -- relationships in the literature.
# reg.MODIS <- list(
#   lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * sstDay_mn + (1|location)",
#   lenSt_to_wtFr.lm=paste("logWtFrond ~ logLenStipe + PAR_atDepth + sstDay_mn + fetch +",
#                          "logLenStipe:PAR_atDepth + logLenStipe:fetch + (1|location)"),
#   wtFr_to_arFr.lm=paste("logAreaFrond ~ logWtFrond + PAR_atDepth + sstDay_mn + fetch +",
#                         "logWtFrond:PAR_atDepth + fetch:PAR_atDepth"),
#   arFr_to_wtFr.lm=paste("logWtFrond ~ logAreaFrond + PAR_atDepth + sstDay_mn + fetch +",
#                         "logAreaFrond:PAR_atDepth + fetch:PAR_atDepth"),
#   canopyHeight.lm="maxStipeLen ~ sstDay_mn",
#   FAI.lm="FAI ~ PAR_atDepth + sstDay_mn",
#   N_canopy.lm=c("N ~ PAR_atDepth + fetch + (1|location)", 
#                 "~ PAR_atDepth * sstDay_mn + (1|location)")
# )
# reg.POWER <- list(
#   lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * sstDay_mn + (1|location)",
#   lenSt_to_wtFr.lm=paste("logWtFrond ~ logLenStipe + PAR_atDepth + fetch +",
#                          "logLenStipe:PAR_atDepth + logLenStipe:fetch + (1|location)"),
#   wtFr_to_arFr.lm=paste("logAreaFrond ~ logWtFrond + PAR_atDepth + sstDay_mn + fetch +",
#                         "logWtFrond:PAR_atDepth"),
#   arFr_to_wtFr.lm=paste("logWtFrond ~ logAreaFrond + PAR_atDepth + sstDay_mn + fetch +",
#                         "logAreaFrond:PAR_atDepth"),
#   canopyHeight.lm="maxStipeLen ~ sstDay_mn",
#   FAI.lm="FAI ~ PAR_atDepth + sstDay_mn + fetch",
#   N_canopy.lm=c("N ~ PAR_atDepth + (1|location)", 
#                 "~ PAR_atDepth * sstDay_mn + (1|location)")
# )
# if(PAR_datasource=="MODIS") reg.best <- reg.MODIS
# if(PAR_datasource=="POWER") reg.best <- reg.POWER

reg.best <- list(
  lenSt_to_wtSt.lm="logWtStipe ~ logLenStipe * PAR_atDepth + fetch + ( 1 | location )",
  lenSt_to_wtFr.lm="logWtFrond ~ logLenStipe + sstDay_mn + ( 1 | location )",
  wtFr_to_arFr.lm="logAreaFrond ~ logWtFrond * PAR_atDepth",
  arFr_to_wtFr.lm="logWtFrond ~ logAreaFrond * PAR_atDepth",
  canopyHeight.lm="maxStipeLen ~ sstDay_mn",
  FAI.lm="FAI ~ PAR_atDepth + sstDay_mn + fetch",
  N_canopy.lm=c("N ~ PAR_atDepth * sstDay_mn + fetch + (1|location)",
                "~ PAR_atDepth * sstDay_mn + (1|location)")
)

priors <- c(prior(normal(0, 10), class=b),
            prior(normal(0, 10), class=Intercept),
            prior(cauchy(0, 2), class=sigma),
            prior(cauchy(0, 2), class=sd))
# Prior is still too strong on canopyHeight




########
##-- Munge datasets

reg.dfs <- vector("list", length(reg.full)) %>% setNames(names(reg.full))

reg.dfs$lenSt_to_wtSt.lm <- data.ls$lengthStipe_weightStipe %>%
  mutate(logWtStipe=log(weightStipe), 
         logLenStipe=log(lengthStipe)) %>%
  select(any_of(str_split(reg.full$lenSt_to_wtSt.lm, " ")[[1]]), location)

reg.dfs$lenSt_to_wtFr.lm <- data.ls$lengthStipe_weightFrond %>%
  mutate(logLenStipe=log(lengthStipe), 
         logWtFrond=log(weightFrond)) %>%
  select(any_of(str_split(reg.full$lenSt_to_wtFr.lm, " ")[[1]]), location)

reg.dfs$wtFr_to_arFr.lm <- data.ls$weightFrond_areaFrond %>%
  slice_sample(prop=1) %>% # randomize row order
  group_by(reference, location, depth) %>%
  slice_head(n=20) %>% ungroup %>%
  mutate(logWtFrond=log(weightFrond),
         logAreaFrond=log(areaFrond/1e4)) %>%
  select(any_of(str_split(reg.full$wtFr_to_arFr.lm, " ")[[1]]), location)

reg.dfs$arFr_to_wtFr.lm <- reg.dfs$wtFr_to_arFr.lm %>%
  select(any_of(str_split(reg.full$arFr_to_wtFr.lm, " ")[[1]]), location)

reg.dfs$canopyHeight.lm <- data.ls$depth_maxStipeLen %>% 
  mutate(location="a") %>%
  select(any_of(str_split(reg.full$canopyHeight.lm, " ")[[1]]), location)

reg.dfs$FAI.lm <- data.ls$depth_FAI %>%
  select(any_of(str_split(reg.full$FAI.lm, " ")[[1]]), location)

reg.dfs$N_canopy.lm <- data.ls$lengthStipe_NperSqM %>%
  left_join(., data.ls$lengthStipe_NperSqM %>% 
              group_by(location, PAR_atDepth) %>% 
              summarise(stipeMin=min(lengthStipe), stipeMax=max(lengthStipe),
                        top33=(stipeMax-stipeMin)*2/3 + stipeMin), 
            by=c("PAR_atDepth", "location")) %>%
  mutate(stage=c("canopy", "subcanopy")[1 + (lengthStipe < top33)]) %>%
  group_by(location, reference, habitat, sstDay_mn, PAR_atDepth, fetch, stage) %>%
  summarise(NperSqM=sum(NperSqM)) %>%
  bind_rows(data.ls$depth_NperSqM %>% 
              select(2:4, location, sstDay_mn, PAR_atDepth, fetch) %>%
              rename(canopy=NperSqM, subcanopy=N_subcanopy, recruits=N_recruits) %>%
              pivot_longer(1:3, names_to="stage", values_to="NperSqM")) %>%
  mutate(N=round(NperSqM)) %>%
  ungroup %>%
  filter(stage=="canopy") %>%
  select(any_of(str_split(reg.full$N_canopy.lm, " ")[[1]]), location, N)

reg.dfs <- map(reg.dfs, ~filter(.x, complete.cases(.x)))
reg.dfs_mn_sd <- map(reg.dfs, 
                     ~.x %>% summarise(across(where(is.numeric), 
                                              list(ctr=mean, scl=sd))) %>%
                       pivot_longer(everything(), names_to="Par", values_to="val") %>%
                       mutate(Metric=str_sub(Par, -3, -1),
                              Par=str_sub(Par, 1, -5)) %>%
                       pivot_wider(names_from="Metric", values_from="val"))
dfs.scaled <- map(reg.dfs,
                  ~.x %>% mutate(across(where(is.numeric), 
                                        ~(.x - mean(.x))/sd(.x))))
dfs.scaled$N_canopy.lm$N <- reg.dfs$N_canopy.lm$N



########
##-- Fit regressions

# Model comparison
compareModels <- F

if(compareModels) {
  reg.fit <- vector("list", length(reg.full)) %>% setNames(names(reg.full))
  
  options(na.action="na.fail")
  library(MuMIn)
  reg.fit$lenSt_to_wtSt.lm <- lmer(reg.full$lenSt_to_wtSt.lm, dfs.scaled$lenSt_to_wtSt.lm)
  reg.fit$lenSt_to_wtFr.lm <- lmer(reg.full$lenSt_to_wtFr.lm, dfs.scaled$lenSt_to_wtFr.lm)
  reg.fit$wtFr_to_arFr.lm <- lm(reg.full$wtFr_to_arFr.lm, dfs.scaled$wtFr_to_arFr.lm)
  reg.fit$arFr_to_wtFr.lm <- lm(reg.full$arFr_to_wtFr.lm, dfs.scaled$arFr_to_wtFr.lm)
  reg.fit$canopyHeight.lm <- lm(reg.full$canopyHeight.lm, dfs.scaled$canopyHeight.lm)
  reg.fit$FAI.lm <- lm(reg.full$FAI.lm, dfs.scaled$FAI.lm)
  reg.fit$N_canopy.lm <- glmmTMB(as.formula(reg.best$N_canopy.lm[1]),
                                 ziformula=as.formula(reg.best$N_canopy.lm[2]),
                                 data=dfs.scaled$N_canopy.lm, family=nbinom2)
  reg.compare <- map(reg.fit, dredge)
  map(reg.compare, ~t(cbind(.x[1,])))
  options(na.action="na.omit") 
}




reg.fit <- vector("list", length(reg.best)) %>% setNames(names(reg.best))
if(lmType=="lm") {
  reg.fit$lenSt_to_wtSt.lm <- lmer(reg.best$lenSt_to_wtSt.lm, dfs.scaled$lenSt_to_wtSt.lm)
  reg.fit$lenSt_to_wtFr.lm <- lmer(reg.best$lenSt_to_wtFr.lm, dfs.scaled$lenSt_to_wtFr.lm)
  reg.fit$wtFr_to_arFr.lm <- lm(reg.best$wtFr_to_arFr.lm, dfs.scaled$wtFr_to_arFr.lm)
  reg.fit$arFr_to_wtFr.lm <- lm(reg.best$arFr_to_wtFr.lm, dfs.scaled$arFr_to_wtFr.lm)
  reg.fit$canopyHeight.lm <- lm(reg.best$canopyHeight.lm, dfs.scaled$canopyHeight.lm)
  reg.fit$FAI.lm <- lm(reg.best$FAI.lm, dfs.scaled$FAI.lm)
  reg.fit$N_canopy.lm <- glmmTMB(as.formula(reg.best$N_canopy.lm[1]),
                                 ziformula=as.formula(reg.best$N_canopy.lm[2]),
                                 data=dfs.scaled$N_canopy.lm, family=nbinom2)
} else if(lmType=="brms") {
  reg.fit$lenSt_to_wtSt.lm <- brm(reg.best$lenSt_to_wtSt.lm, dfs.scaled$lenSt_to_wtSt.lm, 
                                  prior=priors, family=gaussian())
  reg.fit$lenSt_to_wtFr.lm <- brm(reg.best$lenSt_to_wtFr.lm, dfs.scaled$lenSt_to_wtFr.lm, 
                                  prior=priors[-4,], family=gaussian())
  reg.fit$wtFr_to_arFr.lm <- brm(reg.best$wtFr_to_arFr.lm, dfs.scaled$wtFr_to_arFr.lm, 
                                 prior=priors[-4,], family=gaussian())
  reg.fit$arFr_to_wtFr.lm <- brm(reg.best$arFr_to_wtFr.lm, dfs.scaled$arFr_to_wtFr.lm, 
                                 prior=priors[-4,], family=gaussian())
  reg.fit$canopyHeight.lm <- brm(reg.best$canopyHeight.lm, dfs.scaled$canopyHeight.lm, 
                                 prior=priors[-4,], 
                                 family=gaussian())
  reg.fit$FAI.lm <- brm(reg.best$FAI.lm, dfs.scaled$FAI.lm, 
                        prior=priors[-4,], family=gaussian())
  reg.fit$N_canopy.lm <- brm(bf(reg.best$N_canopy.lm[1],
                                paste("zi", reg.best$N_canopy.lm[2])), 
                             data=dfs.scaled$N_canopy.lm, 
                             prior=priors[-3,], family=zero_inflated_negbinomial())
}




########
##-- Save output

saveRDS(reg.fit, glue::glue("{lm.dir}_{PAR_datasource}_fits.rds"))
saveRDS(dfs.scaled, glue::glue("data\\dfs_scaled_{PAR_datasource}.rds"))
saveRDS(reg.dfs_mn_sd, glue::glue("data\\dfs_mn_sd_{PAR_datasource}.rds"))

