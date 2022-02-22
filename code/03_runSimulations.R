# KELPER
# Run simulations
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "lme4", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())

# directories
data.dir <- "..\\data\\digitized\\"
supp.f <- "..\\data\\collab\\collab_all.xlsx"

# switches & settings
lmType <- c("lm", "brms")[2]
PAR_datasource <- c("MODIS", "POWER")[1]
gridRes <- 0.25
dynamicLandscape <- F
depths <- c(2, 5, 10, 15)
tmax <- 30
nSim <- 2
options(mc.cores=12)

# load files
grid.sf <- st_read(glue("data\\grid_{gridRes}_{PAR_datasource}.gpkg"))
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn)
data.ls <- compileDatasets(data.dir, supp.f)
covars.ls <- loadCovariates(loadFile="data\\covar_ls.rds")
lm.fit <- readRDS(glue("data\\{lmType}_{PAR_datasource}_fits.rds"))
lm.mnsd <- readRDS(glue("data\\dfs_mn_sd_{PAR_datasource}.rds"))
surv.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageTo=="dead") %>%
  mutate(survRate=1-rate,
         exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))
fecund.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageFrom=="canopy" & stageTo=="recruits") %>%
  mutate(exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))

# simulate landscapes if needed
set.seed(789)
if(dynamicLandscape) {
  covars.full <- loadCovariates_full(loadFile="data\\covarFull_ls.rds")
  grid.sim <- grid.sf %>% select(id) %>%
    simulateLandscape(covars.full$sstDayGrow, tmax, "SST") %>%
    full_join(., 
              simulateLandscape(grid.sf %>% select(id), 
                                covars.full$KD_mn, tmax, "KD") %>%
                st_drop_geometry()) %>%
    full_join(., 
              simulateLandscape(grid.sf %>% select(id), 
                                covars.full$PAR_mn, tmax, "PAR") %>%
                st_drop_geometry()) %>%
    left_join(., grid.sf %>% st_drop_geometry() %>% select(id, fetch, fetchCat))
  rm(covars.full)
}




########
##-- run simulations


library(doSNOW); library(foreach)
cl <- makeCluster(12, methods=F, outfile="temp\\sim_out.txt")
registerDoSNOW(cl)
obj.exclude <- c("covars.ls", "data.ls", "grid.sf")
obj.include <- ls()
out.ls <- foreach(i=1:nrow(grid.i), .combine="c", 
                  .export=obj.include[-match(obj.exclude, obj.include)],
                  .noexport=obj.exclude) %dopar% {
  library(dbplyr); library(brms); library(tidyverse); library(glue)                  
  walk(dir("code", "^00.*R", full.names=T), source)
  if(dynamicLandscape) {
    cell.env <- grid.sim %>% st_drop_geometry() %>% filter(id==grid.i$id[i])
  } else {
    cell.env <- grid.i[i,]
  }
  pop.ls <- mass.ls <- vector("list", length(depths))
  for(j in 1:length(depths)) {
    par_i <- setParameters(
      tmax=tmax, 
      survRate=filter(surv.df, exposure==cell.env$fetchCat[1])$survRate^(1/2),
      settlementRateBg=filter(fecund.df, exposure==cell.env$fetchCat[1])$rate,
      extraPars=list(
        depth=depths[j],
        env=cell.env,
        lenSt_to_wtSt=lm.fit$lenSt_to_wtSt.lm,
        lenSt_to_wtFr=lm.fit$lenSt_to_wtFr.lm,
        wtFr_to_arFr=lm.fit$wtFr_to_arFr.lm,
        arFr_to_wtFr=lm.fit$arFr_to_wtFr.lm,
        N_canopy=lm.fit$N_canopy.lm,
        FAI=lm.fit$FAI.lm,
        canopyHeight=lm.fit$canopyHeight.lm,
        sc.df=lm.mnsd)
    )
    out <- map(1:nSim, ~simulatePopulation(par_i, lmType=lmType, ndraws=50))
    pop.ls[[j]] <- imap_dfr(out, 
                            ~tibble(sim=.y, 
                                    year=rep(1:par_i$tmax, 3),
                                    month=rep(c(1,6,7), each=par_i$tmax),
                                    date=ymd(glue("{year}-{month}-01")),
                                    PAR_atDepth=rep(.x$PAR, 3),
                                    K_N=rep(.x$K_N, 3),
                                    K_FAI=rep(.x$K_FAI, 3),
                                    FAI=c(.x$FAI[3,,]),
                                    N.recruits=c(.x$N[1,,]),
                                    N.subcanopy=c(.x$N[2,,]),
                                    N.canopy=c(.x$N[3,,]),
                                    kappa_FAI=c(.x$kappa[,,1]),
                                    kappa_N=c(.x$kappa[,,2])) %>%
                              pivot_longer(contains("N."), names_to="stage", values_to="N") %>%
                              mutate(stage=factor(str_sub(stage, 3, -1), 
                                                  levels=c("canopy", "subcanopy", "recruits")))) %>%
      mutate(id=i) %>%
      left_join(., par_i$env) %>%
      mutate(depth=par_i$depth)
    mass.ls[[j]] <- imap_dfr(out,
                             ~tibble(sim=.y,
                                     year=rep(1:par_i$tmax, 3),
                                     month=rep(c(1,6,7), each=par_i$tmax),
                                     date=ymd(glue("{year}-{month}-01")),
                                     PAR_atDepth=rep(.x$PAR, 3),
                                     FAI=c(.x$FAI[3,,]),
                                     biomass=c(.x$biomass),
                                     kappa_FAI=c(.x$kappa[,,1]),
                                     kappa_N=c(.x$kappa[,,2]))) %>%
      mutate(id=i) %>%
      left_join(., par_i$env) %>%
      mutate(depth=par_i$depth)
  }
  sim.info <- glue("{str_pad(i,3,'left','0')}_{gridRes}_{lmType}",
                   "_{ifelse(dynamicLandscape,'dynamic','static')}")
  saveRDS(do.call(rbind, pop.ls), glue("out\\pop_{sim.info}.rds"))
  saveRDS(do.call(rbind, mass.ls), glue("out\\mass_{sim.info}.rds"))
  i
}
stopCluster(cl)


