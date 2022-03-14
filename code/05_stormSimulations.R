# KELPER
# Sensitivity analysis
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths with a
# survival and/or loss rate dependent on annual storm intensity




########
##-- set up

# libraries and local functions
pkgs <- c("glue", "tidyverse", "sf", "brms", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())


# directories & files
data.dir <- "data\\raw\\digitized\\"
supp.f <- "data\\raw\\collab\\collab_all.xlsx"
out.dir <- "out\\storms\\"

# switches & settings
gridRes <- c(0.1, 0.25)[2]
options(mc.cores=12)
pars.sim <- list(depths=c(2, 5, 15), 
                 tskip=0,
                 nSim=1, 
                 stochParams=F,
                 landscape="static")

# load files
grid.sf <- st_read(glue("data\\grid_{gridRes}_MODIS.gpkg"))
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn) %>%
  select(id, SST, KD, PAR, fetch, fetchCat)
data.ls <- compileDatasets(data.dir, supp.f)
lm.fit <- readRDS(glue("data\\fits_{gridRes}.rds"))
lm.mnsd <- readRDS(glue("data\\dfs_mn_sd_{gridRes}.rds"))
surv.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageTo=="dead") %>%
  mutate(survRate=1-rate_mn,
         exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))
fecund.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageFrom=="canopy" & stageTo=="recruits") %>%
  mutate(exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))

# set sim length, storm intensity
pars.sim$tmax <- nrow(data.ls$year_stormIndex)
pars.sim$storms <- data.ls$year_stormIndex$stormIndex


# simulate landscapes if needed
set.seed(789)
if(pars.sim$landscape=="dynamic") {
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
} else {
  grid.sim <- NULL
}






########
##-- run simulations


cl <- makeCluster(2, outfile="temp\\sim_out.txt")
obj.exclude <- c("data.ls", "grid.sf")
obj.include <- ls()

cat("Exporting cluster.")
clusterExport(cl, obj.include[-match(obj.exclude, obj.include)])
Sys.sleep(5)
cat("Exported. Starting parallel runs.")
out.ls <- parLapply(cl, X=1:nrow(grid.i), fun=simDepthsWithinCell,
                    grid.i=grid.i, grid.sim=grid.sim, gridRes=gridRes, 
                    pars.sim=pars.sim, surv.df=surv.df, fecund.df=fecund.df, 
                    lm.fit=lm.fit, lm.mnsd=lm.mnsd)
stopCluster(cl)
