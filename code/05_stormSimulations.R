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
data.dir <- "..\\data\\digitized\\"
supp.f <- "..\\data\\collab\\collab_all.xlsx"
out.dir <- "..\\out\\storms\\"

# switches & settings
gridRes <- c(0.1, 0.25)[2]
options(mc.cores=12)
pars.sim <- list(depths=c(2, 15), 
                 tmax=50,
                 tskip=10,
                 nSim=1, 
                 stochParams=F,
                 landscape="static")
pars.sim$storms <- rexp(pars.sim$tmax, 3)

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




