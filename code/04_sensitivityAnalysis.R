# KELPER
# Sensitivity analysis
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths for
# combinations of parameter values




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())


# directories & files
data.dir <- "..\\data\\digitized\\"
supp.f <- "..\\data\\collab\\collab_all.xlsx"
sens.dir <- "out\\sensitivity\\"
parSets.f <- "000_paramSets_exposure_"

# switches & settings
rerun <- T
reanalyse <- T
gridRes <- c(0.1, 0.25)[2]
options(mc.cores=12)
pars.sens <- list(nParDraws=1e2,
                  depths=c(2),
                  tmax=20,
                  tskip=10,
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







########
##-- assign parameter values

# parameter ranges: approximately mean +- 2 sd
loss_mnPrec <- c(0.2340987, 90.1752) # mn, prec for beta distribution
rate_stipe <- cbind(c(194, 195, 58), c(23, 21, 10)) # mn, sd by stage
rate_frond <- cbind(c(1787, 2299, 3979), c(179, 397, 417))/1e4
par.rng <- list(surv=surv.df %>%
                  mutate(stage=stageFrom,
                         valMin=1-pmin(1, pmax(0, rate_mn+2*rate_sd)),
                         valMax=1-pmin(1, pmax(0, rate_mn-2*rate_sd))) %>%
                  select(stage, valMin, valMax, exposure),
                settlement=fecund.df %>%
                  mutate(valMin=pmax(0, rate_mn-2*rate_sd),
                         valMax=pmax(0, rate_mn+2*rate_sd)) %>%
                  select(valMin, valMax, exposure),
                growStipe=tibble(stage=c("recruits", "subcanopy", "canopy"),
                                 valMin=apply(rate_stipe, 1, function(x) x[1]-2*x[2]),
                                 valMax=apply(rate_stipe, 1, function(x) x[1]+2*x[2])),
                growFrond=tibble(stage=c("recruits", "subcanopy", "canopy"),
                                 valMin=apply(rate_frond, 1, function(x) x[1]-2*x[2]),
                                 valMax=apply(rate_frond, 1, function(x) x[1]+2*x[2])),
                loss=tibble(valMin=qbeta(0.025, prod(loss_mnPrec), (1-loss_mnPrec[1])*loss_mnPrec[2]),
                            valMax=qbeta(0.975, prod(loss_mnPrec), (1-loss_mnPrec[1])*loss_mnPrec[2]))) %>%
  imap_dfr(., ~.x %>% mutate(param=.y))

if(rerun) {
  parSets <- map_dfr(1:pars.sens$nParDraws, 
                     ~par.rng %>% mutate(parDraw=.x, 
                                         val=runif(n(), valMin, valMax))) %>%
    select(parDraw, param, stage, exposure, val) %>%
    group_by(exposure) %>% group_split()
  walk2(parSets, c(1:3, "NA"), ~write_csv(.x, glue("{sens.dir}{parSets.f}{.y}.csv")))
}









########
##-- run simulations

if(rerun) {
  library(parallel)
  cl <- makeCluster(4, outfile="temp\\sensitivity_out.txt")
  obj.exclude <- c("data.ls", "grid.sf", "surv.df", "fecund.df")
  obj.include <- ls()
  
  cat("Exporting cluster.")
  clusterExport(cl, obj.include[-match(obj.exclude, obj.include)])
  Sys.sleep(5)
  cat("Exported. Starting parallel runs.")
  out.ls <- parLapply(cl, X=1:4, fun=simSensitivityDepthsWithinCell,
                      grid.i=grid.i, gridRes=gridRes, pars.sens=pars.sens, 
                      lm.fit=lm.fit, lm.mnsd=lm.mnsd, parSets=parSets)
  stopCluster(cl)
}








########
##-- analyse output

if(reanalyse) {
  meta.cols <- c("sim", "month", "stage", "id", "parDraw", "depth",
                 "PAR_atDepth", "K_N", "K_FAI")
  pop.f <- dir(sens.dir, glue("pop_...._{gridRes}.rds"), full.names=T)
  mass.f <- dir(sens.dir, glue("mass_...._{gridRes}.rds"), full.names=T)
  parSets <- map(dir(sens.dir, "000_.*csv", full.names=T), read_csv) %>%
    map(~.x %>% 
          mutate(param=str_remove(paste(param, stage, sep="."), ".NA")) %>%
          select(parDraw, param, val) %>% 
          pivot_wider(values_from="val", names_from="param"))
  params <- map(parSets, names) %>% unlist %>% unique
  params <- params[params!="parDraw"]
  for(i in 1:nrow(grid.i)) {
    siminfo <- str_remove(str_split_fixed(pop.f[i], "pop_", 2)[,2], ".rds")
    pop.i <- readRDS(pop.f[i]) %>% 
      left_join(., parSets[[4]]) %>%
      left_join(., parSets[[filter(grid.i, id==.$id[1])$fetchCat]]) %>%
      filter(stage=="canopy")
    resp <- names(pop.i)
    resp <- resp[!resp %in% c(params, meta.cols)]
    resp <- grep("_mn|_sd", resp, value=T)
    for(j in 1:length(resp)) {
      emulate_sensitivity(pop.i %>% filter(month==7), params, 
                          resp[j], paste0(sens.dir, "BRTs\\"), glue("{siminfo}_july"))
      emulation_summary(resp[j], paste0(sens.dir, "BRTs\\"), glue("{siminfo}_july"))
    }
    
    mass.i <- readRDS(mass.f[i]) %>% 
      left_join(., parSets[[4]]) %>%
      left_join(., parSets[[filter(grid.i, id==.$id[1])$fetchCat]])
    resp <- names(mass.i)
    resp <- resp[!resp %in% c(params, meta.cols)]
    resp <- grep("_mn|_sd|p_k", resp, value=T)
    for(j in 1:length(resp)) {
      emulate_sensitivity(mass.i %>% filter(month==7), params, 
                          resp[j], paste0(sens.dir, "BRTs\\"), glue("{siminfo}_july"))
      emulation_summary(resp[j], paste0(sens.dir, "BRTs\\"), glue("{siminfo}_july"))
    }
  }
}


















