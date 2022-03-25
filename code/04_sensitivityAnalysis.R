# KELPER
# Sensitivity analysis
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths for
# combinations of parameter values




# set up ------------------------------------------------------------------


# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")


# directories & files
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
sens.dir <- glue("out{sep}sensitivity{sep}")
parSets.f <- glue("000_paramSets_exposure_")

# switches & settings
rerun <- T
reanalyse <- T
gridRes <- c(0.1, 0.25)[1]
nCores <- 50
options(mc.cores=nCores)
pars.sens <- list(nParDraws=1e3,
                  depths=c(2, 15),
                  tmax=30,
                  tskip=10,
                  nSim=1,
                  stochParams=F, 
                  landscape="static")

# load files
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg"))
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn) %>%
  select(id, SST, KD, PAR, fetch, fetchCat)
data.ls <- compileDatasets(data.dir, supp.f)
lm.fit <- readRDS(glue("data{sep}fits_{gridRes}.rds"))
lm.mnsd <- readRDS(glue("data{sep}dfs_mn_sd_{gridRes}.rds"))
surv.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageTo=="dead") %>%
  mutate(survRate=1-rate_mn,
         exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))
fecund.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageFrom=="canopy" & stageTo=="recruits") %>%
  mutate(exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))








# assign parameter values -------------------------------------------------


# parameter ranges: approximately mean +- 2 sd
loss_mnPrec <- c(0.2340987, 90.1752) # mn, prec for beta distribution
rate_stipe <- cbind(c(194, 195, 58), c(23, 21, 10)) # mn, sd by stage
rate_frond <- cbind(c(1787, 2299, 3979), c(179, 397, 417))/1e4
par.rng <- list(surv=surv.df %>%
                  mutate(stage=stageFrom,
                         valMean=1-rate_mn,
                         valSD=rate_sd,
                         valMin=1-pmin(1, pmax(0, qnorm(0.975, rate_mn, rate_sd))),
                         valMax=1-pmin(1, pmax(0, qnorm(0.025, rate_mn, rate_sd)))) %>%
                  select(stage, valMean, valSD, valMin, valMax, exposure),
                settlement=fecund.df %>%
                  mutate(valMean=rate_mn,
                         valSD=rate_sd,
                         valMin=pmax(0, qnorm(0.025, rate_mn, rate_sd)),
                         valMax=pmax(0, qnorm(0.975, rate_mn, rate_sd))) %>%
                  select(valMean, valSD, valMin, valMax, exposure),
                growStipe=tibble(stage=c("recruits", "subcanopy", "canopy"),
                                 valMean=rate_stipe[,1],
                                 valSD=rate_stipe[,2],
                                 valMin=apply(rate_stipe, 1, function(x) qnorm(0.025, x[1], x[2])),
                                 valMax=apply(rate_stipe, 1, function(x) qnorm(0.975, x[1], x[2]))),
                growFrond=tibble(stage=c("recruits", "subcanopy", "canopy"),
                                 valMean=rate_frond[,1],
                                 valSD=rate_frond[,2],
                                 valMin=apply(rate_frond, 1, function(x) qnorm(0.025, x[1], x[2])),
                                 valMax=apply(rate_frond, 1, function(x) qnorm(0.975, x[1], x[2]))),
                loss=tibble(valMean=loss_mnPrec[1],
                            valMin=qbeta(0.025, prod(loss_mnPrec), (1-loss_mnPrec[1])*loss_mnPrec[2]),
                            valMax=qbeta(0.975, prod(loss_mnPrec), (1-loss_mnPrec[1])*loss_mnPrec[2])),
                densityEffShape=tibble(valMean=1, valMin=0.5, valMax=1.5)) %>%
  imap_dfr(., ~.x %>% mutate(param=.y))
saveRDS(par.rng, glue("{sens.dir}parameter_ranges.rds"))

if(rerun) {
  parSets <- map_dfr(1:pars.sens$nParDraws, 
                     ~par.rng %>% mutate(parDraw=.x, 
                                         val=runif(n(), valMin, valMax))) %>%
    select(parDraw, param, stage, exposure, val) %>%
    group_by(exposure) %>% group_split()
  walk2(parSets, c(1:3, "NA"), ~write_csv(.x, glue("{sens.dir}{parSets.f}{.y}.csv")))
}










# run simulations ---------------------------------------------------------


if(rerun) {
  cl <- makeCluster(nCores, outfile=glue("temp{sep}sensitivity_out.txt"))
  obj.exclude <- c("data.ls", "grid.sf", "surv.df", "fecund.df")
  obj.include <- ls()
  
  cat("Exporting cluster.\n")
  clusterExport(cl, obj.include[-match(obj.exclude, obj.include)])
  Sys.sleep(5)
  cat("Exported. Starting parallel runs.\n")
  out.ls <- parLapply(cl, X=1:nrow(grid.i), fun=simSensitivityDepthsWithinCell,
                      grid.i=grid.i, gridRes=gridRes, pars.sens=pars.sens, 
                      lm.fit=lm.fit, lm.mnsd=lm.mnsd, parSets=parSets)
  stopCluster(cl)
  cat("Finished simulations.\n")
}










# run BRTs ----------------------------------------------------------------


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
  
  cl <- makeCluster(nCores, outfile=glue("temp{sep}brt_out.txt"))
  obj.exclude <- c("data.ls", "grid.sf", "surv.df", "fecund.df", "lm.fit", "lm.mnsd")
  obj.include <- ls()
  
  cat("Exporting cluster.\n")
  clusterExport(cl, obj.include[-match(obj.exclude, obj.include, nomatch=0)])
  Sys.sleep(5)
  cat("Exported. Starting BRTs.\n")
  out.ls <- parLapply(cl, X=1:nrow(grid.i), fun=runBRTs,
                      pop.f=pop.f, mass.f=mass.f, parSets=parSets,
                      grid.i=grid.i, meta.cols=meta.cols, params=params,
                      brt.dir=glue("{sens.dir}{sep}BRTs{sep}"))
  stopCluster(cl)
  cat("Finished BRTs.\n")
}








# aggregate RIs -----------------------------------------------------------


dir(glue("{sens.dir}{sep}BRTs{sep}summaries{sep}"), "_ri_", full.names=T) %>%
  map_dfr(~read_csv(.x, show_col_types=F)) %>% 
  write_csv(glue("{sens.dir}{sep}RelInf_{gridRes}.csv"))













