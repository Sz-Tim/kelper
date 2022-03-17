# KELPER
# Sensitivity analysis
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths with a
# survival and/or loss rate dependent on annual storm intensity




# set up ------------------------------------------------------------------


# libraries and local functions
pkgs <- c("glue", "tidyverse", "sf", "brms", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories & files
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")

# switches & settings
gridRes <- c(0.1, 0.25)[2]
gridFigs <- F
gridRemake <- F
nCores <- 20
set.seed(789)
options(mc.cores=nCores)
pars.sim <- list(depths=c(2, 5, 10, 15, 20), 
                 tskip=0,
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

# set sim length, storm intensity
pars.sim$tmax <- nrow(data.ls$year_stormIndex)
pars.sim$storms <- data.ls$year_stormIndex$stormIndex


# simulate landscapes if needed
if(pars.sim$landscape=="dynamic") {
  if(gridRemake) {
    covars.full <- loadCovariates_full(loadFile=glue("data{sep}covarFull_ls.rds"))
    grid.sim <- grid.sf %>% select(id) %>%
      simulateLandscape(covars.full$sstDayGrow, pars.sim$tmax, "SST") %>%
      full_join(., 
                simulateLandscape(grid.sf %>% select(id), 
                                  covars.full$KD_mn, pars.sim$tmax, "KD") %>%
                  st_drop_geometry()) %>%
      full_join(., 
                simulateLandscape(grid.sf %>% select(id), 
                                  covars.full$PAR_mn, pars.sim$tmax, "PAR") %>%
                  st_drop_geometry()) %>%
      left_join(., grid.sf %>% st_drop_geometry() %>% select(id, fetch, fetchCat))
    rm(covars.full)
    st_write(grid.sim, glue("data{sep}gridSim_{gridRes}_MODIS.gpkg"), append=FALSE)
  } else {
   grid.sim <- st_read(glue("data{sep}gridSim_{gridRes}_MODIS.gpkg")) 
  }
  
  if(gridFigs) {
    w <- 18; h <- 8; dpi <- 200
    basePlot <- grid.sim %>% 
      mutate(PAR2=PAR*exp(-KD*2), PAR5=PAR*exp(-KD*5), PAR10=PAR*exp(-KD*10),
             PAR15=PAR*exp(-KD*15), PAR20=PAR*exp(-KD*20)) %>%
      ggplot() + facet_wrap(~year, nrow=4) + 
      theme(axis.text=element_blank(), legend.position="bottom")
    
    p <- basePlot + geom_sf(aes(fill=SST), colour=NA) + 
      scale_fill_viridis_c("", option="C") + ggtitle("SST")
    ggsave(glue("figs{sep}sim_SST_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (0m)")
    ggsave(glue("figs{sep}sim_PAR_00m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR2), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (2m)")
    ggsave(glue("figs{sep}sim_PAR_02m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR5), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (5m)")
    ggsave(glue("figs{sep}sim_PAR_05m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR10), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (10m)")
    ggsave(glue("figs{sep}sim_PAR_10m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR15), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (15m)")
    ggsave(glue("figs{sep}sim_PAR_15m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR20), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (20m)")
    ggsave(glue("figs{sep}sim_PAR_20m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    
    basePlotResid <- grid.sim %>% 
      mutate(PAR2=PAR*exp(-KD*2), PAR5=PAR*exp(-KD*5), PAR10=PAR*exp(-KD*10),
             PAR15=PAR*exp(-KD*15), PAR20=PAR*exp(-KD*20)) %>%
      group_by(id) %>%
      mutate(across(c("SST", starts_with("PAR")), list("resid"=~.x-mean(.x)))) %>%
      ggplot() + facet_wrap(~year, nrow=4) + scale_fill_gradient2("") +
      theme(axis.text=element_blank(), legend.position="bottom")
    
    p <- basePlotResid + geom_sf(aes(fill=SST_resid), colour=NA) + ggtitle("SST")
    ggsave(glue("figs{sep}resid_SST_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlotResid + geom_sf(aes(fill=PAR_resid), colour=NA) + ggtitle("PAR (0m)")
    ggsave(glue("figs{sep}resid_PAR_00m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlotResid + geom_sf(aes(fill=PAR2_resid), colour=NA) + ggtitle("PAR (2m)")
    ggsave(glue("figs{sep}resid_PAR_02m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlotResid + geom_sf(aes(fill=PAR5_resid), colour=NA) + ggtitle("PAR (5m)")
    ggsave(glue("figs{sep}resid_PAR_05m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlotResid + geom_sf(aes(fill=PAR10_resid), colour=NA) + ggtitle("PAR (10m)")
    ggsave(glue("figs{sep}resid_PAR_10m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlotResid + geom_sf(aes(fill=PAR15_resid), colour=NA) + ggtitle("PAR (15m)")
    ggsave(glue("figs{sep}resid_PAR_15m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlotResid + geom_sf(aes(fill=PAR20_resid), colour=NA) + ggtitle("PAR (20m)")
    ggsave(glue("figs{sep}resid_PAR_20m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    w <- 6; h <- 8; dpi <- 200
    p <- ggplot(grid.sf) + geom_sf(aes(fill=fetch), colour=NA) + 
      scale_fill_viridis_c("", option="B", end=0.9) + ggtitle("Fetch") +
      theme(axis.text=element_blank(), legend.position="bottom")
    ggsave(glue("figs{sep}mn_fetch_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- ggplot(grid.sf) + geom_sf(aes(fill=fetchCat), colour=NA) + 
      scale_fill_viridis_c("", option="B", end=0.9) + ggtitle("Fetch category") +
      theme(axis.text=element_blank(), legend.position="bottom")
    ggsave(glue("figs{sep}mn_fetchCat_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    rm(basePlot); rm(basePlotResid); rm(p)
  }
} else {
  grid.sim <- NULL
  
  if(gridFigs) {
    w <- 6; h <- 8; dpi <- 200
    basePlot <- grid.sf %>% 
      rename(PAR=PAR_mn, SST=sstDay_mn, KD=KD_mn) %>%
      mutate(PAR2=PAR*exp(-KD*2), PAR5=PAR*exp(-KD*5), PAR10=PAR*exp(-KD*10),
             PAR15=PAR*exp(-KD*15), PAR20=PAR*exp(-KD*20)) %>%
      ggplot() + theme(axis.text=element_blank(), legend.position="bottom")
    
    p <- basePlot + geom_sf(aes(fill=fetch), colour=NA) + 
      scale_fill_viridis_c("", option="B", end=0.9) + ggtitle("Fetch")
    ggsave(glue("figs{sep}mn_fetch_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=fetchCat), colour=NA) + 
      scale_fill_viridis_c("", option="B", end=0.9) + ggtitle("Fetch category")
    ggsave(glue("figs{sep}mn_fetchCat_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=SST), colour=NA) + 
      scale_fill_viridis_c("", option="C") + ggtitle("SST")
    ggsave(glue("figs{sep}mn_SST_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (0m)")
    ggsave(glue("figs{sep}mn_PAR_00m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR2), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (2m)")
    ggsave(glue("figs{sep}mn_PAR_02m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR5), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (5m)")
    ggsave(glue("figs{sep}mn_PAR_05m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR10), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (10m)")
    ggsave(glue("figs{sep}mn_PAR_10m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR15), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (15m)")
    ggsave(glue("figs{sep}mn_PAR_15m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    p <- basePlot + geom_sf(aes(fill=PAR20), colour=NA) + 
      scale_fill_viridis_c("") + ggtitle("PAR (20m)")
    ggsave(glue("figs{sep}mn_PAR_20m_{gridRes}.png"), p, width=w, height=h, dpi=dpi)
    
    rm(basePlot); rm(p)
  }
}








# run simulations ---------------------------------------------------------


cl <- makeCluster(nCores, outfile=glue("temp{sep}sim_out.txt"))
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
