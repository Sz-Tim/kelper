




pkgs <- c("raster", "lubridate", "tidyverse", "sf", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw()); options(mc.cores=12)

gis.dir <- "..\\..\\00_gis\\"
UK_bbox <- st_bbox(c(xmin=-11, xmax=2, ymin=49, ymax=61), crs=st_crs(4326))
coast <- read_sf(paste0(gis.dir, "coastline\\ne_10m_land.shp")) %>%
  add_row(read_sf(paste0(gis.dir, "coastline\\ne_10m_minor_islands.shp"))) %>%
  st_crop(UK_bbox)
data.dir <- "..\\data\\digitized\\"
supp.f <- "..\\data\\collab\\collab_all.xlsx"
data.ls <- compileDatasets(data.dir, supp.f)

site.i <- read_csv(paste0(data.dir, "sitesDigitized.csv"),
                   col_select=1:3, show_col_types=F) %>%
  bind_rows(readxl::read_xlsx(supp.f, "siteLocations", skip=1) %>%
              select(Location, lat, lon) %>% rename(location=Location)) %>%
  group_by(lat, lon) %>% summarise(location=first(location)) %>% ungroup


lmType <- c("lm", "brms")[2]
PAR_datasource <- c("MODIS", "POWER")[1]
gridRes <- 0.5




covars.ls <- loadCovariates(gis.dir, UK_bbox, loadFile="data\\covar_ls.rds")
data.ls <- extractCovarsToDatasets(data.ls, covars.ls, PAR_datasource)
site.sf <- extractCovarsToPts(site.i, covars.ls, PAR_datasource)
grid.domain <- st_read(paste0("data\\grid_", gridRes, ".gpkg"))




lm.fit <- readRDS(paste0("data\\", lmType, "_fits.rds"))




surv.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageTo=="dead") %>%
  mutate(survRate=1-rate,
         exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))
fecund.df <- data.ls$stageFrom_stageTo %>% 
  filter(stageFrom=="canopy" & stageTo=="recruits") %>%
  mutate(exposure=as.numeric(factor(exposure, levels=c("low", "medium", "high"))))




params <- list(tmax=100,
               depth=5,
               prFullHarvest=0, 
               freqHarvest=300
)

sites_pool <- site.sf %>% filter(!is.na(KD_mn)) %>% st_drop_geometry() %>%
  mutate(id=row_number())
grid.i <- grid.domain %>% 






library(foreach)
library(tictoc)
tic()
out.ls <- foreach(i=1:nrow(grid.i),
                  .packages=c("brms", "tidyverse", "lubridate")) %do% {
  par_i <- setParameters(
    tmax=params$tmax, 
    prFullHarvest=params$prFullHarvest, 
    freqHarvest=params$freqHarvest, 
    harvestTarget=params$harvestTarget,
    survRate=filter(surv.df, exposure==grid.i$fetchCat[i])$survRate^(1/2),
    settlementRateBg=filter(fecund.df, exposure==grid.i$fetchCat[i])$rate,
    extraPars=list(
      depth=params$depth,
      env=grid.i[i,],
      lenSt_to_wtSt=lm.fit$lenSt_to_wtSt.lm,
      lenSt_to_wtFr=lm.fit$lenSt_to_wtFr.lm,
      wtFr_to_arFr=lm.fit$wtFr_to_arFr.lm,
      arFr_to_wtFr=lm.fit$arFr_to_wtFr.lm,
      N_canopy=lm.fit$N_canopy.lm,
      FAI=lm.fit$FAI.lm,
      canopyHeight=lm.fit$canopyHeight.lm)
  )
  out <- map(1:2, ~simulatePopulation(par_i, lmType=lmType, ndraws=20))
  list(out=imap_dfr(out, 
                    ~tibble(sim=.y, 
                            year=rep(1:par_i$tmax, 3),
                            month=rep(c(1,6,7), each=par_i$tmax),
                            date=ymd(paste0(year, "-", month, "-01")),
                            PAR_depth=.x$PAR,
                            FAI=c(.x$FAI[3,,]),
                            N.recruits=c(.x$N[1,,]),
                            N.subcanopy=c(.x$N[2,,]),
                            N.canopy=c(.x$N[3,,])) %>%
                      pivot_longer(contains("N."), names_to="stage", values_to="N") %>%
                      mutate(stage=factor(str_sub(stage, 3, -1), 
                                          levels=c("canopy", "subcanopy", "recruits")))) %>%
         add_column(par_i$env) %>%
         mutate(depth=par_i$depth),
       harvest=imap_dfr(out,
                        ~tibble(sim=.y,
                                year=1:par_i$tmax,
                                date=ymd(paste0(2000+year, "-06-01")),
                                PAR_depth=.x$PAR,
                                FAI.June=c(.x$FAI[3,,2]),
                                N.recruits.June=c(.x$N[1,,2]),
                                N.subcanopy.June=c(.x$N[2,,2]),
                                N.canopy.June=c(.x$N[3,,2]),
                                harvest=.x$harvest/1e3) %>%
                          mutate(cumulHarvest=cumsum(harvest))) %>%
         add_column(par_i$env) %>%
         mutate(depth=par_i$depth),
       biomass=imap_dfr(out,
                        ~tibble(sim=.y,
                                year=rep(1:par_i$tmax, 3),
                                month=rep(c(1,6,7), each=par_i$tmax),
                                date=ymd(paste0(year, "-", month, "-01")),
                                PAR_depth=.x$PAR,
                                biomass=c(.x$biomass))) %>%
         add_column(par_i$env) %>%
         mutate(depth=par_i$depth))
}
toc()




out.df <- do.call(rbind, map(out.ls, ~.x$out))
out.sum <- out.df %>% group_by(id, year, month, date, stage) %>%
  select(!c(sim)) %>%
  summarise(#N_q10=quantile(N, 0.1),
    # N_q90=quantile(N, 0.9),
    # FAI_q10=quantile(FAI, 0.1),
    # FAI_q90=quantile(FAI, 0.9),
    across(everything(), mean)) %>%
  ungroup
out.sum.sf <- full_join(grid.domain %>% select(id), out.sum, by="id")

harvest.df <- do.call(rbind, map(out.ls, ~.x$harvest))
harvest.sum <- harvest.df %>% group_by(id, year, date) %>%
  select(!c(sim)) %>%
  summarise(#cumulHarvest_q10=quantile(cumulHarvest, 0.1),
    # cumulHarvest_q90=quantile(cumulHarvest, 0.9),
    across(everything(), mean)) %>%
  ungroup
harvest.sum.sf <- full_join(grid.domain %>% select(id), harvest.sum, by="id")

mass.df <- do.call(rbind, map(out.ls, ~.x$biomass))
mass.sum <- mass.df %>% group_by(id, year, month, date) %>%
  select(!c(sim)) %>%
  summarise(#biomass_q10=quantile(biomass, 0.1),
    # biomass_q90=quantile(biomass, 0.9),
    across(everything(), mean)) %>%
  ungroup
mass.sum.sf <- full_join(grid.domain %>% select(id), mass.sum, by="id")


library(gganimate)
anim <- out.sum.sf %>% filter(stage=="canopy") %>%
  filter(month==6) %>% filter(year>10) %>%
  ggplot() + geom_sf(aes(fill=N), colour=NA) + 
  transition_time(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy abundance at ", 
                 params$depth, "m: Year {frame_time}")) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(paste0("figs\\canopy_N_", params$depth, "m.gif"), 
          anim, nframes=params$tmax)

anim <- out.sum.sf %>% filter(stage=="canopy") %>%
  filter(month==6) %>% filter(year>10) %>%
  ggplot() + geom_sf(aes(fill=FAI), colour=NA) + 
  transition_time(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July FAI at ", 
                 params$depth, "m: Year {frame_time}")) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(paste0("figs\\canopy_FAI_", params$depth, "m.gif"), 
          anim, nframes=params$tmax)

anim <- mass.sum.sf %>% 
  filter(month==6) %>% filter(year>10) %>%
  ggplot() + geom_sf(aes(fill=biomass/1e3), colour=NA) + 
  transition_time(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July biomass at ", 
                 params$depth, "m: Year {frame_time}")) +
  scale_fill_viridis_c("Biomass\nkg/m2") + theme(axis.text=element_blank())
anim_save(paste0("figs\\biomass_", params$depth, ".gif"),
          anim, nframes=params$tmax)
