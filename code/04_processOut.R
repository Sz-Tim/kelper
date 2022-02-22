# KELPER
# Process output
# Tim Szewczyk

# This script processes the output produces in 03_runSimulations.R




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "lme4", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())

# directories
out.dir <- "out\\"
data.dir <- "..\\data\\digitized\\"
supp.f <- "..\\data\\collab\\collab_all.xlsx"

# switches & settings
lmType <- c("lm", "brms")[1]
PAR_datasource <- c("MODIS", "POWER")[1]
gridRes <- 0.25
dynamicLandscape <- F

# load files
grid.sf <- st_read(glue("data\\grid_{gridRes}_{PAR_datasource}.gpkg")) %>%
  select(id, geom, PAR_surface)
sim.info <- glue("..._{gridRes}_{lmType}",
                 "_{ifelse(dynamicLandscape,'dynamic','static')}")
pop.f <- dir(out.dir, glue("pop_{sim.info}.rds"), full.names=T)
pop.df <- map_dfr(pop.f, readRDS)
mass.f <- dir(out.dir, glue("mass_{sim.info}.rds"), full.names=T)
mass.df <- map_dfr(mass.f, readRDS)

y_vars <- c("FAI", "N", "biomass", "logN", "logBiomass", "kappa_N", "kappa_FAI")
x_vars <- c("K_N", "K_FAI", "SST", "KD", "PAR", "PAR_atDepth", "fetch", "fetchCat")
sim.title <- glue("{gridRes} arc-sec grid, {lmType} regr, {PAR_datasource} PAR")
pop.sum <- pop.df %>%
  select(-date, -year) %>%
  mutate(logN=log(N+1)) %>%
  group_by(id, sim, month, stage, depth, landscape) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max)))
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- mass.df %>%
  select(-date, -year, -month) %>%
  mutate(logBiomass=log(biomass+1)) %>%
  group_by(id, sim, depth, landscape) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max)))
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")


saveRDS(pop.sum, "temp\\pop_sum.rds")
saveRDS(mass.sum, "temp\\mass_sum.rds")


# Abundance by PAR
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_md, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_max, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_sd, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, logN_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)


# Abundance maps
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_max)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)

pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_mn>1)) + geom_sf(colour=NA) + 
  facet_grid(landscape~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)

mass.sum.sf %>% 
  ggplot(aes(fill=biomass_md)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=biomass_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)



# Landscape maps
mass.sum.sf %>% 
  ggplot(aes(fill=PAR_atDepth_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=PAR_atDepth_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=PAR_atDepth_sd/PAR_atDepth_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c("PAR CV") + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=PAR_atDepth_max)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=PAR_atDepth_max-PAR_atDepth_min)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c("PAR range") + 
  facet_grid(landscape~depth)

mass.sum.sf %>% 
  ggplot(aes(fill=SST_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=SST_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)


