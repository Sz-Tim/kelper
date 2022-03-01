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
gridRes <- 0.5

# load files
grid.sf <- st_read(glue("data\\grid_{gridRes}_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
sim.info <- glue("...._{gridRes}")
pop.f <- dir(out.dir, glue("pop_{sim.info}"), full.names=T)
pop.df <- map_dfr(pop.f, readRDS)
mass.f <- dir(out.dir, glue("mass_{sim.info}"), full.names=T)
mass.df <- map_dfr(mass.f, readRDS)
obs.ls <- readRDS("data\\dfs_{gridRes}_MODIS.rds")

y_vars <- c("FAI", "N", "biomass", "logN", "logBiomass", "kappa_N", "kappa_FAI")
x_vars <- c("K_N", "K_FAI", "SST", "KD", "PAR", "PAR_atDepth", "fetch", "fetchCat")
sim.title <- glue("{gridRes} arc-sec grid")
pop.sum <- pop.df %>% 
  filter(year > 10) %>%
  select(-date, -year) %>%
  mutate(logN=log(N+1)) %>%
  group_by(id, sim, month, stage, depth, landscape, lmType) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max, pOcc=~sum(.x>1)/n())))
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- mass.df %>%
  filter(year > 10) %>%
  select(-date, -year, -month) %>%
  mutate(logBiomass=log(biomass+1)) %>%
  group_by(id, sim, depth, landscape, lmType) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max)))
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")


saveRDS(pop.sum, glue("temp\\pop_sum_{gridRes}.rds"))
saveRDS(mass.sum, glue("temp\\mass_sum_{gridRes}.rds"))





sim.title <- glue("0.5 arc-sec grid")
grid.sf <- st_read(glue("data\\grid_0.5_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
pop.sum <- readRDS("temp\\pop_sum_0.5.rds") %>% filter(lmType=="brms")
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- readRDS("temp\\mass_sum_0.5.rds") %>% filter(lmType=="brms") 
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")

sim.title <- glue("0.25 arc-sec grid")
grid.sf <- st_read(glue("data\\grid_0.25_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
pop.sum <- readRDS("temp\\pop_sum_0.25.rds")
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- readRDS("temp\\mass_sum_0.25.rds")
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")

sim.title <- glue("0.1 arc-sec grid")
grid.sf <- st_read(glue("data\\grid_0.1_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
pop.sum <- readRDS("temp\\pop_sum_0.1.rds")
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- readRDS("temp\\mass_sum_0.1.rds")
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")


# Abundance by PAR
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_md, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>% 
  ggplot(aes(fetch_mn, N_md, colour=PAR_atDepth_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_max, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_pOcc, colour=SST_mn)) + geom_point(alpha=0.5) +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_sd, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, kappa_FAI_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, kappa_N_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)



# Abundance maps
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_md)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_pOcc)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_mn>1)) + geom_sf(colour=NA) + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_max)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)

pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_mn>1)) + geom_sf(colour=NA) + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=kappa_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=kappa_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_FAI_mn/K_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=kappa_FAI_mn-kappa_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_gradient2() + 
  facet_grid(lmType~depth)

mass.sum.sf %>% 
  ggplot(aes(fill=biomass_md)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=biomass_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=biomass_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(lmType~depth)



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
mass.df %>% 
  filter(landscape=="dynamic" & depth==5 & month==6) %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=PAR_atDepth)) + geom_sf(colour=NA) +
  scale_fill_viridis_c("PAR (5m)") + 
  facet_wrap(~year, nrow=3)

mass.sum.sf %>% 
  ggplot(aes(fill=fetch_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)

mass.sum.sf %>% 
  ggplot(aes(fill=SST_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% 
  ggplot(aes(fill=SST_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.df %>% 
  filter(landscape=="dynamic" & depth==2 & month==6) %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=SST)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(option="B") + 
  facet_wrap(~year, nrow=3)




mass.df %>% 
  filter(landscape=="static" & depth==2 & month==6) %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=biomass)) + geom_sf(colour=NA) +
  scale_fill_viridis_c() + 
  facet_wrap(~year, nrow=3)



samp_id <- sample(grid.sf$id, 10)

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  ggplot(aes(date, kappa_N, group=id, colour=depth)) + geom_line() +
  facet_grid(stage~depth, scales="free_y") + ylim(0, 1)

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  ggplot(aes(date, kappa_FAI, group=id, colour=depth)) + geom_line() +
  facet_grid(stage~depth, scales="free_y") + ylim(0, 1)

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  mutate(kappa=pmax(kappa_N, kappa_FAI)) %>%
  ggplot(aes(date, kappa, group=id, colour=depth)) + geom_line() +
  facet_grid(stage~depth, scales="free_y") + ylim(0, 1)

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  ggplot(aes(date, N, group=id, colour=depth)) + geom_line() +
  facet_grid(stage~depth, scales="free_y")

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  ggplot(aes(date, FAI, group=id, colour=depth)) + geom_line() +
  facet_grid(stage~depth, scales="free_y")

mass.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  ggplot(aes(date, biomass, group=id, colour=depth)) + geom_line() +
  facet_grid(.~depth, scales="free_y")

