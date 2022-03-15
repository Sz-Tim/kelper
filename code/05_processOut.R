# KELPER
# Process output
# Tim Szewczyk

# This script processes the output produces in 03_runSimulations.R




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")

# switches & settings
gridRes <- 0.25

# load files
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
sim.info <- glue("...._{gridRes}")
pop.f <- dir(out.dir, glue("pop_{sim.info}"), full.names=T)
pop.df <- map_dfr(pop.f, readRDS)
mass.f <- dir(out.dir, glue("mass_{sim.info}"), full.names=T)
mass.df <- map_dfr(mass.f, readRDS)
obs.ls <- readRDS(glue("data{sep}dfs_{gridRes}.rds"))

y_vars <- c("FAI", "N", "biomass", "logN", "logBiomass", "kappa_N", "kappa_FAI")
x_vars <- c("K_N", "K_FAI", "SST", "KD", "PAR", "PAR_atDepth", "fetch", "fetchCat")
sim.title <- glue("{gridRes} arc-sec grid")
pop.sum <- pop.df %>% 
  # filter(year > 10) %>%
  select(-date, -year) %>%
  mutate(logN=log(N+1)) %>%
  group_by(id, sim, month, stage, depth, landscape, stochParams) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max, 
                            pOcc=~sum(.x>1)/n(), q90=~quantile(.x, probs=0.9))))
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- mass.df %>%
  # filter(year > 10) %>%
  select(-date, -year, -month) %>%
  mutate(logBiomass=log(biomass+1)) %>%
  group_by(id, sim, depth, landscape, stochParams) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max, q90=~quantile(.x, probs=0.9))))
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")


saveRDS(pop.df, glue("temp{sep}pop_df_{gridRes}.rds"))
saveRDS(mass.df, glue("temp{sep}mass_df_{gridRes}.rds"))
saveRDS(pop.sum, glue("temp{sep}pop_sum_{gridRes}.rds"))
saveRDS(mass.sum, glue("temp{sep}mass_sum_{gridRes}.rds"))






sim.title <- glue("{gridRes} arc-sec grid")
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
pop.sum <- readRDS(glue("temp{sep}pop_sum_{gridRes}.rds"))
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.sum <- readRDS(glue("temp{sep}mass_sum_{gridRes}.rds"))
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")


# Abundance by PAR
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_mn, colour=SST_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title) 
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_q90, colour=SST_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title) 
mass.sum %>% 
  ggplot(aes(PAR_atDepth_mn, biomass_q90, colour=SST_mn)) + geom_point(alpha=0.5) +
  facet_grid(landscape~stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title) 
pop.sum %>% filter(month==6) %>% 
  ggplot(aes(fetch_mn, N_mn, colour=PAR_atDepth_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(fetch, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_max, colour=fetch_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_pOcc, colour=SST_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(PAR_atDepth) %>% 
               summarise(PAR_atDepth=mean(PAR_atDepth), N_pOcc=sum(N>0)/n()) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N_pOcc), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_sd, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, N_sd/N_mn, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(PAR_atDepth_mn, kappa_FAI_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(PAR_atDepth_mn, kappa_N_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(PAR_atDepth_mn, SST_mn, colour=kappa_N_mn-kappa_FAI_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_gradient2() +
  labs(title=sim.title)

pop.sum %>% filter(month==6) %>%
  ggplot(aes(PAR_atDepth_mn, FAI_max, colour=fetch_mn)) + geom_point() +
  geom_point(data=obs.ls$FAI.lm,
             aes(PAR_atDepth, FAI), colour="red", alpha=0.5, shape=1) +
  facet_grid(.~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)



# Abundance maps
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>% 
  group_by(depth, landscape, stochParams, id) %>% summarise(N_mn=mean(N_mn)) %>%
  ggplot(aes(fill=N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==1) %>%
  ggplot(aes(fill=N_pOcc)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=N_md>1)) + geom_sf(colour=NA) + 
  facet_grid(landscape*stochParams~depth)

pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_N_md>1)) + geom_sf(colour=NA) + 
  facet_grid(landscape*stochParams~depth)

pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(depth~landscape*stochParams)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_FAI_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(depth~landscape*stochParams)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=K_FAI_md>1)) + geom_sf(colour=NA) + 
  facet_grid(depth~landscape*stochParams)

pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=kappa_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=kappa_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy" & month==6) %>%
  ggplot(aes(fill=kappa_N_mn-kappa_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_gradient2(limits=c(-1, 1)) +
  facet_grid(landscape*stochParams~depth)



mass.sum.sf %>%
  ggplot(aes(fill=biomass_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
mass.sum.sf %>%
  ggplot(aes(fill=biomass_q90)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
mass.sum.sf %>%
  ggplot(aes(fill=biomass_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
mass.sum.sf %>%
  ggplot(aes(fill=biomass_sd/biomass_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c("biomass CV") +
  facet_grid(landscape*stochParams~depth)



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



samp_id <- sample(unique(pop.df$id), 4)

pop.df %>% filter(id %in% samp_id) %>% filter(month==1) %>%
  ggplot(aes(year, N, group=paste(id, fetchCat), colour=fetchCat)) +
  geom_line(alpha=0.75) + facet_grid(stage~depth, scales="free_y") + 
  scale_colour_viridis_c(end=0.9, limits=c(1,3))

mass.df %>% filter(id %in% samp_id) %>% 
  ggplot(aes(year, biomass, group=paste(id, depth), colour=depth)) +
  geom_line() + facet_grid(.~month, scales="free_y")



pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, N, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + 
  geom_line() +
  facet_grid(id~depth) 
mass.df %>% filter(id %in% samp_id) %>% filter(month==6) %>% 
  ggplot(aes(year, biomass, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + 
  geom_line() +
  facet_grid(id~depth) 

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, kappa_FAI, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth) + ylim(0, 1) 

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>% filter(stage=="canopy") %>%
  mutate(kappa=pmax(kappa_N, kappa_FAI)) %>%
  ggplot(aes(year, kappa, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth) + ylim(0, 1)

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, N, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")

pop.df %>% filter(id %in% samp_id) %>% filter(month==6) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, FAI, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")

mass.df %>% filter(id %in% samp_id) %>% filter(month==6) %>%
  ggplot(aes(year, biomass, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")




pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  ggplot(aes(mnLogLambda, colour=depth, group=depth)) + geom_density() + 
  scale_colour_viridis_c(direction=-1) + 
  facet_grid(stochParams~landscape, scales="free")

pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=sign(mnLogLambda) * (abs(mnLogLambda) > 0.05))) + geom_sf(colour=NA) +
  scale_fill_gradient2(low="red", high="blue", mid="grey70") +
  facet_grid(depth~stochParams*landscape)
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=mnLogLambda > -0.05)) + geom_sf(colour=NA) +
  facet_grid(depth~stochParams*landscape)
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=mnLogLambda>0)) + geom_sf(colour=NA) +
  facet_grid(depth~stochParams*landscape)

pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=mnLogLambda)) + geom_sf(colour=NA) +
  scale_fill_gradient2() +
  facet_grid(depth~stochParams*landscape)
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(sdLogLambda=sd(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=sdLogLambda)) + geom_sf(colour=NA) +
  scale_fill_viridis_c() +
  facet_grid(depth~stochParams*landscape)







library(gganimate)
pop.sf <- full_join(grid.sf, 
                    pop.df %>% filter(stage=="canopy") %>%
                      filter(month==6), 
                    by="id")
anim <- pop.sf %>% 
  ggplot() + geom_sf(aes(fill=N), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy abundance: Year {closest_state}")) +
  facet_grid(.~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_N_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

anim <- pop.sf %>% 
  ggplot() + geom_sf(aes(fill=FAI), colour=NA) + 
  transition_states(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy FAI: Year {closest_state}")) +
  facet_grid(.~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_FAI_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

anim <- mass.sum.sf %>% 
  filter(month==6) %>% filter(year>10) %>%
  ggplot() + geom_sf(aes(fill=biomass), colour=NA) + 
  transition_states(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July biomass at ", 
                 params$depth, "m: Year {closest_state}")) +
  scale_fill_viridis_c("Biomass\nkg/m2") + theme(axis.text=element_blank())
anim_save(glue("figs{sep}biomass_{params$depth}m.gif"),
          anim, nframes=params$tmax-10)

anim <- out.sum.sf %>% filter(stage=="canopy") %>%
  filter(month==6) %>% filter(year>10) %>%
  ggplot() + geom_sf(aes(fill=SST), colour=NA) + 
  transition_states(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("Mean SST (Jan-Jun): Year {closest_state}")) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}SST.gif"), 
          anim, nframes=params$tmax-10)

anim <- out.sum.sf %>% filter(stage=="canopy") %>%
  filter(month==6) %>% filter(year>10) %>%
  ggplot() + geom_sf(aes(fill=PAR_atDepth), colour=NA) + 
  transition_states(year) + ease_aes('cubic-in-out') +
  ggtitle(paste0("Mean PAR (Jan-Jun) at ", 
                 params$depth, "m: Year {closest_state}")) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}PAR_{params$depth}m.gif"), 
          anim, nframes=params$tmax-10)
