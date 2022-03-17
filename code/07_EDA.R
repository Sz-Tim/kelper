# KELPER
# EDA
# Tim Szewczyk

# This script is for exploring the processed output




# set up ------------------------------------------------------------------

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
gridRes <- 0.1















# storm simulations -------------------------------------------------------

sim.title <- glue("{gridRes} arc-sec grid")
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
pop.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}.rds")) %>% filter(month!=6)
pop.sum <- readRDS(glue("summaries{sep}pop_sum_{gridRes}.rds")) %>% filter(month!=6)
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.df <- readRDS(glue("summaries{sep}mass_df_{gridRes}.rds")) %>% filter(month!=6)
mass.sum <- readRDS(glue("summaries{sep}mass_sum_{gridRes}.rds")) %>% filter(month!=6)
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")
obs.ls <- readRDS(glue("data{sep}dfs_{gridRes}.rds"))


# Abundance by PAR
pop.sum %>% filter(month==7) %>%
  ggplot(aes(PAR_atDepth_mn, N_mn, colour=SST_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title) 
pop.sum %>% filter(month==1) %>%
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
  facet_grid(landscape~month) + 
  scale_colour_viridis_c() +
  labs(title=sim.title) 
pop.sum %>% filter(month==7) %>% 
  ggplot(aes(fetch_mn, N_mn, colour=PAR_atDepth_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(fetch, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>%
  ggplot(aes(PAR_atDepth_mn, N_max, colour=fetch_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(location) %>% 
               summarise(across(everything(), mean)) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>%
  ggplot(aes(PAR_atDepth_mn, N_pOcc, colour=SST_mn)) + geom_point(alpha=0.5) +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(PAR_atDepth) %>% 
               summarise(PAR_atDepth=mean(PAR_atDepth), N_pOcc=sum(N>0)/n()) %>%
               mutate(stage="canopy"),
             aes(PAR_atDepth, N_pOcc), colour="red", shape=1) +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>%
  ggplot(aes(PAR_atDepth_mn, N_sd, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>%
  ggplot(aes(PAR_atDepth_mn, N_sd/N_mn, colour=fetch_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(PAR_atDepth_mn, kappa_FAI_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(PAR_atDepth_mn, kappa_N_mn, colour=SST_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)
pop.sum %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(PAR_atDepth_mn, SST_mn, colour=kappa_N_mn-kappa_FAI_mn)) + geom_point() +
  facet_grid(stage~landscape*stochParams, scales="free_y") + 
  scale_colour_gradient2() +
  labs(title=sim.title)

pop.sum %>% filter(month==7) %>%
  ggplot(aes(PAR_atDepth_mn, FAI_max, colour=fetch_mn)) + geom_point() +
  geom_point(data=obs.ls$FAI.lm,
             aes(PAR_atDepth, FAI), colour="red", alpha=0.5, shape=1) +
  facet_grid(.~landscape, scales="free_y") + 
  scale_colour_viridis_c() +
  labs(title=sim.title)



# Abundance maps
pop.sum.sf %>% filter(stage=="canopy") %>% 
  ggplot(aes(fill=N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(month*landscape~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% 
  ggplot(aes(fill=N_pOcc)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(month*landscape~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=N_sd/N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(month*landscape~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=N_md>1)) + geom_sf(colour=NA) + 
  facet_grid(month*landscape~depth)

pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=K_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=K_N_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=K_N_md>1)) + geom_sf(colour=NA) + 
  facet_grid(landscape*stochParams~depth)

pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=K_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=K_FAI_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=K_FAI_md>1)) + geom_sf(colour=NA) +  
  facet_grid(landscape*stochParams~depth)

pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=kappa_N_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>%filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=kappa_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(landscape*stochParams~depth)
pop.sum.sf %>% filter(stage=="canopy") %>% filter(month!=6) %>%
  ggplot(aes(fill=kappa_N_mn-kappa_FAI_mn)) + geom_sf(colour=NA) + 
  scale_fill_gradient2(limits=c(-1, 1)) +
  facet_grid(landscape*stochParams~depth)



mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=biomass_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=biomass_q90)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=biomass_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=biomass_sd/biomass_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c("biomass CV") +
  facet_grid(landscape*stochParams~depth)



# Landscape maps
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=PAR_atDepth_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=PAR_atDepth_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=PAR_atDepth_sd/PAR_atDepth_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c("PAR CV") + 
  facet_grid(landscape~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=PAR_atDepth_max)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=PAR_atDepth_max-PAR_atDepth_min)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c("PAR range") + 
  facet_grid(landscape~depth)
mass.df %>% filter(month==7) %>%
  filter(landscape=="dynamic" & depth==5 & month==7) %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=PAR_atDepth)) + geom_sf(colour=NA) +
  scale_fill_viridis_c("PAR (5m)") + 
  facet_wrap(~year, nrow=3)

mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=fetch_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)

mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=SST_mn)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.sum.sf %>% filter(month==7) %>%
  ggplot(aes(fill=SST_sd)) + geom_sf(colour=NA) + 
  scale_fill_viridis_c() + 
  facet_grid(landscape~depth)
mass.df %>% filter(month==7) %>%
  filter(landscape=="dynamic" & depth==2 & month==7) %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=SST)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(option="B") + 
  facet_wrap(~year, nrow=3)




mass.df %>% 
  filter(landscape=="static" & depth==2 & month==7) %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=biomass)) + geom_sf(colour=NA) +
  scale_fill_viridis_c() + 
  facet_wrap(~year, nrow=3)



samp_id <- sample(unique(pop.df$id), 10)

pop.df %>% filter(id %in% samp_id) %>% filter(stage=="canopy") %>% filter(month==1) %>%
  ggplot(aes(year, N, group=paste(id, landscape), colour=landscape)) +
  geom_line() + facet_grid(id~depth, scales="free_y") 

mass.df %>% filter(id %in% samp_id) %>% filter(month==7) %>%
  ggplot(aes(year, biomass, group=paste(landscape, depth), colour=landscape)) +
  geom_line() + facet_grid(id~depth, scales="free_y")



pop.df %>% filter(id %in% samp_id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, N, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + 
  geom_line() +
  facet_grid(id~depth) 
mass.df %>% filter(id %in% samp_id) %>% filter(month==7) %>% 
  ggplot(aes(year, biomass, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + 
  geom_line() +
  facet_grid(id~depth) 

pop.df %>% filter(id %in% samp_id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, kappa_FAI, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth) + ylim(0, 1) 

pop.df %>% filter(id %in% samp_id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  mutate(kappa=pmax(kappa_N, kappa_FAI)) %>%
  ggplot(aes(year, kappa, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth) + ylim(0, 1)

pop.df %>% filter(id %in% samp_id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, N, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")

pop.df %>% filter(id %in% samp_id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, FAI, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")

mass.df %>% filter(id %in% samp_id) %>% filter(month==7) %>%
  ggplot(aes(year, biomass, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")





mass.sum %>% ungroup %>%
  select(id, month, depth, landscape, logBiomass_mn) %>%
  pivot_wider(names_from="landscape", values_from="logBiomass_mn") %>%
  full_join(grid.sf, ., by="id") %>%
  ggplot(aes(fill=dynamic-static)) + geom_sf(colour=NA) + 
  scale_fill_gradient2() +
  facet_grid(month~depth)









# I think this doesn't really make sense because of the open settlement. Bummer.
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  filter(id %in% samp_id) %>%
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  ggplot(aes(log(lambda), colour=depth, group=depth)) + geom_density() + 
  scale_colour_viridis_c(direction=-1) + 
  facet_grid(id~landscape, scales="free")

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
  filter(lambda < 2.5) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=sign(mnLogLambda) * (abs(mnLogLambda) > 0.05))) + geom_sf(colour=NA) +
  scale_fill_gradient2(low="red", high="blue", mid="grey70") +
  facet_grid(landscape*stochParams~depth)
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  filter(lambda < 2.5) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=mnLogLambda > -0.05)) + geom_sf(colour=NA) +
  facet_grid(landscape*stochParams~depth)
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  filter(lambda < 2.5) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=mnLogLambda>0)) + geom_sf(colour=NA) +
  facet_grid(landscape*stochParams~depth)

pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  filter(lambda < 2.5) %>%
  summarise(mnLogLambda=mean(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=mnLogLambda)) + geom_sf(colour=NA) +
  scale_fill_gradient2() +
  facet_grid(landscape*stochParams~depth)
pop.df %>% filter(stage=="canopy") %>%
  filter(month==1) %>% 
  group_by(id, depth, landscape, stochParams, sim) %>%
  mutate(lambda=N/lag(N, 1)) %>%
  summarise(sdLogLambda=sd(log(lambda), na.rm=T)) %>%
  full_join(grid.sf, .) %>%
  ggplot(aes(fill=sdLogLambda)) + geom_sf(colour=NA) +
  scale_fill_viridis_c() +
  facet_grid(landscape*stochParams~depth)







library(gganimate)
pop.sf <- full_join(grid.sf, 
                    pop.df %>% filter(stage=="canopy") %>%
                      filter(month==7) %>% filter(depth < 15), 
                    by="id")
anim <- pop.sf %>% 
  ggplot() + geom_sf(aes(fill=N), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy abundance: Year {closest_state}")) +
  facet_grid(landscape~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_N_Jul_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

anim <- pop.sf %>% 
  ggplot() + geom_sf(aes(fill=FAI), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy FAI: Year {closest_state}")) +
  facet_grid(landscape~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_FAI_Jul_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

pop.sf <- full_join(grid.sf, 
                    pop.df %>% filter(stage=="canopy") %>%
                      filter(month==1) %>% filter(depth < 15), 
                    by="id")
anim <- pop.sf %>% 
  ggplot() + geom_sf(aes(fill=N), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("January canopy abundance: Year {closest_state}")) +
  facet_grid(landscape~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_N_Jan_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

anim <- pop.sf %>% 
  ggplot() + geom_sf(aes(fill=FAI), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("January canopy FAI: Year {closest_state}")) +
  facet_grid(.~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_FAI_Jan_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

mass.sf <- full_join(grid.sf, mass.df %>% filter(month==7), by="id")
anim <- mass.sf %>% 
  ggplot() + geom_sf(aes(fill=biomass), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy biomass: Year {closest_state}")) +
  facet_grid(.~depth) +
  scale_fill_viridis_c("Biomass\nkg/m2") + theme(axis.text=element_blank())
anim_save(glue("figs{sep}biomass_Jul_{gridRes}.gif"),
          anim, nframes=max(mass.sf$year))

mass.sf <- full_join(grid.sf, mass.df %>% filter(month==1), by="id")
anim <- mass.sf %>% 
  ggplot() + geom_sf(aes(fill=biomass), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy biomass: Year {closest_state}")) +
  facet_grid(.~depth) +
  scale_fill_viridis_c("Biomass\nkg/m2") + theme(axis.text=element_blank())
anim_save(glue("figs{sep}biomass_Jan_{gridRes}.gif"),
          anim, nframes=max(mass.sf$year))




