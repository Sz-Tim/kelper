# KELPER
# EDA
# Tim Szewczyk

# This script is for exploring the processed output




# set up ------------------------------------------------------------------

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "ggpubr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")
sens.dir <- glue("out{sep}sensitivity{sep}")

# switches & settings
gridRes <- 0.1

# plot themes, etc
map_base.gg <- ggplot() +
  theme(legend.position="bottom", 
        panel.grid=element_blank()) +
  scale_x_continuous(breaks=c(-8,0)) +
  scale_y_continuous(breaks=c(52, 58)) 









# regressions -------------------------------------------------------------

lm.fit <- readRDS(glue("data{sep}fits_{gridRes}.rds"))
out.eff <- map(lm.fit, conditional_effects)

ggpubr::ggarrange(plotlist=plot(out.eff[[1]], ask=F, points=T, plot=F, 
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1, widths=c(1,1,1.2))
ggsave(glue("figs{sep}regr{sep}condEff_1_{gridRes}.png"), width=12, height=4)

ggpubr::ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[2]], 
                                    effects=c("logLenStipe", "lPAR_atDepth", "SST")), 
                ask=F, points=T, plot=F, 
                point_args=list(shape=1, size=0.75, alpha=0.8)), 
  nrow=1, widths=c(1,1,1))
ggsave(glue("figs{sep}regr{sep}condEff_2a_{gridRes}.png"), width=12, height=4)
ggpubr::ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[2]], 
                                    effects=c("logLenStipe:lPAR_atDepth", "lPAR_atDepth:SST")), 
                ask=F, points=T, plot=F, 
                point_args=list(shape=1, size=0.75, alpha=0.8)), 
  nrow=1, widths=c(1,1))
ggsave(glue("figs{sep}regr{sep}condEff_2b_{gridRes}.png"), width=9, height=4)

ggpubr::ggarrange(plotlist=plot(out.eff[[3]], ask=F, points=T, plot=F, 
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1, widths=c(1,1,1.2))
ggsave(glue("figs{sep}regr{sep}condEff_4_{gridRes}.png"), width=12, height=4)

ggpubr::ggarrange(plotlist=plot(out.eff[[4]], ask=F, points=T, plot=F, 
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1, widths=c(1,1,1.2))
ggsave(glue("figs{sep}regr{sep}condEff_3_{gridRes}.png"), width=12, height=4)

ggpubr::ggarrange(plotlist=plot(out.eff[[5]], ask=F, points=T, plot=F, 
                                point_args=list(shape=1, size=0.75, alpha=0.8)))
ggsave(glue("figs{sep}regr{sep}condEff_5_{gridRes}.png"), width=4, height=4)

ggpubr::ggarrange(plotlist=plot(out.eff[[6]], ask=F, points=T, plot=F, 
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1)
ggsave(glue("figs{sep}regr{sep}condEff_6_{gridRes}.png"), width=8, height=4)

ggpubr::ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[7]], 
                                    effects=c("PAR_atDepth", "SST", "fetch")), ask=F,
                points=T, plot=F, 
                point_args=list(shape=1, size=0.75, alpha=0.8)),
  nrow=1)
ggsave(glue("figs{sep}regr{sep}condEff_7a_{gridRes}.png"), width=12, height=4)
ggpubr::ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[7]], 
                                    effects=c("PAR_atDepth:SST", "PAR_atDepth:fetch")), ask=F, 
                points=T, plot=F, 
                point_args=list(shape=1, size=0.75, alpha=0.8)),
  nrow=1)
ggsave(glue("figs{sep}regr{sep}condEff_7b_{gridRes}.png"), width=9, height=4)







# storm simulations -------------------------------------------------------

sim.title <- glue("{gridRes} arc-sec grid")
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2))
pop.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}.rds")) %>% filter(month!=6) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2))
pop.sum <- readRDS(glue("summaries{sep}pop_sum_{gridRes}.rds")) %>% filter(month!=6) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2))
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")
mass.df <- readRDS(glue("summaries{sep}mass_df_{gridRes}.rds")) %>% filter(month!=6) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2))
mass.sum <- readRDS(glue("summaries{sep}mass_sum_{gridRes}.rds")) %>% filter(month!=6) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2))
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")
obs.ls <- readRDS(glue("data{sep}dfs_{gridRes}.rds"))





# environmental scatterplots ----------------------------------------------

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






# maps: kelp --------------------------------------------------------------

pop.sum.sf_canopy <- pop.sum.sf %>% filter(stage=="canopy") %>%
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])
mass.sum.sf_canopy <- mass.sum.sf %>% 
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])
map_base.gg <- ggplot() +
  theme(legend.position="bottom", 
        panel.grid=element_blank()) +
  scale_x_continuous(breaks=c(-8,0)) +
  scale_y_continuous(breaks=c(52, 58)) 


# * N: all months, landscapes, depths  ------------------------------------

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy, 
          colour=NA, aes(fill=N_mn)) +
  scale_fill_gradient("N", low="grey90", high="green4", limits=c(0, NA)) +
  facet_grid(month*landscape~depth) +
  ggtitle("Mean canopy density")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_mean_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy, 
          colour=NA, aes(fill=N_md>1)) +
  scale_fill_manual("N > 1", values=c("grey90", "green4")) +
  facet_grid(month*landscape~depth) +
  ggtitle("Canopy presence (median N > 1)")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_Pres_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy, 
          colour=NA, aes(fill=N_pOcc)) +
  scale_fill_gradient("N", low="grey90", high="green4") +
  facet_grid(month*landscape~depth) +
  ggtitle("Proportion of years with canopy presence")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_pOcc_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy, 
          colour=NA, aes(fill=logN_sd)) +
  scale_fill_viridis_c("sd") +
  facet_grid(month*landscape~depth) +
  ggtitle("Canopy variability among years (sd(ln(N))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_sdln_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy, 
          colour=NA, aes(fill=logN_sd/logN_mn)) +
  scale_fill_viridis_c("CV") +
  facet_grid(month*landscape~depth) +
  ggtitle("Canopy variability among years (CV(ln(N))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_CVln_MonLndDep.png"), p, width=6, height=8, dpi=200)





# * N: select scenarios ---------------------------------------------------

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=N_mn)) +
  scale_fill_gradient("N", low="grey90", high="green4", limits=c(0, NA)) +
  facet_grid(.~depth) +
  ggtitle("July mean canopy density")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_mean_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=N_md>1)) +
  scale_fill_manual("N > 1", values=c("grey90", "green4")) +
  facet_grid(.~depth) +
  ggtitle("July canopy presence (median N > 1)")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_Pres_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=N_pOcc)) +
  scale_fill_gradient("N", low="grey90", high="green4") +
  facet_grid(.~depth) +
  ggtitle("Proportion of years with July canopy presence")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_pOcc_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=logN_sd)) +
  scale_fill_viridis_c("sd") +
  facet_grid(.~depth) +
  ggtitle("July canopy variability among years (sd(ln(N))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_sdln_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=logN_sd/logN_mn)) +
  scale_fill_viridis_c("CV") +
  facet_grid(.~depth) +
  ggtitle("July canopy variability among years (CV(ln(N))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_CVln_JulStaDep.png"), p, width=8, height=4, dpi=200)
 






# * N: resid landscape ----------------------------------------------------
pop.resid.sf <- pop.sum.sf_canopy %>% 
  select(id, month, depth, landscape, N_md, N_mn, N_sd, N_pOcc, logN_mn, logN_sd) %>%
  mutate(N_pres=as.numeric(N_md>1),
         logN_CV=logN_sd/logN_mn) %>%
  group_by(id, month, depth) %>% 
  arrange(id, month, depth, landscape) %>%
  mutate(across(where(is.numeric), ~first(.x)-last(.x))) %>%
  slice_head(n=1)

p <- map_base.gg +
  geom_sf(data=pop.resid.sf, 
          colour=NA, aes(fill=N_mn)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Mean canopy abundance: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_mean_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.resid.sf, 
          colour=NA, aes(fill=N_pres)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Canopy presence: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_Pres_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.resid.sf, 
          colour=NA, aes(fill=N_pOcc)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Proportion of years with canopy presence: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_pOcc_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.resid.sf, 
          colour=NA, aes(fill=logN_sd)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Canopy variability among years: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_sdln_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=pop.resid.sf, 
          colour=NA, aes(fill=logN_CV)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Canopy variability among years: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_N_CVln_MonDep_Resid.png"), p, width=8, height=6, dpi=200)




# * biomass: all months, landscapes, depths -------------------------------

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy, 
          colour=NA, aes(fill=biomass_mn)) +
  scale_fill_gradient("Biomass (kg)", low="grey90", high="green4", limits=c(0, NA)) +
  facet_grid(month*landscape~depth) +
  ggtitle("Mean canopy biomass")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_mean_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy, 
          colour=NA, aes(fill=logBiomass_mn)) +
  scale_fill_gradient("ln(biomass) (kg)", low="grey90", high="green4") +
  facet_grid(month*landscape~depth) +
  ggtitle("Mean canopy biomass")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_lnmean_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy, 
          colour=NA, aes(fill=logBiomass_sd)) +
  scale_fill_viridis_c("sd") +
  facet_grid(month*landscape~depth) +
  ggtitle("Canopy biomass variability among years (sd(ln(biomass))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_sdln_MonLndDep.png"), p, width=6, height=8, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy, 
          colour=NA, aes(fill=logBiomass_sd/logBiomass_mn)) +
  scale_fill_viridis_c("CV") +
  facet_grid(month*landscape~depth) +
  ggtitle("Canopy variability among years (CV(log(biomass))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_CVln_MonLndDep.png"), p, width=6, height=8, dpi=200)





# * biomass: select scenarios ---------------------------------------------

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=biomass_mn)) +
  scale_fill_gradient("Biomass (kg)", low="grey90", high="green4", limits=c(0, NA)) +
  facet_grid(.~depth) +
  ggtitle("July mean canopy biomass")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_mean_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=logBiomass_mn)) +
  scale_fill_gradient("ln(biomass) (kg)", low="grey90", high="green4") +
  facet_grid(.~depth) +
  ggtitle("July mean canopy biomass")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_lnmean_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=logBiomass_sd)) +
  scale_fill_viridis_c("sd") +
  facet_grid(.~depth) +
  ggtitle("July canopy biomass variability among years (sd(ln(biomass))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_sdln_JulStaDep.png"), p, width=8, height=4, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & landscape=="static"), 
          colour=NA, aes(fill=logBiomass_sd/logBiomass_mn)) +
  scale_fill_viridis_c("CV") +
  facet_grid(.~depth) +
  ggtitle("July canopy variability among years (CV(log(biomass))")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_CVln_JulStaDep.png"), p, width=8, height=4, dpi=200)







# * biomass: resid landscape ----------------------------------------------
mass.resid.sf <- mass.sum.sf_canopy %>% 
  select(id, month, depth, landscape, biomass_md, biomass_mn, biomass_sd, logBiomass_mn, logBiomass_sd) %>%
  mutate(logBiomass_CV=logBiomass_sd/logBiomass_mn) %>%
  group_by(id, month, depth) %>% 
  arrange(id, month, depth, landscape) %>%
  mutate(across(where(is.numeric), ~first(.x)-last(.x))) %>%
  slice_head(n=1)

p <- map_base.gg +
  geom_sf(data=mass.resid.sf, 
          colour=NA, aes(fill=biomass_mn)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Mean canopy biomass: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_mean_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.resid.sf, 
          colour=NA, aes(fill=logBiomass_mn)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Mean canopy biomass: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_lnmean_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.resid.sf, 
          colour=NA, aes(fill=logBiomass_sd)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Canopy biomass variability among years: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_sdln_MonDep_Resid.png"), p, width=8, height=6, dpi=200)

p <- map_base.gg +
  geom_sf(data=mass.resid.sf, 
          colour=NA, aes(fill=logBiomass_CV)) +
  scale_fill_gradient2("dynamic - static", mid="grey90") +
  facet_grid(month~depth) +
  ggtitle("Canopy variability among years: Environmental stochasticity")
ggsave(glue("figs{sep}storm{sep}eda{sep}map_biomass_CVln_MonDep_Resid.png"), p, width=8, height=6, dpi=200)








# * pub figs --------------------------------------------------------------

fig2a <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & landscape=="dynamic"), 
          colour=NA, aes(fill=logBiomass_mn)) +
  scale_fill_gradient(expression("ln kg/m"^2), low="grey90", high="green4", limits=c(0,NA)) +
  facet_grid(.~depth) + theme(legend.position="right") +
  ggtitle("July canopy biomass mean") + theme_classic()
fig2b <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & landscape=="dynamic"), 
          colour=NA, aes(fill=logBiomass_sd)) +
  scale_fill_viridis_c(expression("ln kg/m"^2), limits=c(0, NA)) +
  facet_grid(.~depth) + theme(legend.position="right") +
  ggtitle("July canopy biomass standard deviation among years") + theme_classic()
fig2 <- ggarrange(fig2a, fig2b, ncol=1, labels="auto")
ggsave(glue("figs{sep}storm{sep}pub{sep}biomass_mn_sd.png"), fig2,
       width=8, height=5.75, dpi=300)







# timeseries --------------------------------------------------------------

samp <- mass.df %>% 
  filter(landscape=="static" & depth==2 & month==1 & year==1) %>%
  group_by(fetchCat) %>%
  sample_n(2)
map_base.gg + 
  geom_sf(data=grid.sf, colour=NA, fill="grey90") +
  geom_sf(data=filter(grid.sf, id %in% samp$id), colour=NA, aes(fill=as.character(id))) +
  scale_fill_brewer(type="qual", palette=2) +
  theme(legend.position="right")
mass.df %>% filter(id %in% samp$id) %>% 
  filter(month==1) %>%
  filter(depth %in% c(2,5,10)) %>% 
  filter(landscape=="static") %>%
  ggplot(aes(year, biomass, group=id, colour=as.character(id))) +
  geom_line() + scale_colour_brewer(type="qual", palette=2) +
  facet_grid(fetchCat~depth, scales="free_y") 




mass.df %>% filter(id %in% samp$id) %>% filter(month==7) %>%
  ggplot(aes(year, biomass, group=paste(landscape, depth), colour=landscape)) +
  geom_line() + facet_grid(id~depth, scales="free_y")



pop.df %>% filter(id %in% samp$id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, N, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + 
  geom_line() +
  facet_grid(id~depth) 
mass.df %>% filter(id %in% samp$id) %>% filter(month==7) %>% 
  ggplot(aes(year, biomass, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + 
  geom_line() +
  facet_grid(id~depth) 

pop.df %>% filter(id %in% samp$id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, kappa_FAI, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth) + ylim(0, 1) 

pop.df %>% filter(id %in% samp$id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  mutate(kappa=pmax(kappa_N, kappa_FAI)) %>%
  ggplot(aes(year, kappa, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth) + ylim(0, 1)

pop.df %>% filter(id %in% samp$id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, N, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")

pop.df %>% filter(id %in% samp$id) %>% filter(month==7) %>% filter(stage=="canopy") %>%
  ggplot(aes(year, FAI, group=paste(landscape, stochParams, sim), 
             colour=landscape, linetype=stochParams)) + geom_line() +
  facet_grid(id~depth, scales="free_y")

mass.df %>% filter(id %in% samp$id) %>% filter(month==7) %>%
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
  filter(id %in% samp$id) %>%
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





# animations --------------------------------------------------------------


library(gganimate)
pop.sf <- full_join(grid.sf, 
                    pop.df %>% filter(stage=="canopy") %>%
                      filter(month==7), 
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
                      filter(month==1), 
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
  facet_grid(landscape~depth) +
  scale_fill_viridis_c() + theme(axis.text=element_blank())
anim_save(glue("figs{sep}canopy_FAI_Jan_{gridRes}.gif"), 
          anim, nframes=max(pop.sf$year))

mass.sf <- full_join(grid.sf, mass.df %>% filter(month==7), by="id")
anim <- mass.sf %>% 
  ggplot() + geom_sf(aes(fill=biomass), colour=NA) + 
  transition_states(year+1942) + ease_aes('cubic-in-out') +
  ggtitle(paste0("July canopy biomass: Year {closest_state}")) +
  facet_grid(landscape~depth) +
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




