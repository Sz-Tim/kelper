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

# grid
sim.title <- glue("{gridRes} arc-sec grid")
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  # select(id, geom, PAR_surface) %>%
  filter(id > ifelse(gridRes==0.1, 5, 2))
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn) %>%
  select(id, SST, KD, PAR, fetch, fetchCat)
obs.ls <- readRDS(glue("data{sep}dfs_{gridRes}.rds"))
data.ls <- compileDatasets(data.dir, supp.f)

# plot themes, etc
map_base.gg <- ggplot() +
  theme(legend.position="bottom", 
        panel.grid=element_blank(),
        axis.title=element_blank()) +
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

gridSim.sum <- readRDS(glue("data{sep}gridSim_{gridRes}.rds")) %>%
  select(id, grid.id, year, PAR, SST, KD) %>%
  group_by(id) %>%
  summarise(SST_mn=mean(SST),
            PAR_0=mean(PAR),
            PAR_2=mean(PAR*exp(-KD*2)),
            PAR_5=mean(PAR*exp(-KD*5)),
            PAR_10=mean(PAR*exp(-KD*10)),
            PAR_15=mean(PAR*exp(-KD*15)),
            PAR_20=mean(PAR*exp(-KD*20))) %>%
  ungroup %>%
  pivot_longer(starts_with("PAR_"), names_to="depth", values_to="PAR_atDepth") %>%
  mutate(depth=as.numeric(str_sub(depth, 5, -1))) %>%
  full_join(., grid.i %>% select(id, fetch, fetchCat))

mass.df <- readRDS(glue("summaries{sep}mass_df_{gridRes}.rds")) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2)) %>%
  filter(!stochParams)
mass.sum <- readRDS(glue("summaries{sep}mass_sum_{gridRes}.rds")) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2)) %>%
  left_join(., gridSim.sum)
mass.sum.sf <- full_join(grid.sf, mass.sum, by="id")

Nrcr.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}_recruits.rds")) %>%
  filter(id > ifelse(gridRes==0.1, 5, 2))
pop.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}.rds")) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2))
pop.sum <- readRDS(glue("summaries{sep}pop_sum_{gridRes}.rds")) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2)) 
pop.sum.sf <- full_join(grid.sf, pop.sum, by="id")

pop.sum.sf_canopy <- pop.sum.sf %>% filter(stage=="canopy") %>%
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])
mass.sum.sf_canopy <- mass.sum.sf %>% 
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])





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


# * biomass ---------------------------------------------------------------

mass.sum %>% 
  filter(month==7) %>%
  filter(!stochParams) %>%
  filter(depth==2) %>%
  summary


mass.sum %>%
  filter(month==7) %>% filter(!stochParams) %>% 
  ungroup %>% 
  mutate(depth=paste0(depth, "m"),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20), "m"))) %>%
  ggplot(aes(biomass_mn, biomass_sd, colour=fetch)) + 
  scale_colour_viridis_c(option="B") + 
  geom_point(shape=1, alpha=0.5) + 
  facet_grid(depth~.) +
  labs(x=bquote(Biomass~mean~(kg/m^2)), y=bquote(Biomass~interannual~sd~(kg/m^2)))
ggsave(glue("figs{sep}pub{sep}biomass_mn_sd_byFetch.png"), 
       width=3.5, height=9, dpi=300)
mass.sum %>%
  filter(month==7) %>% filter(!stochParams) %>% 
  ungroup %>% 
  mutate(depth=paste0(depth, "m"),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20), "m"))) %>%
  ggplot(aes(biomass_mn, biomass_sd/biomass_mn, colour=fetch)) + 
  scale_colour_viridis_c(option="B") + 
  geom_point(shape=1, alpha=0.5) + 
  facet_grid(depth~.) +
  labs(x=bquote(Biomass~mean~(kg/m^2)), y="Biomass interannual CV")
ggsave(glue("figs{sep}pub{sep}biomass_mn_CV_byFetch.png"), 
       width=3.5, height=9, dpi=300)

fig_a1 <- mass.sum %>% filter(month==7) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_mn, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day))),
       y="Mean July canopy biomass (kg)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
fig_b1 <- mass.sum %>% filter(month==7) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(y="SD July canopy biomass (kg)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
fig_c1 <- mass.sum %>% filter(month==7) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd/biomass_mn, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day))),
       y="CV(July canopy biomass)")

fig_a2 <- mass.sum %>% filter(month==7) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_mn, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank())
fig_b2 <- mass.sum %>% filter(month==7) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank())
fig_c2 <- mass.sum %>% filter(month==7) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd/biomass_mn, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day)))) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

fig.sst <- ggarrange(plotlist=list(fig_a1, fig_b1, fig_c1), nrow=3, 
                     common.legend=T, legend="bottom", 
                     labels=c("a.", "b.", "c."), label.x=0.85, label.y=0.95)
fig.fetch <- ggarrange(plotlist=list(fig_a2, fig_b2, fig_c2), nrow=3, 
                       common.legend=T, legend="bottom", 
                       labels=c("d.", "e.", "f."), label.x=0.85, label.y=0.95)
fig.sst_fetch <- ggarrange(fig.sst, fig.fetch, widths=c(1, 0.9))
ggsave(glue("figs{sep}pub{sep}PAR_SST_biomass_Jul_stoch.png"), 
       fig.sst_fetch, width=6, height=9, dpi=300)



fig_a1 <- mass.sum %>% filter(month==1) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_mn, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day))),
       y="Mean January canopy biomass (kg)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
fig_b1 <- mass.sum %>% filter(month==1) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(y="SD January canopy biomass (kg)") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
fig_c1 <- mass.sum %>% filter(month==1) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd/biomass_mn, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day))),
       y="CV(January canopy biomass)")

fig_a2 <- mass.sum %>% filter(month==1) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_mn, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank())
fig_b2 <- mass.sum %>% filter(month==1) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank())
fig_c2 <- mass.sum %>% filter(month==1) %>%
  filter(stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd/biomass_mn, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day)))) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

fig.sst <- ggarrange(plotlist=list(fig_a1, fig_b1, fig_c1), nrow=3, 
                     common.legend=T, legend="bottom", 
                     labels=c("a.", "b.", "c."), label.x=0.85, label.y=0.95)
fig.fetch <- ggarrange(plotlist=list(fig_a2, fig_b2, fig_c2), nrow=3, 
                       common.legend=T, legend="bottom", 
                       labels=c("d.", "e.", "f."), label.x=0.85, label.y=0.95)
fig.sst_fetch <- ggarrange(fig.sst, fig.fetch, widths=c(1, 0.9))
ggsave(glue("figs{sep}pub{sep}PAR_SST_biomass_Jan_stoch.png"), 
       fig.sst_fetch, width=6, height=9, dpi=300)




# maps: kelp --------------------------------------------------------------

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
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & stochParams==F), 
          colour=NA, aes(fill=logBiomass_mn)) +
  scale_fill_gradient(expression("ln kg/m"^2), low="grey90", high="green4", limits=c(0,NA)) +
  facet_grid(.~depth) + theme(legend.position="right") +
  ggtitle("July canopy biomass mean") + theme_classic()
fig2b <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy %>% filter(month=="July" & stochParams==F), 
          colour=NA, aes(fill=logBiomass_sd)) +
  scale_fill_viridis_c(expression("ln kg/m"^2), limits=c(0, NA)) +
  facet_grid(.~depth) + theme(legend.position="right") +
  ggtitle("July canopy biomass standard deviation among years") + theme_classic()
fig2 <- ggarrange(fig2a, fig2b, ncol=1, labels="auto")
ggsave(glue("figs{sep}pub{sep}biomass_mn_sd.png"), fig2,
       width=8, height=5.75, dpi=300)







# timeseries --------------------------------------------------------------

mass_simGrid.df <- read_csv(glue("summaries{sep}mass_simGrid_{gridRes}.csv"))
samp_id <- sample(unique(mass_simGrid.df$id), 4)

mass_simGrid.df %>% filter(id %in% samp_id) %>% 
  filter(!stochParams) %>%
  filter(month(date)==7) %>%
  left_join(., grid.i, by="id") %>%
  mutate(covar=glue("SST: {round(SST,1)}, PAR: {round(PAR,1)}, ",
                    "KD: {round(KD,1)}, fetch: {round(fetch,1)}")) %>%
  ggplot(aes(date, biomass_md, group=depth, colour=depth, fill=depth)) + 
  # geom_ribbon(aes(ymin=biomass_10, ymax=biomass_90),
  #             alpha=0.5, colour=NA) +
  geom_ribbon(aes(ymin=biomass_25, ymax=biomass_75),
              alpha=0.5, colour=NA) +
  geom_line() +
  stat_smooth(data=data.ls$year_stormIndex, 
              aes(date(paste0(year, "-01-01")), (stormIndex*10+10)),
              se=F, method="loess", formula=y~x, span=0.5, linetype=2) +
  scale_colour_viridis_c("Depth (m)", direction=-1, end=0.95) +
  scale_fill_viridis_c("Depth (m)", direction=-1, end=0.95) +
  scale_x_date(breaks=date(c("1950-01-01", "1975-01-01", "2000-01-01")),
               date_labels="%Y") +
  facet_wrap(~covar, ncol=1, dir="h") + theme_classic() + 
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm")) +
  labs(x="", y="Canopy biomass (kg)")
ggsave("figs/storm/timeseries_biomass.png", width=3, height=8, dpi=200)

mass_simGrid.df %>% filter(id %in% samp_id) %>% 
  filter(!stochParams) %>%
  mutate(month=month(date)) %>%
  ggplot(aes(date, biomass_md, group=depth, colour=depth)) + 
  geom_line() + 
  stat_smooth(method="loess", formula=y~x, se=F, span=0.25) +
  scale_colour_viridis_c("Depth (m)", direction=-1, end=0.95) +
  scale_x_date(breaks=date(c("1950-01-01", "1975-01-01", "2000-01-01")),
               date_labels="%Y") +
  facet_grid(id~month) + theme_classic() + 
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm")) +
  labs(y="Mean biomass")


mass_simGrid.df %>% filter(id %in% samp_id) %>% 
  filter(!stochParams) %>%
  mutate(month=month(date)) %>%
  ggplot(aes(date, biomass_sd, group=depth, colour=depth)) + 
  geom_line() + 
  stat_smooth(method="loess", formula=y~x, se=F, span=0.25) +
  scale_colour_viridis_c("Depth (m)", direction=-1, end=0.95) +
  scale_x_date(breaks=date(c("1950-01-01", "1975-01-01", "2000-01-01")),
               date_labels="%Y") +
  facet_grid(id~month) + theme_classic() + 
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm")) +
  labs(y="Biomass sd")


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











# sensitivity analysis ----------------------------------------------------


ri.df <- read_csv(glue("{sens.dir}{sep}RelInf_{gridRes}.csv")) %>%
  group_by(response, id, month, depth) %>%
  filter(smp==max(smp), td==max(td))
stages <- c("canopy", "recruits", "subcanopy")
stages_ord <- c("canopy", "subcanopy", "recruits")
sens.param <- tibble(messy=sort(unique(ri.df$var)),
                     full=c("Density shape", 
                            paste("Frond growth:", stages),
                            paste("Stipe growth:", stages),
                            "Erosion", "Settlement", 
                            paste("Survival:", stages)),
                     eq=c("theta", 
                          paste0("omega[", stages, "]"),
                          paste0("gamma[", stages, "]"),
                          "epsilon", "z",
                          paste0("s[", stages, "]"))) %>%
  mutate(eq=paste0("italic(", eq, ")")) %>%
  mutate(eq=factor(eq, 
                   levels=c("italic(epsilon)",
                            "italic(theta)", 
                            "italic(z)",
                            paste0("italic(omega[", rev(stages_ord), "])"),
                            paste0("italic(gamma[", rev(stages_ord), "])"),
                            paste0("italic(s[", rev(stages_ord), "])")
                            )))

ri.df %>% group_by(response, var) %>%
  summarise(mn=round(mean(rel.inf), 3),
            sd=round(sd(rel.inf), 3)) %>% 
  arrange(response, desc(mn)) %>%
  print.AsIs()

ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  group_by(var, month, depth, response) %>%
  summarise(md=median(rel.inf*100),
            q10=quantile(rel.inf*100, probs=0.1),
            q25=quantile(rel.inf*100, probs=0.25),
            q75=quantile(rel.inf*100, probs=0.75),
            q90=quantile(rel.inf*100, probs=0.9)) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  ggplot(aes(param, md, colour=depth, group=as.factor(depth))) +
  geom_hline(yintercept=0, colour="grey90") +
  geom_errorbar(aes(ymin=q10, ymax=q90), position=position_dodge(width=0.75), width=0.5) +
  geom_linerange(aes(ymin=q25, ymax=q75), position=position_dodge(width=0.75), size=1) +
  geom_point(position=position_dodge(width=0.75)) +
  facet_grid(.~response) + coord_flip() +
  scale_colour_viridis_c("Depth (m)", direction=-1, end=0.95) +
  scale_x_discrete(labels=scales::parse_format()) +
  scale_y_continuous(breaks=c(0, 20, 40)) + 
  theme_classic() +
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(y="Relative influence (%)", x="", 
       title="Canopy biomass sensitivity")
ggsave(glue("figs{sep}pub{sep}sensitivity_pts_mdAmongCells.png"), width=7, height=5, dpi=300)

ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  group_by(var, month, depth, response) %>%
  summarise(mn=mean(rel.inf*100),
            se=sd(rel.inf*100)/sqrt(n())) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  ggplot(aes(param, mn, colour=depth, group=as.factor(depth))) +
  geom_hline(yintercept=0, colour="grey90") +
  geom_errorbar(aes(ymin=mn-2*se, ymax=mn+2*se), position=position_dodge(width=0.75), width=0.5) +
  geom_point(position=position_dodge(width=0.75)) +
  facet_grid(.~response) + coord_flip() +
  scale_colour_viridis_c("Depth (m)", direction=-1, end=0.95) +
  scale_x_discrete(labels=scales::parse_format()) +
  scale_y_continuous(breaks=c(0, 20, 40)) + 
  theme_classic() +
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(y="Relative influence (%)", x="", 
       title="Canopy biomass sensitivity")
ggsave(glue("figs{sep}pub{sep}sensitivity_pts_mnAmongCells.png"), width=7, height=5, dpi=300)


ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  ggplot(aes(rel.inf*100, param, fill=depth)) +
  geom_vline(xintercept=0, colour="grey90") +
  ggridges::geom_density_ridges(alpha=0.5, scale=1, colour="grey30", size=0.25) +
  facet_grid(.~response) +
  scale_fill_viridis_d("Depth (m)", direction=-1, end=0.95) +
  scale_y_discrete(labels=scales::parse_format()) +
  scale_x_continuous(breaks=c(0, 25, 50)) +
  theme_classic() +
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(x="Relative influence", y="", 
       title="Canopy biomass sensitivity")


ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  right_join(grid.i %>% select(id, fetchCat), by="id") %>%
  mutate(rel.inf=rel.inf*100) %>%
  group_by(var, month, depth, response, fetchCat) %>%
  summarise(md=median(rel.inf),
            q10=quantile(rel.inf, probs=0.1),
            q25=quantile(rel.inf, probs=0.25),
            q75=quantile(rel.inf, probs=0.75),
            q90=quantile(rel.inf, probs=0.9)) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)],
         fetchCat=factor(fetchCat, labels=c("Low", "High"))) %>% 
  ggplot(aes(depth, md, colour=fetchCat, fill=fetchCat)) +
  geom_hline(yintercept=0, colour="grey90") +
  geom_ribbon(aes(ymin=q25, ymax=q75), colour=NA, alpha=0.4) + 
  geom_ribbon(aes(ymin=q10, ymax=q90), colour=NA, alpha=0.3) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_line() +
  facet_grid(response~param, 
             labeller=labeller(response=label_value, depth=label_value, param=label_parsed)) + 
  scale_colour_viridis_d("Exposure", option="B", end=0.75) +
  scale_fill_viridis_d("Exposure", option="B", end=0.75) +
  scale_x_continuous(breaks=c(0, 5, 10, 15), limits=c(0,16)) +
  scale_y_continuous(breaks=c(0, 25, 50)) + 
  theme_classic() +
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(y="Relative influence (%)", x="Depth (m)", 
       title="Canopy biomass sensitivity")
ggsave(glue("figs{sep}pub{sep}sensitivity_ribbon_mnAmongCells.png"), width=11, height=4, dpi=300)

ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  right_join(grid.i %>% select(id, fetchCat), by="id") %>%
  mutate(rel.inf=log10(rel.inf*100)) %>%
  group_by(var, month, depth, response, fetchCat) %>%
  summarise(md=median(rel.inf),
            q10=quantile(rel.inf, probs=0.1),
            q25=quantile(rel.inf, probs=0.25),
            q75=quantile(rel.inf, probs=0.75),
            q90=quantile(rel.inf, probs=0.9)) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)],
         fetchCat=factor(fetchCat, labels=c("Low", "High"))) %>% 
  ggplot(aes(depth, md, colour=fetchCat, fill=fetchCat)) +
  geom_ribbon(aes(ymin=q25, ymax=q75), colour=NA, alpha=0.4) + 
  geom_ribbon(aes(ymin=q10, ymax=q90), colour=NA, alpha=0.3) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_line() +
  facet_grid(response~param, 
             labeller=labeller(response=label_value, depth=label_value, param=label_parsed)) + 
  scale_colour_viridis_d("Exposure", option="B", end=0.75) +
  scale_fill_viridis_d("Exposure", option="B", end=0.75) +
  scale_x_continuous(breaks=c(0, 5, 10, 15), limits=c(0,16)) +
  # scale_y_continuous(breaks=c(0, 25, 50)) + 
  theme_classic() +
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(y="log10(Relative influence)", x="Depth (m)", 
       title="Canopy biomass sensitivity")
ggsave(glue("figs{sep}pub{sep}sensitivity_ribbon_mnAmongCells_logRI.png"), width=11, height=4, dpi=300)



ri.df %>%
  filter(month=="july") %>%
  group_by(var, month, depth, response) %>%
  summarise(rel.inf=sd(rel.inf*100)) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  ggplot(aes(param, rel.inf, fill=depth, group=as.factor(depth))) +
  geom_bar(stat="identity", position="dodge", colour="grey30", size=0.25) +
  facet_grid(.~response) + coord_flip() +
  scale_fill_viridis_c("Depth (m)", direction=-1) +
  scale_x_discrete(labels=scales::parse_format()) +
  scale_y_continuous(breaks=c(0, 10, 20)) + 
  theme_classic() + 
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(y="Standard deviation in relative influence among cells", x="",
       title="Canopy biomass sensitivity")
ggsave(glue("figs{sep}pub{sep}sensitivity_bar_sdAmongCells.png"), width=7, height=4, dpi=300)


ri.df %>% filter(month=='july') %>% 
  filter(depth %in% paste0(c(2, 10), "m")) %>%
  filter(response=="biomass_mn") %>%
  group_by(var) %>%
  mutate(max_ri=max(rel.inf)) %>%
  ungroup %>%
  filter(max_ri > 0.1) %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20),"m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  left_join(grid.sf, .) %>%
  ggplot(aes(fill=log10(rel.inf*100))) +
  geom_sf(colour=NA) +
  scale_fill_viridis_c(expression(log[10](Rel.~Inf.)), option="cividis", limits=c(NA,2)) +
  theme_classic() + 
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        panel.grid=element_blank(),
        axis.title=element_blank(), axis.text=element_blank()) +
  scale_x_continuous(breaks=c(-8,0)) +
  scale_y_continuous(breaks=c(52, 58)) +
  facet_grid(depth~param, 
             labeller=labeller(response=label_value, depth=label_value, param=label_parsed))
ggsave(glue("figs{sep}pub{sep}sensitivity_map_logRI.png"), width=6, height=3.75, dpi=300)

ri.df %>% filter(month=='july') %>% 
  filter(depth %in% paste0(c(2, 10), "m")) %>%
  filter(response=="biomass_mn") %>%
  group_by(var) %>%
  mutate(max_ri=max(rel.inf)) %>%
  ungroup %>%
  filter(max_ri > 0.1) %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20),"m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  left_join(grid.sf, .) %>%
  ggplot(aes(fill=rel.inf*100)) +
  geom_sf(colour=NA) +
  scale_fill_viridis_c("Rel. Inf. (%)", option="cividis") +
  theme_classic() + 
  theme(axis.text=element_blank(),
        legend.position="bottom", legend.key.height=unit(0.2, "cm")) +
  facet_grid(depth~param, 
             labeller=labeller(response=label_value, depth=label_value, param=label_parsed))
ggsave(glue("figs{sep}pub{sep}sensitivity_map_RI.png"), width=6, height=3.75, dpi=300)




ri.df %>% filter(response=="biomass_mn") %>%
  right_join(grid.i, .) %>%
  ggplot(aes(rel.inf, colour=as.factor(fetchCat))) +
  geom_density() +
  facet_wrap(~depth*month*var, scales="free_y", ncol=n_distinct(ri.df$var))


ri.df %>%
  right_join(grid.i, .) %>%
  ggplot(aes(fetch, rel.inf, colour=paste(month, depth))) + geom_point(shape=1) +
  facet_grid(response~var)










# lagged effects ----------------------------------------------------------

nLags <- 15

stormIndexLag <- data.ls$year_stormIndex %>% select(year) %>%
  bind_cols(imap_dfc(setNames(1:nLags, glue("lag_{1:nLags}_strm")), 
                     ~lag(data.ls$year_stormIndex$stormIndex, .x)))

gridSim.df <- readRDS(glue("data{sep}gridSim_{gridRes}.rds")) %>%
  select(id, grid.id, year, PAR, SST, KD)



# * lag recruits ----------------------------------------------------------

N.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}_recruits.rds")) %>%
  filter(id > ifelse(gridRes==0.1, 5, 2))
N.df_jul <- N.df %>%
  filter(month==7) %>%
  select(id, grid.id, depth, year, date, N) %>%
  left_join(., gridSim.df) %>%
  mutate(year=year(date),
         PAR=PAR * exp(-depth * KD)) %>%
  group_by(id, grid.id, depth) %>%
  mutate(N_rcr=c(scale(log(N+1))),
         PAR=c(scale(PAR)),
         SST=c(scale(SST))) %>%
  ungroup
rm(N.df); gc()
lag.df_jul <- N.df_jul %>% 
  group_by(id, grid.id, depth) %>%
  multijetlag(PAR, SST, n=nLags) %>%
  ungroup %>%
  full_join(stormIndexLag, by="year")
rm(N.df_jul); gc()
b.df_jul <- lag.df_jul %>% group_by(id, grid.id, depth) %>%
  summarise(across(starts_with("lag_"), ~coef(lm(N~.x))[2])) %>%
  pivot_longer(starts_with("lag_"), names_to="predictor", values_to="b") %>%
  mutate(lag=as.numeric(str_split_fixed(predictor, "_", 3)[,2]),
         covar=str_split_fixed(predictor, "_", 3)[,3])
saveRDS(b.df_jul, "temp\\b_lagEffects_recruits.rds")



# * lag subcanopy ---------------------------------------------------------

N.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}_subcanopy.rds")) %>%
  filter(id > ifelse(gridRes==0.1, 5, 2))
N.df_jul <- N.df %>%
  filter(month==7) %>%
  select(id, grid.id, depth, year, date, N) %>%
  left_join(., gridSim.df) %>%
  mutate(year=year(date),
         PAR=PAR * exp(-depth * KD)) %>%
  group_by(id, grid.id, depth) %>%
  mutate(N=c(scale(log(N+1))),
         PAR=c(scale(PAR)),
         SST=c(scale(SST))) %>%
  ungroup
rm(N.df); gc()
lag.df_jul <- N.df_jul %>% 
  group_by(id, grid.id, depth) %>%
  multijetlag(PAR, SST, n=nLags) %>%
  ungroup %>%
  full_join(stormIndexLag, by="year")
rm(N.df_jul); gc()
b.df_jul <- lag.df_jul %>% group_by(id, grid.id, depth) %>%
  summarise(across(starts_with("lag_"), ~coef(lm(N~.x))[2])) %>%
  pivot_longer(starts_with("lag_"), names_to="predictor", values_to="b") %>%
  mutate(lag=as.numeric(str_split_fixed(predictor, "_", 3)[,2]),
         covar=str_split_fixed(predictor, "_", 3)[,3])
saveRDS(b.df_jul, "temp\\b_lagEffects_subcanopy.rds")



# * lag canopy ------------------------------------------------------------

N.df <- readRDS(glue("summaries{sep}pop_df_{gridRes}_canopy.rds")) %>%
  filter(id > ifelse(gridRes==0.1, 5, 2)) %>% filter(id < 20)
N.df_jul <- N.df %>%
  filter(month==7) %>%
  select(id, grid.id, depth, year, date, N) %>%
  left_join(., gridSim.df) %>%
  mutate(year=year(date),
         PAR=PAR * exp(-depth * KD)) %>%
  group_by(id, grid.id, depth) %>%
  mutate(N=c(scale(log(N+1))),
         PAR=c(scale(PAR)),
         SST=c(scale(SST))) %>%
  ungroup
rm(N.df); gc()
lag.df_jul <- N.df_jul %>% 
  group_by(id, grid.id, depth) %>%
  multijetlag(PAR, SST, n=nLags) %>%
  ungroup %>%
  full_join(stormIndexLag, by="year")
rm(N.df_jul); gc()
b.df_jul <- lag.df_jul %>% group_by(id, grid.id, depth) %>%
  summarise(across(starts_with("lag_"), ~coef(lm(N~.x))[2])) %>%
  pivot_longer(starts_with("lag_"), names_to="predictor", values_to="b") %>%
  mutate(lag=as.numeric(str_split_fixed(predictor, "_", 3)[,2]),
         covar=str_split_fixed(predictor, "_", 3)[,3])
saveRDS(b.df_jul, "temp\\b_lagEffects_canopy.rds")


# * lag mass --------------------------------------------------------------

mass.df_jul <- mass.df %>% 
  # filter(id < 100) %>%
  filter(month==7) %>% 
  filter(grid.id < 100) %>%
  select(id, grid.id, depth, year, date, biomass) %>%
  left_join(., gridSim.df) %>%
  mutate(year=year(date),
         PAR=PAR * exp(-depth * KD)) %>%
  group_by(id, grid.id, depth) %>%
  mutate(biomass=c(scale(biomass)),
         PAR=c(scale(PAR)),
         SST=c(scale(SST))) %>%
  ungroup
rm(mass.df); rm(gridSim.df); gc()
lag.df_jul <- mass.df_jul %>% 
  group_by(id, grid.id, depth) %>%
  multijetlag(PAR, SST, n=nLags) %>%
  ungroup %>%
  full_join(stormIndexLag, by="year")
rm(mass.df_jul); gc()
b.df_jul <- lag.df_jul %>% group_by(id, grid.id, depth) %>%
  summarise(across(starts_with("lag_"), ~coef(lm(biomass~.x))[2])) %>%
  pivot_longer(starts_with("lag_"), names_to="predictor", values_to="b") %>%
  mutate(lag=as.numeric(str_split_fixed(predictor, "_", 3)[,2]),
         covar=str_split_fixed(predictor, "_", 3)[,3])
saveRDS(b.df_jul, "temp\\b_lagEffects_biomass.rds")




# * lag plots -------------------------------------------------------------

lag_foc <- c("recruits", "subcanopy", "canopy", "biomass")[4]
b.df_jul <- readRDS(glue("temp{sep}b_lagEffects_{lag_foc}.rds"))
b.sum_jul <- b.df_jul %>% group_by(lag, depth, covar) %>%
  summarise(b_mn=mean(b), b_md=median(b),
            b_q10=quantile(b, probs=0.1),
            b_q90=quantile(b, probs=0.9))

b.sum_jul %>% filter(depth %in% c(2,5,10,15)) %>%
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         covar=factor(covar, levels=c("strm", "PAR", "SST"),
                      labels=c("Storms", "PAR", "SST"))) %>%
  ggplot(aes(lag, b_md, ymin=b_q10, ymax=b_q90, colour=covar, fill=covar)) +
  geom_hline(yintercept=0, colour="grey80") +
  geom_ribbon(alpha=0.3, colour=NA) +
  geom_line() +
  scale_colour_manual("", values=c("grey40", "goldenrod", "red4")) +
  scale_fill_manual("", values=c("grey40", "goldenrod", "red4")) +
  facet_grid(depth~.) + 
  scale_x_continuous(breaks=c(5, 10, 15)) +
  scale_y_continuous(breaks=c(-0.5, 0, 0.5)) +
  theme_classic() + theme(legend.position="bottom") +
  labs(x="Lag (years)", 
       y=glue("Effect on {lag_foc}\n(standardized slope: median + middle 80%)"))
ggsave(glue("figs{sep}pub{sep}varLagEffects_ribbon_{lag_foc}.png"), 
       width=3, height=7, dpi=300)

b.sum_jul_fetch <- b.df_jul %>% 
  right_join(., grid.i %>% select(id, fetchCat)) %>%
  group_by(lag, depth, fetchCat, covar) %>%
  summarise(b_mn=mean(b), b_md=median(b),
            b_q10=quantile(b, probs=0.1),
            b_q90=quantile(b, probs=0.9))
b.sum_jul_fetch %>% filter(depth %in% c(2,5,10,15)) %>%
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         covar=factor(covar, levels=c("strm", "PAR", "SST"),
                      labels=c("Storms", "PAR", "SST"))) %>%
  ggplot(aes(lag, b_md, ymin=b_q10, ymax=b_q90, colour=covar, fill=covar)) +
  geom_hline(yintercept=0, colour="grey80") +
  geom_ribbon(alpha=0.3, colour=NA) +
  geom_line() +
  scale_colour_manual("", values=c("grey40", "goldenrod", "red4")) +
  scale_fill_manual("", values=c("grey40", "goldenrod", "red4")) +
  facet_grid(depth~fetchCat) + 
  scale_x_continuous(breaks=c(5, 10, 15)) +
  theme_classic() + theme(legend.position="bottom") +
  labs(x="Lag (years)", 
       y=glue("Effect on {lag_foc}\n(standardized slope: median + middle 80%)"))
ggsave(glue("figs{sep}pub{sep}varLagEffects_fetch_ribbon_{lag_foc}.png"),
       width=4, height=7, dpi=300)

b.sum.id_jul <- b.df_jul %>% group_by(id, depth, covar, lag) %>%
  summarise(b_mn=mean(b), b_md=median(b),
            b_q10=quantile(b, probs=0.1),
            b_q90=quantile(b, probs=0.9)) 
  
b_rng <- range(b.sum.id_jul$b_md)

b.amplitude <- b.sum.id_jul %>%
  filter(depth < 15) %>%
  group_by(covar, id, depth) %>%
  summarise(b_amplitude=max(b_mn)-min(b_mn)) %>%
  mutate(depth=factor(paste0(depth, "m"), 
                      levels=paste0(c(2,5,10,15,20), "m")),
         covar=factor(covar, levels=c("strm", "PAR", "SST"),
                      labels=c("Storms", "PAR", "SST"))) %>%
  full_join(grid.sf, ., by="id")

map_base.gg +
  geom_sf(data=b.amplitude, aes(fill=b_amplitude), colour=NA) +
  scale_fill_viridis_c("Effect amplitude", option="inferno", limits=c(0, NA),
                       guide=guide_colorbar(title.position="top")) + 
  facet_grid(depth~covar) +
  theme_classic() +
  theme(legend.position="bottom", 
        legend.key.height=unit(0.2, "cm"),
        legend.key.width=unit(1.2, "cm"),
        title=element_text(size=10)) +
  ggtitle(glue("Interannual fluctuation: {lag_foc}"))
ggsave(glue("figs{sep}pub{sep}covar_effect_amplitude_{lag_foc}.png"), 
       width=4, height=6, dpi=300)

b.amplitude %>% st_drop_geometry() %>%
  mutate(fetchCat=factor(fetchCat, labels=c("low", "high"))) %>%
  ggplot(aes(fetch, b_amplitude)) + 
  geom_point(shape=1, alpha=0.2, size=0.75) + 
  stat_smooth(method="lm", aes(group=fetchCat, colour=fetchCat)) +
  scale_colour_viridis_d("Exposure\nlevel", option="B", end=0.75) +
  facet_grid(depth~covar) + 
  theme_classic() +
  labs(x=bquote(log[10](fetch)), y="Effect amplitude") +
  theme(title=element_text(size=10)) +
  ggtitle(glue("Interannual fluctuation: {lag_foc}"))
ggsave(glue("figs{sep}pub{sep}covar_effect_fetch_amplitude_{lag_foc}.png"),
       width=6, height=6, dpi=300)






# covariates --------------------------------------------------------------



fetch_thresh <- 4.02 # bad practice, but based on Pedersen sites
grid.i %>%
  mutate(fetchCat=factor(fetchCat, labels=c("Low", "High"))) %>%
  ggplot(aes(fetch, fill=fetchCat)) + 
  geom_vline(xintercept=fetch_thresh, linetype=2, colour="grey30") +
  geom_histogram(bins=80) +
  scale_fill_viridis_d("Exposure\ncategory", option="B", end=0.75) +
  theme_classic() +
  labs(x=expression(log[10](fetch)))
ggsave(glue("figs{sep}pub{sep}fetchDistribution.png"), height=4, width=6, dpi=300)


ggplot(data.ls$year_stormIndex, aes(year, stormIndex)) + 
  stat_smooth(method="loess", colour="grey", se=F, span=0.5, size=0.5) +
  geom_point() + geom_line(size=0.25) +
  theme_classic() +
  labs(x="Year", y="Storm index")
ggsave(glue("figs{sep}pub{sep}stormIndex.png"), height=4, width=6, dpi=300)


grid.long <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  filter(id > ifelse(gridRes==0.1, 5, 2)) %>%
  rename(PAR=PAR_mn, SST=sstDay_mn, KD=KD_mn) %>%
  mutate(PAR2=PAR*exp(-KD*2), PAR5=PAR*exp(-KD*5), PAR10=PAR*exp(-KD*10),
         PAR15=PAR*exp(-KD*15), PAR20=PAR*exp(-KD*20)) %>%
  select(id, SST, fetch, fetchCat, PAR, PAR2, PAR5, PAR10, PAR15, PAR20) %>%
  pivot_longer(2:10, names_to="predictor", values_to="value") %>%
  mutate(predictor=factor(predictor, 
                          levels=c("SST", "fetch", "fetchCat", "PAR", "PAR2",
                                   "PAR5", "PAR10", "PAR15", "PAR20"),
                          labels=c("SST (growing season)", "ln(fetch)", 
                                   "Exposure category",
                                   "PAR (surface)", "PAR (2m)", "PAR (5m)",
                                   "PAR (10m)", "PAR (15m)", "PAR (20m)")))

sst.p <- map_base.gg +
  geom_sf(data=filter(grid.long, predictor=="SST (growing season)"), 
          colour=NA, aes(fill=value)) +
  theme_classic() + scale_fill_viridis_c("°C", option="C") + 
  ggtitle("SST (growing season)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
fetch.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="ln(fetch)"), 
          colour=NA, aes(fill=value)) +
  theme_classic() + scale_fill_viridis_c("", option="B", end=0.9) + 
  ggtitle(expression(log[10](wave~fetch))) + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
fetchCat.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="Exposure category"), 
          colour=NA, aes(fill=value)) +
  theme_classic() + 
  scale_fill_viridis_c("", option="B", end=0.75, guide=guide_legend(), breaks=1:2, labels=c("low", "high")) + 
  ggtitle("Exposure category") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR0.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (surface)"), 
          colour=NA, aes(fill=value)) +
  # theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (surface)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR2.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (2m)"), 
          colour=NA, aes(fill=value)) +
  # theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (2m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR5.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (5m)"), 
          colour=NA, aes(fill=value)) +
  # theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (5m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR10.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (10m)"), 
          colour=NA, aes(fill=value)) +
  # theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (10m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR15.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (15m)"), 
          colour=NA, aes(fill=value)) +
  # theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (15m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR20.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (20m)"), 
          colour=NA, aes(fill=value)) +
  # theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  theme_classic() + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (20m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
fig_cov_a <- ggarrange(sst.p, fetch.p, fetchCat.p, nrow=1)
fig_cov_b <- ggarrange(PAR0.p, PAR2.p, PAR5.p, PAR10.p, PAR15.p, PAR20.p, nrow=2, ncol=3,
                       common.legend=T, legend="bottom")
fig_cov_b_alt <- ggarrange(PAR0.p, PAR2.p, PAR5.p, PAR10.p, PAR15.p, PAR20.p, nrow=2, ncol=3,
                       common.legend=F, legend="bottom")
ggsave(glue("figs{sep}pub{sep}env_a.png"), fig_cov_a, width=7, height=4, dpi=300)
ggsave(glue("figs{sep}pub{sep}env_b.png"), fig_cov_b, width=7.5, height=8, dpi=300)
ggsave(glue("figs{sep}pub{sep}env_b_alt.png"), fig_cov_b_alt, width=8.5, height=8, dpi=300)



# WARNING! THIS HAS NOT BEEN UPDATED AND WILL NOT WORK!
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












# parameter distributions -------------------------------------------------

# storm effects on survival rates
par.rng <- readRDS(glue("{sens.dir}000_parameter_ranges.rds")) %>%
  mutate(exposure=case_when(exposure=="low"~1, exposure=="high"~2)) 
surv.rng <- par.rng %>% filter(param=="surv") %>%
  group_by(exposure) %>% group_split()
loss.rng <- par.rng %>% filter(param=="loss")


storm_effect.df <- data.ls$year_stormIndex %>%
  select(year, stormIndex) %>%
  mutate(surv_canopy_1=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[1]], stage=="canopy")$shp1,
                             filter(surv.rng[[1]], stage=="canopy")$shp2),
         surv_canopy_2=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[2]], stage=="canopy")$shp1,
                             filter(surv.rng[[2]], stage=="canopy")$shp2),
         surv_subcanopy_1=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[1]], stage=="subcanopy")$shp1,
                             filter(surv.rng[[1]], stage=="subcanopy")$shp2),
         surv_subcanopy_2=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[2]], stage=="subcanopy")$shp1,
                             filter(surv.rng[[2]], stage=="subcanopy")$shp2),
         surv_recruit_1=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[1]], stage=="recruits")$shp1,
                             filter(surv.rng[[1]], stage=="recruits")$shp2),
         surv_recruit_2=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[2]], stage=="recruits")$shp1,
                             filter(surv.rng[[2]], stage=="recruits")$shp2)) %>%
  pivot_longer(starts_with("surv_"), names_to="param", values_to="value") %>%
  mutate(stage=str_split_fixed(param, "_", 3)[,2] %>%
           factor(levels=c("recruit", "subcanopy", "canopy")),
         exposure=str_split_fixed(param, "_", 3)[,3] %>% 
           factor(labels=c("Low", "High"))) %>%
  mutate(param=paste0("italic(s[", stage, "])"),
         value=sqrt(value)) %>%
  bind_rows(data.ls$year_stormIndex %>%
              select(year, stormIndex) %>%
              mutate(param="italic(epsilon)",
                     exposure="Low",
                     value=qbeta(pnorm(stormIndex, 0, 1),
                                   loss.rng$shp1, loss.rng$shp2))) %>%
  bind_rows(data.ls$year_stormIndex %>%
              select(year, stormIndex) %>%
              mutate(param="italic(epsilon)",
                     exposure="High",
                     value=qbeta(pnorm(stormIndex, 0, 1),
                                 loss.rng$shp1, loss.rng$shp2))) %>%
  mutate(param=factor(param, levels=rev(unique(param))),
         exposure=factor(exposure, levels=c("Low", "High")))

ggplot(storm_effect.df, aes(year, value, colour=param)) + 
  geom_line() + 
  theme_classic() +
  scale_x_continuous(breaks=c(1950, 1975, 2000)) +
  scale_y_continuous("Rate (non-growing season)", breaks=c(0, 0.5, 1)) +
  scale_colour_manual("Parameter", 
                      values=c("#bf812d", "#80cdc1", "#35978f", "#003c30"), 
                      labels=scales::parse_format()) +
  theme(axis.title.x=element_blank(),
        legend.text.align=0) +
  facet_grid(.~exposure) 
ggsave(glue("figs{sep}pub{sep}stormEffect_parameters.png"), width=7, height=2.5, dpi=300)



fetch_cols <- c("grey70", "black")
stage_cols <- c(recruits="#80cdc1", subcanopy="#35978f", canopy="#003c30")
par.rng <- par.rng 

# s_recruits
s_r <- list(xlim=c(0,1), main="Recruit survival rate",
            xlab=expression(italic(s[recruits])), ylab="density")
s_r$x.full <- seq(s_r$xlim[1], 0.9, length.out=1e3)
s_r$c.full <- map(1:2, 
                  ~with(filter(par.rng, param=="surv" & stage=="recruits" & exposureNum==.x),
                        dbeta(s_r$x.full, shp1, shp2)))
s_r$x.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="recruits" & exposureNum==.x),
                           c(valMin, valMax)))
s_r$c.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="recruits" & exposureNum==.x),
                           dbeta(s_r$x.95[[.x]], shp1, shp2)))
s_r$ylim <- c(0, max(unlist(s_r$c.full)[is.finite(unlist(s_r$c.full))]))

# s_subcanopy
s_s <- list(xlim=c(0,1),main="Subcanopy survival rate",
            xlab=expression(italic(s[subcanopy])), ylab="density")
s_s$x.full <- seq(s_s$xlim[1], s_s$xlim[2], length.out=1e3)
s_s$c.full <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="subcanopy" & exposureNum==.x),
                             dbeta(s_s$x.full, shp1, shp2)))
s_s$x.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="subcanopy" & exposureNum==.x),
                           c(valMin, valMax)))
s_s$c.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="subcanopy" & exposureNum==.x),
                           dbeta(s_s$x.95[[.x]], shp1, shp2)))
s_s$ylim <- c(0, max(unlist(s_s$c.full)[is.finite(unlist(s_s$c.full))]))

# s_canopy
s_c <- list(xlim=c(0,1), main="Canopy survival rate",
            xlab=expression(italic(s[canopy])), ylab="density")
s_c$x.full <- seq(s_c$xlim[1], s_c$xlim[2], length.out=5e2)
s_c$c.full <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="canopy" & exposureNum==.x),
                             dbeta(s_c$x.full, shp1, shp2)))
s_c$x.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="canopy" & exposureNum==.x),
                           c(valMin, valMax)))
s_c$c.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="canopy" & exposureNum==.x),
                           dbeta(s_c$x.95[[.x]], shp1, shp2)))
s_c$ylim <- c(0, max(unlist(s_c$c.full)[is.finite(unlist(s_c$c.full))]))

# z
z <- list(xlim=c(0,5e3), main="Settlement rate",
          xlab=expression(italic(z)~(N/m^2)), ylab="density")
z$x.full <- seq(z$xlim[1], z$xlim[2], length.out=5e2)
z$c.full <- map(1:2, ~with(filter(par.rng, param=="settlement" & exposureNum==.x),
                             dnorm(z$x.full, valMean, valSD)))
z$x.95 <- map(1:2, ~with(filter(par.rng, param=="settlement" & exposureNum==.x),
                           c(valMin, valMax)))
z$c.95 <- map(1:2, ~with(filter(par.rng, param=="settlement" & exposureNum==.x),
                           dnorm(z$x.95[[.x]], valMean, valSD)))
z$ylim <- c(0, max(unlist(z$c.full)))

# growStipe
gamma <- list(xlim=c(0,300), main="Stipe growth rate",
              xlab=expression(italic(gamma[x])~(mm/year)), ylab="density")
gamma$x.full <- seq(gamma$xlim[1], gamma$xlim[2], length.out=5e2)
gamma$c.full <- map(1:3, ~with(filter(par.rng, param=="growStipe" & stage==names(stage_cols)[.x]),
                           dnorm(gamma$x.full, valMean, valSD)))
gamma$x.95 <- map(1:3, ~with(filter(par.rng, param=="growStipe" & stage==names(stage_cols)[.x]),
                         c(valMin, valMax)))
gamma$c.95 <- map(1:3, ~with(filter(par.rng, param=="growStipe" & stage==names(stage_cols)[.x]),
                         dnorm(gamma$x.95[[.x]], valMean, valSD)))
gamma$ylim <- c(0, max(unlist(gamma$c.full)))

# growFrond
omega <- list(xlim=c(0,0.6), main="Frond growth rate",
              xlab=expression(italic(omega[x])~(m^2/year)), ylab="density")
omega$x.full <- seq(omega$xlim[1], omega$xlim[2], length.out=5e2)
omega$c.full <- map(1:3, ~with(filter(par.rng, param=="growFrond" & stage==names(stage_cols)[.x]),
                               dnorm(omega$x.full, valMean, valSD)))
omega$x.95 <- map(1:3, ~with(filter(par.rng, param=="growFrond" & stage==names(stage_cols)[.x]),
                             c(valMin, valMax)))
omega$c.95 <- map(1:3, ~with(filter(par.rng, param=="growFrond" & stage==names(stage_cols)[.x]),
                             dnorm(omega$x.95[[.x]], valMean, valSD)))
omega$ylim <- c(0, max(unlist(omega$c.full)))

# erosion
epsilon <- list(xlim=c(0,1), main="Frond erosion rate",
                xlab=expression(italic(epsilon)~(prop.~frond/non-growing~season)), ylab="density")
epsilon$x.full <- seq(epsilon$xlim[1], epsilon$xlim[2], length.out=5e2)
epsilon$c.full <- with(filter(par.rng, param=="loss"), 
                       dbeta(epsilon$x.full, shp1, shp2))
epsilon$x.95 <- with(filter(par.rng, param=="loss"), c(valMin, valMax))
epsilon$c.95 <- with(filter(par.rng, param=="loss"),
                             dbeta(epsilon$x.95, shp1, shp2))
epsilon$ylim <- c(0, max(unlist(epsilon$c.full)))

# density effect shape
theta <- list(xlim=c(0,2), main="Density effect shape",
              xlab=expression(italic(theta)), ylab="density")
theta$x.full <- seq(theta$xlim[1], theta$xlim[2], length.out=5e2)
theta$c.full <- dnorm(theta$x.full, 1, 0.25)
theta$x.95 <- with(filter(par.rng, param=="densityEffShape"), c(valMin, 1.5))
theta$c.95 <- with(filter(par.rng, param=="densityEffShape"), dnorm(theta$x.95, 1, 0.25))
theta$ylim <- c(0, max(unlist(theta$c.full)))




png(glue("figs{sep}pub{sep}parameter_distributions.png"), width=9, height=4.5, res=300, units="in")
{
  par(mfrow=c(2,4), mar=c(5,4,2,1)+0.1)

  # s_recruits
  plot(NA, NA, xlim=s_r$xlim, ylim=s_r$ylim, main=s_r$main, xlab=s_r$xlab, ylab=s_r$ylab)
  for(i in 1:2) {
    segments(x0=s_r$x.95[[i]], y0=0, y1=s_r$c.95[[i]], col=fetch_cols[i])
    lines(s_r$x.full, s_r$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low",  "High"), bty="n", title="Wave fetch")
  
  # s_subcanopy
  plot(NA, NA, xlim=s_s$xlim, ylim=s_s$ylim, main=s_s$main, xlab=s_s$xlab, ylab=s_s$ylab)
  for(i in 1:2) {
    segments(x0=s_s$x.95[[i]], y0=0, y1=s_s$c.95[[i]], col=fetch_cols[i])
    lines(s_s$x.full, s_s$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low", "High"), bty="n", title="Wave fetch")
  
  # s_canopy
  plot(NA, NA, xlim=s_c$xlim, ylim=s_c$ylim, main=s_c$main, xlab=s_c$xlab, ylab=s_c$ylab)
  for(i in 1:2) {
    segments(x0=s_c$x.95[[i]], y0=0, y1=s_c$c.95[[i]], col=fetch_cols[i])
    lines(s_c$x.full, s_c$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low", "High"), bty="n", title="Wave fetch")
  
  # z
  plot(NA, NA, xlim=z$xlim, ylim=z$ylim, main=z$main, xlab=z$xlab, ylab=z$ylab)
  for(i in 1:2) {
    segments(x0=z$x.95[[i]], y0=0, y1=z$c.95[[i]], col=fetch_cols[i])
    lines(z$x.full, z$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low", "High"), bty="n", title="Wave fetch")
  
  # gamma
  plot(NA, NA, xlim=gamma$xlim, ylim=gamma$ylim, main=gamma$main, xlab=gamma$xlab, ylab=gamma$ylab)
  for(i in 1:3) {
    segments(x0=gamma$x.95[[i]], y0=0, y1=gamma$c.95[[i]], col=stage_cols[i])
    lines(gamma$x.full, gamma$c.full[[i]], col=stage_cols[i])    
  }
  legend("topright", col=stage_cols, lty=1, names(stage_cols), bty="n", title="Stage")
  
  # omega
  plot(NA, NA, xlim=omega$xlim, ylim=omega$ylim, main=omega$main, xlab=omega$xlab, ylab=omega$ylab)
  for(i in 1:3) {
    segments(x0=omega$x.95[[i]], y0=0, y1=omega$c.95[[i]], col=stage_cols[i])
    lines(omega$x.full, omega$c.full[[i]], col=stage_cols[i])    
  }
  legend("topright", col=stage_cols, lty=1, names(stage_cols), bty="n", title="Stage")
  
  # epsilon
  plot(NA, NA, xlim=epsilon$xlim, ylim=epsilon$ylim, main=epsilon$main, xlab=epsilon$xlab, ylab=epsilon$ylab)
  segments(x0=epsilon$x.95, y0=0, y1=epsilon$c.95)
  lines(epsilon$x.full, epsilon$c.full)    
  
  # theta
  plot(NA, NA, xlim=theta$xlim, ylim=theta$ylim, main=theta$main, xlab=theta$xlab, ylab=theta$ylab)
  segments(x0=theta$x.95, y0=0, y1=theta$c.95)
  lines(theta$x.full, theta$c.full)    
}
dev.off()








# map: site locations -----------------------------------------------------

sites.sf <- bind_rows(
  read_csv(glue("data{sep}raw{sep}collab{sep}collab_sites.csv"), skip=1) %>%
    select(lat, lon) %>% mutate(source="collab"),
  read_csv(glue("data{sep}raw{sep}digitized{sep}sitesDigitized.csv")) %>%
    select(lat, lon) %>% mutate(source="digitized")
  ) %>% 
  filter(complete.cases(.)) %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326)

map_base.gg +
  geom_sf(data=grid.sf, colour=NA, fill="#80cdc1") + 
  geom_sf(data=sites.sf %>% st_crop(grid.sf), colour="#8c510a") +
  theme_classic() +
  ggtitle("UK data source locations")
ggsave(glue("figs{sep}pub{sep}map_data_locations.png"), width=3, height=4.5, dpi=300)










# observed vs modelled ----------------------------------------------------

mass.sum.sf_canopy %>%
  filter(id %in% obs.ls$N_canopy.lm$id) %>%
  ggplot(aes(PAR_atDepth_mn, biomass_mn)) + 
  geom_point() + stat_smooth(method="lm") +
  geom_point(data=obs.ls$N_canopy.lm %>% group_by(id) %>%
               summarise(PAR_atDepth=mean(PAR_atDepth), 
                         N=median(N)), 
             aes(PAR_atDepth, N), colour="red") +
  stat_smooth(data=obs.ls$N_canopy.lm, method="lm",
             aes(PAR_atDepth, N*0.5), colour="red")


pop.sum.sf_canopy %>%
  filter(id %in% obs.ls$N_canopy.lm$id) %>%
  ggplot(aes(PAR_atDepth_mn, N_max)) + 
  geom_point() + stat_smooth(method="lm") +
  # geom_point(data=obs.ls$N_canopy.lm %>% group_by(id) %>%
  #              summarise(PAR_atDepth=mean(PAR_atDepth), 
  #                        N=median(N)), 
  #            aes(PAR_atDepth, N), colour="red") +
  geom_point(data=obs.ls$N_canopy.lm, 
             aes(PAR_atDepth, N), colour="red") +
  stat_smooth(data=obs.ls$N_canopy.lm, method="lm",
              aes(PAR_atDepth, N), colour="red")

