# Tim Szewczyk

# This script makes the figures for publication




# set up ------------------------------------------------------------------

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "ggpubr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00[0-9]_.*R", full.names=T), source)
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
out.eff <- map(lm.fit, ~conditional_effects(.x, surface=T))

ggpubr::ggarrange(plotlist=plot(out.eff[[1]], ask=F, points=T, plot=F, stype="raster",
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1, widths=c(1,1,1.2))
ggsave(glue("figs{sep}regr{sep}condEff_1_surface_{gridRes}.png"), width=12, height=4)

ggpubr::ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[2]], 
                                    effects=c("logLenStipe", "lPAR_atDepth", "SST")), 
                ask=F, points=T, plot=F, stype="raster",
                point_args=list(shape=1, size=0.75, alpha=0.8)), 
  nrow=1, widths=c(1,1,1))
ggsave(glue("figs{sep}regr{sep}condEff_2a_surface_{gridRes}.png"), width=12, height=4)
ggpubr::ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[2]], surface=T,
                                    effects=c("logLenStipe:lPAR_atDepth", "lPAR_atDepth:SST")), 
                ask=F, points=T, plot=F, stype="raster",
                point_args=list(shape=1, size=0.75, alpha=0.8)), 
  nrow=1, widths=c(1,1))
ggsave(glue("figs{sep}regr{sep}condEff_2b_surface_{gridRes}.png"), width=9, height=4)

ggpubr::ggarrange(plotlist=plot(out.eff[[3]], ask=F, points=T, plot=F, stype="raster",
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1, widths=c(1,1,1.2))
ggsave(glue("figs{sep}regr{sep}condEff_4_surface_{gridRes}.png"), width=12, height=4)

ggpubr::ggarrange(plotlist=plot(out.eff[[4]], ask=F, points=T, plot=F, stype="raster",
                                point_args=list(shape=1, size=0.75, alpha=0.8)), 
                  nrow=1, widths=c(1,1,1.2))
ggsave(glue("figs{sep}regr{sep}condEff_3_surface_{gridRes}.png"), width=12, height=4)

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
  plotlist=plot(conditional_effects(lm.fit[[7]], surface=T,
                                    effects=c("PAR_atDepth:SST", "PAR_atDepth:fetch")), ask=F, 
                points=T, plot=F, stype="raster",
                point_args=list(shape=1, size=0.75, alpha=0.8)),
  nrow=1)
ggsave(glue("figs{sep}regr{sep}condEff_7b_surface_{gridRes}.png"), width=9, height=4)








# biomass broad patterns --------------------------------------------------

mass.sum <- readRDS(glue("summaries{sep}mass_sum_{gridRes}.rds")) %>% 
  filter(id > ifelse(gridRes==0.1, 5, 2)) %>%
  left_join(., gridSim.sum)
mass.sum.sf_canopy <- mass.sum.sf %>% 
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])

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
