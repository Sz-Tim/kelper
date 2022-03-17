# KELPER
# Process output
# Tim Szewczyk

# This script processes the output produces in 03_runSimulations.R





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
gridRes <- 0.25

# load files
grid.sf <- st_read(glue("data{sep}grid_{gridRes}_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
sim.info <- glue("...._{gridRes}")
pop.f <- dir(out.dir, glue("pop_{sim.info}"), full.names=T)
pop.df <- map_dfr(pop.f, readRDS)
mass.f <- dir(out.dir, glue("mass_{sim.info}"), full.names=T)
mass.df <- map_dfr(mass.f, readRDS)





# summarise output --------------------------------------------------------

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
mass.sum <- mass.df %>%
  # filter(year > 10) %>%
  select(-date, -year) %>%
  mutate(logBiomass=log(biomass+1)) %>%
  group_by(id, sim, month, depth, landscape, stochParams) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max, q90=~quantile(.x, probs=0.9))))





# save output -------------------------------------------------------------

saveRDS(pop.df, glue("summaries{sep}pop_df_{gridRes}.rds"))
saveRDS(mass.df, glue("summaries{sep}mass_df_{gridRes}.rds"))
saveRDS(pop.sum, glue("summaries{sep}pop_sum_{gridRes}.rds"))
saveRDS(mass.sum, glue("summaries{sep}mass_sum_{gridRes}.rds"))



