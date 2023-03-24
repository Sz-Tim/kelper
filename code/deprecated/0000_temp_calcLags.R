


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
sens.dir <- glue("out{sep}sensitivity{sep}")

# switches & settings
gridRes <- 0.1

# grid
data.ls_strm <- compileDatasets(data.dir, supp.f)$year_stormIndex






# lagged effects ----------------------------------------------------------

nLags <- 15

stormIndexLag <- data.ls_strm %>% select(year) %>%
  bind_cols(imap_dfc(setNames(1:nLags, glue("lag_{1:nLags}_strm")), 
                     ~lag(data.ls_strm$stormIndex, .x)))

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
saveRDS(b.df_jul, glue("temp{sep}b_lagEffects_canopy.rds"))


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

