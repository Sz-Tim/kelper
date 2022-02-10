# KELPER
# Make grids
# Tim Szewczyk

# This script constructs the grids to use in the simulations and extracts the
# mean environmental conditions within each cell.




########
##-- set up

# libraries and local functions
pkgs <- c("raster", "lubridate", "tidyverse", "sf")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
options(mc.cores=4)

# switches
PAR_datasource <- c("MODIS", "POWER")[1]
gridRes <- 1 # currently in arc-seconds

# directories
gis.dir <- "..\\..\\00_gis\\"

# maximum extent for covariates
UK_bbox <- st_bbox(c(xmin=-11, xmax=2, ymin=49, ymax=61), crs=st_crs(4326))

# datasets
covars.ls <- loadCovariates(gis.dir, UK_bbox, zscore=T, saveFile="data\\covarZ_ls.rds")

# make grid
grid.domain <- UK_bbox %>%
  st_make_grid(cellsize=c(gridRes, gridRes)) %>%
  st_sf(id=1:length(.)) %>%
  extractCovarsToGrid(., covars.ls, PAR_datasource) %>%
  filter(!is.na(KD_mn) & !is.na(PAR_surface) & !is.na(fetch)) %>%
  mutate(id=row_number())

# save as shp
st_write(grid.domain, paste0("data\\grid_", gridRes, ".gpkg"), append=F)
