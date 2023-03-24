

library(tidyverse); library(stars)
gis.dir <- "..\\..\\00_gis\\"


UK_bbox <- st_bbox(c(xmin=-11, xmax=2, ymin=49, ymax=61), crs=st_crs(4326))
coast <- read_sf(paste0(gis.dir, "coastline\\ne_10m_land.shp")) %>%
  add_row(read_sf(paste0(gis.dir, "coastline\\ne_10m_minor_islands.shp"))) %>%
  st_crop(UK_bbox)

# mean wind by season across years... could be useful for background erosion...
wind.st <- read_stars(dir(paste0(gis.dir, "wind\\HadUKgrid"), full.names=T)[-1], var="sfcWind")
plot(wind.st)

map(dir(paste0(gis.dir, "wind"), full.names=T)[-1], 
    ~read_ncdf(.x, var="sfcWind") %>% plot)


# super coarse. Not useful
precExt.st <- read_stars(dir(paste0(gis.dir, "precipExtreme"), "r20", full.names=T))
plot(precExt.st[,,,1400])


# predictions of 5 models, all extremely different....
surge.st <- read_stars(dir(paste0(gis.dir, "stormSurge"), "skew", full.names=T))
plot(surge.st)

# cool, but only really for England and mostly inland
wind.locs <- read_csv(dir(paste0(gis.dir, "wind\\TEMPEST"), "locations", full.names=T), 
                      skip=40, col_select=1:5) %>% 
  setNames(c("place", "lat", "lon", "country", "locType"))
wind.df <- read_csv(dir(paste0(gis.dir, "wind\\TEMPEST"), "events", full.names=T), 
                    skip=134, col_select=c(22, 23, 16, 17)) %>%
  setNames(c("place", "eventType", "month", "year")) %>%
  left_join(., wind.locs, by="place") %>%
  filter(!is.na(lat)) 
wind.sf <- st_as_sf(wind.df, coords=c("lon", "lat"), crs=4326)

ggplot(wind.sf) + geom_sf(data=coast) + geom_sf(alpha=0.5)

