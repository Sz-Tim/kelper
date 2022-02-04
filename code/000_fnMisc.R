# KELPER
# Helper functions
# Tim Szewczyk


# Miscellaneous helper functions associated with the KELPER project

compileDatasets <- function(data.dir) {
  
  library(tidyverse)
  
  site.i <- read_csv(paste0(data.dir, "sitesDigitized.csv"), col_select=1:3, show_col_types=F)
  digitizedMetadata <- read_csv(paste0(data.dir, "metadata.csv"), show_col_types=F) %>%
    left_join(., site.i, by="location")
  fig.f <- unique(str_remove(dir(data.dir, "_fig.*csv"), "[a-z]?.csv"))
  tab.f <- unique(str_remove(dir(data.dir, "_table.*csv"), "[a-z]?.csv"))
  
  digitized <- c(fig.f, tab.f) %>% 
    map(~dir(data.dir, paste0(.x, "[a-z]?.csv")) %>%
          map_dfr(~read_csv(paste0(data.dir, .x), show_col_types=F) %>%
                    mutate(reference=str_remove(.x, ".csv"))) %>%
          left_join(., digitizedMetadata, by="reference", suffix=c("", ".meta")) %>%
          select(-contains(".meta"))) %>%
    setNames(c(fig.f, tab.f))
  
  fig_var <- imap_dfr(digitized, 
                      ~tibble(V1=names(.x)[1], V2=names(.x)[2], reference=.y)) %>%
    arrange(V1, V2) %>%
    mutate(vars=paste0(V1, "_", V2))
  
  unique_comparisons <- unique(fig_var$vars)
  varCombined.ls <- vector("list", length(unique_comparisons)) %>% 
    setNames(unique_comparisons)
  for(i in seq_along(varCombined.ls)) {
    refs <- fig_var$reference[fig_var$vars==unique_comparisons[i]]
    varCombined.ls[[i]] <- do.call('rbind', digitized[refs])
  }
  return(varCombined.ls)
}







#' Generate list of parameters for simulations
#' 
#' All rates must already be scaled to the timeframe of the transition matrix (e.g., growth rate for the season rather than daily).
#'
#' @param path Path to .rds file with stored parameters; if \code{NULL}, creates list instead
#' @param override If `path` is not NULL, named list of parameters to override (e.g., \code{list(N_stages=2)})
#' @param N_stages Number of stages
#' @param N_seasons Number of seasons per year
#' @param tmax Number of years
#' @param growthRateStipeMax Vector of maximum stipe growth rates during spring
#' @param growthRateFrond Vector of maximum frond area growth rates during spring 
#' @param frondAreaMax Carrying capacity for Frond Area Index
#' @param growthRateStipeDensityShape Shape parameter for stipe growth density dependence
#' @param sizeClassLimits $\delta$ Size limits between stages
#' @param survRate $\s$ Survival rate
#' @param settlementRateBg Background recruit settlement rate
#' @param lossRate Vector of frond loss rates during winter
#' @param prFullHarvest Proportion of population to be harvested in full
#' @param freqHarvest Harvest frequency (years)
#' @param harvestTarget Kelp part for biomass calculation (all, stipe, frond)
#' @param extraPars Any additional parameters to include (e.g., regressions)
#'
#' @return
#' @export
#'
#' @examples
setParameters <- function(path=NULL,
                          override=NULL,
                          tmax=20,
                          growthRateStipeMax=c(194, 195, 58),
                          growthRateFrond=c(1787, 2299, 3979)/1e4, # m^2/plant/spring
                          frondAreaMax=5500/1e4,
                          growthRateDensityShape=1.5,
                          sizeClassLimits=(1000 * (0:6)/6)[c(1,2,5,7)],
                          survRate=c(0.4, 0.7, 0.9),
                          settlementRateBg=100,
                          lossRate=c(0.1, 0.2, 0.3),
                          prFullHarvest=0.4,
                          freqHarvest=1,
                          harvestTarget=c("all", "stipe", "frond")[1],
                          extraPars=NULL
                          ) {
  
  if(!is.null(path)) {
    pars <- readRDS(path)
    if(!is.null(override)) {
      for(i in seq_along(override)) {
        pars[names(override)[i]] <- override[i]
      }
    }
  } else {
    pars <- list(
      
      # simulation 
      tmax=tmax,
      
      # growth
      growthRateStipeMax=growthRateStipeMax,
      growthRateFrond=growthRateFrond,
      frondAreaMax=frondAreaMax,
      growthRateDensityShape=growthRateDensityShape,
      sizeClassLimits=sizeClassLimits,
      sizeClassMdpts=(sizeClassLimits+lag(sizeClassLimits))[-1]/2,
      
      # survival
      survRate=survRate,
      
      # settlement
      settlementRateBg=settlementRateBg,
      
      # loss
      lossRate=lossRate,
      
      # harvest
      prFullHarvest=prFullHarvest,
      freqHarvest=freqHarvest,
      harvestTarget=harvestTarget
    ) 
  }
  
  if(!is.null(extraPars)) {
    for(i in seq_along(extraPars)) {
      pars[names(extraPars)[i]] <- extraPars[i]
    }
  }
  
  return(pars)
}






loadCovariates <- function(gis.dir=NULL, bbox=NULL, loadFile=NULL, saveFile=NULL) {
  if(is.null(loadFile)) {
    covars.ls <- list(
      fetch=dir(paste0(gis.dir, "dataGeog"), "^FetchUK_200m", full.names=T) %>% 
        read_csv(show_col_types=F) %>% st_as_sf(coords=c("OSEast", "OSNorth")) %>% 
        st_set_crs(27700) %>% st_transform(4326),
      sst_day=dir(paste0(gis.dir, "dataEnv\\climate"), 
                  "MODISA_L3m_SST_Monthly", full.names=T) %>%
        raster() %>% crop(bbox),
      sst_night=dir(paste0(gis.dir, "dataEnv\\climate"), 
                    "MODISA_L3m_NSST_Monthly", full.names=T) %>%
        raster() %>% crop(bbox),
      chla=dir(paste0(gis.dir, "dataEnv\\chla"), 
               "MODISA_L3m_CHL", full.names=T) %>%
        raster() %>% crop(bbox),
      KD_mn=dir(paste0(gis.dir, "dataEnv\\attenuation"), 
                "MODISA_L3m_KD.*11W_49N", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      KD_sd=dir(paste0(gis.dir, "dataEnv\\attenuation"), 
                "MODISA_L3m_KD.*11W_49N", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      PAR_mn=dir(paste0(gis.dir, "dataEnv\\light"), 
                 "MODISA_L3m_PAR", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      PAR_sd=dir(paste0(gis.dir, "dataEnv\\light"), 
                 "MODISA_L3m_PAR", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      irrad_longterm=dir(paste0(gis.dir, "dataEnv\\light"), 
                         "POWER_Regional_Climatology", full.names=T) %>%
        read_csv(skip=10, show_col_types=F) %>% 
        st_as_sf(coords=c("LON", "LAT"), crs=4326),
      irrad_monthly=dir(paste0(gis.dir, "dataEnv\\light"), 
                        "POWER_Regional_monthly", full.names=T) %>%
        read_csv(skip=10, show_col_types=F) %>% 
        filter(JAN != -999) %>% select(-ANN) %>%
        group_by(PARAMETER, YEAR) %>% mutate(id=row_number()) %>% ungroup %>%
        pivot_longer(5:16, names_to="MONTH", values_to="PAR") %>%
        st_as_sf(coords=c("LON", "LAT"), crs=4326),
      irrad_daily=dir(paste0(gis.dir, "dataEnv\\light"), 
                      "POWER_Regional_Daily", full.names=T) %>%
        map_dfr(., ~read_csv(.x, skip=10, show_col_types=F)) %>% 
        mutate(DATE=case_when(is.na(DOY) ~ ymd(paste(YEAR, MO, DY, sep="-")),
                              !is.na(DOY) ~ as_date(DOY, 
                                                    origin=paste0(YEAR,"-01-01")))) %>% 
        group_by(DATE) %>% mutate(id=row_number()) %>% ungroup %>%
        pivot_longer(contains("_PAR_"), names_to="PARAMETER", values_to="PAR") %>%
        st_as_sf(coords=c("LON", "LAT"), crs=4326)
    )
    
    irrad.grid <- (st_bbox(covars.ls$irrad_longterm) + .25*c(-1,-1,1,1)) %>%
      st_make_grid(cellsize=c(0.5, 0.5)) %>% st_sf(id=1:length(.))
    
    covars.ls$irrad_growing.month <- covars.ls$irrad_monthly %>%
      filter(MONTH %in% c("JAN", "FEB", "MAR", "APR", "MAY", "JUN")) %>%
      group_by(PARAMETER, YEAR, id) %>%
      summarise(PAR=mean(PAR)) %>%
      group_by(PARAMETER, id) %>%
      summarise(mnPAR=mean(PAR), sdPAR=sd(PAR)) %>%
      st_join(irrad.grid, .) %>%
      filter(grepl("ALLSKY", PARAMETER))
    if(!is.null(saveFile)) {
      saveRDS(covars.ls, saveFile)
    }
  } else {
    covars.ls <- readRDS(loadFile)
  }
  return(covars.ls)
}


extractCovarsToDatasets <- function(data.ls, covars.ls, PAR_datasource) {
  fetch_thirds <- quantile(covars.ls$fetch$fetchsum, probs=(1:2)/3)
  for(i in seq_along(data.ls)) {
    if(any(!is.na(data.ls[[i]]$lat))) {
      i_sf <- data.ls[[i]] %>% 
        filter(!is.na(lat) & !is.na(lon)) %>%
        st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
        mutate(sst_day=raster::extract(covars.ls$sst_day, ., 
                                       fun=mean, small=T, method="bilinear"),
               sst_night=raster::extract(covars.ls$sst_night, ., 
                                         fun=mean, small=T, method="bilinear"),
               chla=raster::extract(covars.ls$chla, ., 
                                    fun=mean, small=T, method="bilinear"),
               KD_mn=raster::extract(covars.ls$KD_mn, ., 
                                     fun=mean, small=T, method="bilinear"),
               KD_sd=raster::extract(covars.ls$KD_sd, ., 
                                     fun=mean, small=T, method="bilinear"),
               PAR_mn=raster::extract(covars.ls$PAR_mn, ., 
                                      fun=mean, small=T, method="bilinear"),
               PAR_sd=raster::extract(covars.ls$PAR_sd, ., 
                                      fun=mean, small=T, method="bilinear"),
               PAR_POWER=st_join(., 
                                 covars.ls$irrad_growing.month, 
                                 left=T, join=st_nearest_feature)$mnPAR,
               fetch=st_join(., 
                             covars.ls$fetch, 
                             left=T, join=st_nearest_feature)$fetchsum,
               fetchCat=case_when(fetch < fetch_thirds[1] ~ 1,
                                  between(fetch, fetch_thirds[1], fetch_thirds[2]) ~ 2,
                                  fetch > fetch_thirds[2] ~ 3))
      data.ls[[i]] <- left_join(data.ls[[i]], st_drop_geometry(i_sf))
      if(any(!is.na(data.ls[[i]]$depth))) {
        data.ls[[i]] <- data.ls[[i]] %>%
          mutate(PAR_atDepth_POWER=PAR_POWER * exp(-KD_mn * depth),
                 PAR_atDepth_MODIS=PAR_mn * exp(-KD_mn * depth))
        data.ls[[i]]$PAR_atDepth <- data.ls[[i]][[paste0("PAR_atDepth_", PAR_datasource)]]
      }
    }
  }
  return(data.ls)
}

extractCovarsToPts <- function(site.i, covars.ls, PAR_datasource) {
  fetch_thirds <- quantile(covars.ls$fetch$fetchsum, probs=(1:2)/3)
  site.sf <- site.i %>% 
    filter(!is.na(lat) & !is.na(lon)) %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    mutate(sst_day=raster::extract(covars.ls$sst_day, ., 
                                   fun=mean, small=T, method="bilinear"),
           sst_night=raster::extract(covars.ls$sst_night, ., 
                                     fun=mean, small=T, method="bilinear"),
           chla=raster::extract(covars.ls$chla, ., 
                                fun=mean, small=T, method="bilinear"),
           KD_mn=raster::extract(covars.ls$KD_mn, ., 
                                 fun=mean, small=T, method="bilinear"),
           KD_sd=raster::extract(covars.ls$KD_sd, ., 
                                 fun=mean, small=T, method="bilinear"),
           PAR_mn=raster::extract(covars.ls$PAR_mn, ., 
                                  fun=mean, small=T, method="bilinear"),
           PAR_sd=raster::extract(covars.ls$PAR_sd, ., 
                                  fun=mean, small=T, method="bilinear"),
           PAR_POWER=st_join(., 
                             covars.ls$irrad_growing.month, 
                             left=T, join=st_nearest_feature)$mnPAR,
           fetch=st_join(., 
                         covars.ls$fetch, 
                         left=T, join=st_nearest_feature)$fetchsum,
           fetchCat=case_when(fetch < fetch_thirds[1] ~ 1,
                              between(fetch, fetch_thirds[1], fetch_thirds[2]) ~ 2,
                              fetch > fetch_thirds[2] ~ 3))
  site.sf$PAR_surface <- site.sf[[ifelse(PAR_datasource=="MODIS", "PAR_mn", "PAR_POWER")]]
  return(site.sf)
}

extractCovarsToGrid <- function(grid.domain, covars.ls, PAR_datasource) {
  fetch_thirds <- quantile(covars.ls$fetch$fetchsum, probs=(1:2)/3)
  grid.domain <- grid.domain %>% 
    mutate(sst_day=raster::extract(covars.ls$sst_day, grid.domain, fun=mean, na.rm=T),
           sst_night=raster::extract(covars.ls$sst_night, grid.domain, fun=mean, na.rm=T),
           chla=raster::extract(covars.ls$chla, grid.domain, fun=mean, na.rm=T),
           KD_mn=raster::extract(covars.ls$KD_mn, grid.domain, fun=mean, na.rm=T),
           KD_sd=raster::extract(covars.ls$KD_sd, grid.domain, fun=mean, na.rm=T),
           PAR_mn=raster::extract(covars.ls$PAR_mn, grid.domain, fun=mean, na.rm=T),
           PAR_sd=raster::extract(covars.ls$PAR_sd, grid.domain, fun=mean, na.rm=T)) %>%
    st_join(., 
            covars.ls$irrad_growing.month %>% select(mnPAR) %>% rename(PAR_POWER=mnPAR), 
            left=T) %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    st_join(.,
            covars.ls$fetch %>% select(fetchsum) %>% rename(fetch=fetchsum), 
            left=T) %>% 
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    mutate(fetchCat=case_when(fetch < fetch_thirds[1] ~ 1,
                              between(fetch, fetch_thirds[1], fetch_thirds[2]) ~ 2,
                              fetch > fetch_thirds[2] ~ 3))
  
  grid.domain$PAR_surface <- grid.domain[[ifelse(PAR_datasource=="MODIS", "PAR_mn", "PAR_POWER")]]
  return(grid.domain)
}

  



getPrediction <- function(mod, lmType, ndraws, new.df) {
  if(lmType=="brms") {
    pred <- colMeans(posterior_epred(mod, newdata=new.df, ndraws=ndraws))
  } else {
    pred <- predict(mod, new.df)
  }
  return(pred)
}
