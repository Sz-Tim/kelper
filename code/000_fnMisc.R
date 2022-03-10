# KELPER
# Helper functions
# Tim Szewczyk


# Miscellaneous helper functions associated with the KELPER project

compileDatasets <- function(data.dir, supp.f=NULL) {
  
  library(tidyverse)
  
  if(!is.null(supp.f)) {
    # read in data from Dan Smale -- requires different formatting
    library(readxl)
    supp.sheets <- excel_sheets(supp.f)[-1]
    # Add this chunk somewhere else
    # NOTE: I should assume 4 year olds = subcanopy for the growth rates
    # ggplot(supp.raw$allometry, aes(factor(Age), weightFrond, fill=factor(depth))) + 
    #   geom_boxplot() + scale_fill_viridis_d()
    # BUT facet by region......?
    supp.sites <- read_xlsx(supp.f, "siteLocations", skip=1) %>%
      select(Location, lat, lon) %>%
      rename(location=Location)
    supp.raw <- map(supp.sheets, ~read_xlsx(supp.f, .x, skip=1)) %>% 
      setNames(supp.sheets)
    supp.ls <- list()
    supp.raw$allometry <- supp.raw$allometry %>%
      mutate(across(contains("length"), ~.x*10),
             across(contains("width"), ~.x*10)) %>%
      rename(depth=Depth, location=Location)
    supp.ls$lengthStipe_weightStipe <- supp.raw$allometry %>% 
      select(lengthStipe, weightStipe, depth, location)
    supp.ls$lengthStipe_weightFrond <- supp.raw$allometry %>% 
      select(lengthStipe, weightFrond, depth, location)
    supp.ls$weightStipe_weightFrond <- supp.raw$allometry %>% 
      select(weightStipe, weightFrond, depth, location)
    supp.ls$lengthStipe_weightTotal <- supp.raw$allometry %>% 
      mutate(weightTotal=weightHoldfast + weightStipe + weightFrond) %>%
      select(lengthStipe, weightTotal, depth, location)
    supp.ls$depth_N <- supp.raw$abundStage %>%
      rename(depth=Depth, location=Location, N_canopy=CF_N, N_subcanopy=SC_N) %>%
      mutate(N_recruits=DJ_N+UJ_N) %>%
      select(depth, starts_with("N_"), location) %>%
      bind_rows(supp.raw$densityCanopy %>% 
                  rename(depth=Depth, location=Location, N_canopy=CanopyDensity) %>%
                  select(location, N_canopy, depth)) %>%
      rename(NperSqM=N_canopy) %>%
      mutate(data_id=row_number())
    supp.ls$depth_prCover <- supp.raw$abundStage %>%
      rename(depth=Depth, location=Location, pr_canopy=CF_pct, pr_subcanopy=SC_pct) %>%
      mutate(pr_recruits=DJ_N+UJ_N) %>%
      select(depth, "pr_canopy", location, Rep) %>%
      bind_rows(supp.raw$densityCanopy %>% 
                  rename(depth=Depth, location=Location, 
                         pr_canopy=CanopyPctCover,
                         Rep=QuadratRep) %>%
                  select(location, Rep, pr_canopy, depth)) %>%
      mutate(across(starts_with("pr_"), ~.x/100),
             data_id=row_number())
    supp.ls <- map(supp.ls, ~left_join(.x, supp.sites, by="location")) %>%
      setNames(., paste0("Smale_Moore_", 1:length(.))) %>%
      imap(~.x %>% mutate(reference=.y))

  }
  
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
  
  
  
  fig_var <- imap_dfr(c(digitized, supp.ls), 
                      ~tibble(V1=names(.x)[1], V2=names(.x)[2], reference=.y)) %>%
    arrange(V1, V2) %>%
    mutate(vars=paste0(V1, "_", V2))
  
  unique_comparisons <- unique(fig_var$vars)
  varCombined.ls <- vector("list", length(unique_comparisons)) %>% 
    setNames(unique_comparisons)
  for(i in seq_along(varCombined.ls)) {
    refs <- fig_var$reference[fig_var$vars==unique_comparisons[i]]
    varCombined.ls[[i]] <- do.call('bind_rows', c(digitized, supp.ls)[refs])
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
                          growthRateStipeMax=cbind(c(194, 195, 58),
                                                   c(23, 21, 10)), # Kain 1976
                          # growthRateStipeMax=c(74, 74, 29), # Rinde & Sjotun 2005 southern pop
                          growthRateFrond=cbind(c(1787, 2299, 3979)/1e4,
                                                c(179, 397, 417)/1e4), # m^2/plant/growing season
                          frondAreaMax=5500/1e4,
                          growthRateDensityShape=0.5,
                          sizeClassLimits=(1000 * (0:6)/6)[c(1,2,5,7)],
                          survRate=c(0.4, 0.7, 0.9),
                          settlementRate=400,
                          lossRate=0.2345,
                          prFullHarvest=0,
                          freqHarvest=1e5,
                          harvestTarget=c("all", "stipe", "frond")[1],
                          extraPars=NULL,
                          stochParams=FALSE,
                          stormIntensity=NULL
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
      stochParams=stochParams,
      stormIntensity=stormIntensity,
      
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
      settlementRate=settlementRate,
      
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
    rast.proj <- dir(paste0(gis.dir, "climate"), "MODISA_L3m_SST.*11W", full.names=T)[1] %>%
      raster() %>% crs
    covars.ls <- list(
      fetch=dir(paste0(gis.dir, "fetch"), "log10_UK200m_depth40_0na.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      logSlope=dir(paste0(gis.dir, "bathymetry"), "UK_log10slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      slope=dir(paste0(gis.dir, "bathymetry"), "UK_slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      sstDayGrow_mn=dir(paste0(gis.dir, "climate"), 
                        "MODISA_L3m_SST.*11W.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      sstDayGrow_sd=dir(paste0(gis.dir, "climate"), 
                        "MODISA_L3m_SST.*11W_49N.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      KD_mn=dir(paste0(gis.dir, "attenuation"), 
                "MODISA_L3m_KD.*11W_49N.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      KD_sd=dir(paste0(gis.dir, "attenuation"), 
                "MODISA_L3m_KD.*11W_49N.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      PAR_mn=dir(paste0(gis.dir, "light"), 
                 "MODISA_L3m_PAR.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      PAR_sd=dir(paste0(gis.dir, "light"), 
                 "MODISA_L3m_PAR.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      irrad_monthly=dir(paste0(gis.dir, "light"), 
                        "POWER_Regional_monthly", full.names=T) %>%
        map_dfr(~read_csv(.x, skip=9, show_col_types=F)) %>% 
        filter(JAN != -999) %>% select(-ANN) %>%
        group_by(PARAMETER, YEAR, LAT, LON) %>%
        summarise(across(everything(), mean)) %>%
        group_by(PARAMETER, YEAR) %>% mutate(id=row_number()) %>% ungroup %>%
        pivot_longer(5:16, names_to="MONTH", values_to="PAR") %>%
        st_as_sf(coords=c("LON", "LAT"), crs=4326)
    )
    
    irrad.grid <- (st_bbox(covars.ls$irrad_monthly) + .25*c(-1,-1,1,1)) %>%
      st_make_grid(cellsize=c(0.5, 0.5)) %>% st_sf(id=1:length(.))
    
    covars.ls$irrad_growing.month <- covars.ls$irrad_monthly %>%
      filter(MONTH %in% c("JAN", "FEB", "MAR", "APR", "MAY", "JUN")) %>%
      group_by(YEAR, id) %>%
      summarise(PAR=mean(PAR)) %>%
      group_by(id) %>%
      summarise(mnPAR=mean(PAR), sdPAR=sd(PAR)) %>%
      ungroup() %>% select(-id) %>%
      st_join(irrad.grid %>% select(-id), .)
    covars.ls <- covars.ls[-which(names(covars.ls)=="irrad_monthly")]
    
    if(!is.null(saveFile)) {
      saveRDS(covars.ls, saveFile)
    }
  } else {
    covars.ls <- readRDS(loadFile)
  }
  return(covars.ls)
}



loadCovariates_full <- function(gis.dir=NULL, bbox=NULL, loadFile=NULL, saveFile=NULL) {
  library(stars)
  if(is.null(loadFile)) {
    rast.proj <- dir(paste0(gis.dir, "climate"), "MODISA_L3m_SST.*11W", full.names=T)[1] %>%
      raster() %>% crs
    covars.ls <- list(
      fetch=dir(paste0(gis.dir, "fetch"), "log10_UK200m_depth40_0na.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      logSlope=dir(paste0(gis.dir, "bathymetry"), "UK_log10slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      slope=dir(paste0(gis.dir, "bathymetry"), "UK_slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      sstDayGrow=dir(paste0(gis.dir, "climate"), 
                     "MODISA_L3m_SST.*11W_49N.*tif$", full.names=T) %>%
        setNames(., str_sub(., 84, 87)) %>%
        map(raster) %>% stack() %>% crop(bbox) %>% 
        st_as_stars() %>% st_as_sf() %>% 
        pivot_longer(starts_with("X"), names_to="YEAR", values_to="SST") %>%
        mutate(YEAR=str_sub(YEAR, 2, -1)),
      KD_mn=dir(paste0(gis.dir, "attenuation"), 
                "MODISA_L3m_KD.*11W_49N.*tif$", full.names=T) %>%
        setNames(., str_sub(., 84, 87)) %>%
        map(raster) %>% stack() %>% crop(bbox) %>% 
        st_as_stars() %>% st_as_sf() %>% 
        pivot_longer(starts_with("X"), names_to="YEAR", values_to="KD") %>%
        mutate(YEAR=str_sub(YEAR, 2, -1)),
      PAR_mn=dir(paste0(gis.dir, "light"), 
                 "MODISA_L3m_PAR.*tif$", full.names=T)  %>%
        setNames(., str_sub(., 76, 79)) %>%
        map(raster) %>% stack() %>% crop(bbox) %>% 
        st_as_stars() %>% st_as_sf() %>% 
        pivot_longer(starts_with("X"), names_to="YEAR", values_to="PAR") %>%
        mutate(YEAR=str_sub(YEAR, 2, -1)),
      irrad_monthly=dir(paste0(gis.dir, "light"), 
                        "POWER_Regional_monthly", full.names=T) %>%
        map_dfr(~read_csv(.x, skip=9, show_col_types=F)) %>% 
        filter(JAN != -999) %>% select(-ANN) %>%
        group_by(PARAMETER, YEAR, LAT, LON) %>%
        summarise(across(everything(), mean)) %>%
        group_by(PARAMETER, YEAR) %>% mutate(id=row_number()) %>% ungroup %>%
        pivot_longer(5:16, names_to="MONTH", values_to="PAR") %>%
        st_as_sf(coords=c("LON", "LAT"), crs=4326)
    )
    
    irrad.grid <- (st_bbox(covars.ls$irrad_monthly) + .25*c(-1,-1,1,1)) %>%
      st_make_grid(cellsize=c(0.5, 0.5)) %>% st_sf(id=1:length(.))
    
    covars.ls$irrad_growing.month <- covars.ls$irrad_monthly %>%
      filter(MONTH %in% c("JAN", "FEB", "MAR", "APR", "MAY", "JUN")) %>%
      group_by(YEAR, id) %>%
      summarise(PAR=mean(PAR)) %>%
      ungroup() %>% select(-id) %>%
      st_join(irrad.grid %>% select(-id), .) %>%
      filter(!is.na(YEAR))
    covars.ls <- covars.ls[-which(names(covars.ls)=="irrad_monthly")]
    
    if(!is.null(saveFile)) {
      saveRDS(covars.ls, saveFile)
    }
  } else {
    covars.ls <- readRDS(loadFile)
  }
  return(covars.ls)
}







extractCovarsToDatasets <- function(data.ls, covars.ls=NULL, PAR_datasource, grid.sf=NULL) {
  if(is.null(grid.sf)) {
    # fetch_thirds <- quantile(covars.ls$fetch@data@values, probs=(1:2)/3, na.rm=T)
    fetch_thirds <- c(4.02, 4.33) # bad practice, but based on Pedersen sites
    for(i in seq_along(data.ls)) {
      if(any(!is.na(data.ls[[i]]$lat))) {
        i_sf <- data.ls[[i]] %>% 
          filter(!is.na(lat) & !is.na(lon)) %>%
          st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
          mutate(sstDay_mn=raster::extract(covars.ls$sstDayGrow_mn, ., buffer=8e3,
                                           fun=mean, small=T),
                 sstDay_sd=raster::extract(covars.ls$sstDayGrow_sd, ., buffer=8e3, 
                                           fun=mean, small=T),
                 KD_mn=raster::extract(covars.ls$KD_mn, ., buffer=8e3, 
                                       fun=mean, small=T),
                 KD_sd=raster::extract(covars.ls$KD_sd, ., buffer=8e3, 
                                       fun=mean, small=T),
                 PAR_mn=raster::extract(covars.ls$PAR_mn, ., buffer=8e3, 
                                        fun=mean, small=T),
                 PAR_sd=raster::extract(covars.ls$PAR_sd, ., buffer=8e3, 
                                        fun=mean, small=T),
                 PAR_POWER=st_join(., 
                                   covars.ls$irrad_growing.month, 
                                   left=T, join=st_nearest_feature)$mnPAR,
                 logSlope=raster::extract(covars.ls$logSlope, ., buffer=8e3, 
                                          fun=mean, small=T),
                 slope=raster::extract(covars.ls$slope, ., buffer=8e3, 
                                       fun=mean, small=T),
                 fetch=raster::extract(covars.ls$fetch, ., buffer=8e3, 
                                       fun=mean, small=T),
                 fetchCat=case_when(fetch < fetch_thirds[1] ~ 1,
                                    between(fetch, fetch_thirds[1], fetch_thirds[2]) ~ 2,
                                    fetch > fetch_thirds[2] ~ 3))
        data.ls[[i]] <- left_join(data.ls[[i]], st_drop_geometry(i_sf))
        if(any(!is.na(data.ls[[i]]$depth) & 
               !is.na(data.ls[[i]]$KD_mn) & 
               !is.na(data.ls[[i]]$PAR_mn))) {
          data.ls[[i]] <- data.ls[[i]] %>%
            mutate(PAR_atDepth_POWER=PAR_POWER * exp(-KD_mn * depth),
                   PAR_atDepth_MODIS=PAR_mn * exp(-KD_mn * depth))
          data.ls[[i]]$PAR_atDepth <- data.ls[[i]][[paste0("PAR_atDepth_", PAR_datasource)]]
        }
      }
    }
  } else {
    for(i in seq_along(data.ls)) {
      if(any(!is.na(data.ls[[i]]$lat))) {
        i_sf <- data.ls[[i]] %>% 
          filter(!is.na(lat) & !is.na(lon)) %>%
          st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
          st_join(., grid.sf, left=T)
        data.ls[[i]] <- left_join(data.ls[[i]], st_drop_geometry(i_sf))
        if(any(!is.na(data.ls[[i]]$depth) & 
               !is.na(data.ls[[i]]$KD) & 
               !is.na(data.ls[[i]]$PAR))) {
          data.ls[[i]] <- data.ls[[i]] %>%
            mutate(PAR_atDepth=PAR * exp(-KD * depth))
        }
      }
    }
  }
  
  return(data.ls)
}




extractCovarsToPts <- function(site.i, covars.ls, PAR_datasource) {
  # fetch_thirds <- quantile(covars.ls$fetch@data@values, probs=(1:2)/3, na.rm=T)
  fetch_thirds <- c(4.02, 4.33) # bad practice, but based on Pedersen sites
  site.sf <- site.i %>% 
    filter(!is.na(lat) & !is.na(lon)) %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    mutate(sstDay_mn=raster::extract(covars.ls$sstDayGrow_mn, ., 
                                     fun=mean, small=T, method="bilinear"),
           sstDay_sd=raster::extract(covars.ls$sstDayGrow_sd, ., 
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
           logSlope=raster::extract(covars.ls$logSlope, ., 
                                    fun=mean, small=T, method="bilinear"),
           slope=raster::extract(covars.ls$slope, ., 
                                 fun=mean, small=T, method="bilinear"),
           fetch=raster::extract(covars.ls$fetch, ., 
                                 fun=mean, small=T, method="bilinear"),
           fetchCat=case_when(fetch < fetch_thirds[1] ~ 1,
                              between(fetch, fetch_thirds[1], fetch_thirds[2]) ~ 2,
                              fetch > fetch_thirds[2] ~ 3))
  site.sf$PAR_surface <- site.sf[[ifelse(PAR_datasource=="MODIS", "PAR_mn", "PAR_POWER")]]
  return(site.sf)
}





extractCovarsToGrid <- function(grid.domain, covars.ls, PAR_datasource) {
  # fetch_thirds <- quantile(covars.ls$fetch@data@values, probs=(1:2)/3, na.rm=T)
  fetch_thirds <- c(4.02, 4.33) # bad practice, but based on Pedersen sites
  grid.domain <- grid.domain %>% 
    mutate(sstDay_mn=raster::extract(covars.ls$sstDayGrow_mn, grid.domain, fun=mean, na.rm=T),
           sstDay_sd=raster::extract(covars.ls$sstDayGrow_sd, grid.domain, fun=mean, na.rm=T),
           KD_mn=raster::extract(covars.ls$KD_mn, grid.domain, fun=mean, na.rm=T),
           KD_sd=raster::extract(covars.ls$KD_sd, grid.domain, fun=mean, na.rm=T),
           PAR_mn=raster::extract(covars.ls$PAR_mn, grid.domain, fun=mean, na.rm=T),
           PAR_sd=raster::extract(covars.ls$PAR_sd, grid.domain, fun=mean, na.rm=T),
           logSlope=raster::extract(covars.ls$logSlope, grid.domain, fun=mean, na.rm=T),
           slope=raster::extract(covars.ls$slope, grid.domain, fun=mean, na.rm=T),
           fetch=raster::extract(covars.ls$fetch, grid.domain, fun=mean, na.rm=T)
    ) %>%
    st_join(., 
            covars.ls$irrad_growing.month %>% select(mnPAR) %>% rename(PAR_POWER=mnPAR), 
            left=T) %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    mutate(fetchCat=case_when(fetch < fetch_thirds[1] ~ 1,
                              between(fetch, fetch_thirds[1], fetch_thirds[2]) ~ 2,
                              fetch > fetch_thirds[2] ~ 3))
  
  grid.domain$PAR_surface <- grid.domain[[ifelse(PAR_datasource=="MODIS", "PAR_mn", "PAR_POWER")]]
  return(grid.domain)
}

  



simulateLandscape <- function(grid.sf, var.sf, nYr, colName) {
  # grid.sf: base grid
  # var.sf: variable
  # colName: column name with variable in var.sf
  # nYr: number of years to simulate
  library(mvtnorm)
  # pair var with grid, then extract to matrix
  var.mx <- grid.sf %>% select(id) %>%
    st_join(., 
            var.sf %>%
              pivot_wider(names_from="YEAR", values_from=colName, names_prefix="y_") %>%
              st_as_sf(),
            left=T) %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    st_drop_geometry() %>% select(-id) %>% as.matrix()
  # generate simulated values
  na.ids <- which(is.na(var.mx), arr.ind=T)
  if(nrow(na.ids)>0) {
    for(i in 1:nrow(na.ids)) {
      var.mx[na.ids[i,1], na.ids[i,2]] <- mean(var.mx[na.ids[i,1],], na.rm=T) +
        rnorm(1, 0, 0.01)
    }
  }
  if(colName=="KD") var.mx <- log(var.mx)
  var.sim <- t(rmvnorm(nYr, rowMeans(var.mx), cov(t(var.mx))))
  colnames(var.sim) <- paste0(colName, "_", 1:ncol(var.sim))
  if(colName=="KD") var.sim <- exp(var.sim)
  var.sim[var.sim < 0] <- min(var.sf[[colName]], na.rm=T)
  # bind to grid
  grid.updated <- grid.sf %>% bind_cols(as_tibble(var.sim)) %>% 
    pivot_longer(starts_with(paste0(colName, "_")), 
                 names_to="year", values_to=colName) %>%
    mutate(year=as.numeric(str_sub(year, nchar(colName)+2, -1)))
  return(grid.updated)
}






getPrediction <- function(mod, ndraws, new.df, scale.df, y_var) {
  # center and scale predictors to match regression inputs
  for(i in 1:ncol(new.df)) {
    if(is.numeric(new.df[[i]])) {
      scale.row <- which(scale.df$Par==names(new.df)[i])
      new.df[[i]] <- (new.df[[i]] - scale.df$ctr[scale.row])/scale.df$scl[scale.row]
    }
  }
  
  pred <- colMeans(posterior_epred(mod, newdata=new.df, ndraws=ndraws, re.form=NA))
  
  # de-center and de-scale predictions to natural scale
  if(y_var != "N") {
    scale.row <- which(scale.df$Par==y_var)
    pred <- pred * scale.df$scl[scale.row] + scale.df$ctr[scale.row] 
  }
  return(pred)
}









simDepthsWithinCell <- function(x, grid.i, grid.sim=NULL, gridRes, pars.sim, 
                                surv.df, fecund.df, lm.fit, lm.mnsd) {
  library(glue); library(lubridate); library(sf); library(brms)
  library(Matrix); library(tidyverse)
  if(pars.sim$landscape=="dynamic") {
    cell.env <- grid.sim %>% st_drop_geometry() %>% filter(id==grid.i$id[x])
  } else {
    cell.env <- grid.i[x,]
  }
  pop.ls <- mass.ls <- vector("list", length(pars.sim$depths))
  for(j in 1:length(pars.sim$depths)) {
    par_i <- setParameters(
      tmax=pars.sim$tmax, 
      survRate=cbind(filter(surv.df, exposure==cell.env$fetchCat[1])$survRate,
                     filter(surv.df, exposure==cell.env$fetchCat[1])$rate_sd),
      settlementRate=c(filter(fecund.df, exposure==cell.env$fetchCat[1])$rate_mn,
                       filter(fecund.df, exposure==cell.env$fetchCat[1])$rate_sd),
      lossRate=c(0.2340987, 90.1752),
      stochParams=pars.sim$stochParams,
      stormIntensity=pars.sim$storms,
      extraPars=list(
        depth=pars.sim$depths[j],
        env=cell.env,
        lenSt_to_wtSt=lm.fit$lenSt_to_wtSt.lm,
        lenSt_to_wtFr=lm.fit$lenSt_to_wtFr.lm,
        wtFr_to_arFr=lm.fit$wtFr_to_arFr.lm,
        arFr_to_wtFr=lm.fit$arFr_to_wtFr.lm,
        N_canopy=lm.fit$N_canopy.lm,
        FAI=lm.fit$FAI.lm,
        canopyHeight=lm.fit$canopyHeight.lm,
        sc.df=lm.mnsd)
    )
    out <- map(1:pars.sim$nSim, 
               ~simulatePopulation(par_i, ndraws=ifelse(par_i$stochParams, 50, 4e3)))
    pop.ls[[j]] <- imap_dfr(out, 
                            ~tibble(sim=.y, 
                                    year=rep(1:par_i$tmax, 3),
                                    month=rep(c(1,6,7), each=par_i$tmax),
                                    date=ymd(glue("{year}-{month}-01")),
                                    PAR_atDepth=rep(.x$PAR, 3),
                                    K_N=rep(.x$K_N, 3),
                                    K_FAI=rep(.x$K_FAI, 3),
                                    FAI=c(.x$FAI[3,,]),
                                    N.recruits=c(.x$N[1,,]),
                                    N.subcanopy=c(.x$N[2,,]),
                                    N.canopy=c(.x$N[3,,]),
                                    kappa_FAI=c(.x$kappa[,,1]),
                                    kappa_N=c(.x$kappa[,,2])) %>%
                              pivot_longer(contains("N."), names_to="stage", values_to="N") %>%
                              mutate(stage=factor(str_sub(stage, 3, -1), 
                                                  levels=c("canopy", "subcanopy", "recruits")))) %>%
      mutate(id=x) %>%
      left_join(., par_i$env) %>%
      mutate(depth=par_i$depth,
             landscape=par_ilandscape,
             stochParams=par_istochParams)
    mass.ls[[j]] <- imap_dfr(out,
                             ~tibble(sim=.y,
                                     year=rep(1:par_i$tmax, 3),
                                     month=rep(c(1,6,7), each=par_i$tmax),
                                     date=ymd(glue("{year}-{month}-01")),
                                     PAR_atDepth=rep(.x$PAR, 3),
                                     FAI=c(.x$FAI[3,,]),
                                     biomass=c(.x$biomass),
                                     kappa_FAI=c(.x$kappa[,,1]),
                                     kappa_N=c(.x$kappa[,,2]))) %>%
      mutate(id=x) %>%
      left_join(., par_i$env) %>%
      mutate(depth=par_i$depth,
             landscape=par_ilandscape,
             stochParams=par_istochParams)
  }
  sim.info <- glue("{str_pad(x,4,'left','0')}_{gridRes}_{landscape}",
                   "_{ifelse(stochParams, 'stochPar', 'optPar')}")
  saveRDS(do.call(rbind, pop.ls), glue("out\\pop_{sim.info}.rds"))
  saveRDS(do.call(rbind, mass.ls), glue("out\\mass_{sim.info}.rds"))
  return(x)
}









simSensitivityDepthsWithinCell <- function(x, grid.i, gridRes, pars.sens, 
                                           lm.fit, lm.mnsd, parSets) {
  library(glue); library(lubridate); library(sf); library(brms)
  library(Matrix); library(tidyverse)
  
  # setup landscape and parameters
  cell.env <- grid.i[x,]
  grStipe.draws <- parSets[[4]] %>% filter(param=="growStipe") %>% group_split(parDraw)
  grFrond.draws <- parSets[[4]] %>% filter(param=="growFrond") %>% group_split(parDraw)
  loss.draws <- parSets[[4]] %>% filter(param=="loss")
  surv.draws <- parSets[[cell.env$fetchCat]] %>% filter(param=="surv") %>% group_split(parDraw)
  settle.draws <- parSets[[cell.env$fetchCat]] %>% filter(param=="settlement")
  
  pop.depth <- mass.depth <- vector("list", length(pars.sens$depths))
  for(j in 1:length(pars.sens$depths)) {
    pop.ls <- mass.ls <- vector("list", pars.sens$nParDraws)
    for(k in 1:pars.sens$nParDraws) {
      par_i <- setParameters(
        tmax=pars.sens$tmax, 
        growthRateStipeMax=cbind(grStipe.draws[[k]]$val),
        growthRateFrond=cbind(grFrond.draws[[k]]$val),
        lossRate=loss.draws$val[k],
        survRate=cbind(surv.draws[[k]]$val),
        settlementRate=settle.draws$val[k],
        stochParams=pars.sens$stochParams,
        extraPars=list(
          depth=pars.sens$depths[j],
          env=cell.env,
          lenSt_to_wtSt=lm.fit$lenSt_to_wtSt.lm,
          lenSt_to_wtFr=lm.fit$lenSt_to_wtFr.lm,
          wtFr_to_arFr=lm.fit$wtFr_to_arFr.lm,
          arFr_to_wtFr=lm.fit$arFr_to_wtFr.lm,
          N_canopy=lm.fit$N_canopy.lm,
          FAI=lm.fit$FAI.lm,
          canopyHeight=lm.fit$canopyHeight.lm,
          sc.df=lm.mnsd)
      )
      out <- map(1:pars.sens$nSim, ~simulatePopulation(par_i, ndraws=4e3))
      pop.ls[[k]] <- imap_dfr(out, 
                              ~tibble(sim=.y, 
                                      year=rep(1:par_i$tmax, 3),
                                      month=rep(c(1,6,7), each=par_i$tmax),
                                      date=ymd(glue("{year}-{month}-01")),
                                      N.recruits=c(.x$N[1,,]),
                                      N.subcanopy=c(.x$N[2,,]),
                                      N.canopy=c(.x$N[3,,])) %>%
                                pivot_longer(contains("N."), names_to="stage", values_to="N") %>%
                                mutate(stage=factor(str_sub(stage, 3, -1), 
                                                    levels=c("canopy", "subcanopy", "recruits")))) %>%
        filter(month != 6 & year > pars.sens$tskip) %>%
        mutate(N=replace(N, N<=0, 1e-4)) %>%
        group_by(sim, month, stage) %>%
        summarise(N_mn=mean(N), N_md=median(N), N_sd=sd(N),
                  lN_mn=mean(log(N+1)), lN_md=median(log(N+1)), lN_sd=sd(log(N+1))) %>%
        mutate(parDraw=k)
      mass.ls[[k]] <- imap_dfr(out,
                               ~tibble(sim=.y,
                                       year=rep(1:par_i$tmax, 3),
                                       month=rep(c(1,6,7), each=par_i$tmax),
                                       date=ymd(glue("{year}-{month}-01")),
                                       PAR_atDepth=rep(.x$PAR, 3),
                                       K_N=rep(.x$K_N, 3),
                                       K_FAI=rep(.x$K_FAI, 3),
                                       FAI=c(.x$FAI[3,,]),
                                       biomass=c(.x$biomass),
                                       kappa_FAI=c(.x$kappa[,,1]),
                                       kappa_N=c(.x$kappa[,,2]))) %>%
        filter(month != 6 & year > pars.sens$tskip) %>%
        mutate(kappa=pmax(kappa_FAI, kappa_N)) %>%
        group_by(sim, month) %>%
        summarise(PAR_atDepth=first(PAR_atDepth), K_N=first(K_N), K_FAI=first(K_FAI),
                  FAI_mn=mean(FAI), FAI_md=median(FAI), FAI_sd=sd(FAI),
                  biomass_mn=mean(biomass), biomass_md=median(biomass), biomass_sd=sd(biomass),
                  p_kFAI=mean(kappa_FAI==kappa), p_kN=mean(kappa_N==kappa),
                  kappa_mn=mean(kappa), kappa_md=median(kappa), kappa_sd=sd(kappa),
                  kFAI_mn=mean(kappa_FAI), kFAI_md=median(kappa_FAI), kFAI_sd=sd(kappa_FAI),
                  kN_mn=mean(kappa_N), kN_md=median(kappa_N), kN_sd=sd(kappa_N)) %>%
        mutate(parDraw=k) 
    }
    pop.depth[[j]] <- do.call(rbind, pop.ls) %>% mutate(depth=pars.sens$depths[j])
    mass.depth[[j]] <- do.call(rbind, mass.ls) %>% mutate(depth=pars.sens$depths[j])
  }
  sim.info <- glue("{str_pad(x,4,'left','0')}_{gridRes}")
  saveRDS(do.call(rbind, pop.depth) %>% mutate(id=x), 
          glue("out\\sensitivity\\pop_{sim.info}.rds"))
  saveRDS(do.call(rbind, mass.depth) %>% mutate(id=x), 
          glue("out\\sensitivity\\mass_{sim.info}.rds"))
  return(x)
}











#' Emulate global sensitivity analysis output
#'
#' Emulate the output from a global sensitivity analysis using boosted
#' regression trees with different interaction depths. Based on function
#' described in Prowse et al 2016. Writes files to brt.dir
#' @param sens.out Dataframe of the parameter sets and simulation summaries;
#'   \code{.$results} from \link{global_sensitivity}
#' @param params Vector of parameter names to include
#' @param n.cores \code{1} Number of cores for fitting subsample BRTs
#' @param n.sub \code{10} Number of subsamples for each emulation
#' @param td \code{c(1,3,5)} Vector of regression tree interaction depths to
#'   test
#' @param resp Which response summary to use (column name from \code{sens.out})
#' @param brt.dir \code{"out/brt/"} Directory to store boosted regression tree
#'   output
#' @return Success message
#' @keywords parameters, sensitivity, save, output
#' @export

emulate_sensitivity <- function(sens.out, params, resp, brt.dir, siminfo, 
                                n.sub=10, td=c(1,3,5)) {
  library(tidyverse)
  if(!dir.exists(brt.dir)) dir.create(brt.dir)
  sub.prop <- seq(0.75, 1, length.out=n.sub)
  
  for(i in 1:n.sub) {
    # subset sensitivity results
    sub.samp <- sample_frac(sens.out, sub.prop[i]) %>% as.data.frame
    n <- nrow(sub.samp)
    
    # fit BRTs of different tree complexities for given response variable
    for(j in 1:length(td)) {
      td_j <- td[j]
      brt.fit <- dismo::gbm.step(sub.samp, gbm.x=params, gbm.y=resp, 
                                 max.trees=200000, n.folds=5, 
                                 family="gaussian", tree.complexity=td_j,
                                 bag.fraction=0.8, silent=T, plot.main=F)
      saveRDS(brt.fit, glue::glue("{brt.dir}\\{resp}_{siminfo}_td-{td_j}-{n}.rds"))
    }
  }
  return(paste("Finished emulations of", resp))
}









#' Summarize BRT emulations
#'
#' Summarize the output from global sensitivity analysis boosted regression
#' trees emulations. Based on function described in Prowse et al 2016. Reads BRT
#' output saved via emulate_sensitivity
#' @param resp Which response summary to use (column name from \code{sens.out})
#' @return List of two dataframes: relative influences, and cross-validation
#'   deviances, each across different subsample sizes and tree complexities.
#' @param brt.dir \code{"out/brt/"} Directory with stored boosted regression
#'   tree output
#' @keywords parameters, sensitivity, save, output
#' @export

emulation_summary <- function(resp, brt.dir, siminfo) {
  library(gbm); library(tidyverse)
  f <- dir(brt.dir, paste0(resp, "_"))
  f.i <- str_split_fixed(f, "-", 3)
  id <- as.numeric(str_split_fixed(siminfo, "_", 3)[,1])
  cvDev.df <- tibble(response=resp,
                     td=f.i[,2],
                     smp=as.numeric(str_remove(f.i[,3], ".rds")),
                     id=id)
  cvDev.df$Dev <- NA
  ri.ls <- vector("list", length(f))
  for(i in seq_along(f)) {
    brt <- readRDS(paste0(brt.dir, f[i]))
    # cross validation deviance
    cvDev.df$Dev[i] <- brt$cv.statistics$deviance.mean
    # relative influence
    ri.ls[[i]] <- as_tibble(brt$contributions) %>%
      mutate(td=cvDev.df$td[i], 
             smp=cvDev.df$smp[i], 
             response=cvDev.df$response[i])
  }
  ri.df <- bind_rows(ri.ls) %>% 
    group_by(response, td, smp, var) %>%
    summarise(rel.inf=sum(rel.inf)) %>% 
    ungroup %>% group_by(response, td, smp) %>%
    mutate(rel.inf=rel.inf/sum(rel.inf)) %>%
    ungroup %>% arrange(desc(td), desc(smp), desc(rel.inf)) %>%
    mutate(id=id)
  
  write_csv(cvDev.df, glue::glue("{brt.dir}\\summaries\\{resp}_cv_{siminfo}.csv"))
  write_csv(ri.df, glue::glue("{brt.dir}\\summaries\\{resp}_ri_{siminfo}.csv"))
  return()
}








runBRTs <- function(x, pop.f, mass.f, parSets, grid.i, meta.cols, params, brt.dir) {
  library(glue); library(lubridate); library(sf); library(brms)
  library(lme4); library(glmmTMB); library(tidyverse)
  siminfo <- str_remove(str_split_fixed(pop.f[x], "pop_", 2)[,2], ".rds")
  pop.i <- readRDS(pop.f[x]) %>% 
    left_join(., parSets[[4]]) %>%
    left_join(., parSets[[filter(grid.i, id==.$id[1])$fetchCat]]) %>%
    filter(stage=="canopy")
  resp <- names(pop.i)
  resp <- resp[!resp %in% c(params, meta.cols)]
  resp <- grep("_mn|_sd", resp, value=T)
  for(j in 1:length(resp)) {
    emulate_sensitivity(pop.i %>% filter(month==7), params, n.sub=2, td=5,
                        resp[j], brt.dir, glue("{siminfo}_july"))
    emulation_summary(resp[j], brt.dir, glue("{siminfo}_july"))
  }
  
  mass.i <- readRDS(mass.f[x]) %>%
    left_join(., parSets[[4]]) %>%
    left_join(., parSets[[filter(grid.i, id==.$id[1])$fetchCat]])
  resp <- names(mass.i)
  resp <- resp[!resp %in% c(params, meta.cols)]
  resp <- grep("_mn|_sd|p_k", resp, value=T)
  for(j in 1:length(resp)) {
    emulate_sensitivity(mass.i %>% filter(month==7), params,
                        resp[j], brt.dir, glue("{siminfo}_july"))
    emulation_summary(resp[j], brt.dir, glue("{siminfo}_july"))
  }
}


