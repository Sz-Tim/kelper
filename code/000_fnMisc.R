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
                          N_stages=3,
                          N_seasons=2,
                          tmax=20,
                          growthRateStipeMax=c(194, 195, 58),
                          growthRateFrond=c(1787, 2299, 3979)/1e4, # m^2/plant/spring
                          frondAreaMax=5500/1e4,
                          growthRateStipeDensityShape=1,
                          sizeClassLimits=(2000 * (0:6)/6)[c(1,2,5,7)],
                          survRate=c(0.4, 0.7, 0.9),
                          settlementRateBg=100,
                          lossRate=c(0.4, 0.4, 0.4),
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
      N_stages=N_stages,
      N_seasons=N_seasons,
      tmax=tmax,
      
      # growth
      growthRateStipeMax=growthRateStipeMax,
      growthRateFrond=growthRateFrond,
      frondAreaMax=frondAreaMax,
      growthRateDensityStipeShape=growthRateStipeDensityShape,
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

