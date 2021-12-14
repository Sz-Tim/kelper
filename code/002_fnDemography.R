# KELPER
# Demography functions
# Tim Szewczyk


# Demographic functions used by the KELPER population models



####----------------------------------------------------------------------------
# Growth, breakage, & herbivory
####----------------------------------------------------------------------------


#' Calculate growth rate based on size 
#' \code{Tom's code}
#'
#' @param size Vector of holdfast diameters (lower bound for each class)
#' @param growthRateM Mean daily growth rate (m)
#' @param growthRateS Size effect (-: increase; 0: no effect; +: decline)
#'
#' @return
#' 
growSize <- function(size, growthRateM, growthRateS) {
    growthRateM * exp(-growthRateS * size)
  }


#' Calculate growth rate based on density
#' \code{Tom's code}
#'
#' @param propOccupiedSpace Proportion of canopy area occupied (interactType ==
#'   Symmetric: all plants; interactType == Asymmetric: only larger plants)
#' @param g0 Vector (length sizeClasses) of base growth rates
#' @param g1 Density effect (+); from Roughgarden?
#' @param sizeClass Index of size class
#'
#' @return
#' 
growDens <- function(propOccupiedSpace, g0, g1, sizeClass) {
  max(0, g0[sizeClass] * (1 - propOccupiedSpace / g1))
}








####----------------------------------------------------------------------------
# Survival & mortality
####----------------------------------------------------------------------------


#' Calculate mortality rate: linear with density
#' \code{Tom's code}
#'
#' @param propOccupiedSpace  Proportion of canopy area occupied (interactType ==
#'   Symmetric: all plants; interactType == Asymmetric: only larger plants)
#' @param mu0 Base mortality rate
#' @param mu1 Density effect
#'
#' @return
#' 
mortLinDens <- function(propOccupiedSpace, mu0, mu1) {
  max(0, mu0 + mu1 * propOccupiedSpace)
}


#' Calculate mortality rate: quadratic with density (high mortality at extremes)
#' \code{Tom's code}
#'
#' @param propOccupiedSpace  Proportion of canopy area occupied (interactType ==
#'   Symmetric: all plants; interactType == Asymmetric: only larger plants)
#' @param mu0 Base mortality rate
#' @param mu1 Density effect
#'
#' @return
#' 
mortOptDens <- function(propOccupiedSpace, mu0, mu1) {
  max(0, mu0 + mu1 * (1 - 4 * (propOccupiedSpace - 0.5) ^ 2))
}


#' Calculate mortality rate
#' \code{Tom's code}
#'
#' @param propOccupiedSpace  Proportion of canopy area occupied (interactType ==
#'   Symmetric: all plants; interactType == Asymmetric: only larger plants)
#' @param mu0 Base mortality rate
#' @param mu1 Density effect
#' @param linear Boolean: linear density effect or higher mortality at low &
#'   high densities?
#'
#' @return
#' 
mortDens <- function(propOccupiedSpace, mu0, mu1, linear) {
  if(linear==TRUE) {
    max(0, mu0 + mu1 * propOccupiedSpace)
  } else {
    max(0, mu0 + mu1*(1 - 4*(propOccupiedSpace - 0.5)^2))
  }
}

