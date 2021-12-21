# KELPER
# Demography functions
# Tim Szewczyk


# Demographic functions used by the KELPER population models




####----------------------------------------------------------------------------
# Allometry
####----------------------------------------------------------------------------


#' Convert between 1D metrics of kelp size 
#'
#' @param orig Vector of original measurements
#' @param beta_mn Mean intercept + slope(s) for species
#' @param shape Shape of allometric relationship ("linear", "quadratic", "exponential", "log")
#' @param resid_err Residual standard deviation among individuals (default:
#'   NULL). If included, results are stochastic with Norm(0, resid_err) added to
#'   each individual
#' @param beta_sd Standard deviation in beta estimates (default: NULL). If
#'   included, results are stochastic with beta ~Norm(beta_mn, beta_sd)
#'
#' @return
#' 
convertAllometry1D <- function(orig, beta_mn, shape="linear", 
                               resid_err=NULL, beta_sd=NULL) {
  
  # Incorporate parameter uncertainty
  if(!is.null(beta_sd)) {
    beta <- rnorm(length(beta_mn), beta_mn, beta_sd)
  } else {
    beta <- beta_mn
  }
  
  # Incorporate environmental stochasticity
  if(!is.null(resid_err)) {
    err <- rnorm(length(orig), 0, resid_err)
  } else {
    err <- 0
  }
  
  switch(shape,
         linear=beta[1] + beta[2]*orig + err,
         quadratic=beta[1] + beta[2]*orig + beta[3]*orig^2 + err,
         exponential=exp(beta[1] + beta[2]*orig + err),
         log=log(beta[1] + beta[2]*orig + err))
}





####----------------------------------------------------------------------------
# Growth, breakage, & herbivory
####----------------------------------------------------------------------------


#' Calculate growth rate based on size 
#' \code{Tom's code}
#' NOTE: This is not GROWTH, per se, but TRANSITION PROBABILITY
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
#' NOTE: This is not GROWTH, per se, but TRANSITION PROBABILITY
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
  pmax(0, g0[sizeClass] * (1 - propOccupiedSpace / g1))
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
  pmax(0, mu0 + mu1 * propOccupiedSpace)
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
  pmax(0, mu0 + mu1 * (1 - 4 * (propOccupiedSpace - 0.5) ^ 2))
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
    pmax(0, mu0 + mu1 * propOccupiedSpace)
  } else {
    pmax(0, mu0 + mu1*(1 - 4*(propOccupiedSpace - 0.5)^2))
  }
}

