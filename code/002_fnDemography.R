# KELPER
# Demography functions
# Tim Szewczyk


# Demographic functions used by the KELPER population models




####----------------------------------------------------------------------------
# Allometry
####----------------------------------------------------------------------------


#' Convert between metrics of kelp size 
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
convertAllometry <- function(orig, beta_mn, shape="linear", 
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




#' Calculate growth in frond area
#' 
#' Calculates Frond Area Index at the end of the time step
#'
#' @param FAI_orig Vector of initial frond areas
#' @param N_orig Vector of initial population abundances
#' @param A.mx Transition matrix
#' @param kappa FAI_orig/FAI_K
#' @param logAreaFrond.stage Vector of log frond areas calculated allometrically for each size class
#' @param pars parameter list
#'
#' @return
#' @export
#'
#' @examples
growFrondArea <- function(FAI_orig, N_orig, A.mx, kappa, logAreaFrond.stage, pars) {
  # Laminaria hyperborea exhibits a 'May cast', where the previous year's growth 
  # is shed essentially in its entirety. For other less deciduous species, 
  # FAI_new[j] would include FAI_orig[j] * A.mx[j,j] in the sum
  FAI_new <- rep(0, length(N_orig))
  # recruits: (surviving FAI) + (new growth within stage)
  FAI_new[1] <- N_orig[1]*A.mx[1,1]*pars$growthRateFrond[1]*(1-kappa)
  for(j in 2:length(N_orig)) {
    # others: (surviving FAI) + (new growth within stage) + (FAI of newcomers)
    FAI_new[j] <- N_orig[j]*A.mx[j,j]*pars$growthRateFrond[j] +
      N_orig[j-1]*A.mx[j,j-1]*exp(logAreaFrond.stage[j])
  }
  
  return(FAI_new)
}









calcBiomass <- function(N_t, lwtStipe=NULL, lmFit=NULL, lmType=NULL, 
                        ndraws=NULL, new.df=NULL, scale.df=NULL, y_var=NULL) {
  stipeMass <- N_t * exp(lwtStipe)
  frondMass <- N_t * exp(getPrediction(lmFit, lmType, ndraws, 
                                       new.df, scale.df, y_var))
  return(sum(stipeMass, frondMass, na.rm=T)/1e3)
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

