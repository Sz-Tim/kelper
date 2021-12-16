

library(tidyverse)
walk(dir("code", "^00._fn", full.names=T), source)

# Test out the function with various inputs
inputs <- list(
  
  # Domain/model run discretisation parameters
  domainArea=100, # Canopy area in m^2
  Ntimesteps=100, # Time to run model for
  timestep=20, # days
  NsizeClasses=3, # Number of size classes
  
  # Settlement of new plants
  settle=0.0001, # Settlement per unit free space (per m2 per day)
  
  # Baseline parameters apply when occupiedArea=0
  # Growth. g0-g1 give lowest possible growth rate (bounded at 0)
  N0=100, # Scalar or vector specifying initial number for all size classes
  maxSize=0.15, # metres holdfast diameter
  growthRateM=0.0002, # daily growth rate in m
  growthRateS=5, # daily size-dependent growth rate (exponential: 0=flat, >0=decline)
  g1=1, # Growth density effect (0,Inf) = (intense competition, low competition)
  allom1=100, # Scaling of canopy area to holdfast area
  interactType="Symmetric",
  
  # Mortality. Sum (mu0+mu1) = maximum possible mortality rate 
  mu0=0.0001,# Mortality baseline (daily - scaled internally by timestep)
  mu1=-0.5, # Mortality density effect (positive = competition)
  mortalityLinearDecline=TRUE,
  
  # Harvesting regime
  thinDuration=365,
  thinProp=0.1,
  thinLowLimit=0.5,
  thinHighLimit=1
)


# Work out size structured thinning 
classProps <- seq(0,1-1/(inputs$NsizeClasses),1/(inputs$NsizeClasses))
thinTmp <- array(0,dim=inputs$NsizeClasses)
# Set elements relating to each size class to "thinProp"
thinTmp[classProps >= thinLowLimit & classProps < thinHighLimit] <- thinProp
inputs$thin <- thinTmp

inputs$plotPars <- TRUE
inputs$plotOut <- TRUE

# ------------ Additional elements -----------------
# Light level (OR depth + turbidity coefficient)
# Storm frequency + response to storm (cf harvesting type?)

a <- runMatrixModelFiner(inputs)

