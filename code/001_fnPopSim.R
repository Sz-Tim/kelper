# KELPER
# Simulation functions
# Tim Szewczyk


# Population simulation functions

simulatePopulation <- function(pars, N0=NULL, ndraws=4e3) {
  
  library(tidyverse); library(brms)
  #---- setup landscape
  env.df <- pars$env %>% 
    mutate(PAR_atDepth=PAR * exp(-KD * pars$depth),
           location=NA) %>%
    select(PAR_atDepth, SST, fetch, location)
  if(nrow(pars$env) == 1) {
    env.df <- env.df %>% uncount(pars$tmax)
  }
  #---- setup parameters
  if(pars$stochParams) {
    par.yr <- list(loss=rbeta(pars$tmax, prod(pars$lossRate), (1-pars$lossRate[1])*pars$lossRate[2]),
                   settlement=pmax(0, rnorm(pars$tmax, pars$settlementRate[1], pars$settlementRate[2])),
                   surv=apply(pars$survRate, 1, function(x) pmax(0, pmin(1, rnorm(pars$tmax, x[1], x[2])))),
                   growStipeMax=apply(pars$growthRateStipeMax, 1, function(x) rnorm(pars$tmax, x[1], x[2])),
                   growFrond=apply(pars$growthRateFrond, 1, function(x) rnorm(pars$tmax, x[1], x[2])))
  } else {
    par.yr <- list(loss=rep(pars$lossRate[1], pars$tmax),
                   settlement=rep(pars$settlementRate[1], pars$tmax),
                   surv=apply(pars$survRate, 1, function(x) rep(x[1], pars$tmax)),
                   growStipeMax=apply(pars$growthRateStipeMax, 1, function(x) rep(x[1], pars$tmax)),
                   growFrond=apply(pars$growthRateFrond, 1, function(x) rep(x[1], pars$tmax)))
  }
  #---- setup storm effects
  if(is.null(pars$stormIntensity)) {
    par.yr$surv_storm <- par.yr$surv
  } else {
    # vector of storm intensities
    # affects winter survival, loss
    # should depend on depth...
    par.yr$loss <- qbeta(pnorm(pars.sim$storms, 0, 1), 
                         prod(pars$lossRate), 
                         (1-pars$lossRate[1])*pars$lossRate[2])
    par.yr$surv_strm <- apply(pars$survRate, 1, 
                              function(x) pmax(0, 
                                               pmin(1, 
                                                    qnorm(pnorm(-pars.sim$storms, 
                                                                mean(-pars.sim$storms), 1), 
                                                          x[1], x[2]))))
  }
  par.yr$surv <- sqrt(par.yr$surv) # annual rates to 1/2 year rates
  par.yr$surv_strm <- sqrt(par.yr$surv_strm)
  
  #---- global parameters
  ## maxStipeLen = maximum expected canopy height
  ## sizeClassLimits = boundaries between stages
  ## sizeClassMdPts = midpoint per stage
  ## K_N = carrying capacity on abundance / m2
  ## K_FAI = carrying capacity on frond area / m2
  maxStipeLen <- getPrediction(pars$canopyHeight, ndraws, 
                               env.df %>% summarise(across(.fns=mean)) %>% mutate(location=NA_character_),
                               pars$sc.df$canopyHeight.lm, "maxStipeLen")
  sizeClassLimits <- maxStipeLen * c(0, 0.333, 0.75, 1.25)#(0:3)/3
  sizeClassMdpts <- (sizeClassLimits+lag(sizeClassLimits))[-1]/2
  
  K_N <- pmax(1e-2, getPrediction(pars$N_canopy, ndraws, env.df,
                                 pars$sc.df$N_canopy.lm, "N"))
  K_FAI <- pmax(1e-2, getPrediction(pars$FAI, ndraws, env.df,
                                  pars$sc.df$FAI.lm, "FAI"))
  if(is.null(N0)) N0 <- K_N[1] * c(5,2,1)
  
  
  #---- per capita mass, area by stage
  logWtStipe.stage <- log(sizeClassMdpts) %>%
    map(~getPrediction(pars$lenSt_to_wtSt, ndraws, bind_cols(logLenStipe=.x, env.df),
                       pars$sc.df$lenSt_to_wtSt.lm, "logWtStipe")) %>%
    do.call('cbind', .)
  logWtFrond.stage <- log(sizeClassMdpts) %>%
    map(~getPrediction(pars$lenSt_to_wtFr, ndraws, bind_cols(logLenStipe=.x, env.df),
                       pars$sc.df$lenSt_to_wtFr.lm, "logWtFrond")) %>%
    do.call('cbind', .)
  logAreaFrond.stage <- log(sizeClassMdpts) %>%
    map(~getPrediction(pars$wtFr_to_arFr, ndraws, bind_cols(logWtFrond=.x, env.df),
                       pars$sc.df$wtFr_to_arFr, "logAreaFrond")) %>%
    do.call('cbind', .)
  
  
  #---- storage & initialization
  ## N[stage,year,season] = density/m2
  ## FAI[stage,year,season] = frond area / m2
  ## harvest[year] = log(grams harvested) / m2
  ## kappa[year,season,FAI|N] = proportion of K
  ## biomass[year,season] = kg / m2
  N <- FAI <- array(dim=c(3, pars$tmax, 3))
  N[,1,1] <- N0
  FAI[,1,1] <- N[,1,1] * exp(logAreaFrond.stage[1,])
  harvest <- rep(0, pars$tmax) 
  kappa <- array(dim=c(pars$tmax, 3, 2))
  
  
  #---- transition matrices
  ## A[[growing|non-growing]][stageTo,stageFrom]
  A <- map(1:2, ~matrix(0, 3, 3))
  
  
  #---- simulation loop
  for(year in 1:pars$tmax) {
    harvestYear <- year %% pars$freqHarvest == 0
    
    #---- growing season: 
    season <- 1
    kappa[year,season,] <- pmin(1, c(FAI[3,year,season]/K_FAI[year], N[3,year,season]/K_N[year]))

    # growth
    growRate_i <- par.yr$growStipeMax[year,] * (1-kappa[year,season,1]^pars$growthRateDensityShape)
    prGrowToNext <- pmin(1, pmax(0, growRate_i/(sizeClassLimits-lag(sizeClassLimits))[-1]))
    A[[1]][2,1] <- par.yr$surv[year,1]*prGrowToNext[1]
    A[[1]][3,2] <- par.yr$surv[year,2]*prGrowToNext[2]
    # survival
    A[[1]][1,1] <- par.yr$surv[year,1] - A[[1]][2,1]
    A[[1]][2,2] <- par.yr$surv[year,2] - A[[1]][3,2]
    A[[1]][3,3] <- par.yr$surv[year,3]
    # update population
    N[,year,season+1] <- A[[1]] %*% N[,year,season]
    FAI[,year,season+1] <- growFrondArea(FAI[,year,season], N[,year,season], 
                                         A[[1]], kappa[year,season,1],
                                         logAreaFrond.stage[year,], par.yr$growFrond[year,])
    
    
    #---- harvest:
    season <- 2
    kappa[year,season,] <- pmin(1, c(FAI[3,year,season]/K_FAI[year], N[3,year,season]/K_N[year]))
  
    # harvest
    if(harvestYear) {
      # update population
      N[,year,season+1] <- (1-pars$prFullHarvest) * N[,year,season]
      FAI[,year,season+1] <- (1-pars$prFullHarvest) * FAI[,year,season]
    } else {
      N[,year,season+1] <- N[,year,season]
      FAI[,year,season+1] <- FAI[,year,season]
    }
    
    
    #---- non-growing season
    season <- 3
    kappa[year,season,] <- pmin(1, c(FAI[3,year,season]/K_FAI[year], N[3,year,season]/K_N[year]))
    
    if(year < pars$tmax) {
      # survival
      diag(A[[2]]) <- pmax(0, par.yr$surv_strm[year,])
      # update population
      N[,year+1,1] <- A[[2]] %*% N[,year,season]
      # reproduction
      N[1,year+1,1] <- par.yr$settlement[year]*(1-max(kappa[year,season,]))
      FAI[,year+1,1] <- FAI[,year,season] * pmax(0, diag(A[[2]]) - par.yr$loss[year])
    }
    
  }
  
  
  # biomass calculation
  biomass <- calcBiomass(N, FAI, logWtStipe.stage, pars$arFr_to_wtFr,
                         ndraws, env.df, pars$sc.df$arFr_to_wtFr.lm, stages=3)
  
  return(list(N=N, FAI=FAI, harvest=harvest, kappa=kappa, K_FAI=K_FAI, K_N=K_N,
              biomass=biomass, PAR=env.df$PAR_atDepth))
}














#' Run Matrix Model Finer (Tom's handover code)
#'
#' @param inputs Named list of model parameters
#'
#' @return Named list with simulation outputs and key input parameters
#' 
runMatrixModelFiner <- function(inputs){
  
  # Unlist inputs
  domainArea <- inputs$domainArea  # Area in m^2
  Ntimesteps <- inputs$Ntimesteps
  timestep <- inputs$timestep
  NsizeClasses <- inputs$NsizeClasses
  maxSize <- inputs$maxSize
  settle <- inputs$settle * timestep  # Settlement per m^2 free space
  mu0 <- 1 - (1 - inputs$mu0)^timestep # Mortality in open space
  mu1 <- inputs$mu1*mu0
  # Growth (transition rate) in open space
  g1 <- inputs$g1
  allom1 <- inputs$allom1
  # allom2 <- inputs$allom2
  # allom3 <- inputs$allom3
  N0 <- inputs$N0
  k <- inputs$k  # Max canopy area allowing enough light for subcanopy growth (proportion)
  
  plotPars <- inputs$plotPars
  plotOut <- inputs$plotOut
  
  
  classStartSizes <- c(0:(NsizeClasses - 1))*maxSize/NsizeClasses
  classMidPointSizes <- (maxSize/NsizeClasses) * (c(0:(NsizeClasses - 1)) + 0.5)
  classSizes <- classStartSizes
  
  # Work out size dependent growth rate in m/day based on supplied parameters
  growthRateM <- growSize(classSizes, inputs$growthRateM, inputs$growthRateS)
  
  # Derive model transition rates from size growth rates
  g0 <- growthRateM*timestep*NsizeClasses/maxSize
  
  if (max(g0) > 1) {
    msg <- paste("Growth transition rate > 1",
                 "(number of size classes or timestep too high for given growth rate);",
                 "automatically set to 1")
    warning(msg)
    output$warnings <- msg
    g0 <- g0 / max(g0)
  }
  
  
  # Sort out initial population, which can either be a vector with length NsizeClasses, 
  # or single value (assumed all in first size class)
  if (length(N0)==1) {
    N0 <- c(N0, rep(0, NsizeClasses-1))
  }
  
  
  
  sizeClassHoldfastArea <- pi*(classSizes/2)^2 # Species specific allometry
  sizeClassCanopyArea <- allom1*sizeClassHoldfastArea # Species specific allometry
  
  if (plotPars & plotOut) {
    par(mfrow=c(2,3))
  } else {
    par(mfrow=c(1,3))
  }
  
  # If a constant growth rate (over size) has been specified, replace with a vector covering all size classes
  if (length(g0)==1){
    print("Single (constant) growth rate defined")
    g0 <- array(g0, dim=NsizeClasses)  
  }
  
  # Plot the functions to have a look at them
  xrange <- seq(0, 1, by=0.05)
  growthDensExample <- sapply(xrange, growDens, g0[1], g1, 1)
  
  mortLinearExample <- sapply(xrange, mortLinDens, mu0, mu1)
  mortOptimumExample <- sapply(xrange, mortOptDens, mu0, mu1)
  if (plotPars){
    plot(classSizes, growthRateM, type="l",
         main="", xlab='Holdfast diameter (m)', ylab="Growth rate (m/day)",
         ylim=c(0, growthRateM[1]))
    plot(NA, NA, 
         main="", xlab='Proportion canopy filled', ylab="Growth rate",
         xlim=c(0, 1), ylim=c(0, .1))
    for(i in 1:NsizeClasses) {
      lines(xrange, sapply(xrange, growDens, g0[i], g1, 1), col=i)
    }
    if(inputs$mortalityLinearDecline) {
      plot(xrange, mortLinearExample ,type="l",
           main="", xlab='Proportion canopy filled', ylab="Mortality rate",
           ylim=c(0, abs(mu0) + abs(mu1)))
    } else {
      plot(xrange, mortOptimumExample, type="l", col="red",
           main="", xlab='Proportion canopy filled', ylab="Mortality rate",
           ylim=c(0, abs(mu0) + abs(mu1)))
    }
  }
  
  # Set up  matrix to store population stage (rows) structure over time (cols)
  Nt <- matrix(0, ncol=Ntimesteps, nrow=NsizeClasses) 
  Nt[,1] <- N0 
  
  AreaHoldfast <- matrix(0, ncol=Ntimesteps, nrow=NsizeClasses)
  AreaCanopy <- matrix(0, ncol=Ntimesteps, nrow=NsizeClasses)
  yield <- array(0, dim=Ntimesteps)
  
  A.ar <- array(0, dim=c(NsizeClasses, NsizeClasses, Ntimesteps))
  
  for (step in 1:(Ntimesteps-1)) {     
    
    #Lt <- max(sum(AreaCanopy),0) # Calculate canopy area
    #print(paste("AreaCanopy = ",max(sum(AreaCanopy[,step]),0)))
    
    A <- matrix(0, nrow=NsizeClasses, ncol=NsizeClasses)
    
    for (i in 1:(NsizeClasses - 1)) {
      # Proportion of occupied space: asymmetric = plants larger than i
      if (inputs$interactType == "Asymmetric") {
        propOccupiedSpace <- sum(AreaCanopy[(i+1):NsizeClasses,step])/domainArea
      } else {
        propOccupiedSpace <- sum(AreaCanopy[,step])/domainArea
      }
      
      mortalityRate <- mortDens(propOccupiedSpace, mu0, mu1, 
                                inputs$mortalityLinearDecline)
      
      # Sub-diagonal entries: proportion that move to next class
      A[i+1,i] <- growDens(propOccupiedSpace, g0, g1, i)*(1 - mortalityRate)
      # Diagonal entries: proportion that stay in the same size class
      A[i,i] <- max(0, 1 - mortalityRate - A[i+1,i])
    }
    # Canopy (last stage) mortality
    A[NsizeClasses,NsizeClasses] <- max(0, 1 - mu0)
    A.ar[,,step] <- A
    
    Nt[,step+1]  <- A %*% Nt[,step] # Apply mortality and development
    
    Ft <- max(0, domainArea - sum(AreaHoldfast[,step])) # Calculate free space
    
    # Number of new recruits - assumed independent of population at present
    Nt[1,step+1] <- Nt[1,step+1] + settle*Ft
    
    # Thin stand of plants once every inputs$thinDuration
    if ((step*timestep) %% inputs$thinDuration < timestep) {
      Nt[,step+1] <- Nt[,step+1]*(1 - inputs$thin)
      yield[step + 1] <- yield[step] + sum(Nt[,step+1]*inputs$thin*sizeClassCanopyArea)
    } else {
      yield[step + 1] <- yield[step]
    }
    
    AreaHoldfast[,step+1] <- sizeClassHoldfastArea*Nt[,step+1]
    AreaCanopy[,step+1] <- sizeClassCanopyArea*Nt[,step+1]
  }
  # Add the final step areas to the matrices
  AreaHoldfast[,Ntimesteps] <- sizeClassHoldfastArea*Nt[,Ntimesteps]
  AreaCanopy[,Ntimesteps] <- sizeClassCanopyArea*Nt[,Ntimesteps]
  
  # Transpose Nt for plotting
  TNt <- t(Nt)
  if (plotOut) {
    image(log10(TNt + 1), main = "Number per size class",
          xlab = expression(paste("Time (years)")), ylab = "Holdfast diameter (m)",
          axes = F)
    axis(1, at = seq(0, 1, length = Ntimesteps),
         labels = sprintf("%.2f", c(1:Ntimesteps)*timestep/365),
         lwd = 0, pos = 0
    )
    axis(2, at = seq(0, 1, length = NsizeClasses),
         labels = sprintf("%.2f", c(1:NsizeClasses)*maxSize/NsizeClasses),
         lwd = 0, pos = 0
    )
  }
  TAh <- t(AreaHoldfast)
  TAc <- t(AreaCanopy)
  if (plotOut) {
    x <- c(1:Ntimesteps)*timestep / 365
    matplot(x, cbind(TAc, rowSums(TAc))/domainArea,
            main = "Canopy area (size classes)",
            xlab = "Time", ylab = expression(paste("Area (m" ^ "2", ")")),
            ylim = c(0, 1), type = "l"
    )
    lines(x, rowSums(TAc)/domainArea, lwd = 2)
    
    
  }
  if (plotOut == TRUE) {
    x <- c(1:Ntimesteps)*timestep/365
    plot(x, yield/domainArea,
         main = "Yield",
         xlab = "Time", ylab = expression(paste("Area (m" ^ "2", ")")),
         type = "l"
    )
  }
  
  output <- list(
    TNt=TNt,
    TAh=TAh,
    TAc=TAc,
    yield=yield,
    transitionMatrix=A,
    densrange=xrange,
    growVsSize=growthRateM,
    growVsDens=growthDensExample,
    classSizes=classSizes,
    sizeClassHoldfastArea=sizeClassHoldfastArea,
    sizeClassCanopyArea=sizeClassCanopyArea,
    domainArea=domainArea,
    NsizeClasses=NsizeClasses,
    Ntimesteps=Ntimesteps,
    maxSize=maxSize,
    timestep=timestep,
    mu0=mu0,
    mu1=mu1
  )
  
  if (inputs$mortalityLinearDecline){
    output$mortVsDens <- mortLinearExample
  } else {
    output$mortVsDens <- mortOptimumExample
  }
  
  return(output)
}