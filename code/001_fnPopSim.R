# KELPER
# Simulation functions
# Tim Szewczyk


# Population simulation functions


simulatePopulation_3Seasons <- function(pars, N0=c(1, 1, 1)) {
  
  # calculate global parameters
  logWtStipe.stage <- predict(pars$lenStipe_wtStipe,
                              tibble(logLenStipe=log(pars$sizeClassMdpts)), ndraws=5)[,1]
  K_FAI <- exp(predict(pars$depth_FAI, tibble(depth=pars$depth), ndraws=5)[1,1])
  K_N <- exp(predict(pars$depth_N, tibble(depth=pars$depth, stage="canopy"), ndraws=5)[1,1])
  
  # storage
  N <- FAI <- array(dim=c(pars$N_stages, pars$tmax, 3))
  N[,1,1] <- N0
  logWtFrond.stage <- predict(pars$wtStipe_wtFrond, 
                              newdata=tibble(logWtStipe=logWtStipe.stage, 
                                             propClear=1), ndraws=5)[,1]
  FAI[,1,1] <- exp(predict(pars$wtFrond_areaFrond,
                           tibble(depth=pars$depth, 
                                  logWtFrond=log(N[,1,1]*exp(logWtFrond.stage))), ndraws=5)[,1])
  harvest <- rep(0, pars$tmax) # log(grams of frond)
  kappa <- array(dim=c(pars$tmax, 3, 2))
  
  # transition matrices
  A <- map(1:3, ~matrix(0, pars$N_stages, pars$N_stages)) %>%
    setNames(c("spring", "fall", "winter"))
  
  
  # simulation loop
  for(i in 1:pars$tmax) {

    #-- first season: growth
    # parameters
    kappa[i,1,] <- c(min(1, FAI[3,i,1]/K_FAI), min(1, N[3,i,1]/K_N))
    growRate_i <- pars$growthRateStipeMax * (1-kappa[i,1,1])^pars$growthRateDensityStipeShape
    prGrowToNext <- growRate_i/(pars$sizeClassLimits-lag(pars$sizeClassLimits))[-1]
    # growth
    A$spring[2,1] <- pars$survRate[1]*prGrowToNext[1]
    A$spring[3,2] <- pars$survRate[2]*prGrowToNext[2]
    # survival
    A$spring[1,1] <- pars$survRate[1] - A$spring[2,1]
    A$spring[2,2] <- pars$survRate[2] - A$spring[3,2]
    A$spring[3,3] <- pars$survRate[3]
    # update population
    logWtFrond.stage <- predict(pars$wtStipe_wtFrond, 
                                newdata=tibble(logWtStipe=logWtStipe.stage, 
                                               propClear=1-kappa[i,1,1]), ndraws=5)[,1]
    logAreaFrond.stage <- predict(pars$wtFrond_areaFrond,
                                  tibble(depth=pars$depth,
                                         logWtFrond=logWtFrond.stage), ndraws=5)[,1]
    N[,i,2] <- A$spring %*% N[,i,1]
    FAI[,i,2] <- growFrondArea(FAI[,i,1], N[,i,1], A$spring, kappa[i,1,1],
                               logAreaFrond.stage, pars)
    
    
    #-- harvest
    # parameters
    kappa[i,2,] <- c(min(1, FAI[3,i,2]/K_FAI), min(1, N[3,i,2]/K_N))
    # survival & harvest
    diag(A$fall) <- pars$survRate
    if(i%%pars$freqHarvest==0) {
      A$fall[3,3] <- max(0, A$fall[3,3]-pars$prFullHarvest)
      stipeBiomass <- ifelse(pars$harvestTarget %in% c("all", "stipe"),
                             sum(exp(logWtStipe.stage) * N[,i,2]),
                             0)
      frondBiomass <- ifelse(pars$harvestTarget %in% c("all", "frond"),
                             sum(exp(predict(pars$areaFrond_wtFrond,
                                             tibble(depth=pars$depth, 
                                                    logAreaFrond=log(FAI[,i,2])), 
                                             ndraws=5)[,1])),
                             0)
      harvest[i] <- pars$prFullHarvest * (stipeBiomass + frondBiomass)

    } 
      
    # update population
    N[,i,3] <- A$fall %*% N[,i,2]
    FAI[,i,3] <- diag(A$fall) * FAI[,i,2]
    
    
    #-- second season: loss + reproduction
    if(i < pars$tmax) {
      # parameters
      
      kappa[i,3,] <- c(min(1, FAI[3,i,3]/K_FAI), min(1, N[3,i,3]/K_N))
      # survival
      diag(A$winter) <- pars$survRate
      # reproduction
      A$winter[1,3] <- pars$settlementRateBg*(1-max(kappa[i,3,]))
      # update population
      N[,i+1,1] <- A$winter %*% N[,i,3]
      FAI[,i+1,1] <- FAI[,i,3] * pmax(0, (diag(A$winter) - pars$lossRate))
    }
  
  }
  
  return(list(N=N, FAI=FAI, harvest=harvest, kappa=kappa))
}







simulatePopulation_2Seasons <- function(pars, N0=c(1, 1, 1)) {
  
  # calculate global parameters
  PAR <- pars$env$PAR_surface * exp(-pars$env$KD_mn * pars$depth)
  logWtStipe.stage <- predict(pars$lenStipe_wtStipe,
                              tibble(logLenStipe=log(pars$sizeClassMdpts)), ndraws=5)[,1]
  K_FAI <- exp(predict(pars$depth_FAI, tibble(logPAR=log(PAR)), ndraws=5)[1,1])
  K_N <- exp(predict(pars$depth_N, tibble(PAR_atDepth=PAR, stage="canopy"), ndraws=5)[1,1])
  
  # storage
  N <- FAI <- array(dim=c(pars$N_stages, pars$tmax, 3))
  N[,1,1] <- N0
  logWtFrond.stage <- predict(pars$wtStipe_wtFrond, 
                              newdata=tibble(logWtStipe=logWtStipe.stage, 
                                             propClear=1), ndraws=5)[,1]
  FAI[,1,1] <- exp(predict(pars$wtFrond_areaFrond,
                           tibble(PAR_atDepth=PAR, 
                                  logWtFrond=log(N[,1,1]*exp(logWtFrond.stage))), ndraws=5)[,1])
  harvest <- rep(0, pars$tmax) # log(grams of frond)
  kappa <- array(dim=c(pars$tmax, 3, 2))
  
  # transition matrices
  A <- map(1:2, ~matrix(0, pars$N_stages, pars$N_stages)) %>%
    setNames(c("growing", "nongrowing"))
  
  
  # simulation loop
  for(i in 1:pars$tmax) {
    
    #-- first season: 
    # parameters
    kappa[i,1,] <- c(min(1, FAI[3,i,1]/K_FAI), min(1, N[3,i,1]/K_N))
    growRate_i <- pars$growthRateStipeMax * (1-kappa[i,1,1])^pars$growthRateDensityStipeShape
    prGrowToNext <- growRate_i/(pars$sizeClassLimits-lag(pars$sizeClassLimits))[-1]
    # growth
    A[[1]][2,1] <- pars$survRate[1]*prGrowToNext[1]
    A[[1]][3,2] <- pars$survRate[2]*prGrowToNext[2]
    # survival
    A[[1]][1,1] <- pars$survRate[1] - A[[1]][2,1]
    A[[1]][2,2] <- pars$survRate[2] - A[[1]][3,2]
    A[[1]][3,3] <- pars$survRate[3]
    # update population
    logWtFrond.stage <- predict(pars$wtStipe_wtFrond, 
                                newdata=tibble(logWtStipe=logWtStipe.stage, 
                                               propClear=1-kappa[i,1,1]), ndraws=5)[,1]
    logAreaFrond.stage <- predict(pars$wtFrond_areaFrond,
                                  tibble(PAR_atDepth=PAR,
                                         logWtFrond=logWtFrond.stage), ndraws=5)[,1]
    N[,i,2] <- A[[1]] %*% N[,i,1]
    FAI[,i,2] <- growFrondArea(FAI[,i,1], N[,i,1], A[[1]], kappa[i,1,1],
                               logAreaFrond.stage, pars)
    
    
    #-- harvest
    # parameters
    kappa[i,2,] <- c(min(1, FAI[3,i,2]/K_FAI), min(1, N[3,i,2]/K_N))
    # harvest
    if(i%%pars$freqHarvest==0) {
      stipeBiomass <- ifelse(pars$harvestTarget %in% c("all", "stipe"),
                             sum(exp(logWtStipe.stage) * N[,i,2]),
                             0)
      frondBiomass <- ifelse(pars$harvestTarget %in% c("all", "frond"),
                             sum(exp(predict(pars$areaFrond_wtFrond,
                                             tibble(PAR_atDepth=PAR, 
                                                    logAreaFrond=log(FAI[,i,2])), 
                                             ndraws=5)[,1])),
                             0)
      harvest[i] <- pars$prFullHarvest * (stipeBiomass + frondBiomass)
      
      
      # update population
      N[,i,3] <- diag(1-pars$prFullHarvest, 3, 3) %*% N[,i,2]
      FAI[,i,3] <- (1-pars$prFullHarvest) * FAI[,i,2]
    } else {
      N[,i,3] <- N[,i,2]
      FAI[,i,3] <- FAI[,i,2]
    }
    
    
    #-- non-growing season
    if(i < pars$tmax) {
      # parameters
      kappa[i,3,] <- c(min(1, FAI[3,i,3]/K_FAI), min(1, N[3,i,3]/K_N))
      # survival
      diag(A[[2]]) <- pars$survRate
      # reproduction
      A[[2]][1,3] <- pars$settlementRateBg*(1-max(kappa[i,3,]))
      # update population
      N[,i+1,1] <- A[[2]] %*% N[,i,3]
      FAI[,i+1,1] <- FAI[,i,3] * pmax(0, (diag(A[[2]]) - pars$lossRate))
    }
    
  }
  
  return(list(N=N, FAI=FAI, harvest=harvest, kappa=kappa,
              K_FAI=K_FAI, K_N=K_N, PAR=PAR))
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