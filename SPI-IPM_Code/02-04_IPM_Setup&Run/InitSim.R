## Function to sample initial values
PFC_IPM.inits <- function(IPM.data, IPM.constants, sampleRE = FALSE){
  
  # Function inputs:
  # - IPM.data = data to be passed to IPM
  # - IPM.constants = constants to be passed to IPM
  # - sampleRE = whether or not to sample random effects levels for simulation
  
  # Additional notes:
  # - Covariate effects are set to 0 in the inital value simulation (beta = 0)
  #   for simplicity, but this can be changed. 
  
  
  #----------#
  # 0. SETUP #
  #----------#
  
  ## 0.1 Set number of age classes and time-steps
  A <- IPM.constants$A
  Tmax <- IPM.constants$Tmax
  
  ## 0.2 Prepare matrices/vectors for vital rates, population projection, and covariates
  sJ <- sA <- pNS <- rep(NA, Tmax)
  CS <- sN <- ImmAgeProp <- matrix(NA, nrow = A, ncol = Tmax)
  pB <- matrix(NA, nrow = A, ncol = Tmax+1)
  N <- B <- Juv <- localN <- Imm <- ImmB <- ImmCat <- matrix(NA, nrow = A, ncol = Tmax)
  EggNo.ex <- Fledged.ex <- matrix(NA, nrow = A, ncol = Tmax)

  cov1 <- cov2 <- cov3 <- rep(NA, Tmax)

  #---------------------------------#
  # 1. MISSING COVARIATE SIMULATION #
  #---------------------------------#
  
  for(t in which(is.na(IPM.data$cov1))){
    cov1[t] <- rnorm(1, mean = 0, sd = 1)
  }
  
  for(t in which(is.na(IPM.data$cov2))){
    cov2[t] <- rnorm(1, mean = 0, sd = 1)
  }
  
  for(t in which(is.na(IPM.data$cov3))){
    cov3[t] <- rnorm(1, mean = 0, sd = 1)
  }
  
  # NOTE: We assume that covariates  are a) normally distributed and 
  #       b) standardized (scaled and centered, such that mean = 0 and sd = 1).
  #       Other distributions of covariates can be accounted for by changing
  #       this part. 
  
  #--------------------------#
  # 2. VITAL RATE SIMULATION #
  #--------------------------#
  
  ## 2.1 Simulate intercepts
  Mu.sJ <-  runif(1, 0.05, 0.4)
  Mu.sA <-  runif(1, 0.25, 0.75)
  Mu.pB <-  c(runif(1, 0.1, 0.4), runif(1, 0.5, 0.9))
  Mu.CS <-  runif(2, 4, 7)
  Mu.sN <-  runif(2, 0.75, 0.95)
  Mu.pNS <-  runif(1, 0.6, 0.9)
  Mu.ND <-  runif(1, 0.05, 0.3)
  
  ## 2.2 Simulate covariate effect slopes
  #      (initialized at 0 for simplicity)
  beta1.pNS <- 0
  beta1.sN <- 0
  beta2.pNS <- 0
  beta2.sN <- 0
  beta3.sJ <- 0
  
  ## 2.3 Simulate random effect standard deviations
  sigma.sJ <- runif(1, 0.01, 0.5)
  sigma.sA <- runif(1, 0.01, 0.5)
  sigma.pB <- runif(1, 0.01, 0.5)
  sigma.CS <- runif(1, 0.01, 0.5)
  sigma.sN <- runif(1, 0.01, 0.5)
  sigma.pNS <- runif(1, 0.01, 0.5)
  sigma.ND <- runif(1, 0.01, 0.5)
  
  ## 2.4 Simulate random effect levels 
  #      (if sampleRE = TRUE, else set to 0)
  if(sampleRE == TRUE){
    epsilon.sJ <- runif(Tmax, 0, sigma.sJ)
    epsilon.sA <- runif(Tmax, 0, sigma.sA)
    epsilon.pB <- runif(Tmax+1, 0, sigma.pB)
    epsilon.CS <- runif(Tmax, 0, sigma.CS)
    epsilon.sN <- runif(Tmax, 0, sigma.sN)
    epsilon.ND  <- runif(Tmax, 0, sigma.ND)
    epsilon.pNS <- runif(Tmax, 0, sigma.pNS)
  }else{
    epsilon.sJ <- epsilon.sA <- epsilon.CS <- epsilon.sN <- epsilon.ND <- epsilon.pNS <- rep(0, Tmax)
    epsilon.pB <- rep(0, Tmax+1)
  }
    
  ## 2.5 Simulate age- and time-dependent vital rates
  sJ[1:Tmax] <- plogis(qlogis(Mu.sJ) + epsilon.sJ[1:Tmax])
  
  sA[1:Tmax] <- plogis(qlogis(Mu.sA) + epsilon.sA[1:Tmax])  
  
  ImmAgeProp[1,1:Tmax] <- plogis(qlogis(Mu.ND) + epsilon.ND[1:Tmax]) 
  ImmAgeProp[2,1:Tmax] <- 1 - ImmAgeProp[1,1:Tmax]
  
  
  for(a in 1:A){
    
    pB[a,1:(Tmax+1)] <- plogis(qlogis(Mu.pB[a]) + epsilon.pB[1:(Tmax+1)])
    
    CS[a,1:Tmax] <- exp(log(Mu.CS[a]) + epsilon.CS[1:Tmax])
    
    sN[a,1:Tmax] <- plogis(qlogis(Mu.sN[a]) + epsilon.sN[1:Tmax])
    
    pNS[1:Tmax] <- plogis(qlogis(Mu.pNS) + epsilon.pNS[1:Tmax])
  }
  
  #--------------------------#
  # 3. POPULATION SIMULATION #
  #--------------------------#
  
  ## Extract "average" age distribution of breeding population
  AvgAgeProps <- c(1-mean(IPM.constants$CS_FAge-1), mean(IPM.constants$CS_FAge-1))
  
  ## 3.2 Set initial population sizes
  B[1:A,1] <- round(mean(IPM.data$NestCount)*AvgAgeProps[1:A])
  N[1:A,1] <- round(B[1:A,1]/pB[1:A,1])
  localN[1:A,1] <- N[1:A,1]
  Imm[1:A,1] <- ImmB[1:A,1] <- ImmCat[1:A,1] <- 0
  
  ## 3.3 Project population forward in time
  for(t in 1:Tmax){
    
    ### Breeding decision
    if(t > 1){
      B[1:A,t] <- rbinom(n = A, size =  N[1:A,t], prob = pB[1:A,t])
    }
    
    ### Egg laying
    EggNo.ex[1:A,t] <- B[1:A,t]*CS[1:A,t]
    
    ### Nest success/failure & survival to fledging
    Fledged.ex[1:A,t] <- EggNo.ex[1:A,t]*pNS[t]*sN[1:A,t]
    
    ### Sex ratio
    Juv[1:A,t] <- rpois(n = A, lambda = Fledged.ex[1:A,t]*0.5)
    
    if(t < Tmax){
      ### Annual survival of local birds
      localN[1,t+1] <- rbinom(n = 1, size = sum(Juv[1:A,t]), prob = sJ[t])
      localN[2,t+1] <- rbinom(n = 1, size = sum(N[1:A,t]), prob = sA[t])
      
      ### Immigrant numbers
      AvgImmNo <- round(mean(IPM.data$ImmNoObs))
      ImmB[1:A,t+1] <- rmulti(n = 1, size = AvgImmNo, prob = ImmAgeProp[1:A,t+1])
      Imm[1:A,t+1] <- ifelse(
        round(ImmB[1:A,t+1]/pB[1:A,t+1]) < IPM.constants$AvgImm.limit, 
        round(ImmB[1:A,t+1]/pB[1:A,t+1]), 
        IPM.constants$AvgImm.limit
        ) 
      # NOTE: The ifelse statement above ensures that the simulated number will
      #       not exceed the upper limit used in the model
      ImmCat[1:A,t+1] <- Imm[1:A,t+1] + 1
      
      ### Immigration 
      N[1:A,t+1] <- localN[1:A,t+1] + Imm[1:A,t+1]
    }
  }
  
  ### Immigration summary
  AvgImm <- sigma.Imm <- rep(NA, A)
  for(a in 1:A){
    AvgImm[a] <- mean(Imm[a, 2:Tmax])
    sigma.Imm[a] <- sd(Imm[a, 2:Tmax])
  }
  
  ## 3.4 Calculate summary quantities
  Ntot <- colSums(N)
  Btot <- colSums(B)
  
  ## 3.5 Issue extinction warning (if an extinction was simulated)
  if(any(colSums(localN) == 0)){
    t_ext <- min(which(colSums(localN) == 0))
    warning(paste0("Simulation predicted local extinction at t = ", t_ext))
  }
  
  
  #-----------------------------------#
  # 4. DETECTION PARAMETER SIMULATION #
  #-----------------------------------#
  
  pImmDetect <- IPM.data$PropImmDetect
  pImmDetect[which(is.na(IPM.data$PropImmDetect))] <- runif(length(which(is.na(pImmDetect))), 0.1, 0.5)
  pImmDetect[which(!is.na(IPM.data$PropImmDetect))] <- NA
  
  
  #--------------------------#
  # 5. INITIAL VALUES OUTPUT #
  #--------------------------#
  
  ## 5.1 Arrange initial values in a list
  InitsList <- list(
    
    # Vital rate intercepts
    Mu.sJ = Mu.sJ,
    Mu.sA = Mu.sA,
    Mu.pB = Mu.pB,
    Mu.CS = Mu.CS,
    Mu.sN = Mu.sN,
    Mu.pNS = Mu.pNS,
    
    # Vital rate environmental effects
    beta1.pNS = beta1.pNS, 
    beta1.sN = beta1.sN,
    beta2.pNS = beta2.pNS,
    beta2.sN = beta2.sN,
    beta3.sJ = beta3.sJ,
    
    # Vital rate variances
    sigma.sJ = sigma.sJ,
    sigma.sA = sigma.sA,
    sigma.pB = sigma.pB,
    sigma.CS = sigma.CS,
    sigma.sN = sigma.sN,
    sigma.pNS = sigma.pNS,
    
    # Vital rate random effects
    epsilon.sJ = epsilon.sJ,
    epsilon.sA = epsilon.sA,
    epsilon.pB = epsilon.pB,
    epsilon.CS = epsilon.CS,
    epsilon.sN = epsilon.sN,
    epsilon.pNS = epsilon.pNS,
    
    # Age- and time-dependent vital rates
    sJ = sJ,
    sA = sA,
    pB = pB,
    CS = CS,
    sN = sN,
    pNS = pNS,
    
    # Population-level quantities
    N = N,
    B = B,
    Juv = Juv,
    localN = localN,
    Imm = Imm,
    ImmB = ImmB,
    ImmCont = Imm,
    AvgImm = AvgImm,
    sigma.Imm = sigma.Imm,
    EggNo.ex = EggNo.ex,
    Fledged.ex = Fledged.ex,
    
    # Detection parameters
    PropImmDetect = pImmDetect,
    
    # Missing covariate values
    cov1 = cov1,
    cov2 = cov2,
    cov3 = cov3
  )
  
  ## 5.2 Return initial values
  return(InitsList)
  
}
