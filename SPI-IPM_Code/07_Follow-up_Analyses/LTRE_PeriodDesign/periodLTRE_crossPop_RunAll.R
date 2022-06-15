library(coda)
library(matrixStats)
library(reshape)
library(ggplot2)
library(crayon)

## Set data path
UserName <- 'chloe.nater'
DataPath <- paste0('/Users/', UserName, '/Dropbox/PiedFlycatcher_IPM/IPM_Code/210819_FlycatcherIPM_PostSamples_MS1/')

## Make a list of all populations
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')

## Set the study years for all populations
StudyYearsList <- list(
  DIN = 1980:2017,
  EDM = 1955:2020,
  KAT = 1974:2011,
  NAG = 1981:2019,
  NWA = 1986:2019,
  OKE = 1974:2019,
  TEI = 1981:2020
)

## Define the period for cross-population comparison
CompPeriod <- 1986:2011

## Define the indeces for the comparison period for each population
YearIndexList <- list(
  DIN = which(StudyYearsList$DIN %in% CompPeriod),
  EDM = which(StudyYearsList$EDM %in% CompPeriod),
  KAT = which(StudyYearsList$KAT %in% CompPeriod),
  NAG = which(StudyYearsList$NAG %in% CompPeriod),
  NWA = which(StudyYearsList$NWA %in% CompPeriod),
  OKE = which(StudyYearsList$OKE %in% CompPeriod),
  TEI = which(StudyYearsList$TEI %in% CompPeriod)
) 

## Determine how to treat immigration rates estimated at 0
# NOTE: Occasionally, there may be posterior samples that estimate immigration
#       rates at exactly 0 for entire time-periods. This is problematic for the
#       period-design LTRE, as calculations involve the log mean and log sd
#       difference of vital rates, and log(0) = Inf. 
#       Resulting contributions of immigration rates will default to NA, but 
#       contributions of other components for this posterior sample will still
#       have numerical values. This can potentially bias the posteriors of
#       relative contributions.
#       In the following, we present 2 options for avoiding such bias:

# Option 1:
FixLog0MeanRates <- FALSE

# Removes any samples which contain period-mean estimates of immigration rates
# equaling 0 completely prior to analysis.

#######################
#### DATA ASSEMBLY ####
#######################

## Set study years to period covered by all populations
StudyYears <- CompPeriod

## Load MCMC samples (matrix format) for all populations
PopID.MCMC <- list(
  DIN = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_DIN.rds'))),
  EDM = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_EDM.rds'))),
  KAT = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_KAT.rds'))),
  NAG = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_NAG.rds'))),
  NWA = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_NWA.rds'))),
  OKE = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_OKE.rds'))),
  TEI = as.matrix(readRDS(paste0(DataPath,'SPI-IPM_TEI.rds')))
)

## Set sample and year number
nosamples <- dim(PopID.MCMC[[1]])[1]
noyears <- length(StudyYears)

## Identify unusable samples - Part 1
#  NOTE: This includes samples that predicted population extinction (Ntot = 0) at some point
#        (in the example case here, there were none)
NAsamples1 <- list()
for(p in 1:length(PopID_List)){
  NAsamples1[[p]] <- which(PopID.MCMC[[p]][,paste0('Ntot[', YearIndexList[[p]], ']')] == 0, arr.ind = T)[,1]
}

## Identify unusable samples - Part 2

# Assemble time-specific immigration rates
imm_Y <- imm_A <- array(NA, dim = c(length(PopID_List), nosamples, noyears))
for(p in 1:length(PopID_List)){
  for(t in 1:(noyears-1)){
    imm_Y[p,,t+1] <- PopID.MCMC[[p]][, paste0('Imm[1, ', YearIndexList[[p]][t+1], ']')]/PopID.MCMC[[p]][, paste0('Ntot[', YearIndexList[[p]][t], ']')]
    imm_A[p,,t+1] <- PopID.MCMC[[p]][, paste0('Imm[2, ', YearIndexList[[p]][t+1], ']')]/PopID.MCMC[[p]][, paste0('Ntot[', YearIndexList[[p]][t], ']')]
  }
}


# Identify sample indeces with imm = 0 for the entire time periods
imm_Y_Mean0 <- imm_A_Mean0 <- list()
for(p in 1:length(PopID_List)){
  imm_Y_Mean0[[p]] <- which(rowMeans(imm_Y[p,,2:noyears]) == 0)
  imm_A_Mean0[[p]] <- which(rowMeans(imm_A[p,,2:noyears]) == 0)
}


# Set the number of samples to remove based on period-mean immigration rates of 0 (fix option 1)
if(!FixLog0MeanRates){
  NAsamples2 <- unique(c(unlist(imm_Y_Mean0), unlist(imm_A_Mean0)))
}else{
  NAsamples2 <- c()
}

# Remove temporary objects
rm(imm_Y, imm_A)

## Remove unusable samples
NAsamples <- c(unlist(NAsamples1), NAsamples2)
if(length(NAsamples) > 0){
  for(p in 1:length(PopID_List)){
    PopID.MCMC[[p]] <- PopID.MCMC[[p]][-NAsamples,]
  }
}

## Re-calculate sample number
nosamples <- dim(PopID.MCMC[[1]])[1]

###############
#### SETUP ####
###############

message(cyan('Assembling posterior data...'))

## Prepare matrices to rearrange samples - Vital rates & population sizes

# Time-varying vital rates
sJ <- sA <- pNS <- array(NA, dim = c(length(PopID_List), nosamples, noyears))
pB_Y <- CS_Y <- sN_Y <- array(NA, dim = c(length(PopID_List), nosamples, noyears))
pB_A <- CS_A <- sN_A <- array(NA, dim = c(length(PopID_List), nosamples, noyears))

# Time-varying population sizes/proportions and growth rates
N_Y <- n_Y <- array(NA, dim = c(length(PopID_List), nosamples, noyears))
N_A <- n_A <- array(NA, dim = c(length(PopID_List), nosamples, noyears))
N_tot <- array(NA, dim = c(length(PopID_List), nosamples, noyears))
lambda <- array(NA, dim = c(length(PopID_List), nosamples, noyears))

# Time-varying immigrant numbers & proportions
Imm_Y <- Imm_A <- imm_Y <- imm_A <- array(NA, dim = c(length(PopID_List), nosamples, noyears))

## Fill posterior samples into vectors and matrices
for(p in 1:length(PopID_List)){
  
  for(t in 1:noyears){
    
    # Time-varying population sizes & proportions
    N_Y[p,,t] <- PopID.MCMC[[p]][, paste0('N[1, ', YearIndexList[[p]][t], ']')]
    N_A[p,,t] <- PopID.MCMC[[p]][, paste0('N[2, ', YearIndexList[[p]][t], ']')]
    N_tot[p,,t] <- PopID.MCMC[[p]][, paste0('Ntot[', YearIndexList[[p]][t], ']')]
    
    n_Y[p,,t] <- N_Y[p,,t]/N_tot[p,,t]
    n_A[p,,t] <- N_A[p,,t]/N_tot[p,,t]
    
    # Time-varying vital rates   
    pB_Y[p,,t] <- PopID.MCMC[[p]][, paste0('pB[1, ', YearIndexList[[p]][t], ']')]
    pB_A[p,,t] <- PopID.MCMC[[p]][, paste0('pB[2, ', YearIndexList[[p]][t], ']')]
    
    CS_Y[p,,t] <- PopID.MCMC[[p]][, paste0('CS[1, ', YearIndexList[[p]][t], ']')]
    CS_A[p,,t] <- PopID.MCMC[[p]][, paste0('CS[2, ', YearIndexList[[p]][t], ']')]
    
    pNS[p,,t] <- PopID.MCMC[[p]][, paste0('pNS[', YearIndexList[[p]][t], ']')]
    
    sN_Y[p,,t] <- PopID.MCMC[[p]][, paste0('sN[1, ', YearIndexList[[p]][t], ']')]
    sN_A[p,,t] <- PopID.MCMC[[p]][, paste0('sN[2, ', YearIndexList[[p]][t], ']')]
    
    sJ[p,,t] <- PopID.MCMC[[p]][, paste0('sJ[', YearIndexList[[p]][t], ']')]
    sA[p,,t] <- PopID.MCMC[[p]][, paste0('sA[', YearIndexList[[p]][t], ']')]
    
    # Time-varying immigrant numbers/proportions
    Imm_Y[p,,t] <- PopID.MCMC[[p]][, paste0('Imm[1, ', YearIndexList[[p]][t], ']')]
    Imm_A[p,,t] <- PopID.MCMC[[p]][, paste0('Imm[2, ', YearIndexList[[p]][t], ']')]
    
    if(t > 1){
      imm_Y[p,,t] <- Imm_Y[p,,t]/N_tot[p,,t-1]
      imm_A[p,,t] <- Imm_A[p,,t]/N_tot[p,,t-1]
    }
  }
  
  for(t in 1:(noyears-1)){
    
    # Population growth rate
    lambda[p,,t] <- N_tot[p,,t+1]/N_tot[p,,t]
  }
}


###################################################
#### CALCULATION OF GEOMETRIC MEAN GROWTH RATE ####
###################################################

message(cyan('Calculating geometric mean growth rates...'))

## Prepare vectors
loggeolam <- matrix(NA, ncol = nosamples, nrow = length(PopID_List)) # Mean log lambda for each population

## Calculate geometric mean growth rate for each population and sample
for(p in 1:length(PopID_List)){
  for(i in 1:nosamples){
    loggeolam[p,i] <- mean(log(lambda[p,i,1:(noyears-1)]))
  }
}


## Summary statistics
for(p in 1:length(PopID_List)){
  message(paste0('Mean log lambda for Population ', PopID_List[p], ' (', min(CompPeriod), '-', max(CompPeriod), ')'))
  print(quantile(loggeolam[p,], probs = c(0.025, 0.5, 0.975), na.rm = T))
  message('')
}


#############################################################################
#### SIMULATION OF POPULATION DYNAMICS FOR A "MEAN" REFERENCE POPULATION ####
#############################################################################

message(cyan('Simulating dynamics for the reference population...'))

# NOTE: 
# Reference population = hypothetical population whose vital rates and 
# initial structure reflect per time-step averages across all populations. 
# The average vital rates and projected dynamics of the reference population
# are required to evaluate real-time elasticities for the LTRE.

## Define time period length
PeriodLength <- noyears

## Prepare arrays to store reference values
ref_sJ <- matrix(NA, nosamples, PeriodLength)
ref_sA <- matrix(NA, nosamples, PeriodLength)
ref_pB_Y <- matrix(NA, nosamples, PeriodLength)
ref_pB_A <- matrix(NA, nosamples, PeriodLength)
ref_CS_Y <- matrix(NA, nosamples, PeriodLength)
ref_CS_A <- matrix(NA, nosamples, PeriodLength)
ref_pNS <- matrix(NA, nosamples, PeriodLength)
ref_sN_Y <- matrix(NA, nosamples, PeriodLength)
ref_sN_A <- matrix(NA, nosamples, PeriodLength)

ref_n <- array(NA, dim = c(2, 1, nosamples, PeriodLength))
ref_imm <- array(NA, dim = c(2, 1, nosamples, PeriodLength))
# NOTE: These are given as arrays to ensure correct dimensions for matrix operations

lam_realref <- matrix(NA, nosamples, PeriodLength)

## Project reference population
for(i in 1:nosamples){
  
  # Define reference mean starting population structures
  ref_n[1,1,i,1] <- mean(n_Y[,i,1])
  ref_n[2,1,i,1] <- mean(n_A[,i,1])
  
  # Calculate reference mean vital rates
  ref_sJ[i,] <- colMeans(sJ[,i,])
  ref_sA[i,] <- colMeans(sA[,i,])
  ref_pB_Y[i,] <- colMeans(pB_Y[,i,])
  ref_pB_A[i,] <- colMeans(pB_A[,i,])
  ref_CS_Y[i,] <- colMeans(CS_Y[,i,])
  ref_CS_A[i,] <- colMeans(CS_A[,i,])
  ref_pNS[i,] <- colMeans(pNS[,i,])
  ref_sN_Y[i,] <- colMeans(sN_Y[,i,])
  ref_sN_A[i,] <- colMeans(sN_A[,i,])
  
  ref_imm[1,1,i,2:PeriodLength] <- colMeans(imm_Y[,i,2:PeriodLength])
  ref_imm[2,1,i,2:PeriodLength] <- colMeans(imm_A[,i,2:PeriodLength])
  
  for(t in 1:(PeriodLength-1)){
    
    # Formulate projection matrix 
    A <- matrix(NA, nrow = 2, ncol = 2)
    A[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    A[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    A[2,1] <- ref_sA[i,t]
    A[2,2] <- ref_sA[i,t]
    
    # Project reference population
    n_next <- A %*% ref_n[,,i,t] + ref_imm[,,i,t+1] 
    lam_realref[i,t] <- sum(n_next)
    ref_n[,1,i,t+1] <- n_next/sum(n_next)
  }
}

## Calculate temporal means of reference vital rates
ref_sJ.mu <- rowMeans(ref_sJ[,1:(PeriodLength-1)])
ref_sA.mu <- rowMeans(ref_sA[,1:(PeriodLength-1)])
ref_pB_Y.mu <- rowMeans(ref_pB_Y[,1:(PeriodLength-1)])
ref_pB_A.mu <- rowMeans(ref_pB_A[,1:(PeriodLength-1)])
ref_CS_Y.mu <- rowMeans(ref_CS_Y[,1:(PeriodLength-1)])
ref_CS_A.mu <- rowMeans(ref_CS_A[,1:(PeriodLength-1)])
ref_pNS.mu <- rowMeans(ref_pNS[,1:(PeriodLength-1)])
ref_sN_Y.mu <- rowMeans(ref_sN_Y[,1:(PeriodLength-1)])
ref_sN_A.mu <- rowMeans(ref_sN_A[,1:(PeriodLength-1)])

ref_imm_Y.mu <- rowMeans(ref_imm[1,1,,2:PeriodLength])
ref_imm_A.mu <- rowMeans(ref_imm[2,1,,2:PeriodLength])


## Run LTRE for all populations

for(x in 1:length(PopID_List)){
  PopID <- PopID_List[x]
  PopIDIdx <- x
  
  message(green$underline$bold(paste0('Cross-population period design LTRE analyses for ', PopID)))
  source('periodLTRE_crossPop_PopID.R')
}
