#######################
#### DATA ASSEMBLY ####
#######################

message(green$underline$bold(paste0('Period design LTRE analyses for ', PopID)))

## Extract study years for focal population
StudyYears <- eval(parse(text = paste0('StudyYearsList$', PopID)))

## Load MCMC samples
load(paste0(DataPath,'/SPI-IPM_', PopID, '.RData'))
#str(PFC.IPM) # = posterior samples in coda format

## Transform samples into matrix (all chains combined)
out.mat <- as.matrix(PFC.IPM)
#str(out.mat) # = posterior samples in matrix format

## Set sample and year number
nosamples <- dim(out.mat)[1]
noyears <- length(StudyYears)

## Set time-indexes for equal-length periods
Periods <- eval(parse(text = paste0('PeriodsList$', PopID))) - min(StudyYears) + 1

Period1 <- Periods[1,]
Period2 <- Periods[2,]

## Identify unusable samples - Part 1
#  NOTE: This includes samples that predicted population extinction (Ntot = 0) at some point
#        as well as - optionally - samples with period-mean immigration rates of 0
NAsamples1 <- which(out.mat[,paste0('Ntot[', 1:length(StudyYears), ']')] == 0, arr.ind = T)[,1]

## Identify unusable samples - Part 2

# Assemble time-specific immigration rates
imm_Y <- imm_A <- matrix(NA, nrow = nosamples, ncol = noyears)
for(t in 1:(noyears-1)){
  imm_Y[,t+1] <- out.mat[, paste0('Imm[1, ', t+1, ']')]/out.mat[, paste0('Ntot[', t, ']')]
  imm_A[,t+1] <- out.mat[, paste0('Imm[2, ', t+1, ']')]/out.mat[, paste0('Ntot[', t, ']')]
}

# Identify sample indeces with imm = 0 for the entire time periods
imm_Y_Mean0_P1 <- which(rowMeans(imm_Y[,Period1+1]) == 0)
imm_A_Mean0_P1 <- which(rowMeans(imm_A[,Period1+1]) == 0)
imm_Y_Mean0_P2 <- which(rowMeans(imm_Y[,Period2+1]) == 0)
imm_A_Mean0_P2 <- which(rowMeans(imm_A[,Period2+1]) == 0)

# Set the number of samples to remove based on period-mean immigration rates of 0 (fix option 1)
if(!FixLog0MeanRates){
  NAsamples2 <- unique(c(imm_Y_Mean0_P1, imm_A_Mean0_P1,
                         imm_Y_Mean0_P2, imm_A_Mean0_P2))
}else{
  NAsamples2 <- c()
}

# Remove temporary objects
rm(imm_Y, imm_A)

## Remove unusable samples
NAsamples <- c(NAsamples1, NAsamples2)
if(length(NAsamples) > 0){
  out.mat <- out.mat[-NAsamples,]
}

## Re-calculate sample number
nosamples <- dim(out.mat)[1]

###############
#### SETUP ####
###############

message(cyan('Assembling posterior data...'))

## Prepare matrices to rearrange samples - Vital rates & population sizes

# Time-varying vital rates
sJ <- sA <- pNS <- matrix(NA, nrow = nosamples, ncol = noyears)
pB_Y <- CS_Y <- sN_Y <- matrix(NA, nrow = nosamples, ncol = noyears)
pB_A <- CS_A <- sN_A <- matrix(NA, nrow = nosamples, ncol = noyears)

# Time-varying population sizes and growth rates
N_Y <- matrix(NA, nrow = nosamples, ncol = noyears)
N_A <- matrix(NA, nrow = nosamples, ncol = noyears)
N_tot <- matrix(NA, nrow = nosamples, ncol = noyears)
lambda <- matrix(NA, nrow = nosamples, ncol = noyears)

# Time-varying immigrant numbers
Imm_Y <- Imm_A <- matrix(NA, nrow = nosamples, ncol = noyears)

## Fill posterior samples into vectors and matrices
for(t in 1:noyears){
  
  # Time-varying population sizes
  N_Y[,t] <- out.mat[, paste0('N[1, ', t, ']')]
  N_A[,t] <- out.mat[, paste0('N[2, ', t, ']')]
  N_tot[,t] <- out.mat[, paste0('Ntot[', t, ']')]
  
  # Time-varying immigrant numbers
  Imm_Y[,t] <- out.mat[, paste0('Imm[1, ', t, ']')]
  Imm_A[,t] <- out.mat[, paste0('Imm[2, ', t, ']')]
}

for(t in 1:(noyears-1)){
  
  # Time-varying vital rates   
  pB_Y[,t] <- out.mat[, paste0('pB[1, ', t, ']')]
  pB_A[,t] <- out.mat[, paste0('pB[2, ', t, ']')]
  
  CS_Y[,t] <- out.mat[, paste0('CS[1, ', t, ']')]
  CS_A[,t] <- out.mat[, paste0('CS[2, ', t, ']')]
  
  pNS[,t] <- out.mat[, paste0('pNS[', t, ']')]
  
  sN_Y[,t] <- out.mat[, paste0('sN[1, ', t, ']')]
  sN_A[,t] <- out.mat[, paste0('sN[2, ', t, ']')]
  
  sJ[,t] <- out.mat[, paste0('sJ[', t, ']')]
  sA[,t] <- out.mat[, paste0('sA[', t, ']')]
  
  # Population growth rate
  lambda[,t] <- N_tot[,t+1]/N_tot[,t]
}

## Calculate population/immigrant proportions
n_Y <- N_Y/N_tot
n_A <- N_A/N_tot
imm_Y <- cbind(NA,Imm_Y[,2:noyears]/N_tot[,1:(noyears-1)])
imm_A <- cbind(NA,Imm_A[,2:noyears]/N_tot[,1:(noyears-1)])


###################################################
#### CALCULATION OF GEOMETRIC MEAN GROWTH RATE ####
###################################################

message(cyan('Calculating geometric mean growth rates...'))

## Prepare vectors
loggeolam_1 <- rep(NA, nosamples) # Mean log lambda for period 1
loggeolam_2 <- rep(NA, nosamples) # Mean log lambda for period 2
diffgeolam <- rep(NA, nosamples) # Difference of log lambda (period 1 vs. period 2)


## Calculate geometric mean growth rate for each sample
for(i in 1:nosamples){
  loggeolam_1[i] <- mean(log(lambda[i,Period1]))
  loggeolam_2[i] <- mean(log(lambda[i,Period2]))
  diffgeolam[i] <- loggeolam_2[i]-loggeolam_1[i]
}

## Summary statistics
message('Mean log lambda for Period 1:')
print(quantile(loggeolam_1, probs = c(0.025, 0.5, 0.975), na.rm = T))

message('Mean log lambda for Period 2:')
print(quantile(loggeolam_2, probs = c(0.025, 0.5, 0.975), na.rm = T))

message('Difference in mean log lambda between periods:')
print(quantile(diffgeolam, probs = c(0.025, 0.5, 0.975), na.rm = T))

#############################################################################
#### SIMULATION OF POPULATION DYNAMICS FOR A "MEAN" REFERENCE POPULATION ####
#############################################################################

message(cyan('Simulating dynamics for the reference population...'))

# NOTE: 
# Reference population = hypothetical population whose vital rates and 
# initial structure reflect per time-step averages across both time periods. 
# The average vital rates and projected dynamics of the reference population
# are required to evaluate real-time elasticities for the LTRE.

## Define time period length
PeriodLength <- ncol(Periods)

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

ref_n <- array(NA, dim = c(2, 1, nosamples, PeriodLength+1))
ref_imm <- array(NA, dim = c(2, 1, nosamples, PeriodLength+1))
# NOTE: These are given as arrays to ensure correct dimensions for matrix operations

lam_realref <- matrix(NA, nosamples, PeriodLength)


## Project reference population
for(i in 1:nosamples){
  
  # Define reference mean starting population structures
  ref_n[1,1,i,1] <- mean(n_Y[i,c(Period1[1],Period2[1])])
  ref_n[2,1,i,1] <- mean(n_A[i,c(Period1[1],Period2[1])])
  
  # Calculate reference mean vital rates
  ref_sJ[i,] <- (sJ[i,Period1] + sJ[i,Period2])/2
  ref_sA[i,] <- (sA[i,Period1] + sA[i,Period2])/2
  ref_pB_Y[i,] <- (pB_Y[i,Period1] + pB_Y[i,Period2])/2
  ref_pB_A[i,] <- (pB_A[i,Period1] + pB_A[i,Period2])/2
  ref_CS_Y[i,] <- (CS_Y[i,Period1] + CS_Y[i,Period2])/2
  ref_CS_A[i,] <- (CS_A[i,Period1] + CS_A[i,Period2])/2
  ref_pNS[i,] <- (pNS[i,Period1] + pNS[i,Period2])/2
  ref_sN_Y[i,] <- (sN_Y[i,Period1] + sN_Y[i,Period2])/2
  ref_sN_A[i,] <- (sN_A[i,Period1] + sN_A[i,Period2])/2

  ref_imm[1,1,i,2:(PeriodLength+1)] <- (imm_Y[i,Period1+1] + imm_Y[i,Period2+1])/2
  ref_imm[2,1,i,2:(PeriodLength+1)] <- (imm_A[i,Period1+1] + imm_A[i,Period2+1])/2
  
  for(t in 1:PeriodLength){

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
ref_sJ.mu <- rowMeans(ref_sJ)
ref_sA.mu <- rowMeans(ref_sA)
ref_pB_Y.mu <- rowMeans(ref_pB_Y)
ref_pB_A.mu <- rowMeans(ref_pB_A)
ref_CS_Y.mu <- rowMeans(ref_CS_Y)
ref_CS_A.mu <- rowMeans(ref_CS_A)
ref_pNS.mu <- rowMeans(ref_pNS)
ref_sN_Y.mu <- rowMeans(ref_sN_Y)
ref_sN_A.mu <- rowMeans(ref_sN_A)

ref_imm_Y.mu <- rowMeans(ref_imm[1,1,,2:(PeriodLength+1)])
ref_imm_A.mu <- rowMeans(ref_imm[2,1,,2:(PeriodLength+1)])


####################################################################################
#### CALCULATION OF LOG-DIFFERENCES OF VITAL RATE MEANS AND STANDARD DEVIATIONS ####
####################################################################################

message(cyan('Calculating vital rate log-differences...'))

## Prepare vectors to store numbers
logmudiff_sJ <- logmudiff_sA <- rep(NA, nosamples)
logmudiff_pB_Y <- logmudiff_pB_A <- rep(NA, nosamples)
logmudiff_CS_Y <- logmudiff_CS_A <- rep(NA, nosamples)
logmudiff_pNS <- rep(NA, nosamples)
logmudiff_sN_Y <- logmudiff_sN_A <- rep(NA, nosamples)
logmudiff_imm_Y <- logmudiff_imm_A <- rep(NA, nosamples)

logsigdiff_sJ <- logsigdiff_sA <- rep(NA, nosamples)
logsigdiff_pB_Y <- logsigdiff_pB_A <- rep(NA, nosamples)
logsigdiff_CS_Y <- logsigdiff_CS_A <- rep(NA, nosamples)
logsigdiff_pNS <- rep(NA, nosamples)
logsigdiff_sN_Y <- logsigdiff_sN_A <- rep(NA, nosamples)
logsigdiff_imm_Y <- logsigdiff_imm_A <- rep(NA, nosamples)


for(i in 1:nosamples){
  
  ## Calculate log differences of vital rate means
  logmudiff_sJ[i] <- log(mean(sJ[i,Period1])) - log(mean(sJ[i,Period2]))
  logmudiff_sA[i] <- log(mean(sA[i,Period1])) - log(mean(sA[i,Period2]))
  logmudiff_pB_Y[i] <- log(mean(pB_Y[i,Period1])) - log(mean(pB_Y[i,Period2]))
  logmudiff_pB_A[i] <- log(mean(pB_A[i,Period1])) - log(mean(pB_A[i,Period2]))
  logmudiff_CS_Y[i] <- log(mean(CS_Y[i,Period1])) - log(mean(CS_Y[i,Period2]))
  logmudiff_CS_A[i] <- log(mean(CS_A[i,Period1])) - log(mean(CS_A[i,Period2]))
  logmudiff_pNS[i] <- log(mean(pNS[i,Period1])) - log(mean(pNS[i,Period2]))
  logmudiff_sN_Y[i] <- log(mean(sN_Y[i,Period1])) - log(mean(sN_Y[i,Period2]))
  logmudiff_sN_A[i] <- log(mean(sN_A[i,Period1])) - log(mean(sN_A[i,Period2]))
  
  logmudiff_imm_Y[i] <- log(mean(imm_Y[i,Period1+1])) - log(mean(imm_Y[i,Period2+1]))
  logmudiff_imm_A[i] <- log(mean(imm_A[i,Period1+1])) - log(mean(imm_A[i,Period2+1]))
  
  ## Calculate log differences of vital rate standard deviations
  logsigdiff_sJ[i] <- log(sd(sJ[i,Period1])) - log(sd(sJ[i,Period2]))
  logsigdiff_sA[i] <- log(sd(sA[i,Period1])) - log(sd(sA[i,Period2]))
  logsigdiff_pB_Y[i] <- log(sd(pB_Y[i,Period1])) - log(sd(pB_Y[i,Period2]))
  logsigdiff_pB_A[i] <- log(sd(pB_A[i,Period1])) - log(sd(pB_A[i,Period2]))
  logsigdiff_CS_Y[i] <- log(sd(CS_Y[i,Period1])) - log(sd(CS_Y[i,Period2]))
  logsigdiff_CS_A[i] <- log(sd(CS_A[i,Period1])) - log(sd(CS_A[i,Period2]))
  logsigdiff_pNS[i] <- log(sd(pNS[i,Period1])) - log(sd(pNS[i,Period2]))
  logsigdiff_sN_Y[i] <- log(sd(sN_Y[i,Period1])) - log(sd(sN_Y[i,Period2]))
  logsigdiff_sN_A[i] <- log(sd(sN_A[i,Period1])) - log(sd(sN_A[i,Period2]))
  
  logsigdiff_imm_Y[i] <- log(sd(imm_Y[i,Period1+1])) - log(sd(imm_Y[i,Period2+1]))
  logsigdiff_imm_A[i] <- log(sd(imm_A[i,Period1+1])) - log(sd(imm_A[i,Period2+1]))
}

## Optional: Re-calculate immigration rate log-means and log-sd's that are 
#            (-)Inf (fix option 3)
if(FixLog0MeanRates){
  
  # Extract cross-sample minimum non-zero sd
  immY_sd_min <- min(rowSds(imm_Y, na.rm = T)[which(rowSds(imm_Y, na.rm = T) > 0)]) 
  immA_sd_min <- min(rowSds(imm_A, na.rm = T)[which(rowSds(imm_A, na.rm = T) > 0)])
  
  for(i in 1:nosamples){
    
    # Set the period mean to use in calculations (minimum estimated if 0)
    YP1_mean <- ifelse(i %in% imm_Y_Mean0_P1, 
                       min(imm_Y[which(imm_Y > 0)]), mean(imm_Y[i,Period1+1]))
    AP1_mean <- ifelse(i %in% imm_A_Mean0_P1, 
                       min(imm_A[which(imm_A > 0)]), mean(imm_A[i,Period1+1]))
    YP2_mean <- ifelse(i %in% imm_Y_Mean0_P2, 
                       min(imm_Y[which(imm_Y > 0)]), mean(imm_Y[i,Period2+1]))
    AP2_mean <- ifelse(i %in% imm_A_Mean0_P2, 
                       min(imm_A[which(imm_A > 0)]), mean(imm_A[i,Period2+1]))
    
    # Re-calculate difference of log means for immigration rates
    logmudiff_imm_Y[i] <- log(YP1_mean) - log(YP2_mean)
    logmudiff_imm_A[i] <- log(AP1_mean) - log(AP2_mean)
    
    # Set the period sd to use in calculatios (minimum sd estimated if 0)
    YP1_sd <- ifelse(i %in% imm_Y_Mean0_P1, 
                     immY_sd_min, sd(imm_Y[i,Period1+1]))
    AP1_sd <- ifelse(i %in% imm_A_Mean0_P1, 
                     immA_sd_min, sd(imm_A[i,Period1+1]))
    YP2_sd <- ifelse(i %in% imm_Y_Mean0_P2, 
                     immY_sd_min, sd(imm_Y[i,Period2+1]))
    AP2_sd <- ifelse(i %in% imm_A_Mean0_P2, 
                     immA_sd_min, sd(imm_A[i,Period2+1]))
    
    # Re-calculate difference of log sds for immigration rates
    logsigdiff_imm_Y[i] <- log(YP1_sd) - log(YP2_sd)
    logsigdiff_imm_A[i] <- log(AP1_sd) - log(AP2_sd)
    
  }
}

## Write a message detailing the number of samples with immigration rate 
## estimates of 0 for entire time-periods, and how they were treated
message(paste0('Numbers of samples with period-mean immigration rates of 0:
               
  Yearling: ', length(unique(c(imm_Y_Mean0_P1, imm_Y_Mean0_P2))), '
  ',
  'Adult: ', length(unique(c(imm_A_Mean0_P1, imm_A_Mean0_P2))),
  '
  '))


if(!FixLog0MeanRates){
  message('Affected samples were removed entirely prior to analysis.')
}

if(FixLog0MeanRates){
  message('Immigration rate period-mean and period-sd estimates of 0 were 
replaced with the overall minimum non-0 values in all affected samples prior to
analysis.')
}


#############################################################
#### CALCULATION OF REAL-TIME ELASTICITIES - VITAL RATES ####
#############################################################

message(cyan('Calculating real-time elasticities for direct effects...'))

# NOTE: 
# This first set of real-time elasticities are for the direct effects of changes
# in the vital rates and are calculated according to equation S1.5 and S1.6 in 
# Koons et al. 2016, and evaluated at the reference population simulated above. 

## Set dimension of projection matrix
S <- 2

## Prepare arrays for storing elasticity matrices
eA.mu_sJ <- eA.mu_sA <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.mu_pB_Y <- eA.mu_pB_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.mu_CS_Y <- eA.mu_CS_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.mu_pNS <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.mu_sN_Y <- eA.mu_sN_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.mu_imm_Y <- eA.mu_imm_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))

tot_eA.mu_sJ <- tot_eA.mu_sA <- matrix(NA, nosamples, PeriodLength)
tot_eA.mu_pB_Y <- tot_eA.mu_pB_A <- matrix(NA, nosamples, PeriodLength)
tot_eA.mu_CS_Y <- tot_eA.mu_CS_A <- matrix(NA, nosamples, PeriodLength)
tot_eA.mu_pNS <- matrix(NA, nosamples, PeriodLength)
tot_eA.mu_sN_Y <- tot_eA.mu_sN_A <- matrix(NA, nosamples, PeriodLength)
tot_eA.mu_imm_Y <- tot_eA.mu_imm_A <- matrix(NA, nosamples, PeriodLength)

eA.sig_sJ <- eA.sig_sA <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.sig_pB_Y <- eA.sig_pB_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.sig_CS_Y <- eA.sig_CS_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.sig_pNS <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.sig_sN_Y <- eA.sig_sN_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
eA.sig_imm_Y <- eA.sig_imm_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))

tot_eA.sig_sJ <- tot_eA.sig_sA <- matrix(NA, nosamples, PeriodLength)
tot_eA.sig_pB_Y <- tot_eA.sig_pB_A <- matrix(NA, nosamples, PeriodLength)
tot_eA.sig_CS_Y <- tot_eA.sig_CS_A <- matrix(NA, nosamples, PeriodLength)
tot_eA.sig_pNS <- matrix(NA, nosamples, PeriodLength)
tot_eA.sig_sN_Y <- tot_eA.sig_sN_A <- matrix(NA, nosamples, PeriodLength)
tot_eA.sig_imm_Y <- tot_eA.sig_imm_A <- matrix(NA, nosamples, PeriodLength)


## Calculating direct effect real-time elasticities
for(i in 1:nosamples){
  for(t in 1:PeriodLength){
    
    # a) Make derivative of projection matrix with regards to each vital rate
    #    (dA[t]/dVR[t]), in the matrix format.
    #------------------------------------------------------------------------
    
    # sJ
    dA_dsJ <- matrix(0, nrow = S, ncol = S)
    dA_dsJ[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]
    dA_dsJ[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]
    
    # sA
    dA_dsA <- matrix(c(0, 1, 0, 1), nrow = S)
    
    # pB_Y
    dA_dpB_Y <- matrix(0, nrow = S, ncol = S)
    dA_dpB_Y[1,1] <- 0.5*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]

    # pB_A
    dA_dpB_A <- matrix(0, nrow = S, ncol = S)
    dA_dpB_A[1,2] <- 0.5*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]

    # CS_Y
    dA_dCS_Y <- matrix(0, nrow = S, ncol = S)
    dA_dCS_Y[1,1] <- 0.5*ref_pB_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    
    # CS_A
    dA_dCS_A <- matrix(0, nrow = S, ncol = S)
    dA_dCS_A[1,2] <- 0.5*ref_pB_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    
    # pNS
    dA_dpNS <- matrix(0, nrow = S, ncol = S)
    dA_dpNS[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    dA_dpNS[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    
    # sN_Y
    dA_dsN_Y <- matrix(0, nrow = S, ncol = S)
    dA_dsN_Y[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sJ[i,t]
    
    # sN_A
    dA_dsN_A <- matrix(0, nrow = S, ncol = S)
    dA_dsN_A[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sJ[i,t]
    
    for(m in 1:S){
      for(n in 1:S){
        
        # b) Calculate real-time elasticities to changes in the mean (per matrix 
        #    element) by multiplying the sensitivity from above with the 
        #    period-mean VR and dividing by lambda
        # -----------------------------------------------------------------------
        
        eA.mu_sJ[m,n,i,t] <- ref_sJ.mu[i]*dA_dsJ[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_sA[m,n,i,t] <- ref_sA.mu[i]*dA_dsA[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_pB_Y[m,n,i,t] <- ref_pB_Y.mu[i]*dA_dpB_Y[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_pB_A[m,n,i,t] <- ref_pB_A.mu[i]*dA_dpB_A[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_CS_Y[m,n,i,t] <- ref_CS_Y.mu[i]*dA_dCS_Y[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_CS_A[m,n,i,t] <- ref_CS_A.mu[i]*dA_dCS_A[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_pNS[m,n,i,t] <- ref_pNS.mu[i]*dA_dpNS[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_sN_Y[m,n,i,t] <- ref_sN_Y.mu[i]*dA_dsN_Y[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.mu_sN_A[m,n,i,t] <- ref_sN_A.mu[i]*dA_dsN_A[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        
        
        # c) Calculate real-time elasticities to changes standard dev. (per matrix 
        #    element) by multiplying the sensitivity from above with the
        #    period-sd of the VR and dividing by lambda
        # -----------------------------------------------------------------------
        
        eA.sig_sJ[m,n,i,t] <- (ref_sJ[i,t]-ref_sJ.mu[i])*dA_dsJ[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_sA[m,n,i,t] <- (ref_sA[i,t]-ref_sA.mu[i])*dA_dsA[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_pB_Y[m,n,i,t] <- (ref_pB_Y[i,t]-ref_pB_Y.mu[i])*dA_dpB_Y[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_pB_A[m,n,i,t] <- (ref_pB_A[i,t]-ref_pB_A.mu[i])*dA_dpB_A[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_CS_Y[m,n,i,t] <- (ref_CS_Y[i,t]-ref_CS_Y.mu[i])*dA_dCS_Y[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_CS_A[m,n,i,t] <- (ref_CS_A[i,t]-ref_CS_A.mu[i])*dA_dCS_A[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_pNS[m,n,i,t] <- (ref_pNS[i,t]-ref_pNS.mu[i])*dA_dpNS[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_sN_Y[m,n,i,t] <- (ref_sN_Y[i,t]-ref_sN_Y.mu[i])*dA_dsN_Y[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
        eA.sig_sN_A[m,n,i,t] <- (ref_sN_A[i,t]-ref_sN_A.mu[i])*dA_dsN_A[m,n]*ref_n[n,1,i,t]/lam_realref[i,t]
      }
    }

    
    # d) Sum up elasticities to vital rates through different matrix elements
    # -----------------------------------------------------------------------
    
    tot_eA.mu_sJ[i,t] <- sum(eA.mu_sJ[,,i,t])
    tot_eA.mu_sA[i,t] <- sum(eA.mu_sA[,,i,t])
    tot_eA.mu_pB_Y[i,t] <- sum(eA.mu_pB_Y[,,i,t])
    tot_eA.mu_pB_A[i,t] <- sum(eA.mu_pB_A[,,i,t])
    tot_eA.mu_CS_Y[i,t] <- sum(eA.mu_CS_Y[,,i,t])
    tot_eA.mu_CS_A[i,t] <- sum(eA.mu_CS_A[,,i,t])
    tot_eA.mu_pNS[i,t] <- sum(eA.mu_pNS[,,i,t])
    tot_eA.mu_sN_Y[i,t] <- sum(eA.mu_sN_Y[,,i,t])
    tot_eA.mu_sN_A[i,t] <- sum(eA.mu_sN_A[,,i,t])
    
    tot_eA.sig_sJ[i,t] <- sum(eA.sig_sJ[,,i,t])
    tot_eA.sig_sA[i,t] <- sum(eA.sig_sA[,,i,t])
    tot_eA.sig_pB_Y[i,t] <- sum(eA.sig_pB_Y[,,i,t])
    tot_eA.sig_pB_A[i,t] <- sum(eA.sig_pB_A[,,i,t])
    tot_eA.sig_CS_Y[i,t] <- sum(eA.sig_CS_Y[,,i,t])
    tot_eA.sig_CS_A[i,t] <- sum(eA.sig_CS_A[,,i,t])
    tot_eA.sig_pNS[i,t] <- sum(eA.sig_pNS[,,i,t])
    tot_eA.sig_sN_Y[i,t] <- sum(eA.sig_sN_Y[,,i,t])
    tot_eA.sig_sN_A[i,t] <- sum(eA.sig_sN_A[,,i,t])
    
    # e) Calculate the total elasticity to immigration rates
    # ------------------------------------------------------
    
    tot_eA.mu_imm_Y[i,t] <- (ref_imm_Y.mu[i])/lam_realref[i,t]
    tot_eA.mu_imm_A[i,t] <- (ref_imm_A.mu[i])/lam_realref[i,t]

    tot_eA.sig_imm_Y[i,t] <- (ref_imm[1,1,i,t+1]-ref_imm_Y.mu[i])/lam_realref[i,t]
    tot_eA.sig_imm_A[i,t] <- (ref_imm[2,1,i,t+1]-ref_imm_A.mu[i])/lam_realref[i,t]
    # NOTE: At this stage, the t+1 index of immigration rate is translated into
    #       the t index of the elasticity:
    #       eA_imm[t] = elasticity of lambda[t] with respect to imm[t+1] 
    
  }
}

## Summarise average elasticities across years
avg_eA.mu_sJ <- rowMeans(tot_eA.mu_sJ)
avg_eA.mu_sA <- rowMeans(tot_eA.mu_sA)
avg_eA.mu_pB_Y <- rowMeans(tot_eA.mu_pB_Y)
avg_eA.mu_pB_A <- rowMeans(tot_eA.mu_pB_A)
avg_eA.mu_CS_Y <- rowMeans(tot_eA.mu_CS_Y)
avg_eA.mu_CS_A <- rowMeans(tot_eA.mu_CS_A)
avg_eA.mu_pNS <- rowMeans(tot_eA.mu_pNS)
avg_eA.mu_sN_Y <- rowMeans(tot_eA.mu_sN_Y)
avg_eA.mu_sN_A <- rowMeans(tot_eA.mu_sN_A)

avg_eA.mu_imm_Y <- rowMeans(tot_eA.mu_imm_Y)
avg_eA.mu_imm_A <- rowMeans(tot_eA.mu_imm_A)

avg_eA.sig_sJ <- rowMeans(tot_eA.sig_sJ)
avg_eA.sig_sA <- rowMeans(tot_eA.sig_sA)
avg_eA.sig_pB_Y <- rowMeans(tot_eA.sig_pB_Y)
avg_eA.sig_pB_A <- rowMeans(tot_eA.sig_pB_A)
avg_eA.sig_CS_Y <- rowMeans(tot_eA.sig_CS_Y)
avg_eA.sig_CS_A <- rowMeans(tot_eA.sig_CS_A)
avg_eA.sig_pNS <- rowMeans(tot_eA.sig_pNS)
avg_eA.sig_sN_Y <- rowMeans(tot_eA.sig_sN_Y)
avg_eA.sig_sN_A <- rowMeans(tot_eA.sig_sN_A)

avg_eA.sig_imm_Y <- rowMeans(tot_eA.sig_imm_Y)
avg_eA.sig_imm_A <- rowMeans(tot_eA.sig_imm_A)


###################################################################################
#### CALCULATION OF REAL-TIME ELASTICITIES - POPULATION STRUCTURE PERTURBATION ####
###################################################################################

message(cyan('Calculating real-time elasticities for indirect effects...'))

# NOTE: 
# This second set of real-time elasticities are for the indirect effects of changes
# in the vital rates mediated through population structure, are calculated 
# according to equation S1.7 and S1.8 in Koons et al. 2016, and evaluated at the
# reference population simulated above. 

## Prepare arrays for storing elasticity matrices
en.mu_sJ <- en.mu_sA <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.mu_pB_Y <- en.mu_pB_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.mu_CS_Y <- en.mu_CS_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.mu_pNS <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.mu_sN_Y <- en.mu_sN_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.mu_imm_Y <- en.mu_imm_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))

tot_en.mu_sJ <- tot_en.mu_sA <- matrix(NA, nosamples, PeriodLength)
tot_en.mu_pB_Y <- tot_en.mu_pB_A <- matrix(NA, nosamples, PeriodLength)
tot_en.mu_CS_Y <- tot_en.mu_CS_A <- matrix(NA, nosamples, PeriodLength)
tot_en.mu_pNS <- matrix(NA, nosamples, PeriodLength)
tot_en.mu_sN_Y <- tot_en.mu_sN_A <- matrix(NA, nosamples, PeriodLength)
tot_en.mu_imm_Y <- tot_en.mu_imm_A <- matrix(NA, nosamples, PeriodLength)

en.sig_sJ <- en.sig_sA <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.sig_pB_Y <- en.sig_pB_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.sig_CS_Y <- en.sig_CS_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.sig_pNS <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.sig_sN_Y <- en.sig_sN_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))
en.sig_imm_Y <- en.sig_imm_A <- array(NA, dim = c(S, S, nosamples, PeriodLength))

tot_en.sig_sJ <- tot_en.sig_sA <- matrix(NA, nosamples, PeriodLength)
tot_en.sig_pB_Y <- tot_en.sig_pB_A <- matrix(NA, nosamples, PeriodLength)
tot_en.sig_CS_Y <- tot_en.sig_CS_A <- matrix(NA, nosamples, PeriodLength)
tot_en.sig_pNS <- matrix(NA, nosamples, PeriodLength)
tot_en.sig_sN_Y <- tot_en.sig_sN_A <- matrix(NA, nosamples, PeriodLength)
tot_en.sig_imm_Y <- tot_en.sig_imm_A <- matrix(NA, nosamples, PeriodLength)

## Define identity matrix and vector of ones
I <- diag(S)
e <- matrix(1, S, 1)

## Calculating indirect effect real-time elasticities
for(i in 1:nosamples){
  
  # a) Set up arrays or storing perturbation matrices (C2 and C3) and population
  #    structure perturbations (w)
  #------------------------------------------------------------------------
  
  ## Prepare arrays for perturbation matrices (C2) for direct effects
  C2_sJ <- C2_sA <- array(0, dim = c(S^2, S^2, PeriodLength))
  C2_pB_Y <- C2_pB_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  C2_CS_Y <- C2_CS_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  C2_pNS <- array(0, dim = c(S^2, S^2, PeriodLength))
  C2_sN_Y <- C2_sN_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  C2_imm_Y <- C2_imm_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  
  ## Prepare arrays for perturbation matrices (C3) for indirect effects
  C3_sJ <- C3_sA <- array(0, dim = c(S^2, S^2, PeriodLength))
  C3_pB_Y <- C3_pB_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  C3_CS_Y <- C3_CS_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  C3_pNS <- array(0, dim = c(S^2, S^2, PeriodLength))
  C3_sN_Y <- C3_sN_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  C3_imm_Y <- C3_imm_A <- array(0, dim = c(S^2, S^2, PeriodLength))
  
  ## Prepare arrays for population structure perturbations (w) due to changes
  #  in mean VR
  w.mu_sJ <- w.mu_sA <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.mu_pB_Y <- w.mu_pB_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.mu_CS_Y <- w.mu_CS_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.mu_pNS <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.mu_sN_Y <- w.mu_sN_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.mu_imm_Y <- w.mu_imm_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  
  ## Prepare arrays for population structure perturbations (w) due to changes
  #  in VR standard deviation
  w.sig_sJ <- w.sig_sA <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.sig_pB_Y <- w.sig_pB_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.sig_CS_Y <- w.sig_CS_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.sig_pNS <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.sig_sN_Y <- w.sig_sN_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  w.sig_imm_Y <- w.sig_imm_A <- array(0, dim = c(S, S, S, PeriodLength+1))
  
  for(t in 1:PeriodLength){
    
    # b) Build the reference matrix for time step t -> t+1 and make derivative 
    #    of projection matrix with regards to each vital rate (dA[t]/dVR[t])
    #------------------------------------------------------------------------
    
    ## Build the reference projection matrix
    #  NOTE: This is the same matrix as was used above for projecting the 
    #        reference population.
    
    A <- matrix(NA, nrow = S, ncol = S)
    A[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    A[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    A[2,1] <- ref_sA[i,t]
    A[2,2] <- ref_sA[i,t]
    
    ## Calculate derivatives of reference projection matrix
    # sJ
    dA_dsJ <- matrix(0, nrow = S, ncol = S)
    dA_dsJ[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]
    dA_dsJ[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]
    
    # sA
    dA_dsA <- matrix(c(0, 1, 0, 1), nrow = S)
    
    # pB_Y
    dA_dpB_Y <- matrix(0, nrow = S, ncol = S)
    dA_dpB_Y[1,1] <- 0.5*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    
    # pB_A
    dA_dpB_A <- matrix(0, nrow = S, ncol = S)
    dA_dpB_A[1,2] <- 0.5*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    
    # CS_Y
    dA_dCS_Y <- matrix(0, nrow = S, ncol = S)
    dA_dCS_Y[1,1] <- 0.5*ref_pB_Y[i,t]*ref_pNS[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    
    # CS_A
    dA_dCS_A <- matrix(0, nrow = S, ncol = S)
    dA_dCS_A[1,2] <- 0.5*ref_pB_A[i,t]*ref_pNS[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    
    # pNS
    dA_dpNS <- matrix(0, nrow = S, ncol = S)
    dA_dpNS[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_sN_Y[i,t]*ref_sJ[i,t]
    dA_dpNS[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_sN_A[i,t]*ref_sJ[i,t]
    
    # sN_Y
    dA_dsN_Y <- matrix(0, nrow = S, ncol = S)
    dA_dsN_Y[1,1] <- 0.5*ref_pB_Y[i,t]*ref_CS_Y[i,t]*ref_pNS[i,t]*ref_sJ[i,t]
    
    # sN_A
    dA_dsN_A <- matrix(0, nrow = S, ncol = S)
    dA_dsN_A[1,2] <- 0.5*ref_pB_A[i,t]*ref_CS_A[i,t]*ref_pNS[i,t]*ref_sJ[i,t]
    
    
    # c) Calculate real-time elasticities to changes in population structure
    #    resulting from changes in the mean or standard deviation of vital 
    #    rates, using equations S1.7 and S1.8 and the step-wise calculation 
    #    approach from from Koons et al. 2016
    # -----------------------------------------------------------------------
    
    for(m in 1:S){
      for(n in 1:S){
        
        # Step 1) Calculation of the perturbation matrices for time step t
        # NOTE: Into these matrices goes the product of the reference mean/sd VR 
        #       and the sensitivity of matrix element A[m,n] to perturbation in the VR

        C2_sJ[(S+1)*m-S, (S+1)*n-S, t] <- ref_sJ.mu[i]*dA_dsJ[m,n]
        C2_sA[(S+1)*m-S, (S+1)*n-S, t] <- ref_sA.mu[i]*dA_dsA[m,n]
        C2_pB_Y[(S+1)*m-S, (S+1)*n-S, t] <- ref_pB_Y.mu[i]*dA_dpB_Y[m,n]
        C2_pB_A[(S+1)*m-S, (S+1)*n-S, t] <- ref_pB_A.mu[i]*dA_dpB_A[m,n]
        C2_CS_Y[(S+1)*m-S, (S+1)*n-S, t] <- ref_CS_Y.mu[i]*dA_dCS_Y[m,n]
        C2_CS_A[(S+1)*m-S, (S+1)*n-S, t] <- ref_CS_A.mu[i]*dA_dCS_A[m,n]
        C2_pNS[(S+1)*m-S, (S+1)*n-S, t] <- ref_pNS.mu[i]*dA_dpNS[m,n]
        C2_sN_Y[(S+1)*m-S, (S+1)*n-S, t] <- ref_sN_Y.mu[i]*dA_dsN_Y[m,n]
        C2_sN_A[(S+1)*m-S, (S+1)*n-S, t] <- ref_sN_A.mu[i]*dA_dsN_A[m,n]
        
        C3_sJ[(S+1)*m-S, (S+1)*n-S, t] <- (ref_sJ[i,t]-ref_sJ.mu[i])*dA_dsJ[m,n]
        C3_sA[(S+1)*m-S, (S+1)*n-S, t] <- (ref_sA[i,t]-ref_sA.mu[i])*dA_dsA[m,n]
        C3_pB_Y[(S+1)*m-S, (S+1)*n-S, t] <- (ref_pB_Y[i,t]-ref_pB_Y.mu[i])*dA_dpB_Y[m,n]
        C3_pB_A[(S+1)*m-S, (S+1)*n-S, t] <- (ref_pB_A[i,t]-ref_pB_A.mu[i])*dA_dpB_A[m,n]
        C3_CS_Y[(S+1)*m-S, (S+1)*n-S, t] <- (ref_CS_Y[i,t]-ref_CS_Y.mu[i])
        C3_CS_A[(S+1)*m-S, (S+1)*n-S, t] <- (ref_CS_A[i,t]-ref_CS_A.mu[i])
        C3_pNS[(S+1)*m-S, (S+1)*n-S, t] <- (ref_pNS[i,t]-ref_pNS.mu[i])*dA_dpNS[m,n]
        C3_sN_Y[(S+1)*m-S, (S+1)*n-S, t] <- (ref_sN_Y[i,t]-ref_sN_Y.mu[i])*dA_dsN_Y[m,n]
        C3_sN_A[(S+1)*m-S, (S+1)*n-S, t] <- (ref_sN_A[i,t]-ref_sN_A.mu[i])*dA_dsN_A[m,n]
        
        # NOTE: Immigration rates (imm_X) do not get a perturbation matrix
        # because they do not appear within the projection matrix, and are
        # independent of population structure
        
        
        # Step 2.1) Intermediate steps towards solving equation S1.7
        # NOTE: B is one part of w (equation S1.7): (I-n[t+1]e')A[t]/lambda[t]
        
        K <- I - ref_n[,1,i,t+1] %*% t(e)
        B <- K %*% A / lam_realref[i,t]
        
        
        # Step 2.2) Intermediate steps towards solving equation S1.7
        # NOTE: g is another part of w (equation S1.7): (I-n[t+1]e')C[t]n[t]/lambda[t]
        
        # Matrix-internal VRs (mean)
        g.mu_sJ <- K %*% C2_sJ[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_sA <- K %*% C2_sA[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_pB_Y <- K %*% C2_pB_Y[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_pB_A <- K %*% C2_pB_A[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_CS_Y <- K %*% C2_CS_Y[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_CS_A <- K %*% C2_CS_A[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_pNS <- K %*% C2_pNS[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_sN_Y <- K %*% C2_sN_Y[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.mu_sN_A <- K %*% C2_sN_A[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        
        # Matrix-internal VRs (sd)
        g.sig_sJ <- K %*% C3_sJ[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_sA <- K %*% C3_sA[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_pB_Y <- K %*% C3_pB_Y[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_pB_A <- K %*% C3_pB_A[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_CS_Y <- K %*% C3_CS_Y[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_CS_A <- K %*% C3_CS_A[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_pNS <- K %*% C3_pNS[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_sN_Y <- K %*% C3_sN_Y[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        g.sig_sN_A <- K %*% C3_sN_A[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),t] %*% ref_n[,1,i,t] / lam_realref[i,t]
        
        # Immigration rates (mean, sd)
        # NOTE: Deriving an equivalent of equation S1.7 for immigration rates 
        # following Haridas et al. 2009 shows that the term C[t]n[t] is to be
        # replaced with imm_X.mu (mean) or imm_X[t]-imm.mu (sd)
        g.mu_imm_Y <- K %*% c(ref_imm_Y.mu[i],0) / lam_realref[i,t]
        g.mu_imm_A <- K %*% c(0,ref_imm_A.mu[i]) / lam_realref[i,t]
        g.sig_imm_Y <- K %*% c((ref_imm[1,1,i,t+1]-ref_imm_Y.mu[i]),0) / lam_realref[i,t]
        g.sig_imm_A <- K %*% c(0,(ref_imm[2,1,i,t+1]-ref_imm_A.mu[i])) / lam_realref[i,t]
      
        
        # Step 2.3) Completion of equation S1.7 and calculation of w at the next time step
        # B[t]: (I-n[t+1]e')A[t]/lambda[t]
        # g[t]: (I-n[t+1]e')C[t]n[t]/lambda[t]
        # --> (I-n[t+1]e')A[t]w[t]/lambda[t] + (I-n[t+1]e')C[t]n[t]/lambda[t]
        # --> (I-n[t+1]e')(C[t]n[t] + A[t]w[t])/lambda[t], which is equation S1.7
        
        w.mu_sJ[,m,n,t+1] <- B %*% w.mu_sJ[,m,n,t] + g.mu_sJ
        w.mu_sA[,m,n,t+1] <- B %*% w.mu_sA[,m,n,t] + g.mu_sA
        w.mu_pB_Y[,m,n,t+1] <- B %*% w.mu_pB_Y[,m,n,t] + g.mu_pB_Y
        w.mu_pB_A[,m,n,t+1] <- B %*% w.mu_pB_A[,m,n,t] + g.mu_pB_A
        w.mu_CS_Y[,m,n,t+1] <- B %*% w.mu_CS_Y[,m,n,t] + g.mu_CS_Y
        w.mu_CS_A[,m,n,t+1] <- B %*% w.mu_CS_A[,m,n,t] + g.mu_CS_A
        w.mu_pNS[,m,n,t+1] <- B %*% w.mu_pNS[,m,n,t] + g.mu_pNS
        w.mu_sN_Y[,m,n,t+1] <- B %*% w.mu_sN_Y[,m,n,t] + g.mu_sN_Y
        w.mu_sN_A[,m,n,t+1] <- B %*% w.mu_sN_A[,m,n,t] + g.mu_sN_A
        
        w.mu_imm_Y[,m,n,t+1] <- B %*% w.mu_imm_Y[,m,n,t] + g.mu_imm_Y
        w.mu_imm_A[,m,n,t+1] <- B %*% w.mu_imm_A[,m,n,t] + g.mu_imm_A
        
        w.sig_sJ[,m,n,t+1] <- B %*% w.sig_sJ[,m,n,t] + g.sig_sJ
        w.sig_sA[,m,n,t+1] <- B %*% w.sig_sA[,m,n,t] + g.sig_sA
        w.sig_pB_Y[,m,n,t+1] <- B %*% w.sig_pB_Y[,m,n,t] + g.sig_pB_Y
        w.sig_pB_A[,m,n,t+1] <- B %*% w.sig_pB_A[,m,n,t] + g.sig_pB_A
        w.sig_CS_Y[,m,n,t+1] <- B %*% w.sig_CS_Y[,m,n,t] + g.sig_CS_Y
        w.sig_CS_A[,m,n,t+1] <- B %*% w.sig_CS_A[,m,n,t] + g.sig_CS_A
        w.sig_pNS[,m,n,t+1] <- B %*% w.sig_pNS[,m,n,t] + g.sig_pNS
        w.sig_sN_Y[,m,n,t+1] <- B %*% w.sig_sN_Y[,m,n,t] + g.sig_sN_Y
        w.sig_sN_A[,m,n,t+1] <- B %*% w.sig_sN_A[,m,n,t] + g.sig_sN_A
        
        w.sig_imm_Y[,m,n,t+1] <- B %*% w.sig_imm_Y[,m,n,t] + g.sig_imm_Y
        w.sig_imm_A[,m,n,t+1] <- B %*% w.sig_imm_A[,m,n,t] + g.sig_imm_A
        
        
        # Step 3) Calculation of indirect real-time elasticities (of the matrix 
        # element [m,n]) using equation S1.8 (= product of the entire matrix and
        # w divided by lambda)
        # NOTE: The multiplication with t(e) here makes it so that a single 
        # number (and not a vector of length S) is returned

        en.mu_sJ[m,n,i,t] <- t(e) %*% A %*% w.mu_sJ[,m,n,t] / lam_realref[i,t]
        en.mu_sA[m,n,i,t] <- t(e) %*% A %*% w.mu_sA[,m,n,t] / lam_realref[i,t]
        en.mu_pB_Y[m,n,i,t] <- t(e) %*% A %*% w.mu_pB_Y[,m,n,t] / lam_realref[i,t]
        en.mu_pB_A[m,n,i,t] <- t(e) %*% A %*% w.mu_pB_A[,m,n,t] / lam_realref[i,t]
        en.mu_CS_Y[m,n,i,t] <- t(e) %*% A %*% w.mu_CS_Y[,m,n,t] / lam_realref[i,t]
        en.mu_CS_A[m,n,i,t] <- t(e) %*% A %*% w.mu_CS_A[,m,n,t] / lam_realref[i,t]
        en.mu_pNS[m,n,i,t] <- t(e) %*% A %*% w.mu_pNS[,m,n,t] / lam_realref[i,t]
        en.mu_sN_Y[m,n,i,t] <- t(e) %*% A %*% w.mu_sN_Y[,m,n,t] / lam_realref[i,t]
        en.mu_sN_A[m,n,i,t] <- t(e) %*% A %*% w.mu_sN_A[,m,n,t] / lam_realref[i,t]
        
        en.mu_imm_Y[m,n,i,t] <- t(e) %*% A %*% w.mu_imm_Y[,m,n,t] / lam_realref[i,t]
        en.mu_imm_A[m,n,i,t] <- t(e) %*% A %*% w.mu_imm_A[,m,n,t] / lam_realref[i,t]
        
        en.sig_sJ[m,n,i,t] <- t(e) %*% A %*% w.sig_sJ[,m,n,t] / lam_realref[i,t]
        en.sig_sA[m,n,i,t] <- t(e) %*% A %*% w.sig_sA[,m,n,t] / lam_realref[i,t]
        en.sig_pB_Y[m,n,i,t] <- t(e) %*% A %*% w.sig_pB_Y[,m,n,t] / lam_realref[i,t]
        en.sig_pB_A[m,n,i,t] <- t(e) %*% A %*% w.sig_pB_A[,m,n,t] / lam_realref[i,t]
        en.sig_CS_Y[m,n,i,t] <- t(e) %*% A %*% w.sig_CS_Y[,m,n,t] / lam_realref[i,t]
        en.sig_CS_A[m,n,i,t] <- t(e) %*% A %*% w.sig_CS_A[,m,n,t] / lam_realref[i,t]
        en.sig_pNS[m,n,i,t] <- t(e) %*% A %*% w.sig_pNS[,m,n,t] / lam_realref[i,t]
        en.sig_sN_Y[m,n,i,t] <- t(e) %*% A %*% w.sig_sN_Y[,m,n,t] / lam_realref[i,t]
        en.sig_sN_A[m,n,i,t] <- t(e) %*% A %*% w.sig_sN_A[,m,n,t] / lam_realref[i,t]
        
        en.sig_imm_Y[m,n,i,t] <- t(e) %*% A %*% w.sig_imm_Y[,m,n,t] / lam_realref[i,t]
        en.sig_imm_A[m,n,i,t] <- t(e) %*% A %*% w.sig_imm_A[,m,n,t] / lam_realref[i,t]
      }
    }
    
    
    # d) Sum up elasticities to vital rates through different matrix elements
    # -----------------------------------------------------------------------
    
    tot_en.mu_sJ[i,t] <- sum(en.mu_sJ[,,i,t])
    tot_en.mu_sA[i,t] <- sum(en.mu_sA[,,i,t])
    tot_en.mu_pB_Y[i,t] <- sum(en.mu_pB_Y[,,i,t])
    tot_en.mu_pB_A[i,t] <- sum(en.mu_pB_A[,,i,t])
    tot_en.mu_CS_Y[i,t] <- sum(en.mu_CS_Y[,,i,t])
    tot_en.mu_CS_A[i,t] <- sum(en.mu_CS_A[,,i,t])
    tot_en.mu_pNS[i,t] <- sum(en.mu_pNS[,,i,t])
    tot_en.mu_sN_Y[i,t] <- sum(en.mu_sN_Y[,,i,t])
    tot_en.mu_sN_A[i,t] <- sum(en.mu_sN_A[,,i,t])
    
    tot_en.mu_imm_Y[i,t] <- sum(en.mu_imm_Y[,,i,t])
    tot_en.mu_imm_A[i,t] <- sum(en.mu_imm_A[,,i,t])
    
    tot_en.sig_sJ[i,t] <- sum(en.sig_sJ[,,i,t])
    tot_en.sig_sA[i,t] <- sum(en.sig_sA[,,i,t])
    tot_en.sig_pB_Y[i,t] <- sum(en.sig_pB_Y[,,i,t])
    tot_en.sig_pB_A[i,t] <- sum(en.sig_pB_A[,,i,t])
    tot_en.sig_CS_Y[i,t] <- sum(en.sig_CS_Y[,,i,t])
    tot_en.sig_CS_A[i,t] <- sum(en.sig_CS_A[,,i,t])
    tot_en.sig_pNS[i,t] <- sum(en.sig_pNS[,,i,t])
    tot_en.sig_sN_Y[i,t] <- sum(en.sig_sN_Y[,,i,t])
    tot_en.sig_sN_A[i,t] <- sum(en.sig_sN_A[,,i,t])
    
    tot_en.sig_imm_Y[i,t] <- sum(en.sig_imm_Y[,,i,t])
    tot_en.sig_imm_A[i,t] <- sum(en.sig_imm_A[,,i,t])
    
  }
}

## Summarise average elasticities across years
avg_en.mu_sJ <- rowMeans(tot_en.mu_sJ)
avg_en.mu_sA <- rowMeans(tot_en.mu_sA)
avg_en.mu_pB_Y <- rowMeans(tot_en.mu_pB_Y)
avg_en.mu_pB_A <- rowMeans(tot_en.mu_pB_A)
avg_en.mu_CS_Y <- rowMeans(tot_en.mu_CS_Y)
avg_en.mu_CS_A <- rowMeans(tot_en.mu_CS_A)
avg_en.mu_pNS <- rowMeans(tot_en.mu_pNS)
avg_en.mu_sN_Y <- rowMeans(tot_en.mu_sN_Y)
avg_en.mu_sN_A <- rowMeans(tot_en.mu_sN_A)

avg_en.mu_imm_Y <- rowMeans(tot_en.mu_imm_Y)
avg_en.mu_imm_A <- rowMeans(tot_en.mu_imm_A)

avg_en.sig_sJ <- rowMeans(tot_en.sig_sJ)
avg_en.sig_sA <- rowMeans(tot_en.sig_sA)
avg_en.sig_pB_Y <- rowMeans(tot_en.sig_pB_Y)
avg_en.sig_pB_A <- rowMeans(tot_en.sig_pB_A)
avg_en.sig_CS_Y <- rowMeans(tot_en.sig_CS_Y)
avg_en.sig_CS_A <- rowMeans(tot_en.sig_CS_A)
avg_en.sig_pNS <- rowMeans(tot_en.sig_pNS)
avg_en.sig_sN_Y <- rowMeans(tot_en.sig_sN_Y)
avg_en.sig_sN_A <- rowMeans(tot_en.sig_sN_A)

avg_en.sig_imm_Y <- rowMeans(tot_en.sig_imm_Y)
avg_en.sig_imm_A <- rowMeans(tot_en.sig_imm_A)


##################################################
#### CALCULATION OF PERIOD LTRE CONTRIBUTIONS ####
##################################################

message(cyan('Calculating period-design LTRE contributions...'))

## Contributions from direct changes in mean VRs
contA.mu_sJ <- logmudiff_sJ*avg_eA.mu_sJ
contA.mu_sA <- logmudiff_sA*avg_eA.mu_sA
contA.mu_pB_Y <- logmudiff_pB_Y*avg_eA.mu_pB_Y
contA.mu_pB_A <- logmudiff_pB_A*avg_eA.mu_pB_A
contA.mu_CS_Y <- logmudiff_CS_Y*avg_eA.mu_CS_Y
contA.mu_CS_A <- logmudiff_CS_A*avg_eA.mu_CS_A
contA.mu_pNS <- logmudiff_pNS*avg_eA.mu_pNS
contA.mu_sN_Y <- logmudiff_sN_Y*avg_eA.mu_sN_Y
contA.mu_sN_A <- logmudiff_sN_A*avg_eA.mu_sN_A
contA.mu_imm_Y <- logmudiff_imm_Y*avg_eA.mu_imm_Y
contA.mu_imm_A <- logmudiff_imm_A*avg_eA.mu_imm_A

## Contributions from indirect changes in mean VRs
contn.mu_sJ <- logmudiff_sJ*avg_en.mu_sJ
contn.mu_sA <- logmudiff_sA*avg_en.mu_sA
contn.mu_pB_Y <- logmudiff_pB_Y*avg_en.mu_pB_Y
contn.mu_pB_A <- logmudiff_pB_A*avg_en.mu_pB_A
contn.mu_CS_Y <- logmudiff_CS_Y*avg_en.mu_CS_Y
contn.mu_CS_A <- logmudiff_CS_A*avg_en.mu_CS_A
contn.mu_pNS <- logmudiff_pNS*avg_en.mu_pNS
contn.mu_sN_Y <- logmudiff_sN_Y*avg_en.mu_sN_Y
contn.mu_sN_A <- logmudiff_sN_A*avg_en.mu_sN_A
contn.mu_imm_Y <- logmudiff_imm_Y*avg_en.mu_imm_Y
contn.mu_imm_A <- logmudiff_imm_A*avg_en.mu_imm_A

## Contributions from direct changes in VR variability
contA.sig_sJ <- logsigdiff_sJ*avg_eA.sig_sJ
contA.sig_sA <- logsigdiff_sA*avg_eA.sig_sA
contA.sig_pB_Y <- logsigdiff_pB_Y*avg_eA.sig_pB_Y
contA.sig_pB_A <- logsigdiff_pB_A*avg_eA.sig_pB_A
contA.sig_CS_Y <- logsigdiff_CS_Y*avg_eA.sig_CS_Y
contA.sig_CS_A <- logsigdiff_CS_A*avg_eA.sig_CS_A
contA.sig_pNS <- logsigdiff_pNS*avg_eA.sig_pNS
contA.sig_sN_Y <- logsigdiff_sN_Y*avg_eA.sig_sN_Y
contA.sig_sN_A <- logsigdiff_sN_A*avg_eA.sig_sN_A
contA.sig_imm_Y <- logsigdiff_imm_Y*avg_eA.sig_imm_Y
contA.sig_imm_A <- logsigdiff_imm_A*avg_eA.sig_imm_A

## Contributions from indirect changes in VR variability
contn.sig_sJ <- logsigdiff_sJ*avg_en.sig_sJ
contn.sig_sA <- logsigdiff_sA*avg_en.sig_sA
contn.sig_pB_Y <- logsigdiff_pB_Y*avg_en.sig_pB_Y
contn.sig_pB_A <- logsigdiff_pB_A*avg_en.sig_pB_A
contn.sig_CS_Y <- logsigdiff_CS_Y*avg_en.sig_CS_Y
contn.sig_CS_A <- logsigdiff_CS_A*avg_en.sig_CS_A
contn.sig_pNS <- logsigdiff_pNS*avg_en.sig_pNS
contn.sig_sN_Y <- logsigdiff_sN_Y*avg_en.sig_sN_Y
contn.sig_sN_A <- logsigdiff_sN_A*avg_en.sig_sN_A
contn.sig_imm_Y <- logsigdiff_imm_Y*avg_en.sig_imm_Y
contn.sig_imm_A <- logsigdiff_imm_A*avg_en.sig_imm_A

#####################################
#### ASSEMBLING & SAVING RESULTS ####
#####################################

message(cyan('Assembling & saving results...'))

## Assemble contributions from changes in VR means in a list
ContVecs_mu <- list(
  contA_sJ = contA.mu_sJ, contA_sA = contA.mu_sA,
  contA_pB_Y = contA.mu_pB_Y, contA_pB_A = contA.mu_pB_A,
  contA_CS_Y = contA.mu_CS_Y, contA_CS_A = contA.mu_CS_A,
  contA_pNS = contA.mu_pNS,
  contA_sN_Y = contA.mu_sN_Y, contA_sN_A = contA.mu_sN_A,
  contA_imm_Y = contA.mu_imm_Y, contA_imm_A = contA.mu_imm_A,
  
  contn_sJ = contn.mu_sJ, contn_sA = contn.mu_sA,
  contn_pB_Y = contn.mu_pB_Y, contn_pB_A = contn.mu_pB_A,
  contn_CS_Y = contn.mu_CS_Y, contn_CS_A = contn.mu_CS_A,
  contn_pNS = contn.mu_pNS,
  contn_sN_Y = contn.mu_sN_Y, contn_sN_A = contn.mu_sN_A,
  contn_imm_Y = contn.mu_imm_Y, contn_imm_A = contn.mu_imm_A
)

## Assemble contributions from changes in VR variation in a list
ContVecs_sig <- list(
  contA_sJ = contA.sig_sJ, contA_sA = contA.sig_sA,
  contA_pB_Y = contA.sig_pB_Y, contA_pB_A = contA.sig_pB_A,
  contA_CS_Y = contA.sig_CS_Y, contA_CS_A = contA.sig_CS_A,
  contA_pNS = contA.sig_pNS,
  contA_sN_Y = contA.sig_sN_Y, contA_sN_A = contA.sig_sN_A,
  contA_imm_Y = contA.sig_imm_Y, contA_imm_A = contA.sig_imm_A,
  
  contn_sJ = contn.sig_sJ, contn_sA = contn.sig_sA,
  contn_pB_Y = contn.sig_pB_Y, contn_pB_A = contn.sig_pB_A,
  contn_CS_Y = contn.sig_CS_Y, contn_CS_A = contn.sig_CS_A,
  contn_pNS = contn.sig_pNS,
  contn_sN_Y = contn.sig_sN_Y, contn_sN_A = contn.sig_sN_A,
  contn_imm_Y = contn.sig_imm_Y, contn_imm_A = contn.sig_imm_A
)

## Assemble total contributions from VR changes in a list
ContVecs_tot <- list(
  cont_sJ = contA.mu_sJ + contn.mu_sJ + contA.sig_sJ + contn.sig_sJ, 
  cont_sA = contA.mu_sA + contn.mu_sA + contA.sig_sA + contn.sig_sA,
  cont_pB_Y = contA.mu_pB_Y + contn.mu_pB_Y + contA.sig_pB_Y + contn.sig_pB_Y, 
  cont_pB_A = contA.mu_pB_A + contn.mu_pB_A + contA.sig_pB_A + contn.sig_pB_A,
  cont_CS_Y = contA.mu_CS_Y + contn.mu_CS_Y + contA.sig_CS_Y + contn.sig_CS_Y, 
  cont_CS_A = contA.mu_CS_A + contn.mu_CS_A + contA.sig_CS_A + contn.sig_CS_A, 
  cont_pNS = contA.mu_pNS + contn.mu_pNS + contA.sig_pNS + contn.sig_pNS,
  cont_sN_Y = contA.mu_sN_Y + contn.mu_sN_Y + contA.sig_sN_Y + contn.sig_sN_Y,  
  cont_sN_A = contA.mu_sN_A + contn.mu_sN_A + contA.sig_sN_A + contn.sig_sN_A, 
  cont_imm_Y = contA.mu_imm_Y + contn.mu_imm_Y + contA.sig_imm_Y + contn.sig_imm_Y, 
  cont_imm_A = contA.mu_imm_A + contn.mu_imm_A + contA.sig_imm_A + contn.sig_imm_A
)

## Collate contributions (all samples)
ContData <- data.frame(
  PopID = PopID,
  sample = rep(1:nosamples, 11),
  cont = c(ContVecs_tot$cont_sJ, ContVecs_tot$cont_sA, 
           ContVecs_tot$cont_pB_Y, ContVecs_tot$cont_pB_A,
           ContVecs_tot$cont_CS_Y, ContVecs_tot$cont_CS_A, 
           ContVecs_tot$cont_pNS,
           ContVecs_tot$cont_sN_Y, ContVecs_tot$cont_sN_A, 
           ContVecs_tot$cont_imm_Y, ContVecs_tot$cont_imm_A),
  parameter = rep(c('sJ', 'sA', 'pB_Y', 'pB_A',
                    'CS_Y', 'CS_A', 'pNS',
                    'sN_Y', 'sN_A',
                    'imm_Y', 'imm_A'), each = nosamples)
)

## Assembling results in a list
LTRE_Results <- list(
  PopID = PopID,
  Period1 = Period1+min(StudyYears),
  Period2 = Period1+min(StudyYears),
  ContVecs_mu = ContVecs_mu,
  ContVecs_sig = ContVecs_sig,
  ContVecs_tot = ContVecs_tot,
  ContData = ContData
)

## Saving output
if(!FixLog0MeanRates){
  saveRDS(LTRE_Results, file = paste0('periodLTRE_', PopID, '.rds'))
}
if(FixLog0MeanRates){
  saveRDS(LTRE_Results, file = paste0('periodLTRE_', PopID, '_Log0MeanFix.rds'))
}

