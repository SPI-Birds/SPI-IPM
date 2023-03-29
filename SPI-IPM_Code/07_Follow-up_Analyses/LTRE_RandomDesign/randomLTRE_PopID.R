#######################
#### DATA ASSEMBLY ####
#######################

message(crayon::green$underline$bold(paste0('Random design LTRE analyses for ', PopID)))

## Extract study years for focal population
StudyYears <- eval(parse(text = paste0('StudyYearsList$', PopID)))

## Load MCMC samples
PFC.IPM <- readRDS(paste0(DataPath,'/SPI-IPM_', PopID, '.rds'))
#str(PFC.IPM) # = posterior samples in coda format

## Transform samples into matrix (all chains combined)
out.mat <- as.matrix(PFC.IPM)
#str(out.mat) # = posterior samples in matrix format

## Set sample and year number
nosamples <- dim(out.mat)[1]
noyears <- length(StudyYears)

  
###############
#### SETUP ####
###############

message(crayon::cyan('Assembling posterior data...'))

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
for(i in 1:nosamples){
  
  for(t in 1:noyears){
    
    # Time-varying population sizes
    N_Y[i,t] <- out.mat[i, paste0('N[1, ', t, ']')]
    N_A[i,t] <- out.mat[i, paste0('N[2, ', t, ']')]
    N_tot[i,t] <- out.mat[i, paste0('Ntot[', t, ']')]
    
    # Time-varying immigrant numbers
    Imm_Y[i,t] <- out.mat[i, paste0('Imm[1, ', t, ']')]
    Imm_A[i,t] <- out.mat[i, paste0('Imm[2, ', t, ']')]
  }
  
  for(t in 1:(noyears-1)){
    
    # Time-varying vital rates   
    pB_Y[i,t] <- out.mat[i, paste0('pB[1, ', t, ']')]
    pB_A[i,t] <- out.mat[i, paste0('pB[2, ', t, ']')]
    
    CS_Y[i,t] <- out.mat[i, paste0('CS[1, ', t, ']')]
    CS_A[i,t] <- out.mat[i, paste0('CS[2, ', t, ']')]
    
    pNS[i,t] <- out.mat[i, paste0('pNS[', t, ']')]
    
    sN_Y[i,t] <- out.mat[i, paste0('sN[1, ', t, ']')]
    sN_A[i,t] <- out.mat[i, paste0('sN[2, ', t, ']')]
    
    sJ[i,t] <- out.mat[i, paste0('sJ[', t, ']')]
    sA[i,t] <- out.mat[i, paste0('sA[', t, ']')]

    # Population growth rate
    lambda[i,t] <- N_tot[i,t+1]/N_tot[i,t]
  }
}

## Calculate population/immigrant proportions
n_Y <- N_Y/N_tot
n_A <- N_A/N_tot
imm_Y <- cbind(NA,Imm_Y[,2:noyears]/N_tot[,1:(noyears-1)])
imm_A <- cbind(NA,Imm_A[,2:noyears]/N_tot[,1:(noyears-1)])

## Make time-average population sizes/proportions
N_Y_mean <- rowMeans(N_Y[,1:(noyears-1)])
N_A_mean <- rowMeans(N_A[,1:(noyears-1)])
N_tot_mean <- rowMeans(N_tot[,1:(noyears-1)])

n_Y_mean <- rowMeans(n_Y[,1:(noyears-1)])
n_A_mean <- rowMeans(n_A[,1:(noyears-1)])

## Make average immigrant numbers
Imm_Y_mean <- rowMeans(Imm_Y[,2:noyears])
Imm_A_mean <- rowMeans(Imm_A[,2:noyears])

imm_Y_mean <- rowMeans(imm_Y[,2:noyears])
imm_A_mean <- rowMeans(imm_A[,2:noyears])

## Make time-average vital rates
sJ_mean <- rowMeans(sJ[,1:(noyears-1)])
sA_mean <- rowMeans(sA[,1:(noyears-1)])

pB_Y_mean <- rowMeans(pB_Y[,1:(noyears-1)])
pB_A_mean <- rowMeans(pB_A[,1:(noyears-1)])

CS_Y_mean <- rowMeans(CS_Y[,1:(noyears-1)])
CS_A_mean <- rowMeans(CS_A[,1:(noyears-1)])

pNS_mean <- rowMeans(pNS[,1:(noyears-1)])

sN_Y_mean <- rowMeans(sN_Y[,1:(noyears-1)])
sN_A_mean <- rowMeans(sN_A[,1:(noyears-1)])

## Make time-average population growth rate
lambda_mean <- rowMeans(lambda[,1:(noyears-1)])

## Calculate time-variation in population growth rate
lambda_var <- rowVars(lambda[,1:(noyears-1)])

## Identify unusable samples (lambda = Inf)
NAsamples <- which(lambda_mean == Inf)

################################################
#### CALCULATION OF TRANSIENT SENSITIVITIES ####
################################################

message(crayon::cyan('Calculating sensitivities & elasticities...'))

## Make vectors for storing transient sensitivities

sens_sJ <- rep(NA, nosamples)
sens_sA <- rep(NA, nosamples)

sens_pB_Y <- rep(NA, nosamples)
sens_pB_A <- rep(NA, nosamples)

sens_CS_Y <- rep(NA, nosamples)
sens_CS_A <- rep(NA, nosamples)

sens_pNS <- rep(NA, nosamples)

sens_sN_Y <- rep(NA, nosamples)
sens_sN_A <- rep(NA, nosamples)

sens_N_Y <- rep(NA, nosamples)
sens_N_A <- rep(NA, nosamples)
sens_n_Y <- rep(NA, nosamples)
sens_n_A <- rep(NA, nosamples)

sens_Imm_Y <- rep(NA, nosamples)
sens_Imm_A <- rep(NA, nosamples)
sens_imm_Y <- rep(NA, nosamples)
sens_imm_A <- rep(NA, nosamples)


for(i in 1:nosamples){
  
  ## Calculate transient sensitivities for vital rates
  sens_sJ[i] <- 0.5*pB_Y_mean[i]*CS_Y_mean[i]*pNS_mean[i]*sN_Y_mean[i]*n_Y_mean[i] +
    0.5*pB_A_mean[i]*CS_A_mean[i]*pNS_mean[i]*sN_A_mean[i]*n_A_mean[i]
  
  sens_sA[i] <- n_Y_mean[i] + n_A_mean[i]
  
  sens_pB_Y[i] <- 0.5*CS_Y_mean[i]*pNS_mean[i]*sN_Y_mean[i]*sJ_mean[i]*n_Y_mean[i]
  sens_pB_A[i] <- 0.5*CS_A_mean[i]*pNS_mean[i]*sN_A_mean[i]*sJ_mean[i]*n_A_mean[i]
  
  sens_CS_Y[i] <- 0.5*pB_Y_mean[i]*pNS_mean[i]*sN_Y_mean[i]*sJ_mean[i]*n_Y_mean[i]
  sens_CS_A[i] <- 0.5*pB_A_mean[i]*pNS_mean[i]*sN_A_mean[i]*sJ_mean[i]*n_A_mean[i]
  
  sens_pNS[i] <- 0.5*sJ_mean[i]*(pB_Y_mean[i]*CS_Y_mean[i]*sN_Y_mean[i]*n_Y_mean[i] + 
                                    pB_A_mean[i]*CS_A_mean[i]*sN_A_mean[i]*n_A_mean[i])
  
  sens_sN_Y[i] <- 0.5*pB_Y_mean[i]*CS_Y_mean[i]*pNS_mean[i]*sJ_mean[i]*n_Y_mean[i]
  sens_sN_A[i] <- 0.5*pB_A_mean[i]*CS_A_mean[i]*pNS_mean[i]*sJ_mean[i]*n_A_mean[i]
  
  
  ## Calculate transient sensitivities for population structure
	sens_n_Y[i] <- 0.5*pNS_mean[i]*sJ_mean[i]*((pB_Y_mean[i]*CS_Y_mean[i]*sN_Y_mean[i]) - 
	                                           (pB_A_mean[i]*CS_A_mean[i]*sN_A_mean[i]))
	sens_n_A[i] <- 0.5*pNS_mean[i]*sJ_mean[i]*((pB_A_mean[i]*CS_A_mean[i]*sN_A_mean[i]) - 
	                                           (pB_Y_mean[i]*CS_Y_mean[i]*sN_Y_mean[i]))
	
  sens_imm_Y[i] <- sens_imm_A[i] <- 1
 
}

## Set sensitivities to NA for unusable samples
sens_sJ[NAsamples] <- NA
sens_sA[NAsamples] <- NA
sens_pB_Y[NAsamples] <- NA
sens_pB_A[NAsamples] <- NA
sens_CS_Y[NAsamples] <- NA
sens_CS_A[NAsamples] <- NA
sens_pNS[NAsamples] <- NA
sens_sN_Y[NAsamples] <- NA
sens_sN_A[NAsamples] <- NA
sens_n_Y[NAsamples] <- NA
sens_n_A[NAsamples] <- NA
sens_imm_Y[NAsamples] <- NA
sens_imm_A[NAsamples] <- NA


## Assemble sensitivities in a list
SensVecs <- list(
  sens_sJ = sens_sJ, sens_sA = sens_sA,
  sens_pB_Y = sens_pB_Y, sens_pB_A = sens_pB_A,
  sens_CS_Y = sens_CS_Y, sens_CS_A = sens_CS_A,
  sens_pNS = sens_pNS,
  sens_sN_Y = sens_sN_Y, sens_sN_A = sens_sN_A,
  sens_n_Y = sens_n_Y, sens_n_A = sens_n_A,
  sens_imm_Y = sens_imm_Y, sens_imm_A = sens_imm_A
)

###############################################
#### CALCULATION OF TRANSIENT ELASTICITIES ####
###############################################

## Calculate transient elasticities for vital rates
elas_sJ <- sens_sJ*(sJ_mean/lambda_mean)
elas_sA <- sens_sA*(sA_mean/lambda_mean)

elas_pB_Y <- sens_pB_Y*(pB_Y_mean/lambda_mean)
elas_pB_A <- sens_pB_A*(pB_A_mean/lambda_mean)

elas_CS_Y <- sens_CS_Y*(CS_Y_mean/lambda_mean)
elas_CS_A <- sens_CS_A*(CS_A_mean/lambda_mean)

elas_pNS <- sens_pNS*(pNS_mean/lambda_mean)

elas_sN_Y <- sens_sN_Y*(sN_Y_mean/lambda_mean)
elas_sN_A <- sens_sN_A*(sN_A_mean/lambda_mean)


## Calculate transient elasticities for population structure
elas_n_Y <- sens_n_Y*(n_Y_mean/lambda_mean)
elas_n_A <- sens_n_A*(n_A_mean/lambda_mean)

elas_imm_Y <- sens_imm_Y*(imm_Y_mean/lambda_mean)
elas_imm_A <- sens_imm_A*(imm_A_mean/lambda_mean)

## Assemble elasticities in a list
ElasVecs <- list(
  elas_sJ = elas_sJ, elas_sA = elas_sA,
  elas_pB_Y = elas_pB_Y, elas_pB_A = elas_pB_A,
  elas_CS_Y = elas_CS_Y, elas_CS_A = elas_CS_A,
  elas_pNS = elas_pNS,
  elas_sN_Y = elas_sN_Y, elas_sN_A = elas_sN_A,
  elas_n_Y = elas_n_Y, elas_n_A = elas_n_A,
  elas_imm_Y = elas_imm_Y, elas_imm_A = elas_imm_A
)

#############################################################
#### CALCULATE LTRE CONTRIBUTIONS - POPULATION STRUCTURE ####
#############################################################

message(crayon::cyan('Calculating LTRE contributions...'))

## Prepare vectors to store results
cont_sJ <- rep(NA, nosamples)
cont_sA <- rep(NA, nosamples)

cont_pB_Y <- rep(NA, nosamples)
cont_pB_A <- rep(NA, nosamples)

cont_CS_Y <- rep(NA, nosamples)
cont_CS_A <- rep(NA, nosamples)

cont_pNS <- rep(NA, nosamples)

cont_sN_Y <- rep(NA, nosamples)
cont_sN_A <- rep(NA, nosamples)

cont_n_Y <- rep(NA, nosamples)
cont_n_A <- rep(NA, nosamples)

cont_imm_Y <- rep(NA, nosamples)
cont_imm_A <- rep(NA, nosamples)

cont_total <- rep(NA, nosamples)

est_var <- matrix(NA, nrow = 13, ncol = nosamples)
est_covar <- matrix(NA, nrow = 13, ncol = nosamples)


## Calculate LTRE contributions (random design LTRE)

for(i in 1:nosamples){
  
  ## Make matrix of vital rates and population structures
  dp_stoch <- cbind(sJ[i,1:(noyears-1)],
  					        sA[i,1:(noyears-1)],
  					        pB_Y[i,1:(noyears-1)],
  					        pB_A[i,1:(noyears-1)],
  					        CS_Y[i,1:(noyears-1)],
  					        CS_A[i,1:(noyears-1)],
  					        pNS[i,1:(noyears-1)],
  					        sN_Y[i,1:(noyears-1)],
  					        sN_A[i,1:(noyears-1)],
                    n_Y[i,1:(noyears-1)],
                    n_A[i,1:(noyears-1)],
                    imm_Y[i,2:noyears],
  					        imm_A[i,2:noyears])
  
  ## Derive process variances and covariances
  dp_varcov <- var(dp_stoch)
  
  ## Save total estimated (co)variance per parameter
  est_var[,i] <- diag(dp_varcov)
  est_covar[,i] <- rowSums(dp_varcov)
  
  ## Make a vector of sensitivities
  sensvec <- c(sens_sJ[i],
               sens_sA[i],
               sens_pB_Y[i],
               sens_pB_A[i],
               sens_CS_Y[i],
               sens_CS_A[i],
               sens_pNS[i],
               sens_sN_Y[i],
               sens_sN_A[i],
               sens_n_Y[i],
               sens_n_A[i],
               sens_imm_Y[i],
               sens_imm_A[i])
  
  ## Calculate demographic contributions
  # NOTE: Here we multiply sensitivities and (co)variances
  
  cont.mat <- matrix(NA, nrow = length(sensvec), ncol = length(sensvec))
  for(k in 1:length(sensvec)){
    for(m in 1:length(sensvec)){
      cont.mat[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
    }
  }
  
  ## Summarise contributions (sum of variances and covariances)
  cont <- rowSums(cont.mat)
  
  cont_sJ[i] <- cont[1]
  cont_sA[i] <- cont[2]
  cont_pB_Y[i] <- cont[3]
  cont_pB_A[i] <- cont[4]
  cont_CS_Y[i] <- cont[5]
  cont_CS_A[i] <- cont[6]
  cont_pNS[i] <- cont[7]
  cont_sN_Y[i] <- cont[8]
  cont_sN_A[i] <- cont[9]
  cont_n_Y[i] <- cont[10]
  cont_n_A[i] <- cont[11]
  cont_imm_Y[i] <- cont[12]
  cont_imm_A[i] <- cont[13]
  cont_total[i] <- sum(cont)
}

## Assemble contributions in a list
ContVecs_n <- list(
  cont_sJ = cont_sJ, cont_sA = cont_sA,
  cont_pB_Y = cont_pB_Y, cont_pB_A = cont_pB_A,
  cont_CS_Y = cont_CS_Y, cont_CS_A = cont_CS_A,
  cont_pNS = cont_pNS,
  cont_sN_Y = cont_sN_Y, cont_sN_A = cont_sN_A,
  cont_n_Y = cont_n_Y, cont_n_A = cont_n_A,
  cont_imm_Y = cont_imm_Y, cont_imm_A = cont_imm_A
)

## Collate contributions (all samples)
ContData_n <- data.frame(
  PopID = PopID,
  cont = c(cont_sJ, cont_sA, cont_pB_Y, cont_pB_A,
           cont_CS_Y, cont_CS_A, cont_pNS,
           cont_sN_Y, cont_sN_A, cont_n_Y, cont_n_A,
           cont_imm_Y, cont_imm_A),
  parameter = rep(c('sJ', 'sA', 'pB_Y', 'pB_A',
                'CS_Y', 'CS_A', 'pNS',
                'sN_Y', 'sN_A', 'n_Y', 'n_A',
                'imm_Y', 'imm_A'), each = nosamples)
)


## Compare sum of contributions and variation in lambda
message('Variation in lambda:')
print(quantile(lambda_var, probs = c(0.025, 0.5, 0.975), na.rm = T))

message('Sum of contributions (n LTRE):')
print(quantile(cont_total, probs = c(0.025, 0.5, 0.975), na.rm = T))


#####################################
#### ASSEMBLING & SAVING RESULTS ####
#####################################

message(crayon::cyan('Assembling & saving results...'))

## Assembling results in a list
LTRE_Results <- list(
  PopID = PopID,
  SensVecs = SensVecs,
  ElasVecs = ElasVecs,
  ContVecs_n = ContVecs_n,
  ContData_n = ContData_n
)

## Saving output
saveRDS(LTRE_Results, file = paste0('randomLTRE_', PopID, '.rds'))
