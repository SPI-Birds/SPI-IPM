#######################
#### DATA ASSEMBLY ####
#######################

message(crayon::green$underline$bold(paste0('Fixed design LTRE analyses for ', PopID)))

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


########################################
#### FUNCTION FOR FIXED DESIGN LTRE ####
########################################

# NOTES: 
# - The function runs a fixed design LTRE comparing timesteps t1->t2 and t2->t3
# - All calculations are vectorized, meaning they run on all i samples of the posterior distribution at once without the need for a for loop

fixedD.LTRE = function(t1, t2){
	
	#-------------------------------------------------
	# STEP 1) Organise samples for relevant quantities
	#-------------------------------------------------
	
	#message(crayon::cyan('Assembling posterior data...'))
	
	## Extract vital rates for relevant time-steps (all samples)
	pB_Y <- cbind(out.mat[, paste0('pB[1, ', t1, ']')], out.mat[, paste0('pB[1, ', t2, ']')])
				  
  pB_A <- cbind(out.mat[, paste0('pB[2, ', t1, ']')], out.mat[, paste0('pB[2, ', t2, ']')])
  
  CS_Y <- cbind(out.mat[, paste0('CS[1, ', t1, ']')], out.mat[, paste0('CS[1, ', t2, ']')])
  
  CS_A <- cbind(out.mat[, paste0('CS[2, ', t1, ']')], out.mat[, paste0('CS[2, ', t2, ']')])
  
  pNS <- cbind(out.mat[, paste0('pNS[', t1, ']')], out.mat[, paste0('pNS[', t2, ']')])
  
  sN_Y <- cbind(out.mat[, paste0('sN[1, ', t1, ']')], out.mat[, paste0('sN[1, ', t2, ']')])
  
  sN_A <- cbind(out.mat[, paste0('sN[2, ', t1, ']')], out.mat[, paste0('sN[2, ', t2, ']')])

  sJ <- cbind(out.mat[, paste0('sJ[', t1, ']')], out.mat[, paste0('sJ[', t2, ']')])
		
  sA <- cbind(out.mat[, paste0('sA[', t1, ']')], out.mat[, paste0('sA[', t2, ']')])
  
	
	## Extract population sizes for relevant time-steps (all samples)
	# NOTE: For N_tot, I am also taking the additional time-step t2+1 (for calculating lambda at t2)
	N_Y <- cbind(out.mat[, paste0('N[1, ', t1, ']')], out.mat[, paste0('N[1, ', t2, ']')])
  N_A <- cbind(out.mat[, paste0('N[2, ', t1, ']')], out.mat[, paste0('N[2, ', t2, ']')])
  N_tot <- cbind(out.mat[, paste0('Ntot[', t1, ']')], out.mat[, paste0('Ntot[', t2, ']')], out.mat[, paste0('Ntot[', t2+1, ']')])
    
  ## Extract immigrant numbers for relevant time-steps (all samples)
  # NOTE: Here, I am adding +1 to the time index since lambda_t depends on Imm_t+1
  Imm_Y <- cbind(out.mat[, paste0('Imm[1, ', t1+1, ']')], out.mat[, paste0('Imm[1, ', t2+1, ']')])
  Imm_A <- cbind(out.mat[, paste0('Imm[2, ', t1+1, ']')], out.mat[, paste0('Imm[2, ', t2+1, ']')])
    
  ## Calculate population/immigrant proportions and growth rates
	n_Y <- N_Y/N_tot[,1:2]
	n_A <- N_A/N_tot[,1:2]
	imm_Y <- Imm_Y/N_tot[,1:2]
	imm_A <- Imm_A/N_tot[,1:2]
	
	lambda <- cbind(N_tot[,2]/N_tot[,1], N_tot[,3]/N_tot[,2])
	
	lambda_local <- cbind((N_tot[,2]-Imm_Y[,1]-Imm_A[,1])/N_tot[,1], (N_tot[,3]-Imm_Y[,2]-Imm_A[,2])/N_tot[,2])
	
	#--------------------------------------------------
	# STEP 2) Calculate sensitivities (for mean values)
	#--------------------------------------------------
	
	#message(crayon::cyan('Calculating sensitivities...'))

	## Calculate averages between both time-steps
	sJ_mean <- rowMeans(sJ)
	sA_mean <- rowMeans(sA)

	pB_Y_mean <- rowMeans(pB_Y)
	pB_A_mean <- rowMeans(pB_A)

	CS_Y_mean <- rowMeans(CS_Y)
	CS_A_mean <- rowMeans(CS_A)

	pNS_mean <- rowMeans(pNS)

	sN_Y_mean <- rowMeans(sN_Y)
	sN_A_mean <- rowMeans(sN_A)

	n_Y_mean <- rowMeans(n_Y)
	n_A_mean <- rowMeans(n_A)

	imm_Y_mean <- rowMeans(imm_Y)
	imm_A_mean <- rowMeans(imm_A)
	
	N_Y_mean <- rowMeans(N_Y)
	N_A_mean <- rowMeans(N_A)
	N_tot_mean <- rowMeans(N_tot)
	
	Imm_Y_mean <- rowMeans(Imm_Y)
	Imm_A_mean <- rowMeans(Imm_A)

	lambda_mean <- rowMeans(lambda)

	## Identify unusable samples (lambda = Inf)
	NAsamples <- which(lambda_mean == Inf)
	

	## Calculate transient sensitivities for vital rates
  sens_sJ <- 0.5*pB_Y_mean*CS_Y_mean*pNS_mean*sN_Y_mean*n_Y_mean + 
    0.5*pB_A_mean*CS_A_mean*pNS_mean*sN_A_mean*n_A_mean
  
  sens_sA <- n_Y_mean + n_A_mean
  
  sens_pB_Y <- 0.5*CS_Y_mean*pNS_mean*sN_Y_mean*sJ_mean*n_Y_mean
  sens_pB_A <- 0.5*CS_A_mean*pNS_mean*sN_A_mean*sJ_mean*n_A_mean
  
  sens_CS_Y <- 0.5*pB_Y_mean*pNS_mean*sN_Y_mean*sJ_mean*n_Y_mean
  sens_CS_A <- 0.5*pB_A_mean*pNS_mean*sN_A_mean*sJ_mean*n_A_mean
  
  sens_pNS <- 0.5*sJ_mean*(pB_Y_mean*CS_Y_mean*sN_Y_mean*n_Y_mean + 
                             pB_A_mean*CS_A_mean*sN_A_mean*n_A_mean)
  
  sens_sN_Y <- 0.5*pB_Y_mean*CS_Y_mean*pNS_mean*sJ_mean*n_Y_mean
  sens_sN_A <- 0.5*pB_A_mean*CS_A_mean*pNS_mean*sJ_mean*n_A_mean
  
  
  ## Calculate transient sensitivities for population structure
	sens_n_Y <- 0.5*pNS_mean*sJ_mean*((pB_Y_mean*CS_Y_mean*sN_Y_mean) - 
	                                           (pB_A_mean*CS_A_mean*sN_A_mean))
	sens_n_A <- 0.5*pNS_mean*sJ_mean*((pB_A_mean*CS_A_mean*sN_A_mean) - 
	                                           (pB_Y_mean*CS_Y_mean*sN_Y_mean))
	
  sens_imm_Y <- sens_imm_A <- rep(1, nosamples)
	
	
  ## Calculate transient sensitivities for population sizes
	sens_N_Y <- ((0.5*pB_Y_mean*CS_Y_mean*pNS_mean*sN_Y_mean*sJ_mean) +
	                  sA_mean - lambda_mean)/N_tot_mean
	sens_N_A <- ((0.5*pB_A_mean*CS_A_mean*pNS_mean*sN_A_mean*sJ_mean) +
	                  sA_mean - lambda_mean)/N_tot_mean
	
	sens_Imm_Y <- 	sens_Imm_A <- 1 / N_tot_mean

	
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
	sens_N_Y[NAsamples] <- NA
	sens_N_A[NAsamples] <- NA
	sens_n_Y[NAsamples] <- NA
	sens_n_A[NAsamples] <- NA
	sens_Imm_Y[NAsamples] <- NA
	sens_Imm_A[NAsamples] <- NA
	sens_imm_Y[NAsamples] <- NA
	sens_imm_A[NAsamples] <- NA
	

	#----------------------------------------------------
	# STEP 3) Calculate contributions (fixed design LTRE)
	#----------------------------------------------------

	#message(crayon::cyan('Calculating LTRE contributions...'))

	## Calculate LTRE contributions for vital rates
  cont_sJ <- (sJ[,2] - sJ[,1])*sens_sJ
  cont_sA <- (sA[,2] - sA[,1])*sens_sA
  cont_pB_Y <- (pB_Y[,2] - pB_Y[,1])*sens_pB_Y
  cont_pB_A <- (pB_A[,2] - pB_A[,1])*sens_pB_A
  cont_CS_Y <- (CS_Y[,2] - CS_Y[,1])*sens_CS_Y
  cont_CS_A <- (CS_A[,2] - CS_A[,1])*sens_CS_A
  cont_pNS <- (pNS[,2] - pNS[,1])*sens_pNS
  cont_sN_Y <- (sN_Y[,2] - sN_Y[,1])*sens_sN_Y
  cont_sN_A <- (sN_A[,2] - sN_A[,1])*sens_sN_A
  	
  ## Calculate LTRE contributions for population structure
  cont_n_Y <- (n_Y[,2] - n_Y[,1])*sens_n_Y
  cont_n_A <- (n_A[,2] - n_A[,1])*sens_n_A
  cont_imm_Y <- (imm_Y[,2] - imm_Y[,1])*sens_imm_Y
  cont_imm_A <- (imm_A[,2] - imm_A[,1])*sens_imm_A
  	
  ## Calculate LTRE contributions for population sizes
  cont_N_Y <- (N_Y[,2] - N_Y[,1])*sens_N_Y
  cont_N_A <- (N_A[,2] - N_A[,1])*sens_N_A
  cont_Imm_Y <- (Imm_Y[,2] - Imm_Y[,1])*sens_Imm_Y
  cont_Imm_A <- (Imm_A[,2] - Imm_A[,1])*sens_Imm_A
  	
  	
  #-----------------------------------
	# STEP 4) Arrange and output results
	#-----------------------------------
	
  ## Collate contributions (all samples)
	ContData <- data.frame(
  		PopID = PopID,
  		t1 = unname(t1),
  		t2 = unname(t2),
  		lambda1 = lambda[,1], 
  		lambda2 = lambda[,2], 
  		lambda_local1 = lambda_local[,1], 
  		lambda_local2 = lambda_local[,2], 
  		cont = c(cont_sJ, cont_sA, cont_pB_Y, cont_pB_A,
        		 cont_CS_Y, cont_CS_A, cont_pNS,
        		 cont_sN_Y, cont_sN_A, 
        		 cont_n_Y, cont_n_A, cont_imm_Y, cont_imm_A,
        		 cont_N_Y, cont_N_A, cont_Imm_Y, cont_Imm_A),
  		parameter = rep(c('sJ', 'sA', 'pB_Y', 'pB_A',
                	  	  'CS_Y', 'CS_A', 'pNS',
                	  	  'sN_Y', 'sN_A', 
                	  	  'n_Y', 'n_A', 'imm_Y', 'imm_A',
                	  	  'N_Y', 'N_A', 'Imm_Y', 'Imm_A'), each = nosamples)
	)
	
  ## Calculate and add absolute contribution sums
  contSum_n_all <- 
    abs(cont_sJ) + abs(cont_sA) + abs(cont_pB_Y) + abs(cont_pB_A) +
    abs(cont_CS_Y) + abs(cont_CS_A) + abs(cont_pNS) +
    abs(cont_sN_Y) + abs(cont_sN_A) + 
    abs(cont_n_Y) + abs(cont_n_A) + abs(cont_imm_Y) + abs(cont_imm_A)
  
  contSum_n_local <- 
    abs(cont_sJ) + abs(cont_sA) + abs(cont_pB_Y) + abs(cont_pB_A) +
    abs(cont_CS_Y) + abs(cont_CS_A) + abs(cont_pNS) +
    abs(cont_sN_Y) + abs(cont_sN_A) + 
    abs(cont_n_Y) + abs(cont_n_A)
  
  ContData$contSum_n_all <- c(rep(contSum_n_all, 13), rep(NA, 4*length(contSum_n_all)))
  
  ContData$contSum_n_local <- c(rep(contSum_n_local, 11), rep(NA, 6*length(contSum_n_local)))
  
  ## Calculate relative contributions
  ContData$relcont_n_all <- abs(ContData$cont) / ContData$contSum_n_all
  ContData$relcont_n_local <- abs(ContData$cont) / ContData$contSum_n_local
  
	## Return data
	return(ContData)
    
}


##################################################################################
#### RUNNING FIXED DESIGN LTRE FOR ALL SUBSEQUENT TIME-STEPS & SAVING RESULTS ####
##################################################################################

## Make a list of pairs of matrices / time-steps to compare
t1 <- c(1:(length(StudyYears)-2))
t2 <- t1 + 1
years <- cbind(t1, t2)

## Run analysis for all pairs of matrices / time-steps
message(crayon::cyan('Running fixed design LTRE on all subsequent time steps...'))
LTRE_Results <- do.call("rbind", sapply(1:nrow(years), FUN = function(x) fixedD.LTRE(t1 = years[x,1], t2 = years[x,2]), simplify = FALSE))

## Add actual years compared
LTRE_Results$year1 <- LTRE_Results$t1 + min(StudyYears) - 1
LTRE_Results$year2 <- LTRE_Results$t2 + min(StudyYears) - 1

## Saving output
saveRDS(LTRE_Results, file = paste0('fixedLTRE_', PopID, '.rds'))
