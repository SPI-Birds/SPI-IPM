#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)

## Set seed
mySeed <- 0
set.seed(mySeed)

## Make a list of PopID's
PopIDs <- c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT')

## Set the study duration for each population
PopTmax <- list(
  DIN = 38,
  EDM = 66,
  KAT = 38,
  NAG = 39,
  NWA = 34,
  OKE = 46,
  TEI = 40
)

#------------------------------------------------------------#
# FUNCTION FOR GENERATING AND STORING YEAR-SPECIFIC MATRICES #
#------------------------------------------------------------#

# data.Tmax = maximum time-step in the data

generateProjMat <- function(data.Tmax, out.mat){
  
  ## Prepare an array to store projection matrices
  mat.Y <- array(0, dim = c(2, 2, data.Tmax-1))
  
  ## Prepare a matrix to store initial population sizes
  N0.Y <- matrix(NA, nrow = 2, ncol = data.Tmax-1)
  
  ## Prepare a matrix to store immigrant numbers
  Imm.Y <- matrix(NA, nrow = 2, ncol = data.Tmax-1)
  
  ## Extract year-specific matrices and population sizes based on posterior means
  for(t in 1:(data.Tmax-1)){
    
    # Extract vital rates
    sJ <- median(out.mat[, paste("sJ[", t, "]", sep = "")])
    sA <- median(out.mat[, paste("sA[", t, "]", sep = "")])
    
    pB_Y <- median(out.mat[, paste("pB[1, ", t, "]", sep = "")])
    pB_A <- median(out.mat[, paste("pB[2, ", t, "]", sep = "")])
    
    CS_Y <- median(out.mat[, paste("CS[1, ", t, "]", sep = "")])
    CS_A <- median(out.mat[, paste("CS[2, ", t, "]", sep = "")])
    
    pNS <- median(out.mat[, paste("pNS[", t, "]", sep = "")])
    
    sN_Y <- median(out.mat[, paste("sN[1, ", t, "]", sep = "")])
    sN_A <- median(out.mat[, paste("sN[2, ", t, "]", sep = "")])
    
    
    # Write projection matrix
    mat.Y[1, 1, t] <- 0.5*pB_Y*CS_Y*pNS*sN_Y*sJ
    mat.Y[1, 2, t] <- 0.5*pB_A*CS_A*pNS*sN_A*sJ
    
    mat.Y[2, 1, t] <- sA
    mat.Y[2, 2, t] <- sA
    
    
    # Extract starting population size
    N0.Y[1, t] <- median(out.mat[, paste("N[1, ", t, "]", sep = "")])
    N0.Y[2, t] <- median(out.mat[, paste("N[2, ", t, "]", sep = "")])
      
    # Extract immigrant numbers
    Imm.Y[1, t] <- median(out.mat[, paste("Imm[1, ", t+1, "]", sep = "")])
    Imm.Y[2, t] <- median(out.mat[, paste("Imm[2, ", t+1, "]", sep = "")])
  }
  
  ## Combine in list and return
  ProjMat.list <- list(mat.Y = mat.Y, N0.Y = N0.Y, Imm.Y = Imm.Y)
  return(ProjMat.list)
}




#--------------------------------------------------------------------------#
# FUNCTION FOR STOCHASTIC PROJECTION USING A PRE-DEFINED SEQUENCE OF YEARS #
#--------------------------------------------------------------------------#

stoch.Proj = function(YearSeq, ProjMat.list){
	
  # Determine projection length
	Tmax <- length(YearSeq)
	
	# Make population vector
	N <- matrix(NA, nrow = 2, ncol = Tmax+1)
	
	# Set starting population size
	N[,1] <- ProjMat.list$N0.Y[,YearSeq[1]]
	
	# Project population
	for(t in 1:Tmax){
		
		N[,t+1] <- (ProjMat.list$mat.Y[,,YearSeq[t]]%*%N[,t]) + ProjMat.list$Imm.Y[,YearSeq[t]]
	}
	
	return(N)	
}

#--------------------------------------------------------------------------#
# FUNCTION TO RUN A SERIES OF PROJECTIONS AND RETURN RESULTS AS DATA FRAME #
#--------------------------------------------------------------------------#

stoch.Sims = function(i, data.Tmax, sim.Tmax, out.mat){ # i = simulation number
	
  # Assemble estimated projection matrices and population sizes
  ProjMat.list <- generateProjMat(data.Tmax = data.Tmax, out.mat)
    
	# Make a random sequence of years
	YearSeq <- sample(c(1:(data.Tmax-1)), sim.Tmax, replace = T)
	
	# Population projection
	N <- stoch.Proj(YearSeq, ProjMat.list)
	
	# Arrange projection results in a data frame
	output <- data.frame(
	SimNo = i, 
	SimYear = c(1:(length(YearSeq)+1)), 
	Ntot = colSums(N),
	N_Y = N[1,],
	N_A = N[2,],
	p_Y = N[1,]/colSums(N),
	p_A = N[2,]/colSums(N)
	)
	
	# Return results
	return(output)
}


#-------------------------------------------------------------------#
# WRAPPER FUNCTION FOR RUNNING SIMULATIONS ON A SPECIFIC POPULATION #
#-------------------------------------------------------------------#

# PopID = Population identifier
# SimNoMax = Number of simulations to run
# SimYearMax = Number of years in each simulation

run.stochSim = function(PopID, SimNoMax, SimYearMax){
  
  ## Load posterior samples for specified population
  load(paste0('FlycatcherIPM_CovA_Sub_', PopID, '.RData'))
  out.mat <- as.matrix(PFC.IPM)
  
  ## Set the index of the final year with data for the selected population
  data.Tmax <- PopTmax[[which(names(PopTmax) == PopID)]]
  
  ## Running simulations
  sim.results <- do.call("rbind", sapply(1:SimNoMax, FUN = function(i) stoch.Sims(i, data.Tmax = data.Tmax, sim.Tmax = SimYearMax, out.mat = out.mat), simplify = FALSE))
  
  ## Adding PopID
  sim.results$PopID <- PopID
  
  ## Returning results
  return(sim.results)
}

#------------------------------------------#
# RUNNING SIMULATIONS AND PLOTTING RESULTS #
#------------------------------------------#

## Set number of replicates to simulate
SimNoMax <- 200

## Set number of years to simulate for
SimYearMax <- 100

## Running simulations for each population
#sim.test <- run.stochSim(PopID = 'EDM', SimNoMax, SimYearMax)
sim.all <- do.call("rbind", sapply(1:length(PopIDs), FUN = function(x) run.stochSim(PopID = PopIDs[x], SimNoMax, SimYearMax), simplify = FALSE))

## Order PopIDs by latitude
sim.all$PopID <- factor(sim.all$PopID, levels = c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Define custom color code for PopIDs
PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')

## Plot simulations over time

# Total population size
pdf("Plots/ModelAssm/Simulation_PopSize.pdf", width = 8.25, height = 11.75)
ggplot(sim.all, aes(x = SimYear, y = Ntot, group = SimNo)) + 
  geom_line(aes(color = PopID), alpha = 0.3, size = 0.3) + 
  scale_color_manual(values = PFC_ColorCode) +
  ylab('Population Size') + xlab('Simulation Year') + 
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
dev.off()

# Total population size
pdf("Plots/ModelAssm/Simulation_LogPopSize.pdf", width = 8.25, height = 11.75)
ggplot(sim.all, aes(x = SimYear, y = log(Ntot), group = SimNo)) + 
  geom_line(aes(color = PopID), alpha = 0.3, size = 0.3) + 
  scale_color_manual(values = PFC_ColorCode) +
  ylab('Log Population Size') + xlab('Simulation Year') + 
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')
dev.off()


# Population age structure
pdf("Plots/ModelAssm/Simulation_PopStructure.pdf", width = 8.25, height = 11.75)
ggplot(sim.all, aes(x = SimYear, group = SimNo)) + 
  geom_line(aes(y = p_Y), color = viridis(5)[2], alpha = 0.3, size = 0.3) + 
  geom_line(aes(y = p_A), color = viridis(5)[4], alpha = 0.3, size = 0.3) + 
  ylab('Proportions per age class') + xlab('Simulation Year') + 
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

