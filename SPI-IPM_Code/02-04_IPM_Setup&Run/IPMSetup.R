#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

##########################
#### DATA PREPARATION ####
##########################

## Load workspace containing formatted demographic data
load(paste0(PopID, '_IPMData.RData'))

str(PFC.data)

## Data overview (elements contained in EDM.data)

# y[i,t] = CJS capture history of individual i over all time steps t
# CHs[i,t] = multistate capture history of individual i over all time steps t (1 = juvenile, 2 = yearling/adult)
# ageclass[i,t] = age class of individual i during time step t (1 = juvenile, 2 = yearling/adult)
# first[i] = first capture of individual i
# last[i] = last occasion to consider for individual i

# y.sum[i,t] = Unique summarized CJS capture history i over all time steps t
# CHs.sum[i,t] = Unique summarized multistate capture history i over all time steps t (1 = juvenile, 2 = yearling/adult)
# ageclass.sum[i,t] = age class of unique capture history i during time step t (1 = juvenile, 2 = yearling/adult)
# first.sum[i] = first capture of unique capture history i
# last.sum[i] = last occasion to consider for unique capture history i
# CHs.count[i] = number of occurences of unique capture history i

# A.marray[t1,t2] = Number of yearlings/adults released at occasion t1 and recaptured at occasion t2 (m-array)
# J.marray[t1,t2] = Number of fledglings/juveniles released at occasion t1 and recaptured at occasion t2 (m-array)

# NestCount[t] = Total number of active nests in year t

# EggNoTot[t] = Total number of eggs laid in year t (summed over all broods)
# EggNoSP[t] = Proportion of broods that EggNoTot[t] is based on

# ClutchSize[x] = Clutch size for brood x
# CS_year[x] = Year of brood x
# CS_FAge[x] = Age class of the mother of brood x (1 = yearling, 2 = adult)

# FledgedTot[t] = Total number of fledglings produced in year t (summed over all broods)
# FledgedSP[t] = Proportion of broods that EggNoTot[t] is based on

# Fledged[x] = Number of fledglings from brood x
# Laid[x] = Number of eggs laid in brood x
# F_year[x] = Year of brood x
# F_FAge[x] = Age class of the mother of brood x (1 = yearling, 2 = adult)
# anyFledged[x] = Whether any fledglings were produced in brood x
# NoFledged[y] = Number of fledglings from successful brood y
# NoLaid[y] = Number of eggs laid in successful brood y
# NoF_year[y] = Year of successful brood y
# NoF_FAge[y] = Age class of the mother of successful brood y (1 = yearling, 2 = adult)
 
# ImmNoObs[t] = Total number of immigrants (= unmarked adults) captured in year t
# ImmAgeProps[a] = All-time average proportion of known age immigrants in age class a (1 = yearling, 2 = adult) 
# ImmAgeProps_t[a,t] = Proportion of known age immigrants in age class a (1 = yearling, 2 = adult) in year t 

# NS_Data[t] = Whether or not any nest data was collected in year t
# PropCapBrood[t] = the proportion of nests for which the female was captured in year t
# AR_Data[t] = Whether or not any adults were ringed/recaptured in year t

# YearIndeces = Table documenting assignment of time indeces to study years


## Make a data variable for immigrant detection
PropImmDetect <- PFC.data$PropCapBrood
PropImmDetect[which(PFC.data$AR_Data == 1 & PFC.data$PropCapBrood == 0)] <- NA

## Load workspace containing environmental data
load('210413_EnvData.RData')
env.data <- eval(parse(text = paste0('EnvData$', PopID, '$WindowY')))
#env.data <- eval(parse(text = paste0('EnvData$', PopID, '$WindowG')))

## Arrange constants
PFC.IPMconstants <- list(Tmax = length(PFC.data$NestCount), A = 2,
						             ageclass.sum = PFC.data$ageclass.sum, 
						             first.sum = PFC.data$first.sum, 
						             last.sum = PFC.data$last.sum,
						             n.CH = dim(PFC.data$y.sum)[1],
                    	   EggNoSP = PFC.data$EggNoSP, 
                    	   CS_year = PFC.data$CS_year, 
						             CS_FAge = PFC.data$CS_FAge, 
                    	   CS_X = length(PFC.data$ClutchSize),
                    	   FledgedSP = PFC.data$FledgedSP,
                    	   F_year = PFC.data$F_year, 
						             F_FAge = PFC.data$F_FAge, 
                    	   F_X = length(PFC.data$Fledged),
						             NoF_year = PFC.data$NoF_year, 
						             NoF_FAge = PFC.data$NoF_FAge, 
						             NoF_X = length(PFC.data$NoFledged),
                    	   #ImmAgeProp = PFC.data$ImmAgeProp,
						             NS_Data = PFC.data$NS_Data,
						             AR_Data = PFC.data$AR_Data,
                    	   N1.limit = round(max(PFC.data$NestCount)*2),
						             AvgImm.limit = round(max(PFC.data$ImmNoObs)*5))

## Arrange data
PFC.IPMdata <- list(y.sum = PFC.data$y.sum, 
                    CHs.count = PFC.data$CHs.count,
					          PropCapBrood = c(PFC.data$PropCapBrood, 1), # Adding 1 for Tmax+1 (required in model, but does not affect estimation)
					          NestCount = PFC.data$NestCount,
                    EggNoTot = PFC.data$EggNoTot,
                    ClutchSize = PFC.data$ClutchSize,
                    FledgedTot = PFC.data$FledgedTot,
                    anyFledged = PFC.data$anyFledged, 
					          NoFledged = PFC.data$NoFledged,
					          NoLaid = PFC.data$NoLaid,
                    ImmNoObs = PFC.data$ImmNoObs,
					          PropImmDetect = PropImmDetect,
					          cov1 = as.vector(env.data$temp_n),
					          cov2 = as.vector(env.data$rain_n),
					          cov3 = as.vector(env.data$rain_pf7))

##########################################################
#### SPECIFICATION OF CUSTOM DISTRIBUTION FOR CMR DATA ###
##########################################################

dCJS_vv_sum <- nimbleFunction(
  # It is assumed that the individual has already been captured.
  # Therefore, the first entry in x represents the first possible recapture event.
  # probSurvive[t] represents survival from t-1 to t.
  # probCapture[t] represents capture probability at time t.
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(1),
                 mult = double(0), #! NEWLY ADDED: argument stating number of occurences of same capture history in entire dataset 
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(probSurvive) < length(x) - 1)
      stop("Length of probSurvive must be at least length of data minus 1.")
    if (length(x) != length(probCapture)) stop("Length of probCapture does not match length of data.")
    if (x[1] != 1) stop("dCJS requires specifying first capture. x[1] must equal 1.")
    
    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    if (len == 0) {  ## l<1 should not occur, but just in case:
      len <- length(x)
    }
    for (t in 2:len) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive[t - 1]
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture[t]
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture[t])
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
        logProbData <- logProbData + log(probThisObs) * mult #! NEWLY ADDED: "mult"
      }
    }
    if (log) {
      return(logProbData)
    }
    return(exp(logProbData))
    returnType(double())
  }
)

rCJS_vv_sum <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(1),
                 probCapture = double(1),
                 mult = double(0),
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if (len < 2)
      stop("len must be greater than 1.")
    if(length(probSurvive) != len - 1)
      stop("Length of probSurvive is not the same as len - 1.")
    if(length(probCapture) != len)
      stop("Length of probCapture is not the same as len.")
    ans <- numeric(length = len, init = FALSE)
    ans[1] <- 1
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 2:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive[i - 1])
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture[i])
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)

registerDistributions(list(
  dCJS_vv_sum = list(
    BUGSdist = 'dCJS_vv_sum(probSurvive, probCapture, mult, len)',
    types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double(1)', 'mult = double(0)', 'len = double()'),
    discrete = TRUE
  )
))



####################
#### MODEL CODE ####
####################

PFC.IPMcode <- nimbleCode({
	
	#******************#
	# POPULATION MODEL #
	#******************#
	
	## Parameters:
	
	## Data & constants:
	
	
	#---------------#
	# Process model #
	#---------------#
	
 	for (t in 1:Tmax){
 		
 		## 1) Breeding decision
 		for(a in 1:A){
 			B[a,t] ~ dbin(pB[a,t], N[a,t])
 			
 		## 2) Offspring production
 		Juv[a,t] ~ dpois(B[a,t]*CS[a,t]*pNS[t]*sN[a,t]*0.5)
 		
 		}
 	}	
 		
 		## 3) Annual survival of local birds
		
 	  for (t in 1:(Tmax-1)){
		  # Juveniles -> Yearlings
		  localN[1,t+1] ~ dbin(sJ[t], sum(Juv[1:A,t]))
		
		  # Yearlings/Adults -> adults
 		  localN[2,t+1] ~ dbin(sA[t], sum(N[1:A,t]))
 		
 		  ## 4) Immigration
 		  for(a in 1:A){
 			  N[a,t+1] <- localN[a,t+1] + Imm[a,t+1]
 		  }
 		
 		  ## Sum population size
 		  Ntot[t] <- sum(N[1:A,t])
 		  Btot[t] <- sum(B[1:A,t])
 	  }	
    
    Ntot[Tmax] <- sum(N[1:A,Tmax])
  
    #----------------------#
    # Priors & constraints #
    #----------------------#
    
     ## Initial population sizes

     for(a in 1:A){
     	localN[a,1] ~ dcat(DU.prior[1:N1.limit])
      N[a,1] <- localN[a,1] + Imm[a,1]
     }
     
     DU.prior[1:N1.limit] <- 1/N1.limit
     
	# NOTES: 
	# The above formulation assumes a discrete uniform prior for initial population size,
	# which requires an upper boundary (N1.limit) to be provided as data
	# More informative priors (e.g. rounded normal with data-informed mean) could be considered.

	
	#**********************************#
	# BREEDING POPULATION COUNT MODULE #
	#**********************************#
	
	## Parameters:
	
	## Data & constants:

	
	#------------#
	# Likelihood #
	#------------#
	
 	for(t in 1:Tmax){ 		
     	NestCount[t] ~ dpois(sum(B[1:A,t])*NS_Data[t])
  }
	
	# NOTES: 
	# Alternatively, one could use a truncated normal distribution here instead of the Poisson
	
	#------------------------#
	# Priors and constraints #
	#------------------------#
	

  #********************#
	# IMMIGRATION MODULE #
	#********************#
	
	## Parameters:
	
	## Data and constants:

	
     #------------#
     # Likelihood #
     #------------#
     
     ## Likelihood for the number/age distribution of immigrant females
     for (t in 2:Tmax){
       
       # Latent true number of breeding immigrants (count observation model)
       ImmNoObs[t] ~ dpois(sum(ImmB[1:A,t])*PropImmDetect[t])
       
       # Number of breeding immigrants per age class
       for(a in 1:A){
         ImmB[a,t] ~ dbin(pB[a,t], Imm[a,t])
       }
     } 
     
     ImmNoObs[1] <- 0
     ImmB[1:A,1] <- 0
     
     
     #------------------------#
     # Priors and constraints #
     #------------------------#
     
     ## Latent number of all immigrants (breeding & non-breeding) per age class
     Imm[1:A,1] <- 0
     
     for(t in 2:Tmax){
       for(a in 1:A){

         ImmCont[a,t] ~ T(dnorm(AvgImm[a], sd = sigma.Imm[a]), 0, Inf)
         Imm[a,t] <- round(ImmCont[a,t])
       }
     }
     
    for(a in 1:A){
      AvgImm[a] ~ dunif(0, AvgImm.limit)
      sigma.Imm[a] ~ dunif(0, AvgImm.limit*10)
    }
     
     ## Immigrant detection (= marking) probability
     for(t in 1:Tmax){
       PropImmDetect[t] ~ dunif(0, 1)
     }
	
	#***********************#
	# MARK-RECAPTURE MODULE #
	#***********************#
	
	## Parameters:
	
	## Data & constants:
	
	#----------------------------------#
	# Likelihood (custom distribution) #
	#----------------------------------#
	
	## Define individual sequences of survival and breeding (=recapture) probabilities
	for(i in 1:n.CH){
	  for(t in first.sum[i]:last.sum[i]){
	    phi.CH[i,t] <- phi[ageclass.sum[i,t], t]
	    p.CH[i,t+1] <- pB[ageclass.sum[i,t], t+1]*PropCapBrood[t+1]
	  }
	  p.CH[i,first.sum[i]] <- 1
	}
	p.ind[1:n.CH,1] <- 1
	
	## Likelihood with custom distribution
	for (i in 1:n.CH){
	  
	  y.sum[i, first.sum[i]:last.sum[i]] ~ dCJS_vv_sum(
	    probSurvive = phi.CH[i, first.sum[i]:last.sum[i]],
	    probCapture = p.CH[i, first.sum[i]:last.sum[i]],
	    len = last.sum[i]-first.sum[i]+1,
	    mult = CHs.count[i])
	} #i
	
	
	#------------------------#
	# Priors and constraints #
	#------------------------#
	
	## Age- and time-dependent survival probabilities
	phi[1,1:Tmax] <- sJ[1:Tmax]*0.5 # Adjusted for sex ratio
	phi[2,1:Tmax] <- sA[1:Tmax]
	
	logit(sJ[1:Tmax]) <- logit(Mu.sJ) + beta3.sJ*cov3[1:Tmax] + epsilon.sJ[1:Tmax]
	logit(sA[1:Tmax]) <- logit(Mu.sA) + epsilon.sA[1:Tmax]
	
	Mu.sJ ~ dunif(0, 1)
	Mu.sA ~ dunif(0, 1)
	
	## Age- and time-dependent breeding probabilities
	for (a in 1:A){		
	  logit(pB[a,1:(Tmax+1)]) <- logit(Mu.pB[a]) + epsilon.pB[1:(Tmax+1)]
	  Mu.pB[a] ~ dunif(0, 1)
	}
	
	## Temporal random effects
	for(t in 1:Tmax){
	  epsilon.sJ[t] ~ dnorm(0, sd = sigma.sJ)
	  epsilon.sA[t] ~ dnorm(0, sd = sigma.sA)
	}
	
	for(t in 1:(Tmax+1)){
	  epsilon.pB[t] ~ dnorm(0, sd = sigma.pB)
	}
	
	sigma.sJ ~ dunif(0, 10)
	sigma.sA ~ dunif(0, 10)
	sigma.pB ~ dunif(0, 10)
	                 
  ## Covariate effects
	beta3.sJ ~ dunif(-5, 5)


	#********************#
	# CLUTCH SIZE MODULE #
	#********************#
	
	## Parameters:
	
	## Data & constants:

	
  #------------------------------#
	# Likelihood - age-independent #
	#------------------------------#
	
	for(t in 1:Tmax){
		
		# Expected "true" egg number (by mother age class)
		EggNo.ex[1:A,t] <- B[1:A,t]*CS[1:A,t] 
		
		# Observed egg number (corrected by data availaility)
		EggNoTot[t] ~ dpois(sum(EggNo.ex[1:A,t])*EggNoSP[t])
		
	}
    
	#------------------------------#
	# Likelihood - age-dependent   #
	#------------------------------#
	
	for(x in 1:CS_X){
		ClutchSize[x] ~ dpois(CS[CS_FAge[x], CS_year[x]])
	}

	# NOTES: 
	# The age-dependent likelihood can only be used when we know the age of 
	# the mother of each brood z in ClutchSize
	
	#------------------------#
	# Priors and constraints #
	#------------------------#

	## Age-dependent clutch size   
    for(a in 1:A){	
    	log(CS[a,1:Tmax]) <- log(Mu.CS[a]) + epsilon.CS[1:Tmax]
    	Mu.CS[a] ~ dunif(0, 10)
    }
    
    ## Temporal random effects
    for(t in 1:Tmax){
    	epsilon.CS[t] ~ dnorm(0, sd = sigma.CS)
    }
    
    sigma.CS ~ dunif(0, 10)
    
    
	#*************************#
	# FLEDGING SUCCESS MODULE #
	#*************************#
	
	## Parameters:
	
	## Data & constants:

	
  #------------------------------#
	# Likelihood - age-independent #
	#------------------------------#
	
	for(t in 1:Tmax){
		
		# Expected "true" fledgling number (by mother age class)
		Fledged.ex[1:A,t] <- EggNo.ex[1:A,t]*pNS[t]*sN[1:A,t] 
		
		# Observed egg number (corrected by data availaility)
		FledgedTot[t] ~ dpois(sum(Fledged.ex[1:A,t])*FledgedSP[t])
		
		# Alternative:
		#FledgedTot[t] ~ dbin(sN[t], LaidTot[t])
		
	}
	
	# NOTES:
	# The alternative likelihood is subsantially simpler, but only possible if we assume
	# survival to fledging to be independent of mother age
    
    
	#------------------------------#
	# Likelihood - age-dependent #
	#------------------------------#
	
  for(x in 1:F_X){
    anyFledged[x] ~ dbern(pNS[F_year[x]])
  }
    
	for(x in 1:NoF_X){
		NoFledged[x] ~ dbin(sN[NoF_FAge[x], NoF_year[x]], NoLaid[x])
	}

	# NOTES: 
	# The age-dependent likelihood can only be used when we know the age of 
	# the mother of each brood z in Fledged/Laid
	
	#------------------------#
	# Priors and constraints #
	#------------------------#
    
	## Time- and mother age-dependent fledgling survival  
    for(a in 1:A){	
    	logit(sN[a,1:Tmax]) <- logit(Mu.sN[a]) + beta1.sN*cov1[1:Tmax] + beta2.sN*cov2[1:Tmax] + epsilon.sN[1:Tmax]
    	Mu.sN[a] ~ dunif(0, 1)
    }
    
    logit(pNS[1:Tmax]) <- logit(Mu.pNS) + beta1.pNS*cov1[1:Tmax] + beta2.pNS*cov2[1:Tmax] + epsilon.pNS[1:Tmax]
    Mu.pNS ~ dunif(0, 1)
    
    ## Time-dependent fledgling survival  
    #logit(sN[1:Tmax]) <- logit(Mu.sN) + epsilon.sN[1:Tmax]
    #Mu.sN ~ dunif(0, 1)

    
    ## Temporal random effects
    for(t in 1:Tmax){
    	epsilon.sN[t] ~ dnorm(0, sd = sigma.sN)
      epsilon.pNS[t] ~ dnorm(0, sd = sigma.pNS) 
    }
    
    sigma.sN ~ dunif(0, 10)
    sigma.pNS ~ dunif(0, 10)
    
    ## Covariate effects
    beta1.sN ~ dunif(-5, 5)
    beta2.sN ~ dunif(-5, 5)
    
    beta1.pNS ~ dunif(-5, 5)
    beta2.pNS ~ dunif(-5, 5)
    
    
    #*****************************#
    # COVARIATE IMPUTATION MODULE #
    #*****************************#
    
    for(t in 1:Tmax){
      cov1[t] ~ dnorm(0, sd = 1)
      cov2[t] ~ dnorm(0, sd = 1)
      cov3[t] ~ dnorm(0, sd = 1)
    }
    
    # NOTE: The imputation model assumes that covariates 
    #       are a) normally distributed and b) standardized (scaled and centered).
    #       Other distributions of covariates can be accounted for by changing
    #       this part. 
    
})


########################
#### INITIAL VALUES ####
########################

## Source initial value function
source('InitSim.R')


####################
#### MCMC SETUP ####
####################

## Setting parameters to monitor
parameters <- c(
  'Mu.sJ', 'Mu.sA', 'Mu.pB', 'Mu.CS', 'Mu.sN', 'Mu.pNS', 
  'sJ', 'sA', 'pB', 'CS', 'sN', 'pNS', 
  'sigma.sJ', 'sigma.sA', 'sigma.pB', 'sigma.CS', 'sigma.sN', 'sigma.pNS',  
	'epsilon.sJ', 'epsilon.sA', 'epsilon.pB', 'epsilon.CS', 'epsilon.sN', 'epsilon.pNS',
  'AvgImm', 'sigma.Imm',
	'N', 'B', 'Imm', 'ImmB', 'Juv',
	'Ntot', 'Btot',
	'beta1.pNS', 'beta1.sN', 'beta2.pNS', 'beta2.sN', 'beta3.sJ'
  )
