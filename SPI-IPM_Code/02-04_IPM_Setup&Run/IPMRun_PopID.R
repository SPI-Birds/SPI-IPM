#******************************************************************#
# Code for running SPI-Birds Integrated Population Model (SPI-IPM) #
#            Example: EAST DARTMOOR (Pied Flycatcher)              #
#******************************************************************#

library(coda)
library(nimble)
#library(nimbleEcology)
nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

## Set seed
mySeed <- 1

## Set population code
PopID <- 'EDM'

## Source IPM setup code
source('FlycatcherIPM_Setup_CovA.R')


####################
#### MCMC SETUP ####
####################

## MCMC settings
ni <- 200000
nb <- 50000
nt <- 30
nc <- 4

#ni <- 2
#nb <- 0
#nt <- 1
#nc <- 1
#nc <- 3


## Sample initial values
#Inits <- list(PFC_IPM.inits(IPM.data = PFC.IPMdata, IPM.constants = PFC.IPMconstants, sampleRE = FALSE))
Inits <- list(PFC_IPM.inits(IPM.data = PFC.IPMdata, IPM.constants = PFC.IPMconstants, sampleRE = FALSE),
PFC_IPM.inits(IPM.data = PFC.IPMdata, IPM.constants = PFC.IPMconstants, sampleRE = FALSE),
PFC_IPM.inits(IPM.data = PFC.IPMdata, IPM.constants = PFC.IPMconstants, sampleRE = FALSE),
PFC_IPM.inits(IPM.data = PFC.IPMdata, IPM.constants = PFC.IPMconstants, sampleRE = FALSE))

####################
#### RUN NIMBLE ####
####################

PFC.IPM <- nimbleMCMC(code = PFC.IPMcode, constants = PFC.IPMconstants, data = PFC.IPMdata, inits = Inits, monitors = parameters, niter = ni, nburnin = nb, nchains = nc, thin = nt, setSeed = mySeed, samplesAsCodaMCMC = TRUE)

save(PFC.IPM, file = paste0("SPI-IPM_", PopID, ".RData"))
