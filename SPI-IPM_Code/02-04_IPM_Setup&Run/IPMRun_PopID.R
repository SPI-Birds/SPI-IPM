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
set.seed(mySeed)

## Set population code
PopID <- 'EDM'

## Source IPM setup code
source('IPMSetup.R')


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
#Inits <- list(SPI_IPM.inits(IPM.data = SPI.IPMdata, IPM.constants = SPI.IPMconstants, sampleRE = FALSE))
Inits <- list(SPI_IPM.inits(IPM.data = SPI.IPMdata, IPM.constants = SPI.IPMconstants, sampleRE = FALSE),
SPI_IPM.inits(IPM.data = SPI.IPMdata, IPM.constants = SPI.IPMconstants, sampleRE = FALSE),
SPI_IPM.inits(IPM.data = SPI.IPMdata, IPM.constants = SPI.IPMconstants, sampleRE = FALSE),
SPI_IPM.inits(IPM.data = SPI.IPMdata, IPM.constants = SPI.IPMconstants, sampleRE = FALSE))

####################
#### RUN NIMBLE ####
####################

SPI.IPM <- nimbleMCMC(code = SPI.IPMcode, constants = SPI.IPMconstants, data = SPI.IPMdata, inits = Inits, monitors = parameters, niter = ni, nburnin = nb, nchains = nc, thin = nt, setSeed = mySeed, samplesAsCodaMCMC = TRUE)

save(SPI.IPM, file = paste0("SPI-IPM_", PopID, ".RData"))
