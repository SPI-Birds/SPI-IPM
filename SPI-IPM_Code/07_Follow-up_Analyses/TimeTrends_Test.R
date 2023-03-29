#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)

# Loading posteriors #
#--------------------#

## Load posterior samples from all 7 runs
DIN.IPM <- readRDS('SPI-IPM_DIN.rds')

EDM.IPM <- readRDS('SPI-IPM_EDM.rds')

KAT.IPM <- readRDS('SPI-IPM_KAT.rds')

NAG.IPM <- readRDS('SPI-IPM_NAG.rds')

NWA.IPM <- readRDS('SPI-IPM_NWA.rds')

OKE.IPM <- readRDS('SPI-IPM_OKE.rds')

TEI.IPM <- readRDS('SPI-IPM_TEI.rds')


# Re-organizing data #
#--------------------#

## Collect MCMC chains in a list of matrices
out.mat <- list(
  TEI = as.matrix(TEI.IPM),
  EDM = as.matrix(EDM.IPM),
  OKE = as.matrix(OKE.IPM),
  NAG = as.matrix(NAG.IPM),
  DIN = as.matrix(DIN.IPM),
  NWA = as.matrix(NWA.IPM),
  KAT = as.matrix(KAT.IPM)
)

PopIDs <- c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT')


# Setting up for correlation analysis of time-dependent parameters #
#------------------------------------------------------------------#

## Set the study years for each population
StudyYears <- list(
  TEI = 1981:2020,
  EDM = 1955:2020,
  OKE = 1974:2019,
  NAG = 1981:2019,
  DIN = 1980:2017,
  NWA = 1986:2019,
  KAT = 1974:2011
)

## Extract Tmax for all populations
PopTmax <- c(length(StudyYears$TEI),
             length(StudyYears$EDM),
             length(StudyYears$OKE),
             length(StudyYears$NAG),
             length(StudyYears$DIN),
             length(StudyYears$NWA),
             length(StudyYears$KAT)
             )

## Make an empty list to contain estimates for each population
sum.data <- list(
  TEI = NA,
  EDM = NA,
  OKE = NA,
  NAG = NA,
  DIN = NA,
  NWA = NA,
  KAT = NA
)

## Define a short-cut function for posterior summaries
sam.summary <- function(x, na.rm){
  unname(quantile(x, probs = c(0.025, 0.5, 0.995), na.rm = na.rm))
}

## Set sample number
n.sample <- dim(out.mat[[1]])[1]


# Test for time-trends in population sizes and vital rates for each population #
#------------------------------------------------------------------------------#

for(i in 1:length(PopIDs)){
  
  ## Set sample matrix
  sam.mat <- out.mat[[i]]
  
  ## Prepare a dataframe to store results
  trend.data <- data.frame(
    PopID = PopIDs[i],
    Parameter = c('Ntot', 'Btot', 'Rtot', 'Immtot', 'immtot', 'YAratio', 'Rrate', 'sJ', 'sA', 'pB', 'CS', 'pNS', 'sN'),
    r_mean = NA, 
    r_lCI = NA,
    r_median = NA,
    r_uCI = NA
  )
  
  ## Assemble time-dependent population measures and vital rates (link scale)
  Pop_t <- list(
    Ntot_t = sam.mat[,paste0('Ntot[', 1:PopTmax[i], ']')],
    Btot_t = sam.mat[,paste0('Btot[', 1:(PopTmax[i]-1), ']')],
    Rtot_t = sam.mat[,paste0('Juv[1, ', 1:(PopTmax[i]-1), ']')] + sam.mat[,paste0('Juv[2, ', 1:(PopTmax[i]-1), ']')],
    Immtot_t = sam.mat[,paste0('Imm[1, ', 1:PopTmax[i], ']')] + sam.mat[,paste0('Imm[2, ', 1:PopTmax[i], ']')],
    
    immtot_t = (sam.mat[,paste0('Imm[1, ', 1:PopTmax[i], ']')] + sam.mat[,paste0('Imm[2, ', 1:PopTmax[i], ']')])/(sam.mat[,paste0('Ntot[', 1:PopTmax[i], ']')]),
    
    YAratio_t = sam.mat[,paste0('N[1, ', 1:PopTmax[i], ']')] / sam.mat[,paste0('N[2, ', 1:PopTmax[i], ']')],
    
    Rrate_t = (sam.mat[,paste0('Juv[1, ', 1:(PopTmax[i]-1), ']')] + sam.mat[,paste0('Juv[2, ', 1:(PopTmax[i]-1), ']')])/(sam.mat[,paste0('Btot[', 1:(PopTmax[i]-1), ']')]),
    
    linkvar.sJ_t = qlogis(sam.mat[,paste0('sJ[', 1:PopTmax[i], ']')]) - qlogis(sam.mat[,'Mu.sJ']),
    linkvar.sA_t = qlogis(sam.mat[,paste0('sA[', 1:PopTmax[i], ']')]) - qlogis(sam.mat[,'Mu.sA']),
    linkvar.pB_t = log(sam.mat[,paste0('pB[2, ', 1:PopTmax[i], ']')]) - log(sam.mat[,'Mu.pB[2]']),
    linkvar.CS_t = log(sam.mat[,paste0('CS[2, ', 1:PopTmax[i], ']')]) - log(sam.mat[,'Mu.CS[2]']),
    linkvar.pNS_t = qlogis(sam.mat[,paste0('pNS[', 1:PopTmax[i], ']')]) - qlogis(sam.mat[,'Mu.pNS']),
    linkvar.sN_t = qlogis(sam.mat[,paste0('sN[2, ', 1:PopTmax[i], ']')]) - qlogis(sam.mat[,'Mu.sN[2]'])
    )
  
  
  ## For each quantity: calculate time-trend by sample, then summarize
  for(x in 1:nrow(trend.data)){
    r_sam <- p_sam <- rep(NA, n.sample)
   
    for(n in 1:n.sample){
      r_sam[n] <- cor.test(Pop_t[[x]][n,], 1:length(Pop_t[[x]][n,]))$estimate
    }
    
    trend.data$r_mean[x] <- mean(r_sam, na.rm = T)
    trend.data[x, c('r_lCI', 'r_median', 'r_uCI')] <- sam.summary(r_sam, na.rm = T)
  }
  
  ## Insert summarised data into list
  sum.data[[i]] <- trend.data
}


# Combine results and save as data frame #
#----------------------------------------#

## Collate data from all populations
allPop.data <- dplyr::bind_rows(sum.data)

## Add a column indicating evidence strength (95% CI overlap with 0)
allPop.data$Evidence <- ifelse(sign(allPop.data$r_lCI) == sign(allPop.data$r_uCI), '*', '-')

## Save as csv
write.csv(allPop.data, 'TrendEstimates.csv', row.names = FALSE)

