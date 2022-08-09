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
DIN.IPM <- readRDS('FlycatcherIPM_CovA_Sub_DIN.rds')

EDM.IPM <- readRDS('FlycatcherIPM_CovA_Sub_EDM.rds')

KAT.IPM <- readRDS('FlycatcherIPM_CovA_Sub_KAT.rds')

NAG.IPM <- readRDS('FlycatcherIPM_CovA_Sub_NAG.rds')

NWA.IPM <- readRDS('FlycatcherIPM_CovA_Sub_NWA.rds')

OKE.IPM <- readRDS('FlycatcherIPM_CovA_Sub_OKE.rds')

TEI.IPM <- readRDS('FlycatcherIPM_CovA_Sub_TEI.rds')



# Re-organizing data #
#--------------------#

## Collect MCMC chains in a list of matrices
out.mat <- list(
  DIN = as.matrix(DIN.IPM),
  EDM = as.matrix(EDM.IPM),
  KAT = as.matrix(KAT.IPM),
  NAG = as.matrix(NAG.IPM),
  NWA = as.matrix(NWA.IPM),
  OKE = as.matrix(OKE.IPM),
  TEI = as.matrix(TEI.IPM)
)


# Function for estimating vital rate - covariate relationships #
#--------------------------------------------------------------#

## Make a vector of standardized covariates
SDcov <- seq(-2, 3.25, length.out = 100)

cov.predict <- function(MCMC.mat, effect, RepAgeClass, PopID){
  
  # Prepare data frames
  pNS.pred <- sN.pred <- sJ.pred <- data.frame(PopID = PopID,
                                               SDcov = SDcov, 
                                               median = rep(NA, length(SDcov)), 
                                               lCI = rep(NA, length(SDcov)), 
                                               uCI = rep(NA, length(SDcov)), 
                                               Covariate = effect)
  
  pNS.pred$VitalRate <- 'Nest success prob.'
  sN.pred$VitalRate <- 'Nestling surv.'
  sJ.pred$VitalRate <- 'Juvenile annual surv.'
  
  # Set effect sizes  
  if(effect == 'Temperature'){
    beta.pNS <- MCMC.mat[,"beta1.pNS"]
    beta.sN <- MCMC.mat[,"beta1.sN"]
    beta.sJ <- 0
  }
  
  if(effect == 'Rainfall'){
    beta.pNS <- MCMC.mat[,"beta2.pNS"]
    beta.sN <- MCMC.mat[,"beta2.sN"]
    beta.sJ <- MCMC.mat[,"beta3.sJ"]
  }
  
  
  # Make predictions 
  for(i in 1:length(SDcov)){
    
    # Nest success probability
    pNS <- plogis(qlogis(MCMC.mat[,'Mu.pNS']) + beta.pNS*SDcov[i])
    
    pNS.pred[i,'median'] <- quantile(pNS, 0.5)
    pNS.pred[i,'lCI'] <- quantile(pNS, 0.025)
    pNS.pred[i,'uCI'] <- quantile(pNS, 0.975)
    
    # Nestling survival
    sN <- plogis(qlogis(MCMC.mat[,paste0('Mu.sN[', RepAgeClass, ']')]) + beta.sN*SDcov[i])
    
    sN.pred[i,'median'] <- quantile(sN, 0.5)
    sN.pred[i,'lCI'] <- quantile(sN, 0.025)
    sN.pred[i,'uCI'] <- quantile(sN, 0.975)
    
    # Juvenile annual survival
    sJ <- plogis(qlogis(MCMC.mat[,'Mu.sJ']) + beta.sJ*SDcov[i])
    
    sJ.pred[i,'median'] <- quantile(sJ, 0.5)
    sJ.pred[i,'lCI'] <- quantile(sJ, 0.025)
    sJ.pred[i,'uCI'] <- quantile(sJ, 0.975)
    
  }
  
  results <- rbind(pNS.pred, sN.pred, sJ.pred)
  
  return(results)
}



# Estimating vital rate - covariate relationships for all populations #
#---------------------------------------------------------------------#

## Make a grid with PopID - covariate combinations
PopID_Eff <- expand.grid(PopID = names(out.mat), Cov = c('Temperature', 'Rainfall'))

## Run all combinations for temperature effects (both age classes)
TempEff_yr <- do.call("rbind", sapply(1:7, FUN = function(x) cov.predict(out.mat[[x]], 'Temperature', 1, names(out.mat)[x]), simplify = FALSE))
TempEff_ad <- do.call("rbind", sapply(1:7, FUN = function(x) cov.predict(out.mat[[x]], 'Temperature', 2, names(out.mat)[x]), simplify = FALSE))

TempEff_yr$AgeClass <- 'Yearling'
TempEff_ad$AgeClass <- 'Adult'

## Run all combinations for rainfall effects (both age classes)
RainEff_yr <- do.call("rbind", sapply(1:7, FUN = function(x) cov.predict(out.mat[[x]], 'Rainfall', 1, names(out.mat)[x]), simplify = FALSE))
RainEff_ad <- do.call("rbind", sapply(1:7, FUN = function(x) cov.predict(out.mat[[x]], 'Rainfall', 2, names(out.mat)[x]), simplify = FALSE))

RainEff_yr$AgeClass <- 'Yearling'
RainEff_ad$AgeClass <- 'Adult'

## Combine results
CovEff <- rbind(TempEff_yr, TempEff_ad, RainEff_yr, RainEff_ad)
CovEff$PopID <- factor(CovEff$PopID, levels = c('EDM', 'TEI', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Discard data on sJ (plotting with data not possible)
CovEff <- subset(CovEff, VitalRate != 'Juvenile annual surv.')

## Discard duplicate data for nest success and rename age class
CovEff <- subset(CovEff, !(VitalRate == 'Nest success prob.' & AgeClass == 'Yearling'))
CovEff$AgeClass <- ifelse(CovEff$VitalRate == 'Nest success prob.', 'Combined', CovEff$AgeClass)

## Re-order factor levels
CovEff$AgeClass <- factor(CovEff$AgeClass, levels = c('Combined', 'Yearling', 'Adult'))


# Plotting estimated relationships including raw data #
#-----------------------------------------------------#

# NOTE: Nest success and nestling survival only, as data for juvenile survival 
#       not plottable (latent survival)

## Load in environmental data
load('210413_EnvData.RData')

## Assemble raw data
raw.data <- list(
  DIN = NA,
  EDM = NA,
  KAT = NA,
  NAG = NA,
  NWA = NA,
  OKE = NA,
  TEI = NA
)

PopIDs <- names(raw.data)

for (i in 1:length(PopIDs)){
  
  # Set PopID
  PopID <- PopIDs[i]
  
  # Load input data
  PFC.data <- readRDS(paste0(PopID, '_IPMData.rds'))
  
  # Set max time index
  Tmax <- dim(PFC.data$YearIndeces)[1]
  
  # Collect relevant input data (individual level) in a data frame
  data.ind <- data.frame(
    rbind(
      cbind('Nest success prob.', PFC.data$F_FAge, PFC.data$F_year, PFC.data$anyFledged),
      cbind('Nestling surv.', PFC.data$NoF_FAge, PFC.data$NoF_year, PFC.data$NoFledged/PFC.data$NoLaid)
    )
  )
  colnames(data.ind) <- c('VitalRate', 'AgeClass', 'Year', 'Value')
  data.ind$Value <- as.numeric(data.ind$Value)
  data.ind$AgeClass <- dplyr::case_when(data.ind$VitalRate == 'Nest success prob.' ~ 'Combined',
                                        data.ind$AgeClass == 1 ~ 'Yearling', 
                                        data.ind$AgeClass == 2 ~ 'Adult')
  data.ind$Year <- as.integer(data.ind$Year) + PFC.data$YearIndeces[1,1] - 1  
  data.ind$AgeClass <- factor(data.ind$AgeClass, levels = c('Combined', 'Yearling', 'Adult'))
  
  # Assemble environmental data
  env.data <- eval(parse(text = paste0('EnvData$', PopID, '$WindowY')))
  data.env <- data.frame(Year = PFC.data$YearIndeces[,1], Temp = env.data$temp_n, Rain = env.data$rain_n)
  
  # Merge environmental data
  data.temp <- merge(data.ind, data.env[,c('Year', 'Temp')], by = 'Year')
  colnames(data.temp)[5] <- 'SDcov'
  data.temp$Covariate <- 'Temperature'
  
  data.rain <- merge(data.ind, data.env[,c('Year', 'Rain')], by = 'Year')
  colnames(data.rain)[5] <- 'SDcov'
  data.rain$Covariate <- 'Rainfall'
  
  # Merge and store data
  data.out <- rbind(data.temp, data.rain)
  data.out$PopID <- PopID
  
  raw.data[[i]] <- data.out
}

## Merge data into a single data frame
data.out <- data.frame(Reduce(rbind, raw.data))

## Re-order factor levels
data.out$PopID <- factor(data.out$PopID, levels = c('EDM', 'TEI', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Make summarised version of data (year-averaged)
data.out.sum <- plyr::ddply(data.out, .(VitalRate, AgeClass, Year, SDcov, Covariate, PopID), summarise,
                      AvgValue = mean(Value))

## Add combined factor levels
CovEff$Relationship <- dplyr::case_when(CovEff$AgeClass == 'Combined' ~ CovEff$VitalRate,
                                 TRUE ~ paste0(CovEff$VitalRate, ' (', CovEff$AgeClass, ')'))
data.out$Relationship <- dplyr::case_when(data.out$AgeClass == 'Combined' ~ data.out$VitalRate,
                                        TRUE ~ paste0(data.out$VitalRate, ' (', data.out$AgeClass, ')'))
data.out.sum$Relationship <- dplyr::case_when(data.out.sum$AgeClass == 'Combined' ~ data.out.sum$VitalRate,
                                        TRUE ~ paste0(data.out.sum$VitalRate, ' (', data.out.sum$AgeClass, ')'))

## Plot
pdf('Plots/EnvCov/CovEffects_withData.pdf', height= 10, width = 8)
ggplot(CovEff, aes(x = SDcov, y = median)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Covariate), alpha = 0.2) + 
  geom_jitter(data = data.out, aes(x = SDcov, y = Value, color = Covariate), shape = 1, alpha = 0.1, size = 1) + 
  geom_point(data = data.out.sum, aes(x = SDcov, y = AvgValue, color = Covariate), shape = 8, size = 1) + 
  geom_line(aes(color = Covariate)) + 
  xlab('Standardized covariate value') + 
  scale_color_manual(values = c('#00204DFF', '#A69D75FF')) + 
  scale_fill_manual(values = c('#00204DFF', '#A69D75FF')) + 
  facet_grid(rows = vars(PopID), cols = vars(Relationship), scales = 'free') + 
  theme_bw() + theme(legend.position = 'bottom', panel.grid = element_blank(), axis.title.y = element_blank())
dev.off()