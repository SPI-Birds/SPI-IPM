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
SDcov <- seq(-2.5, 2.5, length.out = 100)

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

## Run all combinations for temperature effects
TempEff <- do.call("rbind", sapply(1:7, FUN = function(x) cov.predict(out.mat[[x]], 'Temperature', 2, names(out.mat)[x]), simplify = FALSE))

## Run all combinations for rainfall effects
RainEff <- do.call("rbind", sapply(1:7, FUN = function(x) cov.predict(out.mat[[x]], 'Rainfall', 2, names(out.mat)[x]), simplify = FALSE))

## Combine results
CovEff <- rbind(TempEff, RainEff)
CovEff$PopID <- factor(CovEff$PopID, levels = c('EDM', 'TEI', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Discard non-existent temp effect on sJ
CovEff <- subset(CovEff, !(Covariate == 'Temperature' & VitalRate == 'Juvenile annual surv.'))


# Plotting effects within populations #
#-------------------------------------#

PopIDs <- names(out.mat)
PlotTitles <- c('Dinas (DIN)', 'East Dartmoor (EDM)', 'Loch Katrine (KAT)',
                'Forest of Dean (NAG)', 'North Wales (NWA)', 'Okehampton (OKE)', 'Teign (TEI)')


for(i in 1:7){
  
  pdf(paste0('Plots/EnvCov/CovEffects_byPop', PopIDs[i], '.pdf'), width = 5, height = 7) 
  
  print(ggplot(subset(CovEff, PopID == PopIDs[i]), aes(x = SDcov, y = median, group = Covariate)) + 
          geom_line(aes(color = Covariate, linetype = Covariate)) + 
          geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Covariate), color = NA, alpha = 0.3) + 
          ylab('Estimate') + xlab('Standardized covariate value') + 
          ggtitle(PlotTitles[i]) + 
          scale_color_manual(values = c('#00204DFF', '#A69D75FF')) + 
          scale_fill_manual(values = c('#00204DFF', '#A69D75FF')) + 
          facet_wrap(~VitalRate, scales = 'free_y', ncol = 1) + 
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)))
  
  dev.off()
}


# Plotting effects across populations #
#-------------------------------------#

PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')


## Nest success probability
pdf('Plots/EnvCov/CovEffects_pNS_crossPop.pdf', width = 5, height = 6) 

ggplot(subset(CovEff, VitalRate == 'Nest success prob.'), aes(x = SDcov, y = median, group = PopID)) + 
      geom_line(aes(color = PopID, linetype = PopID)) + 
      geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), color = NA, alpha = 0.1) + 
      ylab('Estimate') + xlab('Standardized covariate value') + 
      ggtitle('Nest Success Probability') + 
      scale_color_manual(values = PFC_ColorCode) + 
      scale_fill_manual(values = PFC_ColorCode) + 
      scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
      facet_wrap(~Covariate, scales = 'free_y', ncol = 1) + 
      theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

dev.off()

## Nestling survival
pdf('Plots/EnvCov/CovEffects_sN_crossPop.pdf', width = 5, height = 6) 

ggplot(subset(CovEff, VitalRate == 'Nestling surv.'), aes(x = SDcov, y = median, group = PopID)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), color = NA, alpha = 0.1) + 
  ylab('Estimate') + xlab('Standardized covariate value') + 
  ggtitle('Nestling Survival') + 
  scale_color_manual(values = PFC_ColorCode) + 
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  facet_wrap(~Covariate, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

dev.off()


## Nestling survival
pdf('Plots/EnvCov/CovEffects_sJ_crossPop.pdf', width = 5, height = 4) 

ggplot(subset(CovEff, VitalRate == 'Juvenile annual surv.'), aes(x = SDcov, y = median, group = PopID)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), color = NA, alpha = 0.1) + 
  ylab('Estimate') + xlab('Standardized covariate value') + 
  ggtitle('Juvenile Annual Survival') + 
  scale_color_manual(values = PFC_ColorCode) + 
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  facet_wrap(~Covariate, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

dev.off()
