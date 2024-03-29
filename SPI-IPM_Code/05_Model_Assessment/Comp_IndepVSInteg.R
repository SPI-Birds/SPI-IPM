#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)


PopID <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')

#--------------------------------------------------------#
# Covariate effects from different models, by population #
#--------------------------------------------------------#

for(i in 1:7){
  
  # Loading posteriors #
  #--------------------#
  
  ## Integrated analysis
  PFC.IPM <- readRDS(paste0('SPI-IPM_', PopID[i], '.rds'))
  
  out.mat <- as.matrix(PFC.IPM)
  
  sam.A <- melt(out.mat)
  sam.A$Model <- 'Integrated'  

  
  ## Independent analysis
  PFC.IPM <- readRDS(paste0('SPI-Indep_', PopID[i], '.rds'))
  
  out.mat <- as.matrix(PFC.IPM)

  sam.B <- melt(out.mat)
  sam.B$Model <- 'Independent'
  
  
  # Combining data #
  #----------------#
  
  ## Bind data together
  sam.data <- rbind(sam.A, sam.B)
  
  ## Re-name columns
  colnames(sam.data) <- c('Sample', 'Parameter', 'Estimate', 'Model')
  
  
  # Plotting - Vital rates #
  #------------------------#
  
  ## Define VR parameters 
  VR.params <- c('Mu.sJ', 'Mu.sA','sigma.sJ', 'sigma.sA',
                 'Mu.pB[1]', 'Mu.pB[2]', 'sigma.pB',
                 'Mu.CS[1]', 'Mu.CS[2]', 'sigma.CS',
                 'Mu.sN[1]', 'Mu.sN[2]', 'sigma.sN',
                 'Mu.pNS', 'sigma.pNS')
  
  ## Plot VRs
  pdf(paste0('Plots/ModelAssm/VR_Posteriors_IndepVSIntegr_', PopID[i], '.pdf'), width = 11, height = 8)
  
  print(ggplot(subset(sam.data, Parameter %in% VR.params), aes(x = Estimate, group = Model)) +
    geom_density(aes(color = Model, fill = Model), alpha = 0.25) +
    facet_wrap(~Parameter, scales = 'free') +
    ggtitle(paste0('Vital rate parameters (', PopID[i], ')')) +
    scale_fill_viridis(discrete = T) +
    scale_color_viridis(discrete = T) +
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold')))
  
  dev.off()

  rm(sam.data)
}

