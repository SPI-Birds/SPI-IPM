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
  load(paste0('/Users/chloe.nater/Dropbox/PiedFlycatcher_IPM/IPM_Code/210819_FlycatcherIPM_PostSamples_MS1/FlycatcherIPM_CovA_Sub_', PopID[i], '.RData'))
  
  out.mat <- as.matrix(PFC.IPM)
  
  sam.A <- melt(out.mat)
  sam.A$Model <- 'Integrated'  

  
  ## Independent analysis
  load(paste0('/Users/chloe.nater/Dropbox/PiedFlycatcher_IPM/IPM_Code/210803_FlycatcherIPM_ClusterRuns_IndepAnalysis/FlycatcherIPM_Indep_CovA_', PopID[i], '.RData'))
  
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

