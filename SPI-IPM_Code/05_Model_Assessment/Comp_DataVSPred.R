#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(plyr)
library(gridExtra)

# Setup #
#-------#

## Make a list of PopIDs
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')

## Set path for input data and MCMC samples
DataPath <- '/Users/chloe.nater/Dropbox/PiedFlycatcher_IPM/IPM_Code/210819_FlycatcherIPM_PostSamples_MS1/'

## Define a short-cut function for posterior summaries
sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.995), na.rm = T))
}

## Defining parameter factor labels for plotting
paramLabels <- data.frame(
  Parameter = c('Btot', 'ImmBtot', 'pNS', 'CS', 'sN', 'pC_RO'),
  ParameterL = factor(c('Obs. breeding population size', 'Obs. # breeding immigrants', 
                        'Proportion successful nests', 'Clutch size', 'Proportion surviving nestlings',
                        'Per capita reproductive output'), 
                      levels = c('Obs. breeding population size', 'Obs. # breeding immigrants', 
                                 'Clutch size', 'Proportion successful nests', 'Proportion surviving nestlings',
                                 'Per capita reproductive output'))
)


for (i in 1:7){
  
  # Assembling input data #
  #-----------------------#
  
  ## Set PopID
  PopID <- PopID_List[[i]]
  
  ## Load and rename input data
  load(paste0(DataPath, '201210_', PopID, '_IPMData.RData'))
  PFC.data <- eval(parse(text = paste0(PopID, '.data')))
  
  ## Set max time index
  Tmax <- dim(PFC.data$YearIndeces)[1]
  
  ## Reconstruct data variable for immigrant detection
  PropImmDetect <- PFC.data$PropCapBrood
  PropImmDetect[which(PFC.data$AR_Data == 1 & PFC.data$PropCapBrood == 0)] <- NA
  
  ## Collect relevant input data (population level) in a data frame
  data.pop <- data.frame(
    Parameter = rep(c('Btot', 'ImmBtot'), each = Tmax),
    AgeClass = rep(c('Combined', 'Combined'), each = Tmax),
    Year = rep(PFC.data$YearIndeces[,1], 2),
    Value = c(PFC.data$NestCount, c(NA, PFC.data$ImmNoObs[2:Tmax]))
  )
  data.pop <- merge(data.pop, paramLabels, all.x = T)
  data.pop$AgeClass <- factor(data.pop$AgeClass, levels = c('Combined', 'Yearling', 'Adult'))
  
  ## Collect relevant input data (individual level) in a data frame
  data.ind <- data.frame(
    rbind(
      cbind('CS', PFC.data$CS_FAge, PFC.data$CS_year, PFC.data$ClutchSize),
      cbind('pNS', PFC.data$F_FAge, PFC.data$F_year, PFC.data$anyFledged),
      cbind('sN', PFC.data$NoF_FAge, PFC.data$NoF_year, PFC.data$NoFledged/PFC.data$NoLaid),
      cbind('pC_RO', PFC.data$F_FAge, PFC.data$F_year, PFC.data$Fledged)
    )
  )
  colnames(data.ind) <- c('Parameter', 'AgeClass', 'Year', 'Value')
  data.ind$Value <- as.numeric(data.ind$Value)
  data.ind$AgeClass <- dplyr::case_when(data.ind$Parameter == 'pNS' ~ 'Combined',
                                        data.ind$AgeClass == 1 ~ 'Yearling', 
                                        data.ind$AgeClass == 2 ~ 'Adult')
  data.ind$Year <- as.integer(data.ind$Year) + PFC.data$YearIndeces[1,1] - 1  
  data.ind <- merge(data.ind, paramLabels, all.x = T)
  data.ind$AgeClass <- factor(data.ind$AgeClass, levels = c('Combined', 'Yearling', 'Adult'))
  
  ## Summarize individual-level input data
  data.ind.sum <- ddply(data.ind, .(ParameterL, AgeClass, Year), summarise,
                              AvgValue = mean(Value))
  
  
  # Making corresponding posterior predictions #
  #--------------------------------------------#
  
  ## Load posterior samples
  load(paste0(DataPath, 'FlycatcherIPM_CovA_Sub_', PopID, '.RData'))
  sam.mat <- as.matrix(PFC.IPM)
  
  ## Prepare matrices to store results
  Btot <- ImmBtot <- matrix(NA, nrow = Tmax, ncol = 3, dimnames = list(NULL, c('lCI', 'Median', 'uCI')))
  CS_Y <- CS_A <- pNS <- sN_Y <- sN_A <- RO_Y <- RO_A <- Btot
  
  ## Make posterior summaries of relevant quantities
  for(t in 1:Tmax){
    
    if(t < Tmax){
      Btot[t,] <- sam.summary(sam.mat[,paste0('Btot[', t, ']')])*PFC.data$NS_Data[t]
      ImmBtot[t+1,] <- sam.summary(sam.mat[,paste0('ImmB[1, ', t+1, ']')] + sam.mat[,paste0('ImmB[2, ', t+1, ']')])*PropImmDetect[t+1]
    }
    
    CS_Y[t,] <- sam.summary(sam.mat[,paste0('CS[1, ', t, ']')])
    CS_A[t,] <- sam.summary(sam.mat[,paste0('CS[2, ', t, ']')])
    pNS[t,] <- sam.summary(sam.mat[,paste0('pNS[', t, ']')])
    sN_Y[t,] <- sam.summary(sam.mat[,paste0('sN[1, ', t, ']')])
    sN_A[t,] <- sam.summary(sam.mat[,paste0('sN[2, ', t, ']')])
    RO_Y[t,] <- sam.summary(sam.mat[,paste0('Juv[1, ', t, ']')]/sam.mat[,paste0('B[1, ', t, ']')])*2
    RO_A[t,] <- sam.summary(sam.mat[,paste0('Juv[2, ', t, ']')]/sam.mat[,paste0('B[2, ', t, ']')])*2
    
  }
  
  ## Organise summaries in a data frame
  post.pred <- data.frame(
    Parameter = c(rep(c('Btot', 'ImmBtot', 'pNS', 'CS', 'sN', 'pC_RO', 'CS', 'sN', 'pC_RO'), each = Tmax)),
    AgeClass = c(rep('Combined', 3*Tmax), rep(c('Yearling', 'Adult'), each = 3*Tmax)),
    Year = rep(PFC.data$YearIndeces[,1], 9)
  )
  
  post.pred <- cbind(
    post.pred, 
    rbind(Btot, ImmBtot, pNS, CS_Y, sN_Y, RO_Y, CS_A, sN_A, RO_A))
  
  ## Add and re-order factor labels for plotting
  post.pred <- merge(post.pred, paramLabels, all.x = T)
  post.pred$AgeClass <- factor(post.pred$AgeClass, levels = c('Combined', 'Yearling', 'Adult'))
  
  # Plotting: population-level quantities #
  #---------------------------------------#
  
  plot.pop <- 
    ggplot(subset(post.pred, Parameter %in% c('Btot', 'ImmBtot')), aes(x = Year, y = Median)) + 
    geom_line(color = magma(7)[3]) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI), fill = magma(7)[3], alpha = 0.2) + 
    geom_point(data = data.pop, aes(x = Year, y = Value), color = magma(7)[3], shape = 8) + 
    scale_x_continuous(breaks = seq(PFC.data$YearIndeces[1,1], PFC.data$YearIndeces[Tmax,1], by = 5), minor_breaks = PFC.data$YearIndeces[,1]) + 
    facet_wrap(~ParameterL, nrow = 2, scales = 'free_y') + 
    theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

  
  # Plotting: individual-level quantities #
  #---------------------------------------#
  
  plot.ind <- 
    ggplot(subset(post.pred, !(Parameter %in% c('Btot', 'ImmBtot'))), aes(x = Year, y = Median)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = AgeClass), alpha = 0.2) + 
    geom_jitter(data = data.ind, aes(x = Year, y = Value, color = AgeClass), shape = 1, alpha = 0.1) + 
    geom_point(data = data.ind.sum, aes(x = Year, y = AvgValue, color = AgeClass), shape = 8) + 
    geom_line(aes(color = AgeClass)) + 
    scale_x_continuous(breaks = seq(PFC.data$YearIndeces[1,1], PFC.data$YearIndeces[Tmax,1], by = 5), minor_breaks = PFC.data$YearIndeces[,1]) + 
    #scale_color_manual(values = viridis(5)[c(2,1,3)]) + 
    #scale_fill_manual(values = viridis(5)[c(2,1,3)]) + 
    scale_color_manual(values = magma(7)[c(3,4,1)]) + 
    scale_fill_manual(values = magma(7)[c(3,4,1)]) + 
    facet_wrap(~ParameterL, nrow = 2, scales = 'free_y') + 
    theme_bw() + theme(legend.position = 'bottom', panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))
  
  # Plotting: final combined plot #
  #-------------------------------#
  
  pdf(paste0('Data_vs_Predictions_', PopID, '.pdf'), width = 8, height = 10)
    print(grid.arrange(plot.pop, plot.ind, ncol = 1, top = paste0('Data vs. Predictions: ', PopID)))
  dev.off()
  
}
