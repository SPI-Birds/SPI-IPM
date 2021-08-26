library(reshape)
library(ggplot2)
library(viridis)
library(plyr)
library(colortools)

###############
#### SETUP ####
###############

## Make a list of all populations
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')


###########################################################
#### DATA ASSEMBLY AND FORMATTING - LTRE CONTRIBUTIONS ####
###########################################################

## Function for assembling and formatting fixed design LTRE results for a population
assemble_LTREdata <- function(PopID){
  
  # Set path for LTRE results
  DataPath <- paste0('fixedLTRE_', PopID, '.RData')

  # Load LTRE data
  load(DataPath)
  
  ## Make summaries of posteriors (medians + 95% CIs)
  LTRE.sum <- ddply(LTRE_Results, .(PopID, year1, year2, parameter), summarise, 
                    median = median(cont, na.rm = T), 
                    lCI = quantile(cont, probs = 0.025, na.rm = T), 
                    uCI = quantile(cont, probs = 0.975, na.rm = T))
  
  return(LTRE.sum)
}  

## Run for all populations
fLTRE_results <- do.call("rbind", sapply(1:7, FUN = function(i) assemble_LTREdata(PopID_List[i]), simplify = FALSE))

## Re-order factor levels
fLTRE_results$parameter <- factor(fLTRE_results$parameter, 
                               levels = c('pB_Y', 'pB_A', 'pB', 'CS_Y', 'CS_A', 
                                          'pNS', 'sN_Y', 'sN_A', 'sN',
                                          'sJ', 'sA',  
                                          'n_Y', 'n_A', 'N_Y', 'N_A',
                                          'imm_Y', 'imm_A', 'Imm_Y', 'Imm_A'))

fLTRE_results$PopID <- factor(fLTRE_results$PopID, 
                                  levels = c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

#############################################################
#### PLOTTING RESULTS - LTRE CONTRIBUTIONS  - ALL LEVELS ####
#############################################################

## All levels - standardized axes
pdf('ResultsAll_AxisSTD.pdf', width = 8.3, height = 11.7)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A'))), aes(x = year1, y = median, group = parameter)) + 
  geom_bar(aes(fill = parameter, color = parameter), stat = 'identity', position = 'stack') + 
  ylab('Contribution') + xlab('Year') + 
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T)  + 
  geom_hline(yintercept = 0, linetype = 'dotted') +
  facet_wrap(~ PopID, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

pdf('ResultsAll_Abs_AxisSTD.pdf', width = 8.3, height = 11.7)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A'))), aes(x = year1, y = abs(median), group = parameter)) + 
  geom_area(aes(fill = parameter, color = parameter), stat = 'identity', position = 'fill') + 
  ylab('Absolute contribution') + xlab('Year') + 
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T)  + 
  facet_wrap(~ PopID, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


## All levels - population-specific axes
pdf('ResultsAll_AxisSPC.pdf', width = 11.7, height = 8.3)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A'))), aes(x = year1, y = median, group = parameter)) + 
  geom_bar(aes(fill = parameter, color = parameter), stat = 'identity', position = 'stack') + 
  ylab('Contribution') + xlab('Year') + 
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T)  + 
  geom_hline(yintercept = 0, linetype = 'dotted') +
  facet_wrap(~ PopID, scales = 'free', ncol = 2) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


pdf('Plots_General/ResultsAll_Abs_AxisSPC.pdf', width = 11.7, height = 8.3)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A'))), aes(x = year1, y = abs(median), group = parameter)) + 
  geom_area(aes(fill = parameter, color = parameter), stat = 'identity', position = 'fill') + 
  ylab('Absolute contribution') + xlab('Year') + 
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T)  + 
  facet_wrap(~ PopID, scales = 'free', ncol = 2) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


########################################################################
#### PLOTTING RESULTS - LTRE CONTRIBUTIONS  - LOCAL PARAMETERS ONLY ####
########################################################################

## Local dynamics - standardized axes
pdf('Plots_General/ResultsLocal_AxisSTD.pdf', width = 8.3, height = 11.7)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = median, group = parameter)) + 
  geom_bar(aes(fill = parameter, color = parameter), stat = 'identity', position = 'stack') + 
  ylab('Contribution') + xlab('Year') + 
  scale_fill_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  scale_color_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  geom_hline(yintercept = 0, linetype = 'dotted') +
  facet_wrap(~ PopID, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


pdf('Plots_General/ResultsLocal_Abs_AxisSTD.pdf', width = 8.3, height = 11.7)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = abs(median), group = parameter)) + 
  geom_area(aes(fill = parameter, color = parameter), stat = 'identity', position = 'fill') + 
  ylab('Absolute contribution') + xlab('Year') + 
  scale_fill_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  scale_color_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  facet_wrap(~ PopID, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


## Local dynamics - population-specific axes
pdf('Plots_General/ResultsLocal_AxisSPC.pdf', width = 11.7, height = 8.3)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = median, group = parameter)) + 
  geom_bar(aes(fill = parameter, color = parameter), stat = 'identity', position = 'stack') + 
  ylab('Contribution') + xlab('Year') + 
  scale_fill_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  scale_color_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  geom_hline(yintercept = 0, linetype = 'dotted') +
  facet_wrap(~ PopID, scales = 'free', ncol = 2) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


pdf('Plots_General/ResultsLocal_Abs_AxisSPC.pdf', width = 11.7, height = 8.3)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = abs(median), group = parameter)) + 
  geom_area(aes(fill = parameter, color = parameter), stat = 'identity', position = 'fill') + 
  ylab('Absolute contribution') + xlab('Year') + 
  scale_fill_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  scale_color_manual(values = c(plasma(12)[9:1], plasma(12)[11:12])) + 
  facet_wrap(~ PopID, scales = 'free', ncol = 2) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


###################################################################
#### PLOTTING RESULTS - LTRE CONTRIBUTIONS  - VITAL RATES ONLY ####
###################################################################

## Vital rates - standardized axes
greens <- colorRampPalette(c("#00A69D", "white"))
pinks <- colorRampPalette(c("#8C085E", "white"))
plot(rep(1,10),col=pinks(10),pch=19,cex=3)
plot(rep(1,10),col=greens(10),pch=19,cex=3)

pdf('Plots_General/ResultsVR_AxisSTD.pdf', width = 8.3, height = 11.7)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'n_Y', 'n_A', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = median, group = parameter)) + 
  geom_bar(aes(fill = parameter, color = parameter), stat = 'identity', position = 'stack') + 
  ylab('Contribution') + xlab('Year') + 
  scale_fill_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  scale_color_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  geom_hline(yintercept = 0, linetype = 'dotted') +
  facet_wrap(~ PopID, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


pdf('Plots_General/ResultsVR_Abs_AxisSTD.pdf', width = 8.3, height = 11.7)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'n_Y', 'n_A', 'Imm_Y', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = abs(median), group = parameter)) + 
  geom_area(aes(fill = parameter, color = parameter), stat = 'identity', position = 'fill') + 
  ylab('Absolute contribution') + xlab('Year') + 
  scale_fill_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  scale_color_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  facet_wrap(~ PopID, scales = 'free_y', ncol = 1) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

## Vital rates - population-specific axes
pdf('Plots_General/ResultsVR_AxisSPC.pdf', width = 11.7, height = 8.3)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'Imm_Y', 'n_Y', 'n_A', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = median, group = parameter)) + 
  geom_bar(aes(fill = parameter, color = parameter), stat = 'identity', position = 'stack') + 
  ylab('Contribution') + xlab('Year') + 
  scale_fill_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  scale_color_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  geom_hline(yintercept = 0, linetype = 'dotted') +
  facet_wrap(~ PopID, scales = 'free', ncol = 2) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

pdf('Plots_General/ResultsVR_Abs_AxisSPC.pdf', width = 11.7, height = 8.3)
ggplot(subset(fLTRE_results, !(parameter %in% c('N_Y', 'N_A', 'n_Y', 'n_A', 'Imm_Y', 'Imm_A', 'imm_Y', 'imm_A'))), aes(x = year1, y = abs(median), group = parameter)) + 
  geom_area(aes(fill = parameter, color = parameter), stat = 'identity', position = 'fill') + 
  ylab('Absolute contribution') + xlab('Year') + 
  scale_fill_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  scale_color_manual(values = c(greens(10)[1:7], pinks(10)[c(6,2)])) + 
  facet_wrap(~ PopID, scales = 'free', ncol = 2) + 
  theme_bw() + theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

