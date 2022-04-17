library(reshape)
library(ggplot2)
library(ggridges)
library(tidyr)
library(viridis)

###############
#### SETUP ####
###############

## Make a list of all populations
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')


############################################################################
#### TOTAL DYNAMICS - DATA ASSEMBLY AND FORMATTING - LTRE CONTRIBUTIONS ####
############################################################################

## Function for assembling and formatting LTRE results for a population
assemble_periodLTREdata <- function(PopID, ComponentSum){
  
  # Setup #
  #-------#
  
  # Set path and name for relevant LTRE's
  DataPath <- paste0('periodLTRE_', PopID, '.rds') # Full LTRE
  
  # Load LTRE data
  LTRE_Results <- readRDS(DataPath)
  
  if(ComponentSum){
    
    # Collate summed contributions (three summary levels) #
    #-----------------------------------------------------#
    
    # Extract non-summarized contributions
    ContData_sum <- LTRE_Results$ContData
    ContData_sum$SummaryLevel <- 'None'
    ContData_sum$sample <- NULL
    
    # Convert contribution vectors into a matrix
    ContMat_sum <- matrix(
      unlist(LTRE_Results$ContVecs_tot), 
      ncol = length(LTRE_Results$ContVecs_tot),
      dimnames = list(NULL, names(LTRE_Results$ContVecs_tot)))
    
    # Assemble data for per-category contribution summaries
    ContData_sum_cat <- data.frame(
      PopID = PopID,
      cont = c(rowSums(ContMat_sum[,c('cont_sJ','cont_sA')]),
               rowSums(ContMat_sum[,c('cont_pB_Y','cont_pB_A')]),
               rowSums(ContMat_sum[,c('cont_CS_Y','cont_CS_A')]),
               ContMat_sum[,'cont_pNS'],
               rowSums(ContMat_sum[,c('cont_sN_Y','cont_sN_A')]),
               rowSums(ContMat_sum[,c('cont_imm_Y','cont_imm_A')])),
      parameter = rep(c('s', 'pB', 'CS', 'pNS', 'sN', 'imm'), each = dim(ContMat_sum)[1]),
      SummaryLevel = 'Category'
    )
    
    # Assemble data for overall contribution summaries
    ContData_sum_ove <- data.frame(
      PopID = PopID,
      cont = c(rowSums(ContMat_sum[,c('cont_sJ','cont_sA')]),
               rowSums(ContMat_sum[,c('cont_pB_Y','cont_pB_A',
                                      'cont_CS_Y','cont_CS_A', 
                                      'cont_pNS',
                                      'cont_sN_Y','cont_sN_A')]),
               rowSums(ContMat_sum[,c('cont_imm_Y','cont_imm_A')])),
      parameter = rep(c('Survival', 'Reproduction', 'Immigration'), each = dim(ContMat_sum)[1]),
      SummaryLevel = 'Overall'
    )
    
    # Combine and store data in list  
    ContData_out <- rbind(ContData_sum, ContData_sum_cat, ContData_sum_ove)
    
  }else{
    
    # Collate per-component contributions #
    #-------------------------------------#
    
    # Convert mean contribution vectors into a data frame
    ContData_mean <- 
      melt(matrix(
        unlist(LTRE_Results$ContVecs_mu), 
        ncol = length(LTRE_Results$ContVecs_mu),
        dimnames = list(NULL, names(LTRE_Results$ContVecs_mu)))) %>%
      
      dplyr::rename('Sample' = X1, 'FullName' = X2, 'Contribution' = value) %>%
      dplyr::mutate(PopID = PopID,
                    Parameter = stringr::str_remove_all(string = FullName, pattern = 'contA_|contn_'),
                    Type = 'Mean',
                    Pathway = ifelse(grepl('contn', FullName), 'Indirect', 'Direct'))
    
    # Convert sd contribution vectors into a data frame
    ContData_sd <- 
      melt(matrix(
        unlist(LTRE_Results$ContVecs_sig), 
        ncol = length(LTRE_Results$ContVecs_sig),
        dimnames = list(NULL, names(LTRE_Results$ContVecs_sig)))) %>%
      
      dplyr::rename('Sample' = X1, 'FullName' = X2, 'Contribution' = value) %>%
      dplyr::mutate(PopID = PopID,
                    Parameter = stringr::str_remove_all(string = FullName, pattern = 'contA_|contn_'),
                    Type = 'SD',
                    Pathway = ifelse(grepl('contn', FullName), 'Indirect', 'Direct'))
    
    # Combine and store data in list  
    ContData_out <- rbind(ContData_mean, ContData_sd)
    
  }
return(ContData_out)
}  

## Run for all populations - sum of component contributions
LTRE_SumVR <- do.call("rbind", sapply(1:7, FUN = function(i) assemble_periodLTREdata(PopID_List[i], ComponentSum = T), simplify = FALSE))

## Run for all populations - by component contributions
LTRE_CompVR <- do.call("rbind", sapply(1:7, FUN = function(i) assemble_periodLTREdata(PopID_List[i], ComponentSum = F), simplify = FALSE))



######################################################################
#### LONG-TERM DYNAMICS - PLOTTING RESULTS - SUMMED CONTRIBUTIONS ####
######################################################################

## Re-order factor levels
LTRE_SumVR$parameter <- factor(LTRE_SumVR$parameter, 
                               levels = c('pB_Y', 'pB_A', 'pB', 'CS_Y', 'CS_A', 'CS',
                                          'pNS', 'sN_Y', 'sN_A', 'sN', 'Reproduction',
                                          'sJ', 'sA', 's', 'Survival', 
                                          'imm_Y', 'imm_A', 'imm', 'Immigration'))

LTRE_SumVR$PopID <- factor(LTRE_SumVR$PopID, levels = c('EDM', 'TEI', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Define custom color scale
PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')

## Plotting: All levels
pdf('Plots_ComponentSums/SumCont_Ind_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(LTRE_SumVR, SummaryLevel == 'None'), aes(x = cont, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Contribution') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Plotting: Categories
pdf('Plots_ComponentSums/SumCont_Cat_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(LTRE_SumVR, SummaryLevel == 'Category'), aes(x = cont, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Contribution') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Plotting: Overall
pdf('Plots_ComponentSums/SumCont_All_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(LTRE_SumVR, SummaryLevel == 'Overall'), aes(x = cont, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Contribution') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Plotting: Survival vs. reproduction vs. immigration
pdf('Plots_ComponentSums/SumCont_SurvRepImm_AllPops.pdf', width = 8, height = 6)
ggplot(subset(LTRE_SumVR, parameter %in% c('Survival', 'Reproduction', 'Immigration') & cont > -0.25 & cont < 0.6), aes(x = cont, y = PopID)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = parameter, color = parameter), alpha = 0.6, scale = 1) +
  geom_vline(aes(xintercept = 0), color = 'white', linetype = 'dotted', size = 0.4) + 
  xlab('Contribution') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = c('#00A69D', '#8C085E', 'orange')) + 
  scale_color_manual(values = c('#00A69D', '#8C085E', 'orange')) + 
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'top', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())

dev.off()


############################################################################
#### LONG-TERM DYNAMICS - PLOTTING RESULTS - CONTRIBUTIONS BY COMPONENT ####
############################################################################

## Re-order factor levels
LTRE_CompVR$Parameter <- factor(LTRE_CompVR$Parameter, 
                               levels = c('pB_Y', 'pB_A', 'CS_Y', 'CS_A', 
                                          'sN_Y', 'sN_A',
                                          'sJ', 'sA', 
                                          'imm_Y', 'imm_A',
                                          'pNS'))
LTRE_CompVR$Type <- factor(LTRE_CompVR$Type, levels = c('Mean', 'SD'))
LTRE_CompVR$Pathway <- factor(LTRE_CompVR$Pathway, levels = c('Direct', 'Indirect'))


## Plotting: All levels - by population

for(i in 1:7){
  
  pdf(paste0('Plots_ByComponents/CompCont_', PopID_List[[i]], '.pdf'), width = 8, height = 10)
  print(ggplot(subset(LTRE_CompVR, PopID == PopID_List[[i]])) + 
    geom_density(aes(x = Contribution, y = ..scaled.., fill = Type, color = Type, linetype = Pathway), alpha = 0.3) +
    xlab('Density') + ylab('Contribution') + ggtitle(paste0('Contributions to long-term growth rate: ', PopID_List[[i]])) + 
    scale_fill_manual(values = viridis(6)[c(3,5)]) + 
    scale_color_manual(values = viridis(6)[c(3,5)]) + 
    facet_wrap(~Parameter, scales = 'free', ncol = 2) +
    theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()))
  dev.off()
  
  
}
