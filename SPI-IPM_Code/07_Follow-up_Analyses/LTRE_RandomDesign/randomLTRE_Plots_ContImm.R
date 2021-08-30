library(reshape)
library(ggplot2)
library(ggridges)

###############
#### SETUP ####
###############

## Make a list of all populations
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')


############################################################################
#### TOTAL DYNAMICS - DATA ASSEMBLY AND FORMATTING - LTRE CONTRIBUTIONS ####
############################################################################

## Function for assembling and formatting LTRE results for a population
assemble_LTREdata <- function(PopID){
  
  # Set path and name for relevant LTRE's
  DataPath <- paste0('randomLTRE_', PopID, '.RData') # Full LTRE
  
  # Make an empty data frame to store combined data
  ContData_out <- data.frame()
  
  # Load LTRE data
  load(DataPath)
    
  # Extract non-summarized contributions
  #ContData_n <- subset(LTRE_Results$ContData_n, !(parameter %in% c('imm_Y', 'imm_A')))
  ContData_n <- LTRE_Results$ContData_n
  ContData_n$SummaryLevel <- 'None'
    
  # Convert contribution vectors into a matrix
  ContMat_n <- matrix(
    unlist(LTRE_Results$ContVecs_n), 
    ncol = length(LTRE_Results$ContVecs_n),
    dimnames = list(NULL, names(LTRE_Results$ContVecs_n)))
    
  # Assemble data for per-category contribution summaries
  ContData_n_cat <- data.frame(
    PopID = PopID,
    cont = c(rowSums(ContMat_n[,c('cont_sJ','cont_sA')]),
              rowSums(ContMat_n[,c('cont_pB_Y','cont_pB_A')]),
              rowSums(ContMat_n[,c('cont_CS_Y','cont_CS_A')]),
              ContMat_n[,'cont_pNS'],
              rowSums(ContMat_n[,c('cont_sN_Y','cont_sN_A')]),
              rowSums(ContMat_n[,c('cont_n_Y','cont_n_A')]),
             rowSums(ContMat_n[,c('cont_imm_Y','cont_imm_A')])),
    parameter = rep(c('s', 'pB', 'CS', 'pNS', 'sN', 'n', 'imm'), each = dim(ContMat_n)[1]),
    SummaryLevel = 'Category'
  )
    
  # Assemble data for overall contribution summaries
  ContData_n_ove <- data.frame(
    PopID = PopID,
    cont = c(rowSums(ContMat_n[,c('cont_sJ','cont_sA')]),
             rowSums(ContMat_n[,c('cont_pB_Y','cont_pB_A',
                                  'cont_CS_Y','cont_CS_A', 
                                  'cont_pNS',
                                  'cont_sN_Y','cont_sN_A')]),
              rowSums(ContMat_n[,c('cont_n_Y','cont_n_A')]),
             rowSums(ContMat_n[,c('cont_imm_Y','cont_imm_A')])),
    parameter = rep(c('Survival', 'Reproduction', 'Pop. structure', 'Immigration'), each = dim(ContMat_n)[1]),
    SummaryLevel = 'Overall'
  )
    
  ContData_n <- rbind(ContData_n, ContData_n_cat, ContData_n_ove)
  ContData_out <- rbind(ContData_out, ContData_n)
  
return(ContData_out)
}  

## Run for all populations
LTRE_total <- do.call("rbind", sapply(1:7, FUN = function(i) assemble_LTREdata(PopID_List[i]), simplify = FALSE))

## Re-order factor levels
LTRE_total$parameter <- factor(LTRE_total$parameter, 
                               levels = c('pB_Y', 'pB_A', 'pB', 'CS_Y', 'CS_A', 'CS',
                                          'pNS', 'sN_Y', 'sN_A', 'sN', 'Reproduction',
                                          'sJ', 'sA', 's', 'Survival', 
                                          'n_Y', 'n_A', 'n', 'Pop. structure',
                                          'imm_Y', 'imm_A', 'imm', 'Immigration'))

LTRE_total$PopID <- factor(LTRE_total$PopID, levels = c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Define custom color scale
PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')


################################################################
#### TOTAL DYNAMICS - PLOTTING RESULTS - LTRE CONTRIBUTIONS ####
################################################################

## All levels
pdf('Plots_ContImm/ResultsInd_Imm_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(LTRE_total, SummaryLevel == 'None'), aes(x = cont, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Contribution') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 7), axis.ticks.y = element_blank())
dev.off()

## Categories
pdf('Plots_ContImm/ResultsCat_Imm_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(LTRE_total, SummaryLevel == 'Category'), aes(x = cont, y = parameter)) + 
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

## Overall
pdf('Plots_ContImm/ResultsAll_Imm_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(LTRE_total, SummaryLevel == 'Overall'), aes(x = cont, y = parameter)) + 
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

## Survival vs. reproduction vs. immigration
pdf('Plots_ContImm/ResultsSurvRepImm_AllPops.pdf', width = 8, height = 6)
ggplot(subset(LTRE_total, parameter %in% c('Survival', 'Reproduction', 'Immigration') & cont > -0.025 & cont < 1), aes(x = cont, y = PopID)) + 
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
