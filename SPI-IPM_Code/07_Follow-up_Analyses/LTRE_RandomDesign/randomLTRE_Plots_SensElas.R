library(reshape)
library(ggplot2)
library(ggridges)

###############
#### SETUP ####
###############

## Make a list of all populations
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')


#####################################################################################
#### TOTAL DYNAMICS - DATA ASSEMBLY AND FORMATTING - SENSITIVITES & ELASTICITIES ####
#####################################################################################

## Function for assembling and formatting local sensitivity and elasticity results for a population
assemble_SensElas <- function(PopID){
  
  # Set path and name for relevant LTRE's
  DataPath <- paste0('randomLTRE_', PopID, '.RData') # Full LTRE
  
  # Make empty data frames to store combined data
  SensData_out <- data.frame()
  ElasData_out <- data.frame()
  
  # Load LTRE data
  load(DataPath)
  
  # Convert sensitivity/elasticity vectors into a matrix
  SensMat <- matrix(
    unlist(LTRE_Results$SensVecs), 
    ncol = length(LTRE_Results$SensVecs),
    dimnames = list(NULL, names(LTRE_Results$SensVecs)))
    
  ElasMat <- matrix(
    unlist(LTRE_Results$ElasVecs), 
    ncol = length(LTRE_Results$ElasVecs),
    dimnames = list(NULL, names(LTRE_Results$ElasVecs)))
    
  # Assemble data for single sensitivities/elasticities
  SensData_n <- data.frame(
    PopID = PopID,
    Quantity = 'Sensitivity', 
    value = c(SensMat[,c(1:9, 12:13, 16:17)]),
    parameter = rep(c('sJ', 'sA', 'pB_Y', 'pB_A', 'CS_Y', 'CS_A', 'pNS', 'sN_Y', 'sN_A', 'n_Y', 'n_A', 'imm_Y', 'imm_A'), each = dim(SensMat)[1]),
    SummaryLevel = 'None'
  )
    
  ElasData_n <- data.frame(
    PopID = PopID,
    Quantity = 'Elasticity', 
    value = c(SensMat[,c(1:9, 12:13, 16:17)]),
    parameter = rep(c('sJ', 'sA', 'pB_Y', 'pB_A', 'CS_Y', 'CS_A', 'pNS', 'sN_Y', 'sN_A', 'n_Y', 'n_A', 'imm_Y', 'imm_A'), each = dim(SensMat)[1]),
    SummaryLevel = 'None'
  )
    
  # Assemble data for per-category sensitivity/elasticity summaries
  SensData_n_cat <- data.frame(
    PopID = PopID,
    Quantity = 'Sensitivity',
    value = c(rowSums(SensMat[,c('sens_sJ','sens_sA')]),
              rowSums(SensMat[,c('sens_pB_Y','sens_pB_A')]),
              rowSums(SensMat[,c('sens_CS_Y','sens_CS_A')]),
              SensMat[,'sens_pNS'],
              rowSums(SensMat[,c('sens_sN_Y','sens_sN_A')]),
              rowSums(SensMat[,c('sens_n_Y','sens_n_A')]),
              rowSums(SensMat[,c('sens_imm_Y','sens_imm_A')])),
    parameter = rep(c('s', 'pB', 'CS', 'pNS', 'sN', 'n', 'imm'), each = dim(SensMat)[1]),
    SummaryLevel = 'Category'
  )
    
  ElasData_n_cat <- data.frame(
    PopID = PopID,
    Quantity = 'Elasticity',
    value = c(rowSums(ElasMat[,c('elas_sJ','elas_sA')]),
              rowSums(ElasMat[,c('elas_pB_Y','elas_pB_A')]),
              rowSums(ElasMat[,c('elas_CS_Y','elas_CS_A')]),
              ElasMat[,'elas_pNS'],
              rowSums(ElasMat[,c('elas_sN_Y','elas_sN_A')]),
              rowSums(ElasMat[,c('elas_n_Y','elas_n_A')]),
              rowSums(ElasMat[,c('elas_imm_Y','elas_imm_A')])),
    parameter = rep(c('s', 'pB', 'CS', 'pNS', 'sN', 'n', 'imm'), each = dim(ElasMat)[1]),
    SummaryLevel = 'Category'
  )
    
  # Assemble data for overall sensitivity/elasticity summaries
  SensData_n_ove <- data.frame(
    PopID = PopID,
    Quantity = 'Sensitivity',
    value = c(rowSums(SensMat[,c('sens_sJ','sens_sA')]),
              rowSums(SensMat[,c('sens_pB_Y','sens_pB_A',
                                 'sens_CS_Y','sens_CS_A', 
                                 'sens_pNS',
                                 'sens_sN_Y','sens_sN_A')]),
              rowSums(SensMat[,c('sens_n_Y','sens_n_A')]),
              rowSums(SensMat[,c('sens_imm_Y','sens_imm_A')])),
    parameter = rep(c('Survival', 'Reproduction', 'Pop. structure', 'Immigration'), each = dim(SensMat)[1]),
    SummaryLevel = 'Overall'
  )
    
  SensData_n <- rbind(SensData_n, SensData_n_cat, SensData_n_ove)
  SensData_out <- rbind(SensData_out, SensData_n)
    
  ElasData_n_ove <- data.frame(
    PopID = PopID,
    Quantity = 'Elasticity',
    value = c(rowSums(ElasMat[,c('elas_sJ','elas_sA')]),
              rowSums(ElasMat[,c('elas_pB_Y','elas_pB_A',
                                 'elas_CS_Y','elas_CS_A', 
                                 'elas_pNS',
                                 'elas_sN_Y','elas_sN_A')]),
              rowSums(ElasMat[,c('elas_n_Y','elas_n_A')]),
              rowSums(ElasMat[,c('elas_imm_Y','elas_imm_A')])),
    parameter = rep(c('Survival', 'Reproduction', 'Pop. structure', 'Immigration'), each = dim(ElasMat)[1]),
    SummaryLevel = 'Overall'
  )
    
  ElasData_n <- rbind(ElasData_n, ElasData_n_cat, ElasData_n_ove)
  ElasData_out <- rbind(ElasData_out, ElasData_n)

  SensElasData_out <- rbind(SensData_out, ElasData_out)
  
  return(SensElasData_out)
}  

## Run for all populations
SensElas_total <- do.call("rbind", sapply(1:7, FUN = function(i) assemble_SensElas(PopID_List[i]), simplify = FALSE))

## Re-order factor levels
SensElas_total$parameter <- factor(SensElas_total$parameter, 
                                   levels = c('pB_Y', 'pB_A', 'pB', 'CS_Y', 'CS_A', 'CS',
                                              'pNS', 'sN_Y', 'sN_A', 'sN', 'Reproduction',
                                              'sJ', 'sA', 's', 'Survival', 
                                              'n_Y', 'n_A', 'n', 'Pop. structure',
                                              'imm_Y', 'imm_A', 'imm', 'Immigration'))

SensElas_total$PopID <- factor(SensElas_total$PopID, levels = c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Define custom color scale
PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')


###########################################################
#### TOTAL DYNAMICS - PLOTTING RESULTS - SENSITIVITIES ####
###########################################################

## All levels
pdf('Plots_SensElas/ResultsInd_Sens_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(SensElas_total, Quantity == 'Sensitivity' & SummaryLevel == 'None'), aes(x = value, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Sensitivity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 7), axis.ticks.y = element_blank())
dev.off()

## Categories
pdf('Plots_SensElas/ResultsCat_Sens_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(SensElas_total, Quantity == 'Sensitivity' & SummaryLevel == 'Category'), aes(x = value, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Sensitivity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Overall
pdf('Plots_SensElas/ResultsAll_Sens_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(SensElas_total, Quantity == 'Sensitivity' & SummaryLevel == 'Overall'), aes(x = value, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Sensitivity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Survival vs. reproduction vs. immigration
pdf('Plots_SensElas/ResultsSurvRepImm_Sens_AllPops.pdf', width = 8, height = 6)
ggplot(subset(SensElas_total, Quantity == 'Sensitivity' & parameter %in% c('Survival', 'Reproduction', 'Immigration')), aes(x = value, y = PopID)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = parameter, color = parameter), alpha = 0.6, scale = 1) +
  xlab('Sensitivity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = c('#00A69D', '#8C085E', 'orange')) + 
  scale_color_manual(values = c('#00A69D', '#8C085E', 'orange')) + 
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'top', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())

dev.off()


##########################################################
#### LOCAL DYNAMICS - PLOTTING RESULTS - ELASTICITIES ####
##########################################################

## All levels
pdf('Plots_SensElas/ResultsInd_Elas_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(SensElas_total, Quantity == 'Elasticity' & SummaryLevel == 'None'), aes(x = value, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Elasticity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 7), axis.ticks.y = element_blank())
dev.off()

## Categories
pdf('Plots_SensElas/ResultsCat_Elas_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(SensElas_total, Quantity == 'Elasticity' & SummaryLevel == 'Category'), aes(x = value, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Elasticity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Overall
pdf('Plots_SensElas/ResultsAll_Elas_AllPops.pdf', width = 8.3, height = 11.7)
ggplot(subset(SensElas_total, Quantity == 'Elasticity' & SummaryLevel == 'Overall'), aes(x = value, y = parameter)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = PopID, color = PopID), alpha = 0.6, scale = 2) +
  xlab('Elasticity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = PFC_ColorCode) + 
  scale_color_manual(values = PFC_ColorCode) + 
  facet_wrap(~PopID, scales = 'free_x', ncol = 1) +
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'none', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())
dev.off()

## Survival vs. reproduction vs. immigration
pdf('Plots_SensElas/ResultsSurvRepImm_Elas_AllPops.pdf', width = 8, height = 6)
ggplot(subset(SensElas_total, Quantity == 'Elasticity' & parameter %in% c('Survival', 'Reproduction', 'Immigration')), aes(x = value, y = PopID)) + 
  geom_density_ridges(aes(height = ..ndensity.., fill = parameter, color = parameter), alpha = 0.6, scale = 1) +
  xlab('Elasticity') + ylab('Category') +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0), limits = rev) +
  scale_fill_manual(values = c('#00A69D', '#8C085E', 'orange')) + 
  scale_color_manual(values = c('#00A69D', '#8C085E', 'orange')) + 
  theme_ridges(font_size = 10, grid = FALSE) +
  theme(legend.title = element_blank(), legend.position = 'top', axis.text.y = element_text(face = 'bold'), axis.ticks.y = element_blank())

dev.off()
