#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)

# Loading posteriors #
#--------------------#

## Load posterior samples from all 7 runs
load('FlycatcherIPM_CovA_Sub_DIN.RData')
DIN.IPM <- PFC.IPM

load('FlycatcherIPM_CovA_Sub_EDM.RData')
EDM.IPM <- PFC.IPM

load('FlycatcherIPM_CovA_Sub_KAT.RData')
KAT.IPM <- PFC.IPM

load('FlycatcherIPM_CovA_Sub_NAG.RData')
NAG.IPM <- PFC.IPM

load('FlycatcherIPM_CovA_Sub_NWA.RData')
NWA.IPM <- PFC.IPM

load('FlycatcherIPM_CovA_Sub_OKE.RData')
OKE.IPM <- PFC.IPM

load('FlycatcherIPM_CovA_Sub_TEI.RData')
TEI.IPM <- PFC.IPM


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

## Set vital rate parameters to plot
VR.params <- c('Mu.sJ', 'Mu.sA', 'sigma.sJ', 'sigma.sA', 
               'Mu.pB[1]', 'Mu.pB[2]', 'sigma.pB', 
               'Mu.CS[1]', 'Mu.CS[2]', 'sigma.CS', 
               'Mu.sN[1]', 'Mu.sN[2]', 'sigma.sN',
               'Mu.pNS', 'sigma.pNS')

## List PopIDs and associated population labels
PopIDs <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')
PopID_Labels <- c('Dinas (DIN)', 'East Dartmoor (EDM)', 'Loch Katrine (KAT)',
                 'Forest of Dean (NAG)', 'North Wales (NWA)', 'Okehampton (OKE)',
                 'Teign (TEI)')

## Re-organize data for whole posteriors of vital rate parameters
post.data <- data.frame(Sample = NA, Parameter = NA, Estimate = NA, PopID = NA)

for(i in 1:7){
  out.data <- melt(out.mat[[i]][,VR.params])
  colnames(out.data) <- c('Sample', 'Parameter', 'Estimate')
  out.data$PopID <- PopIDs[i]
  post.data <- rbind(post.data, out.data)
}

post.data <- post.data %>%
  dplyr::filter(!is.na(Parameter)) %>%
  dplyr::mutate(AgeClass = dplyr::case_when(stringr::str_detect(Parameter, 'sJ') ~ 'Fledgling',
                                            stringr::str_detect(Parameter, 'sA') ~ 'Yearling/Adult',
                                            stringr::str_detect(Parameter, '1') ~ 'Yearling',
                                            stringr::str_detect(Parameter, '2') ~ 'Adult',
                                            TRUE ~ 'Combined'),
                Parameter2 = stringr::str_remove_all(Parameter, '\\[|\\]|1|2|A|J'))


post.data$AgeClass <- factor(post.data$AgeClass, levels = c('Fledgling', 'Yearling', 'Yearling/Adult', 'Adult', 'Combined'))


# Plotting: Vital rate posteriors, by population #
#------------------------------------------------#

## Subset data and assign vital rate names 
plot.data <- subset(post.data, Parameter2 %in% c('Mu.CS', 'Mu.pB', 'Mu.ND', 'Mu.pNS', 'Mu.s', 'Mu.sN')) %>%
  dplyr::mutate(VR = dplyr::case_when(Parameter2 == 'Mu.CS' ~ 'Clutch size',
                                      Parameter2 == 'Mu.pB' ~ 'Breeding probability',
                                      Parameter2 == 'Mu.ND' ~ 'Prop. yearling immigrants',
                                      Parameter2 == 'Mu.pNS' ~ 'Nest success probability',
                                      Parameter2 == 'Mu.s' ~ 'Annual survival',
                                      Parameter2 == 'Mu.sN' ~ 'Survival to fledging'))
plot.data$VR <- factor(plot.data$VR, levels = c('Breeding probability', 'Clutch size', 'Nest success probability',
                                                'Survival to fledging', 'Annual survival', 'Prop. yearling immigrants'))
## Plot to pdf
for(i in 1:7){
  pdf(paste0('Plots/VR_Posteriors_', PopIDs[i], '.pdf'), width = 7, height = 5)
  
  print(
  ggplot(subset(plot.data, PopID == PopIDs[i]), aes(x = Estimate, group = AgeClass)) + 
    geom_density(aes(color = AgeClass, fill = AgeClass), alpha = 0.5) +
    facet_wrap(~VR, nrow = 2, scales = 'free') +
    ggtitle(PopID_Labels[i]) + 
    scale_fill_viridis(discrete = T) + 
    scale_color_viridis(discrete = T) +
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))
  )
  
  dev.off()
  
}


# Plotting: Vital rate posteriors, across population #
#----------------------------------------------------#

## Order PopIDs by latitude
post.data$PopID <- factor(post.data$PopID, levels = c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Define custom color code for PopIDs
PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')

pdf('Plots/VR_Posteriors_crossPop.pdf', width = 12, height = 8)
ggplot(post.data, aes(x = Estimate, group = PopID)) + 
  geom_density(aes(color = PopID, linetype = PopID), fill = NA) +
  facet_wrap(~Parameter, scales = 'free') +
  scale_color_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))
dev.off()


# Extracting summaries for time-dependent parameters #
#----------------------------------------------------#

## Set the study years for each population
StudyYears <- list(
  DIN = 1980:2017,
  EDM = 1955:2020,
  KAT = 1974:2011,
  NAG = 1981:2019,
  NWA = 1986:2019,
  OKE = 1974:2019,
  TEI = 1981:2020
)

## Extract Tmax for all populations
PopTmax <- c(length(StudyYears$DIN),
             length(StudyYears$EDM),
             length(StudyYears$KAT),
             length(StudyYears$NAG),
             length(StudyYears$NWA),
             length(StudyYears$OKE),
             length(StudyYears$TEI)
             )

## Make an empty list for storing estimates for each population
sum.data <- list(
  DIN = NA,
  EDM = NA,
  KAT = NA,
  NAG = NA,
  NWA = NA,
  OKE = NA,
  TEI = NA
)

## Define a short-cut function for posterior summaries (median, 95% CI)
sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.995)))
}

## For each population, extract relevant measures
for(i in 1:length(PopIDs)){
  
  # Set sample matrix
  sam.mat <- out.mat[[i]]
  
  # Prepare matrices to store results
  Ntot <- Btot <- matrix(NA, nrow = PopTmax[i], ncol = 3, dimnames = list(NULL, c('lCI', 'Median', 'uCI')))
  linkvar.sJ <- linkvar.sA <- linkvar.pB <- linkvar.CS <- linkvar.pNS <- linkvar.sN <- Ntot
  var.immY <- var.immA <- Ntot
  
  # Make posterior summaries of relevant quantities
  for(t in 1:PopTmax[i]){
    
    Ntot[t,] <- sam.summary(sam.mat[,paste0('Ntot[', t, ']')])
    if(t < PopTmax[i]){
      Btot[t,] <- sam.summary(sam.mat[,paste0('Btot[', t, ']')])
      
      var.immY[t,] <- sam.summary(sam.mat[,paste0('Imm[1, ', t+1, ']')] / sam.mat[,paste0('Ntot[', t, ']')])
      var.immA[t,] <- sam.summary(sam.mat[,paste0('Imm[2, ', t+1, ']')] / sam.mat[,paste0('Ntot[', t, ']')])
    }
  
    linkvar.sJ[t,] <- sam.summary(qlogis(sam.mat[,paste0('sJ[', t, ']')]) - qlogis(sam.mat[,'Mu.sJ']))
    linkvar.sA[t,] <- sam.summary(qlogis(sam.mat[,paste0('sA[', t, ']')]) - qlogis(sam.mat[,'Mu.sA']))
    linkvar.pB[t,] <- sam.summary(log(sam.mat[,paste0('pB[2, ', t, ']')]) - log(sam.mat[,'Mu.pB[2]']))
    linkvar.CS[t,] <- sam.summary(log(sam.mat[,paste0('CS[2, ', t, ']')]) - log(sam.mat[,'Mu.CS[2]']))
    linkvar.pNS[t,] <- sam.summary(qlogis(sam.mat[,paste0('pNS[', t, ']')]) - qlogis(sam.mat[,'Mu.pNS']))
    linkvar.sN[t,] <- sam.summary(qlogis(sam.mat[,paste0('sN[2, ', t, ']')]) - qlogis(sam.mat[,'Mu.sN[2]']))
  }
  
  # Organise summaries in a data frame
  pop.sum.data <- data.frame(
    PopID = PopIDs[i],
    Year = StudyYears[[i]],
    Parameter = rep(c('Ntot', 'Btot', 'linkvar.sJ', 'linkvar.sA', 'linkvar.pB', 'linkvar.CS', 'linkvar.pNS', 'linkvar.sN', 'var.immY', 'var.immA'), each = PopTmax[i])
  )
  
  pop.sum.data <- cbind(
    pop.sum.data, 
    rbind(Ntot, Btot, linkvar.sJ, linkvar.sA, linkvar.pB, linkvar.CS, linkvar.pNS, linkvar.sN, var.immY, var.immA))
  
  # Insert summarised data into list
  sum.data[[i]] <- pop.sum.data
}

## Make a combined data frame with all populations
allPop.data <- dplyr::bind_rows(sum.data, .id = "column_label")

allPop.data$PopID <- factor(allPop.data$PopID, levels = c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))


# Plotting: Time-variation in populations sizes and vital rates - across populations #
#------------------------------------------------------------------------------------#

## Population sizes
pdf('Plots/PopSizes_Years_crossPop.pdf', width = 8.25, height = 11.75)

# ggplot(subset(allPop.data, Parameter == 'Ntot'), aes(x = Year, y = Median)) + 
#   geom_line(aes(color = PopID)) + 
#   geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
#   ggtitle('Total Population Size') + 
#   ylab('Estimate') + 
#   scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
#   scale_color_manual(values = PFC_ColorCode) +
#   scale_fill_manual(values = PFC_ColorCode) +
#   facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
#   theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

# ggplot(subset(allPop.data, Parameter == 'Btot'), aes(x = Year, y = Median)) + 
#   geom_line(aes(color = PopID)) + 
#   geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
#   ggtitle('Breeding Population Size') + 
#   ylab('Estimate') + 
#   scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
#   scale_color_manual(values = PFC_ColorCode) +
#   scale_fill_manual(values = PFC_ColorCode) +
#   facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
#   theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter %in% c('Ntot', 'Btot')), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = Parameter)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID, alpha = Parameter)) + 
  ggtitle('Total & Breeding Population Size') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_alpha_discrete(range = c(0.5, 0.25)) + 
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

dev.off()


## Vital rate variation
pdf('Plots/VRVar_Years_crossPop.pdf', width = 8.25, height = 11.75)

ggplot(subset(allPop.data, Parameter == 'linkvar.sJ'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
  ggtitle('Variation in Juvenile Annual Survival') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.sA'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
  ggtitle('Variation in Adult Annual Survival') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.pB'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
  ggtitle('Variation in Breeding Probability') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.CS'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
  ggtitle('Variation in Clutch Size') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.pNS'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
  ggtitle('Variation in Nest Success') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.sN'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = PopID), alpha = 0.5) + 
  ggtitle('Variation in Survival to Fledging') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  facet_wrap(~PopID, ncol = 1, scales = 'free_y') + 
  theme_bw() + theme(legend.position = 'none', panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

dev.off()


# Plotting: Time-variation in populations sizes and vital rates - across populations, with overlay #
#--------------------------------------------------------------------------------------------------#

## Standardize median population sizes
allPop.std <- allPop.data %>%
  dplyr::group_by(PopID, Parameter) %>%
  dplyr::summarise(Median_std = (Median - mean(Median, na.rm = T))/sd(Median, na.rm = T), 
                   Year = Year, .groups = "keep")

## Population sizes
pdf('Plots/PopSizes_Years_crossPop2.pdf', width = 8, height = 3.5)

ggplot(subset(allPop.data, Parameter == 'Ntot'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Total Population Size') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.std, Parameter == 'Ntot'), aes(x = Year, y = Median_std)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Total Population Size (standardized)') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'Btot'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Breeding Population Size') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.std, Parameter == 'Btot'), aes(x = Year, y = Median_std)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Breeding Population Size (standardized)') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

dev.off()


## Vital rate variation
pdf('Plots/VRVar_Years_crossPop2.pdf', width = 8, height = 3.5)

ggplot(subset(allPop.data, Parameter == 'linkvar.sJ'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Juvenile Annual Survival') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.sA'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Adult Annual Survival') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.pB'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Breeding Probability') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.CS'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Clutch Size') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.pNS'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Nest Success') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'linkvar.sN'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Survival to Fledging') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'var.immY'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Yearling Immigration Rate') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

ggplot(subset(allPop.data, Parameter == 'var.immA'), aes(x = Year, y = Median)) + 
  geom_line(aes(color = PopID, linetype = PopID)) + 
  ggtitle('Variation in Adult Immigration Rate') + 
  ylab('Estimate') + 
  scale_x_continuous(breaks = seq(min(StudyYears$EDM), max(StudyYears$EDM), by = 5), minor_breaks = StudyYears$EDM, limits = c(1974, 2020)) + 
  scale_color_manual(values = PFC_ColorCode) +
  scale_fill_manual(values = PFC_ColorCode) +
  scale_linetype_manual(values = c(rep(c('solid', 'dashed'), 3), 'solid')) + 
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5), plot.title = element_text(face = 'bold'))

dev.off()