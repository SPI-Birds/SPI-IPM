#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(ggstance)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)

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

## List PopIDs and associated population labels
PopIDs <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')
PopID_Labels <- c('Dinas (DIN)', 'East Dartmoor (EDM)', 'Loch Katrine (KAT)',
                  'Forest of Dean (NAG)', 'North Wales (NWA)', 'Okehampton (OKE)',
                  'Teign (TEI)')

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

## Define a short-cut function for posterior summaries (median, 95% CI)
sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.995)))
}

## Function to summarise average VRs
sum.avgVR <- function(pop.mat, PopID, Tmax){
  
  # Calculate posterior distribution for average immigration rates
  immY <- immA <- matrix(NA, nrow = dim(pop.mat)[1], Tmax)
  
  for(t in 1:(Tmax-1)){
    immY[,t+1] <- pop.mat[,paste0('Imm[1, ', t+1, ']')] / pop.mat[,paste0('Ntot[', t, ']')]
    immA[,t+1] <- pop.mat[,paste0('Imm[2, ', t+1, ']')] / pop.mat[,paste0('Ntot[', t, ']')]
  }
  
  Mu.immY <- rowMeans(immY, na.rm = T)
  Mu.immA <- rowMeans(immA, na.rm = T)
  
  # Collect posterior summaries of average VRs into a dataframe
  avgVR.data <- rbind(
    sam.summary(Mu.immY), sam.summary(Mu.immA),
    sam.summary(pop.mat[,'Mu.pB[1]']), sam.summary(pop.mat[,'Mu.pB[2]']),
    sam.summary(pop.mat[,'Mu.CS[1]']), sam.summary(pop.mat[,'Mu.CS[2]']),
    sam.summary(pop.mat[,'Mu.pNS']),
    sam.summary(pop.mat[,'Mu.sN[1]']), sam.summary(pop.mat[,'Mu.sN[2]']),
    sam.summary(pop.mat[,'Mu.sJ']), sam.summary(pop.mat[,'Mu.sA'])
  )
  colnames(avgVR.data) <- c('lCI', 'Median', 'uCI')
  avgVR.data <- as.data.frame(avgVR.data)
  
  # Add parameter information
  avgVR.data$Parameter <- c('Mu.immY', 'Mu.immA', 
                            'Mu.pB[1]', 'Mu.pB[2]',
                            'Mu.CS[1]', 'Mu.CS[2]',
                            'Mu.pNS',
                            'Mu.sN[1]', 'Mu.sN[2]',
                            'Mu.sJ', 'Mu.sA')
  avgVR.data$VitalRate <- c(rep(c('Immigration rate', 'Breeding probability',
                                'Clutch size'), each = 2),
                            'Nest success probability',
                            rep(c('Nestling survival', 'Annual survival'), each = 2))
  avgVR.data$AgeClassPlot <- as.factor(c(1,2,1,2,1,2,2,1,2,1,2))
  
  # Add population ID
  avgVR.data$PopID <- PopID
  
  # Return data
  return(avgVR.data)
  
}

## Get summary data for all populations
post.data <- data.frame()

for(i in 1:7){
  pop.data <- sum.avgVR(out.mat[[i]], PopIDs[i], PopTmax[i])
  post.data <- rbind(post.data, pop.data)
}


# Plotting: Vital rate posteriors, across population #
#----------------------------------------------------#

## Add a dummy covariate for stacking populations
dummy.stack <- data.frame(
  PopID = rep(c('EDM', 'TEI', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'), 2),
  AgeClassPlot = rep(c(1,2), each = 7),
  Spacing = c(1-c(0, 1:6*1/6)-0.05, 1-c(0, 1:6*1/6))
)

post.data <- merge(post.data, dummy.stack, by = c('PopID', 'AgeClassPlot'), all.x = T)


## Order PopIDs by latitude
post.data$PopID <- factor(post.data$PopID, levels = c('EDM', 'TEI', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT'))

## Order vital rates
post.data$VitalRate <- factor(post.data$VitalRate, levels = c('Breeding probability', 'Clutch size',
                                                              'Nest success probability', 'Nestling survival', 
                                                              'Annual survival', 'Immigration rate'))

## Define custom color code for PopIDs
PFC_ColorCode <- c('#B43AA5', '#F2309B', '#F23E1D', '#E7AA24', '#A5D85F', '#32A638', '#376BAD')



pdf('Plots/VR_Posteriors_crossPop_Forest1.pdf', width = 8, height = 7)
ggplot(post.data) + 
  geom_pointrangeh(aes(x = Median, y = Spacing, xmin = lCI, xmax = uCI, colour = PopID, shape = AgeClassPlot, linetype = AgeClassPlot), size = 0.3, fatten = 4) + 
  scale_color_manual(values = PFC_ColorCode) + 
  guides(colour = guide_legend(order = 1)) + 
  scale_shape_manual(name = 'Age Class', values = c(1,19), labels = c('Juvenile/yearling', 'Adult/combined')) + 
  scale_linetype_manual(name = 'Age Class', values = c(1,1), labels = c('Juvenile/yearling', 'Adult/combined')) + 
  facet_wrap(~VitalRate, ncol = 2, scales = 'free_x') + 
  xlab('Estimate') + 
  theme_classic() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(),
        strip.background = element_rect(color=NA, fill="grey80", size=1.5),
        legend.position = 'right', panel.spacing.y = unit(1, "lines"))
dev.off()


pdf('Plots/VR_Posteriors_crossPop_Forest2.pdf', width = 6*0.9, height = 7*0.9)
ggplot(post.data) + 
  geom_pointrange(aes(y = Median, x = Spacing, ymin = lCI, ymax = uCI, colour = PopID, shape = AgeClassPlot, linetype = AgeClassPlot), size = 0.275, fatten = 4) + 
  scale_color_manual(values = PFC_ColorCode) + 
  guides(colour = guide_legend(nrow = 1, order = 1)) + 
  scale_x_reverse() + 
  scale_shape_manual(name = 'Age Class', values = c(1,19), labels = c('Juvenile/yearling', 'Adult/combined')) + 
  scale_linetype_manual(name = 'Age Class', values = c(1,1), labels = c('Juvenile/yearling', 'Adult/combined')) +   facet_wrap(~VitalRate, ncol = 2, scales = 'free_y') + 
  xlab('Estimate') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_line(color = 'grey80'), axis.ticks.y = element_line(color = 'grey80'),
        strip.background = element_rect(color=NA, fill="grey80"),
        legend.position = 'bottom', legend.box = 'vertical', legend.title = element_blank())
dev.off()

