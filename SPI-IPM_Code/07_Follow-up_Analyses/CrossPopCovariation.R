#*******************************#
#  Integrated Population Model  #
# Pied Flycatcher (PiedFlyNet)  #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(corrplot)

# Loading posteriors #
#--------------------#

## Load posterior samples from all 7 runs
load('SPI-IPM_DIN.RData')
DIN.IPM <- PFC.IPM

load('SPI-IPM_EDM.RData')
EDM.IPM <- PFC.IPM

load('SPI-IPM_KAT.RData')
KAT.IPM <- PFC.IPM

load('SPI-IPM_NAG.RData')
NAG.IPM <- PFC.IPM

load('SPI-IPM_NWA.RData')
NWA.IPM <- PFC.IPM

load('SPI-IPM_OKE.RData')
OKE.IPM <- PFC.IPM

load('SPI-IPM_TEI.RData')
TEI.IPM <- PFC.IPM


# Re-organizing data #
#--------------------#

## Collect MCMC chains in a list of matrices
out.mat <- list(
  TEI = as.matrix(TEI.IPM),
  EDM = as.matrix(EDM.IPM),
  OKE = as.matrix(OKE.IPM),
  NAG = as.matrix(NAG.IPM),
  DIN = as.matrix(DIN.IPM),
  NWA = as.matrix(NWA.IPM),
  KAT = as.matrix(KAT.IPM)
)

PopIDs <- c('TEI', 'EDM', 'OKE', 'NAG', 'DIN', 'NWA', 'KAT')


# Extracting summaries for time-dependent parameters #
#----------------------------------------------------#

## Set the study years for each population
StudyYears <- list(
  TEI = 1981:2020,
  EDM = 1955:2020,
  OKE = 1974:2019,
  NAG = 1981:2019,
  DIN = 1980:2017,
  NWA = 1986:2019,
  KAT = 1974:2011
)

## Extract Tmax for all populations
PopTmax <- c(length(StudyYears$TEI),
             length(StudyYears$EDM),
             length(StudyYears$OKE),
             length(StudyYears$NAG),
             length(StudyYears$DIN),
             length(StudyYears$NWA),
             length(StudyYears$KAT)
             )

## Make an empty list to contain estimates for each population
sum.data <- list(
  TEI = NA,
  EDM = NA,
  OKE = NA,
  NAG = NA,
  DIN = NA,
  NWA = NA,
  KAT = NA
)

## Define a short-cut function for posterior summaries
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
  immY <- immA <- Ntot
  
  # Make posterior summaries of relevant quantities
  for(t in 1:PopTmax[i]){
    
    Ntot[t,] <- sam.summary(sam.mat[,paste0('Ntot[', t, ']')])
    if(t < PopTmax[i]){
      Btot[t,] <- sam.summary(sam.mat[,paste0('Btot[', t, ']')])
      immY[t+1,] <- sam.summary(sam.mat[,paste0('Imm[1, ', t+1, ']')]/sam.mat[,paste0('Ntot[', t, ']')])
      immA[t+1,] <- sam.summary(sam.mat[,paste0('Imm[2, ', t+1, ']')]/sam.mat[,paste0('Ntot[', t, ']')])
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
    Parameter = rep(c('Ntot', 'Btot', 'linkvar.sJ', 'linkvar.sA', 'linkvar.pB', 'linkvar.CS', 'linkvar.pNS', 'linkvar.sN', 'immY', 'immA'), each = PopTmax[i])
  )
  
  pop.sum.data <- cbind(
    pop.sum.data, 
    rbind(Ntot, Btot, linkvar.sJ, linkvar.sA, linkvar.pB, linkvar.CS, linkvar.pNS, linkvar.sN, immY, immA))
  
  # Reformat
  pop.sum.data <- pop.sum.data[,c('Year', 'Parameter', 'Median')]
  colnames(pop.sum.data)[3] <- PopIDs[i]
  
  
  # Insert summarised data into list
  sum.data[[i]] <- pop.sum.data
}

## Match dataframes by year
allPop.data <- plyr::join_all(sum.data, by = c('Year', 'Parameter'), type = 'full')


# Plotting relationships
#-----------------------

## Set-up for correlation plot matrix
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use = 'complete.obs'), digits=2)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt, cex = 1)
}

upper.panel<-function(x, y){
  points(x,y, pch = 19, col = viridis(4)[2])
}

## Set parameters to plot
PlotParams <- c('Ntot', 'Btot', 'linkvar.sJ', 'linkvar.sA', 'linkvar.pB', 'linkvar.CS', 'linkvar.pNS', 'linkvar.sN', 'immY', 'immA')

for(i in 1:length(PlotParams)){
  
  data <- subset(allPop.data, Parameter == PlotParams[i])[,-2]
  
  pdf(paste0('Plots/Covariation/MedianEstimates_Corr1_', PlotParams[i], '.pdf'), width = 8, height = 8)
  print(pairs(data, 
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        main = PlotParams[i]))
  dev.off()
  
  pdf(paste0('Plots/Covariation/MedianEstimates_Corr2_JointOverlap', PlotParams[i], '.pdf'), width = 6, height = 6)
  # print(corrplot(cor(data, use = 'complete.obs'), type = 'upper', tl.col="black"))
  # print(corrplot(cor(data, use = 'complete.obs'), sig.level = 0.05, insig = 'blank', type = 'upper', p.mat = cor.mtest(data)$p, tl.col="black"))
  print(corrplot(cor(data, use = 'complete.obs'), 
           sig.level = 0.05, insig = "label_sig", pch.col = "white", 
           diag = FALSE,
           type = 'upper', p.mat = cor.mtest(data)$p, 
           tl.col="black", method = 'color'))  
  dev.off()
  
  pdf(paste0('Plots/Covariation/MedianEstimates_Corr2_PairwiseOverlap', PlotParams[i], '.pdf'), width = 6, height = 6)
  # print(corrplot(cor(data, use = 'complete.obs'), type = 'upper', tl.col="black"))
  # print(corrplot(cor(data, use = 'complete.obs'), sig.level = 0.05, insig = 'blank', type = 'upper', p.mat = cor.mtest(data)$p, tl.col="black"))
  print(corrplot(cor(data, use = 'pairwise.complete.obs'), 
                 sig.level = 0.05, insig = "label_sig", pch.col = "white", 
                 diag = FALSE,
                 type = 'upper', p.mat = cor.mtest(data)$p, 
                 tl.col="black", method = 'color'))  
  dev.off()
}
