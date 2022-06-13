library(coda)
library(matrixStats)
library(reshape)
library(ggplot2)
library(crayon)

## Set data path
UserName <- 'chloe.nater'
#UserName <- 'chloern'
DataPath <- paste0('/Users/', UserName, '/Dropbox/PiedFlycatcher_IPM/IPM_Code/210819_FlycatcherIPM_PostSamples_MS1/')

## Make a list of all populations
PopID_List <- c('DIN', 'EDM', 'KAT', 'NAG', 'NWA', 'OKE', 'TEI')

## Set the study years for all populations
StudyYearsList <- list(
  DIN = 1980:2017,
  EDM = 1955:2020,
  KAT = 1974:2011,
  NAG = 1981:2019,
  NWA = 1986:2019,
  OKE = 1974:2019,
  TEI = 1981:2020
)

## Define the period for cross-population comparison
CompPeriod <- 1986:2011

## Define the indeces for the comparison period for each population
YearIndexList <- list(
  DIN = which(StudyYearsList$DIN %in% CompPeriod),
  EDM = which(StudyYearsList$EDM %in% CompPeriod),
  KAT = which(StudyYearsList$KAT %in% CompPeriod),
  NAG = which(StudyYearsList$NAG %in% CompPeriod),
  NWA = which(StudyYearsList$NWA %in% CompPeriod),
  OKE = which(StudyYearsList$OKE %in% CompPeriod),
  TEI = which(StudyYearsList$TEI %in% CompPeriod)
) 

## Determine how to treat immigration rates estimated at 0
# NOTE: Occasionally, there may be posterior samples that estimate immigration
#       rates at exactly 0 for entire time-periods. This is problematic for the
#       period-design LTRE, as calculations involve the log mean and log sd
#       difference of vital rates, and log(0) = Inf. 
#       Resulting contributions of immigration rates will default to NA, but 
#       contributions of other components for this posterior sample will still
#       have numerical values. This can potentially bias the posteriors of
#       relative contributions.
#       In the following, we present 2 options for avoiding such bias:

# Option 1:
FixLog0MeanRates <- FALSE

# Removes any samples which contain period-mean estimates of immigration rates
# equaling 0 completely prior to analysis.


## Run LTRE for all populations

for(x in 1:length(PopID_List)){
  PopID <- PopID_List[x]
  PopIDIdx <- x
  source('periodLTRE_crossPop_PopID.R')
}
