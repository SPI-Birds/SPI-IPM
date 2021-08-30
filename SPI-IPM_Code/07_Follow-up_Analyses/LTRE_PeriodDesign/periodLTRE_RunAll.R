library(coda)
library(matrixStats)
library(reshape)
library(ggplot2)
library(crayon)

## Set data path
UserName <- 'chloenater'
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

## Define the periods for comparison for all populations
PeriodsList <- list(
  DIN = matrix(c(1985:2000,  # decreasing
                 2001:2016), # increasing
               nrow = 2, byrow = T),
  EDM = matrix(c(1990:2004,  # decreasing
                 2005:2019), # increasing
               nrow = 2, byrow = T),
  KAT = matrix(c(1974:1987,  # increasing
                 1988:2001), # decreasing
               nrow = 2, byrow = T),
  NAG = matrix(c(1981:1994,  # increasing to stable
                 1995:2008), # decreasing 
               nrow = 2, byrow = T),
  NWA = matrix(c(1990:2003,  # decreasing
                 2004:2017), # increasing
               nrow = 2, byrow = T),
  OKE = matrix(c(1989:2003,  # decreasing
                 2004:2018), # decreasing (but note: fewer nestboxes)
               nrow = 2, byrow = T),
  TEI = matrix(c(1988:2003,  # stable
                 2004:2019), # increasing 
               nrow = 2, byrow = T)
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

# Option 2: 
FixLog0MeanRates <- TRUE
# Replaces all immigration rate means of 0 with the smallest non-0 individual 
# value estimated across all samples and time-steps prior to analysis
# (following advice by Warton and Hui, Ecology, 92(1), 2011, pp. 3–10).
# The same is done for immigration rate standard deviations of 0, but here the
# smallest non-0 value estimated over the entire study period across all samples
# is used. 

# For the cases tested so far, which option was used hardly affected posterior
# distributions of LTRE contributions at all. This is due to the fact that
# the number of samples with period-mean immigration rates of 0 is very low
# compared to the total number of samples. If such samples are more common,
# there may be non-negligible differences with the two options above.


## Run LTRE for all populations

for(i in 1:length(PopID_List)){
  PopID <- PopID_List[i]
  source('periodLTRE_PopID.R')
}
