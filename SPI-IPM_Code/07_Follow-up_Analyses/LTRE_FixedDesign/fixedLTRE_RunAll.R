library(coda)
library(matrixStats)
library(reshape)
library(ggplot2)
library(crayon)

## Set data path
DataPath <- '...' # Set your data path

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

## Run LTRE for all populations
for(i in 1:length(PopID_List)){
  PopID <- PopID_List[i]
  source('fixedLTRE_PopID.R')
}
