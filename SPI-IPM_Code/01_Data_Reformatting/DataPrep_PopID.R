#*********************************************************************#
# Code for formatting SPI-Birds standard format data for IPM analysis #
#               Example: EAST DARTMOOR (Pied Flycatcher)              #
#*********************************************************************#
library(dplyr)
options(tibble.width = Inf)

## Load data in SPI-birds format
load('SPIBirds_StandardData.RData')
str(standard_data)

## Set target population 
tPopID <- 'EDM'

## Extract relevant data tables - Pied Flycatchers from selected population only
BroodData <- subset(standard_data$Brood_data, PopID == tPopID & Species == 'FICHYP')
CapData <- subset(standard_data$Capture_data, CapturePopID == tPopID & Species == 'FICHYP' & BreedingSeason <= max(BroodData$BreedingSeason))
IndData <- subset(standard_data$Individual_data, PopID == tPopID & Species %in% c('FICHYP', 'CCCCCC'))

if(tPopID == 'KAT'){
  CapData <- subset(CapData, !(CapturePlot %in% c('HRL', 'HRN')))
}
# NOTE: For KAT, capture data contains two locations that were not part of the nest monitoring --> exclude

## Make a list of each individual's age and age class (= life history stage) in each breeding season
table(CapData$Age_calculated, CapData$Age_observed)
# --> This shows the combinations of calculated age (rows) and observed age (columns)
#     There should be no instanced where one age == 1 and the other >= 4
#     (If such instances appear, they have to be checked in the raw data)

AgeTable <- CapData[,c('IndvID','BreedingSeason','Age_observed','Age_calculated')] %>%
  dplyr::distinct(IndvID, BreedingSeason, .keep_all=T) %>%
  dplyr::mutate(AgeClass = dplyr::case_when(Age_observed < 4 | Age_calculated < 4 ~ 'chick',
                                            Age_observed > 5 | Age_calculated > 5 ~ 'adult',  
                                            Age_observed %in% c(4,5) & Age_calculated == 5 ~ 'yearling',
                                            Age_observed %in% c(4,5) & Age_calculated == 4 ~ NA_character_
  ))

table(AgeTable$Age_calculated, AgeTable$Age_observed, AgeTable$AgeClass, useNA = "always")
# --> This shows how the different combinations of age information are assigned to an age class (chick, yearling, adult, NA)

# NOTE:
# The assumptions of the age class assignment are as follows:
# 1) Any individual either observed or calculated as EURING codes 1-3 is a chick
# 2) Any individual either observed or calculated as EURING code > 5 is an adult (= more than 1 year old)
# 3) Any individual observed as EURING code 4 or 5 and whose age is confirmed to be 5 by the calculation is a yearling
# 4) Any individual observed as EURING code 4 or 5 and whose age verified by the calculation (i.e Age_calculated = 4) is inconclusive, and age class therefore NA.



## Merge parent age class into BroodData
F_AgeTable <- AgeTable %>%
	dplyr::select(IndvID, BreedingSeason, AgeClass) %>%
	dplyr::rename(FemaleID = IndvID, FemaleAgeClass = AgeClass)

M_AgeTable <- AgeTable %>%
	dplyr::select(IndvID, BreedingSeason, AgeClass) %>%
	dplyr::rename(MaleID = IndvID, MaleAgeClass = AgeClass)

BroodData <- BroodData %>%
	dplyr::left_join(F_AgeTable, by = c('FemaleID', 'BreedingSeason')) %>%
	dplyr::left_join(M_AgeTable, by = c('MaleID', 'BreedingSeason'))

# NOTE: 
# This is data that - in the raw data - were already included in the nest files. 
# However, I re-match it here for the data in the SPI-birds format, so that this code works more generally.


#--------------#
# A) NEST DATA #
#--------------#

# A.1) Nest count
#----------------

## Check for mismatches and missing values in ClutchType assignment (raw data vs. pipeline)
print(subset(BroodData, ClutchType_observed != ClutchType_calculated)[,c('BroodID', 'ClutchType_observed', 'ClutchType_calculated')], n = 150)
# --> There are 53 entries with mismatches
# --> The majority of these are cases where the pipeline assigns "replacement" to a "first" clutch (likely based on laying date)

print(subset(BroodData, is.na(ClutchType_observed)), n = 150)
print(subset(BroodData, is.na(ClutchType_calculated)), n = 150)
# --> There are no missing values in ClutchType observed, but several for the calculated (usually linked to missing LayingDate and/or FemaleID)

# --> We will use ClutchType_observed (= as recorded in raw data) here

## Subset to exclude 2nd broods (replacement and second clutches)
FirstBroodData <- subset(BroodData, ClutchType_observed == 'first')	

## Count the number of broods = breeding females per year
FirstBroodCount <- FirstBroodData %>%
	dplyr::group_by(BreedingSeason) %>%
	dplyr::summarise(Count = n()) %>%
	dplyr::ungroup() %>%
	dplyr::arrange(BreedingSeason)


## Extract the range of years of study
MinYear <- min(BroodData$BreedingSeason)
MaxYear <- max(BroodData$BreedingSeason)

StudyYears <- c(MinYear:MaxYear)

## Double-check data is available for each year (1955-2020)
StudyYears %in% FirstBroodCount$BreedingSeason
StudyYears[which(!(StudyYears %in% FirstBroodCount$BreedingSeason))]
# --> Data is missing for the year 1971 (index = 17)

## Insert missing year into brood count data
FirstBroodCount <- merge(data.frame(BreedingSeason = StudyYears), FirstBroodCount, all.x = T)

## Extract vector NestCount[t]
NestCount <- FirstBroodCount$Count

## Make a vector to indicate status of nest/brood data collection (yes/no)
NS_Data <- rep(1, length(NestCount)) 
NS_Data[which(is.na(NestCount))] <- 0

## Set NestCount in missing year to 0
NestCount[which(is.na(NestCount))] <- 0


# NOTE:
# We base the nest count on first clutches only since we use it as a proxy for the number of breeding females in the IPM. For reproduction parameters below (clutch size, fledging success), on the other hand, we can use all clutches. 


# A.2) Clutch size - Population level
#------------------------------------

## Count the total number of eggs (and broods with available data)
EggCount <- BroodData %>%
	dplyr::filter(!is.na(ClutchSize_observed)) %>%
	dplyr::group_by(BreedingSeason) %>%
	dplyr::summarise(EggNo = sum(ClutchSize_observed),
					 NestNoData = n()) %>%
	dplyr::ungroup() %>%
	dplyr::arrange(BreedingSeason)

## Count the number of all broods (with and without clutch size data)
AllBroodCount <- BroodData %>%
	dplyr::group_by(BreedingSeason) %>%
	dplyr::summarise(NestNoTotal = n()) %>%
	dplyr::ungroup() %>%
	dplyr::arrange(BreedingSeason)

## Double-check data is available for each year (1955-2020)
StudyYears %in% EggCount$BreedingSeason
StudyYears[which(!(StudyYears %in% EggCount$BreedingSeason))]

StudyYears %in% AllBroodCount$BreedingSeason
StudyYears[which(!(StudyYears %in% AllBroodCount$BreedingSeason))]
# --> Again, 1971 is missing

## Add missing years to data
EggCount <- merge(data.frame(BreedingSeason = StudyYears), EggCount, all.x = T)
AllBroodCount <- merge(data.frame(BreedingSeason = StudyYears), AllBroodCount, all.x = T)

## Extract vectors EggNoTot[t] (= total number of eggs) and EggNoSP[t] (= proportion of all nests that EggNoTot[t] is based on)
EggNoTot <- EggCount$EggNo
EggNoSP <- EggCount$NestNoData / AllBroodCount$NestNoTotal

## Replace NA values in EggNoTot[t] and EggNoSP[t]
EggNoTot[which(is.na(EggNoTot))] <- 0
EggNoSP[which(is.na(EggNoSP))] <- 0


# A.3) Clutch size - Nest level
#------------------------------

## Subset the data to contain complete entries (clutch size and mother age class available)
ClutchSizeData <- subset(BroodData, !is.na(ClutchSize_observed) & !is.na(FemaleAgeClass))

## Extract vectors ClutchSize[x], CS_year[x], CS_FAge[x]
ClutchSize <- ClutchSizeData$ClutchSize_observed
CS_year <- ClutchSizeData$BreedingSeason - min(BroodData$BreedingSeason) + 1
CS_FAge <- ifelse(ClutchSizeData$FemaleAgeClass == 'yearling', 1, 2)


# A.4) Fledging success - Population level
#-----------------------------------------

## Count the total number of fledglings (and broods with available data)
FledgeCount <- BroodData %>%
	dplyr::filter(!is.na(NumberFledged_observed)) %>%
	dplyr::group_by(BreedingSeason) %>%
	dplyr::summarise(FledgeNo = sum(NumberFledged_observed),
					 NestNoData = n()) %>%
	dplyr::ungroup() %>%
	dplyr::arrange(BreedingSeason)


## Double-check data is available for each year (1955-2020)
StudyYears %in% FledgeCount$BreedingSeason
StudyYears[which(!(StudyYears %in% FledgeCount$BreedingSeason))]
# --> Again, 1971 is missing

## Add missing years to data
FledgeCount <- merge(data.frame(BreedingSeason = StudyYears), FledgeCount, all.x = T)

## Extract vectors FledgedTot[t] (= total number of fledglings) and FledgedSP[t] (= proportion of all nests that FledgedTot[t] is based on)
FledgedTot <- FledgeCount$FledgeNo
FledgedSP <- FledgeCount$NestNoData / AllBroodCount$NestNoTotal

## Replace NA values in FledgedTot[t] and FledgedSP[t]
FledgedTot[which(is.na(FledgedTot))] <- 0
FledgedSP[which(is.na(FledgedSP))] <- 0



# A.5) Fledging success - Nest level
#------------------------------------

## Subset the data to contain complete entries (egg number, fledgling number and mother age class available)
FledgeData <- subset(BroodData, !is.na(ClutchSize_observed) & !is.na(NumberFledged_observed) & !is.na(FemaleAgeClass))

## Further subset the data to remove any erroneous entries (i.e. where # Fledged > Clutch size)
FledgeData <- subset(FledgeData, NumberFledged_observed <= ClutchSize_observed)

## Extract vectors Fledged[x], Laid[x], F_year[x], F_FAge[x]
Fledged <- FledgeData$NumberFledged_observed
Laid <- FledgeData$ClutchSize_observed
F_year <- FledgeData$BreedingSeason - min(BroodData$BreedingSeason) + 1
F_FAge <- ifelse(FledgeData$FemaleAgeClass == 'yearling', 1, 2)

## Separate data into nest success (0 vs. 1+ fledglings) & number fledged (given success)
NoFledged <- Fledged[which(Fledged > 0)]
NoLaid <- Laid[which(Fledged > 0)]
NoF_year <- F_year[which(Fledged > 0)]
NoF_FAge <- F_FAge[which(Fledged > 0)]

anyFledged <- Fledged
anyFledged[which(Fledged > 0)] <- 1


# A.6) Proportion uncaptured breeding females
#--------------------------------------------

## Count the number of nests with un-identified (un-captured) breeding females
UKFCount <- BroodData %>%
	dplyr::filter(is.na(FemaleID)) %>%
	dplyr::group_by(BreedingSeason) %>%
	dplyr::summarise(UKFNo = n()) %>%
	dplyr::ungroup() %>%
	dplyr::arrange(BreedingSeason)

## Double-check data is available for each year (1955-2020)
StudyYears %in% UKFCount$BreedingSeason
# --> Several years are missing, meaning there were a) no un-identified females in those years or b) no sampling going on (1971)

## Merge into data on all nests and calculate proportion identified females
AllBroodCount <- AllBroodCount %>%
	dplyr::left_join(UKFCount, by = 'BreedingSeason') %>%
	dplyr::mutate(KFProp = dplyr::case_when(is.na(NestNoTotal) ~ 0, # with no sampling, the proportion of breeding females sampled is 0
											is.na(UKFNo) ~ 1, 
											!is.na(UKFNo) ~ (NestNoTotal-UKFNo)/NestNoTotal))

## Extract vector PropCapBrood[t]
PropCapBrood <- AllBroodCount$KFProp


#------------------------#
# B) MARK-RECAPTURE DATA #
#------------------------#

# B.1) CJS Capture histories
#---------------------------

## Merge 'calculated' (= summarised) sex into capture data
CapData <- dplyr::left_join(CapData, IndData[,c('IndvID', 'Sex_calculated')])

NoSexInd <- unique(subset(CapData, Sex_calculated == 'C' | Sex_observed != Sex_calculated)$IndvID)
NoSexInd
subset(CapData, IndvID %in%NoSexInd)
# --> There are no individuals with conflicting sex information

## Remove individuals with conflicting sex
CapData2 <- subset(CapData, !(IndvID%in%NoSexInd))

## Remove duplicate captures (same individual multiple times in same season)
CapData2 <- CapData2 %>%
	dplyr::distinct(IndvID, BreedingSeason, .keep_all=T) %>%

	## Assign a capture state code (life stage / sex combination)
	# 1 = chick
	# 2 = yearling/adult female
	# 3 = yearling/adult male
	dplyr::mutate(State = dplyr::case_when(Age_observed == 1 ~ 1,
										   Sex_calculated == 'F' ~ 2,
										   Sex_calculated == 'M' ~ 3),
				  ## Assign capture session (year starting from 1)
				  Session = BreedingSeason - min(BreedingSeason) + 1,
				  ## Assign a running number to each individual
				  IndvNo = as.numeric(as.factor(IndvID))) %>%
	
	
  ## Remove individuals/captures for which no state could be assigned
  dplyr::filter(!is.na(State)) %>%
  
	## Sort according to individual number, then session
	dplyr::arrange(IndvNo, Session)

# NOTE: 
# Adjustments have to be made in this section if there were any dead recoveries mixed in with life captures


## Function for transforming longitudinal data into multistate capture histories
MSMR.data = function(data, Session, IndNo, state, trait){
	# data = data frame with longitudinal data
	# Session = character string with name of capture session variable in data
	# IndNo = character string with name of individual running number variable in data
	# state = character string with name of individual state variable in data
	# trait = character string with name of individual trait measurement variable in data 
    
    data <- as.data.frame(data) # Convert data to data frame (if necessary)
    
   	Ncapture <- nrow(data)  # Count the number of capture events
   	Nanimal <- data[Ncapture, IndNo]  # Count the number of individuals
    
   	occasion <- data[, Session] - min(data[, Session]) + 1 # Reformat Session
    
   	CH <- matrix(0, Nanimal, max(occasion)) # Define matrix for storing capture histories
   	#traitH <- matrix(NA, Nanimal, max(occasion)) # Define matrix for storing trait observations
    
   	for (x in 1:Ncapture) {
     	CH[data[x, IndNo], occasion[x]] <- data[x, state]
     	#traitH[data[x, IndNo], occasion[x]] <- data[x, trait]
   	}
    
   	return(list(CH = CH#, 
   			   	#traitH = traitH
   			   	))
}

CH.output <- MSMR.data(CapData2, 'Session', 'IndvNo', 'State', NA)


## Remove all adult male captures
CHs <- CH.output$CH
CHs[which(CHs == 3)] <- 0

## Remove all "empty" capture histories 
## (= CHs that previously only contained captures of adult males or individuals whose state could not be assigned due to missing sex information)
CHs <- CHs[-which(rowSums(CHs) %in% c(0, NA)),]

## Reformat to CJS capture histories for IPM analyses (matrix y[i,t])
y <- CHs
y[which(CHs == 2)] <- 1


## Extract first capture and last occasion to analyse
first <- last <- rep(NA, dim(y)[1])

for(i in 1:dim(y)[1]){
	
	first[i] <- min(which(y[i,] == 1))
	
	last.capture <- max(which(y[i,] == 1))
	last[i] <- ifelse(last.capture+20 > dim(y)[2], dim(y)[2], last.capture+20)
}

## Check for potential recaptures in years when PropCapBrood = 0
y.recap <- y
for(i in 1:dim(y.recap)[1]){
  y.recap[i,first[i]] <- 0
  
}

SamplingInfo <- data.frame(StudyYears, PropCapBrood, NoCaptures = colSums(y.recap))
subset(SamplingInfo, PropCapBrood == 0 & NoCaptures > 0)
# In 1970, 1971, and 1998 there was one capture (per year) despite no adult rings having been recorded in the nest file
# (These will have to be removed to avoid violating IPM assumptions)

## Remove recaptures in years when PropCapBrood = 0
for(i in 1:dim(y)[1]){
  if(first[i] < dim(y)[2]){
    for(t in (first[i]+1):last[i]){
      if(PropCapBrood[t] == 0){
        y[i,t] <- 0
        CHs[i,t] <- 0
      }
    }
  }
}
  
  
## Make associated age class matrix (ageclass[i,t])
ageclass <- CHs
for(i in 1:dim(CHs)[1]){
	
	# ageclass = NA prior to ringing
	if(first[i] > 1){
		ageclass[i,1:(first[i]-1)] <- NA
	}
	
	# ageclass = 2 (yearling/adult) for all years after ringing
	if(first[i] < dim(CHs)[2]){
		ageclass[i,(first[i]+1):dim(CHs)[2]] <- 2	
	}
}

## Remove CH's of individuals ringed at the last occasion 
## (these hold no information for survival analysis)

RingedLast <- which(first == dim(CHs)[2])

y <- y[-RingedLast,]
ageclass <- ageclass[-RingedLast,]
first <- first[-RingedLast]
last <- last[-RingedLast]

CHs <- CHs[-RingedLast,]


# B.2) Summarised CJS capture histories
#--------------------------------------

## Identify and count unique capture histories
library(dplyr)
uniqueCHs <- as.matrix(as.data.frame(CHs) %>% group_by_all %>% count)
CHs.sum <- unname(uniqueCHs[,1:(dim(uniqueCHs)[2]-1)])
CHs.count <- unname(uniqueCHs[,dim(uniqueCHs)[2]])
dim(CHs)
dim(CHs.sum)
# --> The entire 12234 CH's (with 2 stages) can be summarised into just 296 unique capture histories

## Reformat to CJS capture histories for IPM analyses (matrix y.sum[i,t])
y.sum <- CHs.sum
y.sum[which(CHs.sum == 2)] <- 1

## Extract first capture and last occasion to analyse
first.sum <- last.sum <- rep(NA, dim(y.sum)[1])

for(i in 1:dim(y.sum)[1]){
	first.sum[i] <- min(which(y.sum[i,] == 1))
	
	last.capture <- max(which(y.sum[i,] == 1))
	last.sum[i] <- ifelse(last.capture+20 > dim(y.sum)[2], dim(y.sum)[2], last.capture+20)
}

## Make associated age class matrix (ageclass.sum[i,t])
ageclass.sum <- CHs.sum
for(i in 1:dim(CHs.sum)[1]){
	
	# ageclass = NA prior to ringing
	if(first.sum[i] > 1){
		ageclass.sum[i,1:(first.sum[i]-1)] <- NA
	}
	
	# ageclass = 2 (yearling/adult) for all years after ringing year
	if(first.sum[i] < dim(CHs.sum)[2]){
		ageclass.sum[i,(first.sum[i]+1):dim(CHs.sum)[2]] <- 2	
	}
}


# B.3) Age-specific m-arrays
#---------------------------

## Function to create m-arrays based on CJS-type capture histories (from BPA Chapter 7.10)
marray <- function(CH){
  
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0,ncol = n.occasions + 1, nrow = n.occasions)
  
  # Calculate number of released individuals at each time-period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  
  for (i in 1:nind){
    pos <- which(CH[i,] != 0)
    g <- length(pos)
    
    for (z in 1:(g - 1)){
      m.array[pos[z], pos[z+1]] <- m.array[pos[z], pos[z+1]] + 1
    }
  } 
  
  # Calculate number of individuals that are never recaptured
  for (t in 1:n.occasions){
    m.array[t, n.occasions+1] <- m.array[t,1] - sum(m.array[t, 2:n.occasions])
  }
  
  out <- m.array[1:(n.occasions-1), 2:(n.occasions+1)]
  return(out)
}


## Isolate adult captures and transform into m-array
adult.y <- y
adult.y[which(ageclass < 2)] <- 0

A.marray <- marray(adult.y)

## Extract ringing age class and isolate CH's of individuals ringed as chicks
RingAge <- rep(NA, dim(y)[1])
for(i in 1:dim(y)[1]){
  RingAge[i] <- ageclass[i, first[i]]
}

chick.y <- y[which(RingAge == 1), ]

## Isolate CH's of chicks that are never recaptured and transform into m-array
chick.y.N <- chick.y[which(rowSums(chick.y) == 1), ]
J.N.marray <- marray(chick.y.N)

## Isolate CH's of chicks recaptured as adults, and transform into m-array
chick.y.R <- chick.y[which(rowSums(chick.y) > 1), ]

for(i in 1:dim(chick.y.R)[1]){
  
  recap1 <- which(chick.y.R[i,] == 1)[2]
  if(recap1 < dim(chick.y.R)[2]){
    chick.y.R[i, (recap1+1):dim(chick.y.R)[2]] <- 0
  }
}

J.R.marray <- marray(chick.y.R)
J.R.marray[, dim(J.R.marray)[2]] <- 0
# NOTE: The last column (representing chicks never recaptured) is set to 0, 
#       because all of these have been recaptured and re-released as adults

## Summarise juvenile M-array
J.marray <- J.N.marray + J.R.marray



# B.3) Immigrant numbers
#-----------------------

## Identify all females ringed as adults 
ImmData <- CapData %>%
  dplyr::arrange(CapturePopID, IndvID, BreedingSeason, CaptureDate, CaptureTime) %>%
  dplyr::group_by(CapturePopID, IndvID) %>%
  
  dplyr::summarise(
    RingSeason = first(BreedingSeason),
    RingAge = ifelse(any(Age_calculated %in% c(1, 3)), "chick", ifelse(min(Age_calculated) == 2, NA_character_, "adult")), 
    .groups = "keep") %>%
  dplyr::rowwise() %>%
  dplyr::ungroup() %>%
  dplyr::left_join(IndData[,c('IndvID', 'Sex_calculated')], by = 'IndvID') %>%
  dplyr::filter(RingAge == 'adult' & Sex_calculated == 'F')

# NOTE: We are not using the information "RingAge" from individual data because
# that will not deal correctly with individuals marked as chicks in another population.

## Merge in age class information
ImmData <- ImmData %>%
	dplyr::mutate(BreedingSeason = RingSeason) %>%
	dplyr::left_join(AgeTable, by = c('IndvID', 'BreedingSeason'))

table(ImmData$AgeClass, useNA = 'ifany')
# --> 99 out of 706 immigrants could be assigned an age class
# --> Out of these, 42 were yearlings and 57 were adults
table(ImmData$AgeClass, ImmData$RingSeason, useNA = 'ifany')
# --> Known age immigrants are spread across the dataset (except for earlier years)

## Set vector ImmAgeProp[a] using the whole dataset (and assuming no time-dependence)
if(all(c('yearling', 'adult')%in%ImmData$AgeClass)){
  ImmAgeProps <- as.numeric(table(ImmData$AgeClass)[2:1]/sum(table(ImmData$AgeClass)[1:2]))
}else{
  if('adult' %in%ImmData$AgeClass){
    ImmAgeProps <- c(0,as.numeric(table(ImmData$AgeClass)[1]/sum(table(ImmData$AgeClass)[1])))
  }else{
    ImmAgeProps <- c(NA, NA)
  }
}

ImmPropTable <- table(ImmData$AgeClass, ImmData$RingSeason)
ImmAgeProps_t <- matrix(NA, nrow = 2, ncol = dim(ImmPropTable)[2])
for(t in 1:dim(ImmAgeProps_t)[2]){
  if(sum(ImmPropTable[,t]) > 0){
    if(dim(ImmPropTable)[1] > 1){
      ImmAgeProps_t[,t] <- ImmPropTable[,t]/sum(ImmPropTable[,t])
    }else{
      if(dimnames(ImmPropTable)[[1]] == "adult"){
        ImmAgeProps_t[2,t] <- ImmPropTable[,t]/sum(ImmPropTable[,t])
        ImmAgeProps_t[1,t] <- 1 - ImmAgeProps_t[2,t]
      }
      if(dimnames(ImmPropTable)[[1]] == "yearling"){
        ImmAgeProps_t[1,t] <- ImmPropTable[,t]/sum(ImmPropTable[,t])
        ImmAgeProps_t[2,t] <- 1 - ImmAgeProps_t[1,t]
      }
    }
  }
}

## Count the number of immigrants ringed each year
ImmCount <- ImmData %>%
	dplyr::group_by(RingSeason) %>%
	dplyr::summarise(ImmNo = n()) %>%
	dplyr::ungroup() %>%
	dplyr::arrange(RingSeason)

## Double-check data is available for each year (1955-2020)
StudyYears %in% ImmCount$RingSeason
StudyYears[which(!(StudyYears %in% ImmCount$RingSeason))]
# --> There are NAs for several years, and these may be due to a) no immigrants present or b) no ringing at all
# --> The affected years are 1962 (8), 1973 (19), 1979 (25), 1988 (34), 1989 (35), 2003 (49)

## Check in which years there was no ringing of adults
which(rowSums(A.marray) == 0)
cbind(StudyYears, c(rowSums(A.marray),NA), PropCapBrood)
# --> No adults were ringed in years 1973 (19), 1979 (25), 1988 (34), 1989 (35), 2003 (49), and this is consistent with nest data too (PropCapBrood)

CMR_Activity <- data.frame(RingSeason = StudyYears, Activity = 1)
CMR_Activity$Activity[which(rowSums(A.marray) == 0)] <- 0

## Adjust immigrant numbers depending on whether there was ringing activity
ImmCount <- dplyr::right_join(ImmCount, CMR_Activity, by = 'RingSeason') %>%
	dplyr::mutate(ImmNo_adj = dplyr::case_when(!is.na(ImmNo) ~ ImmNo,
												Activity == 1 ~ 0L)) %>%
	dplyr::arrange(RingSeason)
	

## Extract vector ImmNoObs[t]
ImmNoObs <- ImmCount$ImmNo_adj
ImmNoObs[which(is.na(ImmNoObs))] <- 0

## Extract vector AR_Data[t]
AR_Data <- ImmCount$Activity

#-------------------#
# C) DATA COLLATION #
#-------------------#

## Collate study years and year indices
Index <- c(1:length(StudyYears))
YearIndeces <- cbind(StudyYears, Index)
dimnames(YearIndeces)[[2]] <- c('Year', 'Index')

## Organise all data into a list
PFC.data <- list(y = y, CHs = CHs,
                 ageclass = ageclass,
                 first = first,
                 last = last,
                 y.sum = y.sum, CHs.sum = CHs.sum,
                 ageclass.sum = ageclass.sum,
                 first.sum = first.sum,
                 last.sum = last.sum,
                 CHs.count = CHs.count,	
                 A.marray = A.marray,
                 J.marray = J.marray,			 
                 NestCount = NestCount,				 
                 EggNoTot = EggNoTot,
                 EggNoSP = EggNoSP,
                 ClutchSize = ClutchSize,
                 CS_year = CS_year,
                 CS_FAge = CS_FAge,
                 FledgedTot = FledgedTot,
                 FledgedSP = FledgedSP,
                 Fledged = Fledged,
                 Laid = Laid,
                 F_year = F_year,
                 F_FAge = F_FAge,
                 anyFledged = anyFledged,
                 NoFledged = NoFledged,
                 NoLaid = NoLaid,
                 NoF_year = NoF_year,
                 NoF_FAge = NoF_FAge,
                 ImmNoObs = ImmNoObs,
                 ImmAgeProps = ImmAgeProps,
                 ImmAgeProps_t = ImmAgeProps_t,
                 NS_Data = NS_Data,
                 PropCapBrood = PropCapBrood,
                 AR_Data = AR_Data,
                 YearIndeces = YearIndeces)

## Save data
saveRDS(PFC.data, file = paste0(PopID, '_IPMData.rds'))
