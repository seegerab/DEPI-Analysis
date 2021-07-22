####################################################################################
# The goal of this script is to create .csv files of the data with the outliers removed, 
# as well as quantile normalization. These .csv files will be used in later R scripts to 
# create visualizations and calculate epistasis and selection coefficient values. This
# allows me to run this code only once and have concise scripts in the future. 
# 
# IN: MAPK_DEPI_data_070121.txt 
# OUT: Fitness_NormalizedValues_OutliersRemoved_07212021.csv
#      Fitness_NormalizedValues_OutliersIncluded_07212021.csv
####################################################################################

###Load Necessary Packages
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
library(preprocessCore)
library(Routliers)
library(tidyr)
library(ggthemes)
library(extrafont)
library(stringi)
### Read in the data
### Note: this is using the data updated July 1, 2021
data <- read.delim("C:/Users/Owner/Documents/Research/Shiu_Lab/Shiu_Lab_R/Data/MAPK_DEPI_data_070121.txt", header = TRUE)
#Make sure flat, row, column, and genotype are coded as factors
data$Flat<-as.factor(data$Flat)
data$Row<-as.factor(data$Row)
data$Column<-as.factor(data$Column)
data$Genotype <- as.factor(data$Genotype)

####################################################################################
# Outlier removal and summary
####################################################################################

### Create a unique identifier for each plant:
data$ID <- paste(data$Experiment, "_", data$Flat, "_", data$Number, sep = "")
### Find the ID's for the measured values that are already encoded as NA
MeasuredValue_NA_ID <- c(filter(data, is.na(data$SN))$ID,
filter(data, is.na(data$TSC))$ID,
filter(data, is.na(data$SPF))$ID)
### Melissa asked me to replace the NA values for DEPI1_2_92 with 0
MeasuredValue_NA_ID <- MeasuredValue_NA_ID[MeasuredValue_NA_ID != "DEPI1_2_92"]
### Remove this ID's from the data frame:
data <- data%>%
  filter(!(ID %in% MeasuredValue_NA_ID))
### Now, the only remained NA values should be DEPI1_2_92 for SPF and TSC
### Verify this:
sum(is.na(data$SPF))
sum(is.na(data$TSC))
### Replace these remaining NA values with 0 
data$SPF[is.na(data$SPF)] <- 0  
data$TSC[is.na(data$TSC)] <- 0  
### Create a data frame with the outliers removed
data_outliers_removed <- data%>%
  ###Group by experiment and genotype
  group_by(Experiment, Genotype)%>%
  ###For each of the measurements, replace the measured value with NA if the measured value is classified as an outlier
  mutate(SN = replace(SN, outliers_mad(SN, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  mutate(SPF = replace(SPF, outliers_mad(SPF, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  mutate(TSC = replace(TSC, outliers_mad(TSC, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  ungroup()%>%
  ###Rearrange the data by experiment and genotype
  arrange(Experiment, Genotype)
### Here are the unique plants that are classified as outliers
SN_outliers <- filter(data_outliers_removed, is.na(data_outliers_removed$SN))$ID
TSC_outliers <- filter(data_outliers_removed, is.na(data_outliers_removed$TSC))$ID
SPF_outliers <- filter(data_outliers_removed, is.na(data_outliers_removed$SPF))$ID
### Create a text file for this output:
n <- max(length(SN_outliers), length(TSC_outliers), length(SPF_outliers))
length(SN_outliers) <- n
length(TSC_outliers) <- n
length(SPF_outliers) <- n
cbind(SPF_outliers, SN_outliers, TSC_outliers)
write.csv(cbind(SPF_outliers, SN_outliers, TSC_outliers), 
            file = "Outliers_Per_Measurement.csv")
###Summarize the proportion of measured values that are classified as outliers for TSC
data_outliers_removed%>%
  group_by(Experiment)%>%
  summarize(PropOutlier = sum(is.na(TSC)) / length(TSC))
###Summarize the proportion of measured values that are classified as outliers for SPF
data_outliers_removed%>%
  group_by(Experiment)%>%
  summarize(PropOutlier = sum(is.na(SPF)) / length(SPF))
###Summarize the proportion of measured values that are classified as outliers for SN
data_outliers_removed%>%
  group_by(Experiment)%>%
  summarize(PropOutlier = sum(is.na(SN)) / length(SN))

####################################################################################
# Quantile normalization - flat level
####################################################################################
#Get a list of the double mutants in this data frame:
all_double_mutants = list()
for (gen in levels(data$Genotype)) {
  if (str_detect(gen, "_") == T) {
    all_double_mutants = c(all_double_mutants, gen)
  }
}
###This function allows us to create a data frame where the columns are different lengths because there are a different number of plants per flat
###Add NA to the end of each column until each column is the same length
addToDF <- function(df, v){
  nRow <- nrow(df)
  lngth <- length(v)
  if(nRow > lngth){
    length(v) <- nRow
  }else if(nRow < lngth){
    df[(nRow+1):lngth, ] <- NA
  }
  cbind(df,v)
}
###Create a loop for each measurement and experiment
for (i in c("SN", "TSC", "SPF")){
  ###Set tmp.name to be the measurement
  tmp.name <- i
  ###For each of the experiments:
  for (j in c("DEPI1", "DEPI2", "DEPI3")){
    ###Filter to the data for the experiment
    experiment.data <- filter(data_outliers_removed, Experiment == j)
    ###Initialize an empty data frame
    temp.df <- data.frame()
    ###Loop through each flat 
    for (k in 1:length(unique(experiment.data$Flat))){
      ### Create a temporary vector of the measured values for each flat
      flat.measurements <- experiment.data%>%
        ### Filter to flat k
        filter(Flat == k)%>%
        ### Select the column with measurement i
        select(i)
      ### Add the flat.measurements to temp.df
      ### Note as.data.frame is needed because temp.vector is technically a tibble
      temp.df <- addToDF(temp.df, as.data.frame(flat.measurements)[,1])
      }
    ###Normalize across the flat (temp_normalize is a data frame, but temp.df needs to be converted to a matrix to use the normalize.quantiles function)
    temp_normalize <- as.data.frame(normalize.quantiles(as.matrix(temp.df)))%>%
      ###Add columns with the measurement and experiment to verify flat and experiment information
      mutate(measurement = i, experiment= j)
    ###Create a name to give the normalized data based on the measurement and experiment
    temp_name <- paste(tolower(j), "_", tmp.name , "_normalize", sep = "")
    ###Rename the columns, depending on the number of flats per experiment
    ###3 flats
    if(ncol(temp_normalize) == 5){
      temp_normalize <- temp_normalize%>%
        rename(flat_1 = V1, flat_2 = V2, flat_3 = V3)
    ###4 flats
    }else if (ncol(temp_normalize) == 6){
      temp_normalize <- temp_normalize%>%
        rename(flat_1 = V1, flat_2 = V2, flat_3 = V3, flat_4 = V4)
    ###5 flats
    }else{
      temp_normalize <- temp_normalize%>%
        rename(flat_1 = V1, flat_2 = V2, flat_3 = V3, flat_4 = V4, flat_5 = V5)}
    ###Assign the name to the data frame
    assign(temp_name, temp_normalize)}
}

####################################################################################
# Verify Quantile Normalization is Correct!
####################################################################################

### The ultimate goal is to replace the measured values with the quantile normalized measured values.
### First, start with a test so I know this is the correct method to use. 

###Filter to only include Flat 2 and DEPI1 from the data frame with the outliers removed 
test <- subset(data_outliers_removed, Flat == 2 & Experiment == "DEPI1")
###Add a column with the normalized value. These normalized values are the second column, with the string of NA values removed from the column
test$NormalizedSN <- depi1_SN_normalize[,2][1:132]
###Add a column with the count to the data frame
test$Count <- 1:nrow(test)
###Sort first by quantile normalized value, and then by the un-normalized value. Because the order is maintained when we quantile normalized, this should be the same data frame!
order1 <- test%>%
  arrange(SN)
order2 <- test%>%
  arrange(NormalizedSN)
### order1 == order2 is a dataframe of boolean values. If the data frames match, all values will be TRUE, so the sum of each column will be the same - the number of rows in the data frame. 
apply(order1 == order2, 2, sum)
### This looks good! Apply the same approach to the entire data frame. First, divide the data by flat and experiment, do the above approach, and re-combine the data.

####################################################################################
# Create final data frame with quantile normalized values
####################################################################################

###Reshape the data, so instead of having a column for each measurement, we have a column with the measurement type - either SN, TSC, or SPF
data_outliers_removed <- data_outliers_removed%>%
  gather(Measurement, Measured_Value, c("SN", "TSC", "SPF"))%>%
  select(-sSPF)
###For each experiment and flat, remove NA values from the data and add a column with the normalized values to the data
### Initialize an empty data frame for the DEPI1 experiment for the SN measurement
DEPI1_SN <- data.frame()
### For each flat in DEPI1
for (i in c(1:4)){
  ### Filter to the ith flat and the DEPI1 experiment
  ###Remove rows that have NA (remember, these are the NA values appended to the end of the rows in order to combine flats that have a different number of plants)
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
  ### Add a column with the normalized value from the data frame that we created in the loop above
  tempFlat$NormalizedMeasuredValue <- na.omit(depi1_SN_normalize[,i])
  DEPI1_SN <- rbind(DEPI1_SN, tempFlat)
}

DEPI1_TSC <- data.frame()
for (i in 1:4){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "TSC")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi1_TSC_normalize[,i])
  DEPI1_TSC <- rbind(DEPI1_TSC, tempFlat)
}
DEPI1_SPF <- data.frame()
for (i in 1:4){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "SPF")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi1_SPF_normalize[,i])
  DEPI1_SPF <- rbind(DEPI1_SPF, tempFlat)
}
DEPI2_SN <- data.frame()
for (i in 1:3){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI2", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi2_SN_normalize[,i])
  DEPI2_SN <- rbind(DEPI2_SN, tempFlat)
}
DEPI2_TSC <- data.frame()
for (i in 1:3){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI2", Measurement == "TSC")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi2_TSC_normalize[,i])
  DEPI2_TSC <- rbind(DEPI2_TSC, tempFlat)
}
DEPI2_SPF <- data.frame()
for (i in 1:3){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI2", Measurement == "SPF")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi2_SPF_normalize[,i])
  DEPI2_SPF <- rbind(DEPI2_SPF, tempFlat)
}
DEPI3_SN <- data.frame()
for (i in 1:5){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI3", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi3_SN_normalize[,i])
  DEPI3_SN <- rbind(DEPI3_SN, tempFlat)
}
DEPI3_TSC <- data.frame()
for (i in 1:5){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI3", Measurement == "TSC")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi3_TSC_normalize[,i])
  DEPI3_TSC <- rbind(DEPI3_TSC, tempFlat)
}
DEPI3_SPF <- data.frame()
for (i in 1:5){
  tempFlat <- data_outliers_removed%>%
    filter(Flat == i, Experiment == "DEPI3", Measurement == "SPF")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi3_SPF_normalize[,i])
  DEPI3_SPF <- rbind(DEPI3_SPF, tempFlat)
}
###Combine all of these data frames into one with all normalized values
NormalizedData <- rbind(DEPI1_SN,
                        DEPI1_SPF,
                        DEPI1_TSC,
                        DEPI2_SN,
                        DEPI2_SPF,
                        DEPI2_TSC,
                        DEPI3_SN,
                        DEPI3_SPF,
                        DEPI3_TSC
                        )
### SAVE THESE NORMALIZED VALUES TO USE IN LATER SCRIPTS. REMEMBER, THIS IS WITH REMOVING THE OUTLIERS
write.csv(NormalizedData, file = "Fitness_NormalizedValues_OutliersRemoved_07212021.csv")

####################################################################################
####################################################################################
### Repeat the above process, WITHOUT removing the outliers
####################################################################################
####################################################################################

###Create a loop for each measurement and experiment
for (i in c("SN", "TSC", "SPF")){
  ###Set tmp.name to be the measurement
  tmp.name <- i
  ###For each of the experiments:
  for (j in c("DEPI1", "DEPI2", "DEPI3")){
    ###Filter to the data for the experiment
    experiment.data <- filter(data, Experiment == j)
    ###Initialize an empty data frame
    temp.df <- data.frame()
    ###Loop through each flat 
    for (k in 1:length(unique(experiment.data$Flat))){
      ### Create a temporary vector of the measured values for each flat
      flat.measurements <- experiment.data%>%
        ### Filter to flat k
        filter(Flat == k)%>%
        ### Select the column with measurement i
        select(i)
      ### Add the flat.measurements to temp.df
      ### Note as.data.frame is needed because temp.vector is technically a tibble
      temp.df <- addToDF(temp.df, as.data.frame(flat.measurements)[,1])
    }
    ###Normalize across the flat (temp_normalize is a data frame, but temp.df needs to be converted to a matrix to use the normalize.quantiles function)
    temp_normalize <- as.data.frame(normalize.quantiles(as.matrix(temp.df)))%>%
      ###Add columns with the measurement and experiment to verify flat and experiment information
      mutate(measurement = i, experiment= j)
    ###Create a name to give the normalized data based on the measurement and experiment
    temp_name <- paste(tolower(j), "_", tmp.name , "_normalize", sep = "")
    ###Rename the columns, depending on the number of flats per experiment
    ###3 flats
    if(ncol(temp_normalize) == 5){
      temp_normalize <- temp_normalize%>%
        rename(flat_1 = V1, flat_2 = V2, flat_3 = V3)
      ###4 flats
    }else if (ncol(temp_normalize) == 6){
      temp_normalize <- temp_normalize%>%
        rename(flat_1 = V1, flat_2 = V2, flat_3 = V3, flat_4 = V4)
      ###5 flats
    }else{
      temp_normalize <- temp_normalize%>%
        rename(flat_1 = V1, flat_2 = V2, flat_3 = V3, flat_4 = V4, flat_5 = V5)}
    ###Assign the name to the data frame
    assign(temp_name, temp_normalize)}
}

####################################################################################
# Verify Quantile Normalization is Correct!
####################################################################################

### The ultimate goal is to replace the measured values with the quantile normalized measured values.
### First, start with a test so I know this is the correct method to use. 

###Filter to only include Flat 2 and DEPI1 from the data frame without the outliers removed 
test <- subset(data, Flat == 2 & Experiment == "DEPI1")
###Add a column with the normalized value. These normalized values are the second column, with the string of NA values removed from the column
test$NormalizedSN <- depi1_SN_normalize[,2][1:132]
###Add a column with the count to the data frame
test$Count <- 1:nrow(test)
###Sort first by quantile normalized value, and then by the un-normalized value. Because the order is maintained when we quantile normalized, this should be the same data frame!
order1 <- test%>%
  arrange(SN)
order2 <- test%>%
  arrange(NormalizedSN)
### order1 == order2 is a dataframe of boolean values. If the data frames match, all values will be TRUE, so the sum of each column will be the same - the number of rows in the data frame. 
apply(order1 == order2, 2, sum)
### This looks good! Apply the same approach to the entire data frame. First, divide the data by flat and experiment, do the above approach, and re-combine the data.

####################################################################################
# Create final data frame with quantile normalized values
####################################################################################

###Reshape the data, so instead of having a column for each measurement, we have a column with the measurement type - either SN, TSC, or SPF
data <- data%>%
  gather(Measurement, Measured_Value, c("SN", "TSC", "SPF"))%>%
  select(-sSPF)
###For each experiment and flat, remove NA values from the data and add a column with the normalized values to the data
### Initialize an empty data frame for the DEPI1 experiment for the SN measurement
DEPI1_SN <- data.frame()
### For each flat in DEPI1
for (i in c(1:4)){
  ### Filter to the ith flat and the DEPI1 experiment
  ###Remove rows that have NA (remember, these are the NA values appended to the end of the rows in order to combine flats that have a different number of plants)
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
  ### Add a column with the normalized value from the data frame that we created in the loop above
  tempFlat$NormalizedMeasuredValue <- na.omit(depi1_SN_normalize[,i])
  DEPI1_SN <- rbind(DEPI1_SN, tempFlat)
}

DEPI1_TSC <- data.frame()
for (i in 1:4){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "TSC")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi1_TSC_normalize[,i])
  DEPI1_TSC <- rbind(DEPI1_TSC, tempFlat)
}
DEPI1_SPF <- data.frame()
for (i in 1:4){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "SPF")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi1_SPF_normalize[,i])
  DEPI1_SPF <- rbind(DEPI1_SPF, tempFlat)
}
DEPI2_SN <- data.frame()
for (i in 1:3){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI2", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi2_SN_normalize[,i])
  DEPI2_SN <- rbind(DEPI2_SN, tempFlat)
}
DEPI2_TSC <- data.frame()
for (i in 1:3){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI2", Measurement == "TSC")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi2_TSC_normalize[,i])
  DEPI2_TSC <- rbind(DEPI2_TSC, tempFlat)
}
DEPI2_SPF <- data.frame()
for (i in 1:3){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI2", Measurement == "SPF")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi2_SPF_normalize[,i])
  DEPI2_SPF <- rbind(DEPI2_SPF, tempFlat)
}
DEPI3_SN <- data.frame()
for (i in 1:5){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI3", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi3_SN_normalize[,i])
  DEPI3_SN <- rbind(DEPI3_SN, tempFlat)
}
DEPI3_TSC <- data.frame()
for (i in 1:5){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI3", Measurement == "TSC")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi3_TSC_normalize[,i])
  DEPI3_TSC <- rbind(DEPI3_TSC, tempFlat)
}
DEPI3_SPF <- data.frame()
for (i in 1:5){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI3", Measurement == "SPF")
  tempFlat <- na.omit(tempFlat)
  tempFlat$NormalizedMeasuredValue <- na.omit(depi3_SPF_normalize[,i])
  DEPI3_SPF <- rbind(DEPI3_SPF, tempFlat)
}
###Combine all of these data frames into one with all normalized values
NormalizedData <- rbind(DEPI1_SN,
                        DEPI1_SPF,
                        DEPI1_TSC,
                        DEPI2_SN,
                        DEPI2_SPF,
                        DEPI2_TSC,
                        DEPI3_SN,
                        DEPI3_SPF,
                        DEPI3_TSC
)
### SAVE THESE NORMALIZED VALUES TO USE IN LATER SCRIPTS. REMEMBER, THIS IS WITHOUT REMOVING THE OUTLIERS
write.csv(NormalizedData, file = "Fitness_NormalizedValues_OutliersIncluded_07212021.csv")

