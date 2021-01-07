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

data <- read.csv("MAPK_DEPI_fitness_data_110920.csv", header = TRUE)
###The "Experiment" column name reads in strange, with weird characters at the begining
###Rename this column:
colnames(data)[1] <- "Experiment"
#Make sure flat, row, column, and genotype are coded as categorical variables
data$Flat<-as.factor(data$Flat)
data$Row<-as.factor(data$Row)
data$Column<-as.factor(data$Column)
data$Genotype <- as.factor(data$Genotype)


data <- data%>%
  ###Group by experiment and genotype
  group_by(Experiment, Genotype)%>%
  ###For each of the measurements, replace the measured value with NA if the measured value is classified as an outlier
  mutate(SN = replace(SN, outliers_mad(SN, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  mutate(SPF = replace(SPF, outliers_mad(SPF, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  mutate(TSC = replace(TSC, outliers_mad(TSC, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  ungroup()%>%
  ###Rearrange the data by experiment and genotype
  arrange(Experiment, Genotype)

###Here, summarize the proportion of measured values that are classified as outliers for each measurement:
###TSC
data%>%
  group_by(Experiment)%>%
  summarize(PropOutlier = sum(is.na(TSC)) / length(TSC))
###SPF
data%>%
  group_by(Experiment)%>%
  summarize(PropOutlier = sum(is.na(SPF)) / length(SPF))
###SN
data%>%
  group_by(Experiment)%>%
  summarize(PropOutlier = sum(is.na(SN)) / length(SN))

###Generally, DEPI2 has more outliers than the other two experiments
###Leaning towards not showing SPF - focus on SN and TSC
###See which ones are thrown out

#Get a list of the double mutants in this data frame:
all_double_mutants = list()
for (gen in levels(data$Genotype)) {
  if (str_detect(gen, "_") == T) {
    all_double_mutants = c(all_double_mutants, gen)
  }
}

###This function allows us to create a data frame where the columns are different lengths
###This is needed, because there are a different number of plants per flat
###To achieve this, add NA to the end of each column until each column is the same length
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
for (i in 8:10){
  ###Set the name to be either SN, SPF, or TSC, depending on the column
  tmp.name <- colnames(data)[i]
  ###For each of the experiments:
  for (j in c("DEPI1", "DEPI2", "DEPI3")){
    ###Filter to select each specific month
    temp_vector <- filter(data, Experiment == j)
    ###Initialize an empty data frame
    temp_df <- data.frame()
    ###Loop through each flat 
    for (k in 1:length(unique(temp_vector$Flat))){
      ###Create a temporary data frame - each column is the measured values for each flat
      temp <- filter(temp_vector, Flat == k)[,i]
      temp <- as.data.frame(temp)[,1]
      temp_df <- addToDF(temp_df, temp)}
    
    ###Normalize across the flats
    temp_normalize <- as.data.frame(normalize.quantiles(as.matrix(temp_df)))%>%
      ###Add columns with the measurement and experiment to be certain there hasn't been any mix-ups
      mutate(col.number = i, experiment= j)
    ###Create a name to give the normalized data based on the measurement and experiment
    temp_name <- paste(tolower(j), "_", tmp.name , "_normalize", sep = "")
    ###Rename the columns
    
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

###Replace the measured values with the quantile normalized measured values:
###First, start with a test so I know this is the correct method

###Filter to only include Flat 2 and DEPI1
test <- subset(data, Flat == 2 & Experiment == "DEPI1")
###Add a column with the normalized value
###These normalized values are the second column, with the string of NA
###values removed from the column
test$NormalizedSN <- depi1_SN_normalize[,2][1:132]
###Add a column with the count
test$Count <- 1:nrow(test)
###Sort first by quantile normalized value, and then by the un-normalized value
###Because the order is maintained when we quantile normalized, this should be the same data frame!

order1 <- test%>%
  arrange(SN)
order2 <- test%>%
  arrange(NormalizedSN)

###Plan:
###Divide the data by Flat and Experiment
###Do the above approach
###Re-combine the data

###Reshape the data, so instead of having a column for each measurement, we
###have a column with the measurement type - either SN, TSC, or SPF
data <- data%>%
  gather(Measurement, Measured_Value, 8:10)

###For each experiment and flat:

###Remove NA values from the data, and add a column with the normalized
###values to the data

DEPI1_SN <- data.frame()
for (i in 1:4){
  tempFlat <- data%>%
    filter(Flat == i, Experiment == "DEPI1", Measurement == "SN")
  tempFlat <- na.omit(tempFlat)
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

###Visualize the normalizations to make sure the quantile normalized
###flats each have the same distribution for all experiments



for (i in c("DEPI1", "DEPI2", "DEPI3")){
  for (m in c("TSC", "SN", "SPF")){
  plotData <- NormalizedData%>%
    filter(Measurement == m , Experiment == i)
  title1 <- paste(i, m, "- Before Normalization", sep = " ")
  title2 <- paste(i, m, "- After Normalization", sep = " ")
  plot1 <- ggplot(data = plotData, aes(x = Measured_Value, color = Flat))+
    geom_density()+
    stat_density(geom = "line", position = "identity")+
    theme_linedraw()+
    labs(x = m,
         y = "Density",
         title = title1)
  plot2 <- ggplot(data = plotData, aes(x = NormalizedMeasuredValue, color = Flat))+
    geom_density()+
    stat_density(geom = "line", position = "identity")+
    theme_linedraw()+
    labs(x = paste(m ,"Normalized by Flat", sep = " "),
         y = "Density",
         title = title2)
  grid.arrange(plot1, plot2)
  ggsave(file=QuantFileName, arrangeGrob(plot1, plot2), width = 7, height = 7)
  } 
}

###Something is wrong with SN! Talk to Melissa about this...

###Investigate why SN quantile normalization isn't correct

###Epistasis Calculations:
###Initialize an empty data frame to populate with information:
geneticInteractions <- data.frame(DoubleMutant = rep(NA, 0), 
                                     MutantA = rep(NA, 0),
                                     MutantB = rep(NA, 0), 
                                     AdditiveEpistasis = rep(NA, 0),
                                     ProportionalEpistatis = rep(NA, 0),
                                     Experiment = rep(NA, 0),
                                     Measurement = rep(NA, 0))

###Loop through each experiment and measurement:
for(e in c("DEPI1", "DEPI2", "DEPI3")){
  for (m in c("SN", "TSC", "SPF")){
    ###Filter to each specific experiment and measurement
    tempData <- filter(NormalizedData, Experiment == e, Measurement == m)
    ###Create an empty data frame to fill with the information and calcuations:
    geneticInteractionsTmp <- data.frame(DoubleMutant = rep(NA, 25), 
                                      MutantA = rep(NA, 25),
                                      MutantB = rep(NA, 25), 
                                      AdditiveEpistasis = rep(NA, 25),
                                      ProportionalEpistatis = rep(NA, 25),
                                      Experiment = rep(NA, 25),
                                      Measurement = rep(NA, 25))
    ###Initialize a row count to use to populate the data frame
    rowCount <- 1  
    ###For each of the double mutants:
    for (dm in unlist(all_double_mutants)){
      ###Extract the single mutants from the double mutant
      ma <- unlist(strsplit(dm, "_"))[1]
      mb <- paste("mpk", unlist(strsplit(dm, "_"))[2], sep = "")
      ###Calculate the fitness of the dm, ma, mb, and wt
      fdm <- mean(filter(tempData, Genotype == dm)$NormalizedMeasuredValue)
      fwt <- mean(filter(tempData, Genotype == "Col")$NormalizedMeasuredValue)
      fma <- mean(filter(tempData, Genotype == ma)$NormalizedMeasuredValue)
      fmb <- mean(filter(tempData, Genotype == mb)$NormalizedMeasuredValue)
      ###Calculate Additive and Proportional Epistasis
      AddEp <- fdm + fwt - (fma + fmb)
      PropEp <- log((fdm * fwt)/ (fma * fmb))
      ###Populate the data frame with this information:
      geneticInteractionsTmp[rowCount, 1] <- dm
      geneticInteractionsTmp[rowCount, 2] <- ma
      geneticInteractionsTmp[rowCount, 3] <- mb
      geneticInteractionsTmp[rowCount, 4] <- AddEp
      geneticInteractionsTmp[rowCount, 5] <- PropEp
      geneticInteractionsTmp[rowCount, 6] <- e
      geneticInteractionsTmp[rowCount, 7] <- m
      rowCount <- rowCount + 1
    }
    ###Add the rows of the temporary genetic interaction information to the main data frame
    geneticInteractions <- rbind(geneticInteractions, geneticInteractionsTmp)
  }
}

###Now, calculate Selection Coefficient:

###Initialize an empty data frame to populate with information:
selectionCoef<- data.frame(Mutant = rep(NA, 0), 
                                  SelectionCoefficient = rep(NA, 0),
                                  Experiment = rep(NA, 0),
                                  Measurement = rep(NA, 0))

for(e in c("DEPI1", "DEPI2", "DEPI3")){
  for (m in c("SN", "TSC", "SPF")){
    ###Create an empty data frame to fill with the information and calcuations:
    selectionCoefTmp <- data.frame(Mutant = rep(NA, 38), 
                                   SelectionCoefficient = rep(NA, 38),
                                   Experiment = rep(NA, 38), 
                                   Measurement = rep(NA, 38))
    count <- 1
    for (g in unique(NormalizedData$Genotype)){
      fm <- mean(filter(NormalizedData, Genotype == g, Experiment == e, Measurement == m)$NormalizedMeasuredValue)
      fwt <- mean(filter(NormalizedData, Genotype == "Col", Experiment == e, Measurement == m)$NormalizedMeasuredValue)
      TempSelectionCoef <- (fm - fwt) / fwt
      selectionCoefTmp[count, 1] <- g
      selectionCoefTmp[count, 2] <- TempSelectionCoef
      selectionCoefTmp[count, 3] <- e
      selectionCoefTmp[count, 4] <- m
      
      count <- count + 1
    }
    selectionCoef <- rbind(selectionCoef, selectionCoefTmp)
  }
}

selectionCoef <- selectionCoef%>%
  arrange(Measurement, Mutant)

library(stringi)

###This function adds a column that is used to sort the genotype in the plots:
add_number <- function(data_frame){
  ###First, if the genotype is Col0 (only genotype with length 4), assign 0 as number
  ###Else, assign number as genotype with "mpk" removed
  ###Example: mpk1 will be 1, mpk1-17 will be 1-17
  data_frame <- data_frame%>%
    mutate(number = ifelse(Mutant != "Col",(stri_sub(Mutant, 4, length(Mutant))), 0))
  ###Next, for all double mutants, replace "-" with "0"
  ###Example: 1-17 becomes 1017
  data_frame$number <- as.numeric(gsub("_","0", data_frame$number))
  ###Almost there! There's a problem with two single digit double mutants
  ###We need a four digit number to sort correctly 
  ###Example: mpk1-3 -> 1-3 -> 103, but we need it to be 1003 to sort correctly
  data_frame$number[data_frame$number == "103"] <- "1003"
  data_frame$number[data_frame$number == "506"] <- "5006"
  data_frame$number[data_frame$number == "608"] <- "6008"
  data_frame$number[data_frame$number == "609"] <- "6009"
  ###Convert number to a numberic in order to sort
  data_frame$number <- as.numeric(data_frame$number)
  data_frame <- data_frame%>%arrange(number)
  data_frame <- data_frame%>%mutate(number_2 = number)
  data_frame$number_2[nchar(data_frame$number_2) == 4] <- 0
  data_frame$number_2[nchar(data_frame$number_2) == 5] <- 0
  return(data_frame)
}

###This function adds a column that is used to sort the genotype in the plots:
add_number2 <- function(data_frame){
  ###First, if the genotype is Col0 (only genotype with length 4), assign 0 as number
  ###Else, assign number as genotype with "mpk" removed
  ###Example: mpk1 will be 1, mpk1-17 will be 1-17
  data_frame <- data_frame%>%
    mutate(number = ifelse(DoubleMutant != "Col",(stri_sub(DoubleMutant, 4, length(DoubleMutant))), 0))
  ###Next, for all double mutants, replace "-" with "0"
  ###Example: 1-17 becomes 1017
  data_frame$number <- as.numeric(gsub("_","0", data_frame$number))
  ###Almost there! There's a problem with two single digit double mutants
  ###We need a four digit number to sort correctly 
  ###Example: mpk1-3 -> 1-3 -> 103, but we need it to be 1003 to sort correctly
  data_frame$number[data_frame$number == "103"] <- "1003"
  data_frame$number[data_frame$number == "506"] <- "5006"
  data_frame$number[data_frame$number == "608"] <- "6008"
  data_frame$number[data_frame$number == "609"] <- "6009"
  ###Convert number to a numberic in order to sort
  data_frame$number <- as.numeric(data_frame$number)
  data_frame <- data_frame%>%arrange(number)
  data_frame <- data_frame%>%mutate(number_2 = number)
  data_frame$number_2[nchar(data_frame$number_2) == 4] <- 0
  data_frame$number_2[nchar(data_frame$number_2) == 5] <- 0
  return(data_frame)
}

selectionCoef <- add_number(selectionCoef)
selectionCoef$Mutant <- reorder(selectionCoef$Mutant, desc(selectionCoef$number))


SelectionFileNameList <- c()
SelectionPlots <- list()
count <- 1
for (m in c("SN", "TSC", "SPF")){ 
  tmpTitle <- paste("Selection Coefficient - ", m, sep = "")
  SelectionFileNameList[count] <- paste("Selection", m, ".tiff", sep = "")
  
  plot <- ggplot(data = filter(selectionCoef, Measurement == m), aes(x = Experiment, y = Mutant, fill = SelectionCoefficient))+
    geom_tile()+
    labs(x = "",
         y = "Genotype",
         title = tmpTitle)+
    theme_tufte(base_family = "Calibri")+
    scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(-0.55, .5), breaks = c(-0.55, 0, 0.5), labels = c("-0.55", "0", "0.5"))
  
  SelectionPlots[[count]] <- plot
  count <- count + 1
  
  print(plot)
}


for (i in 1:3) {
  file_name = SelectionFileNameList[i]
  tiff(file_name, units = "in", width = 7, height = 7, res = 500)
  print(SelectionPlots[[i]])
  dev.off()
}


geneticInteractions <- add_number2(geneticInteractions)
geneticInteractions$DoubleMutant <- reorder(geneticInteractions$DoubleMutant, desc(geneticInteractions$number))


AddEpFileNameList <- c()
additiveEpistasisPlots <- list()
count <- 1
for (m in c("SN", "TSC", "SPF")){
  tmpTitle <- paste("Additive Epistasis - ", m, sep = "")
  AddEpFileNameList[count] <- paste("AdditiveEpistasis", m, ".tiff", sep = "")
  
  plot <- ggplot(data = filter(geneticInteractions, Measurement == m), aes(x = Experiment, y = DoubleMutant, fill = AdditiveEpistasis))+
    geom_tile()+
    labs(x = "",
         y = "Genotype",
         title = tmpTitle)+
    theme_tufte(base_family = "Calibri")+
    scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(-85, 85), breaks = c(-85, 0, 85), labels = c("-85", "0", "85"))
  print(plot)
  additiveEpistasisPlots[[count]] <- plot
  count <- count + 1
}

for (i in 1:3) {
  file_name = AddEpFileNameList[i]
  tiff(file_name, units = "in", width = 7, height = 7, res = 500)
  print(additiveEpistasisPlots[[i]])
  dev.off()
}

PropEpFileNameList <- c()
proportionalEpistasisPlots <- list()
count <- 1
for (m in c("SN", "TSC", "SPF")){
  tmpTitle <- paste("Proportional Epistasis - ", m, sep = "")
  PropEpFileNameList[count] <- paste("ProportionalEpistasis", m, ".tiff", sep = "")
  
  plot <- ggplot(data = filter(geneticInteractions, Measurement == m), aes(x = Experiment, y = DoubleMutant, fill = ProportionalEpistatis))+
    geom_tile()+
    labs(x = "",
         y = "Genotype",
         title = tmpTitle)+
    theme_tufte(base_family = "Calibri")+
    scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(-0.51, 0.74), breaks = c(-0.51, 0, 0.74), labels = c("-0.51", "0", "0.74"))

  print(plot)
  proportionalEpistasisPlots[[count]] <- plot
  count <- count + 1
}

for (i in 1:3) {
  file_name = PropEpFileNameList[i]
  tiff(file_name, units = "in", width = 7, height = 7, res = 500)
  print(proportionalEpistasisPlots[[i]])
  dev.off()
}

###Write the geneticInteractions and selectionCoef data frames to a .csv file
###Remove the number and number2 columns - these were just used to order the genotypes correctly in the plots
geneticInteractions <- geneticInteractions[, -c(8:9)]
selectionCoef <- selectionCoef[, -c(5:6)]

###Write these files to .csv:
write.csv(geneticInteractions, file = "epistasis.csv", row.names = FALSE)
write.csv(selectionCoef, file = "selectionCoefficients.csv", row.names = FALSE)

###Notes from 11/23 Meeting with Melissa

###Direction of the effect is different 
###In some cases, effect is just less strong
###In others, there is an opposite effect
###Third experiment = less stress = smaller interaction?

###To Do
###Meet with Melissa once these heatmaps are done [triplet and all] [pdf is okay]

###Create a similar set of plots for phi2 for a couple of time points
###Based on the heat maps if there's a time point that looks interesting





