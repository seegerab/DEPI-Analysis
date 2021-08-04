####################################################################################
# The goal of this script is to use the .csv files created in the CreaveCSVFiles_07212021.R script
# to calculuate epistasis values and selection coefficients for the data. I will
# use the .csv files with the outliers removed, but this can easily be changed by reading 
# in the .csv file without the outliers removed. 
#
# Further, this script creates visualizations for the selection coefficients and epistasis
# calculations, and computes p-values for the fitness values (SN, TSC, SPF).
#
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

### Read in the .csv file created in the CreateCSVFiles.R script with quantile normalized values and outliers included 
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/Data/Fitness_NormalizedValues_OutliersIncluded_07212021.csv", header = TRUE)

####################################################################################
### add_number and add_number2 are soley for plotting purposes!
####################################################################################

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

####################################################################################
### Visualize Normalization
####################################################################################

### Convert the flat column to a factor
data$Flat <- as.factor(data$Flat)
###Visualize the normalization to make sure quantile normalized flats have the same distribution for all experiments
for (i in c("DEPI1", "DEPI2", "DEPI3")){
  for (m in c("TSC", "SN", "SPF")){
    plotData <- data%>%
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
  } 
}

### Generally I'm happy with these distributions, although SN seems less in alignment than the other measurements 

####################################################################################
### Epistasis Calculations
####################################################################################

### Here are all the double mutants:
all_double_mutants = list()
for (gen in levels(as.factor(data$Genotype))) {
  if (str_detect(gen, "_") == T) {
    all_double_mutants = c(all_double_mutants, gen)
  }
}

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
    tempData <- filter(data, Experiment == e, Measurement == m)
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

### Add a column with the normalization scheme
geneticInteractions <- geneticInteractions%>%
  mutate(Normalization = "ByFlat")

####################################################################################
### Genetic Interationcs Calculations - Repeat with Normalization by Exp
####################################################################################

### Read in the data that is normalized by experiment
normByExp <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/Fitness_NormalizedValuesByExp_OutliersIncluded_07232021.csv", header = TRUE, row.names = 1)

###Initialize an empty data frame to populate with information:
geneticInteractionsExp <- data.frame(DoubleMutant = rep(NA, 0), 
                                  MutantA = rep(NA, 0),
                                  MutantB = rep(NA, 0), 
                                  AdditiveEpistasis = rep(NA, 0),
                                  ProportionalEpistatis = rep(NA, 0),
                                  Experiment = rep(NA,0),
                                  Measurement = rep(NA, 0))

###Loop through each measurement:
for (m in c("SN", "TSC", "SPF")){
    ###Filter to each specific experiment and measurement
    tempData <- filter(normByExp, Measurement == m)
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
      fdm <- mean(filter(tempData, Genotype == dm)$NormalizedExp)
      fwt <- mean(filter(tempData, Genotype == "Col")$NormalizedExp)
      fma <- mean(filter(tempData, Genotype == ma)$NormalizedExp)
      fmb <- mean(filter(tempData, Genotype == mb)$NormalizedExp)
      ###Calculate Additive and Proportional Epistasis
      AddEp <- fdm + fwt - (fma + fmb)
      PropEp <- log((fdm * fwt)/ (fma * fmb))
      ###Populate the data frame with this information:
      geneticInteractionsTmp[rowCount, 1] <- dm
      geneticInteractionsTmp[rowCount, 2] <- ma
      geneticInteractionsTmp[rowCount, 3] <- mb
      geneticInteractionsTmp[rowCount, 4] <- AddEp
      geneticInteractionsTmp[rowCount, 5] <- PropEp
      geneticInteractionsTmp[rowCount, 6] <- "Doesn't Apply"
      geneticInteractionsTmp[rowCount, 7] <- m
      rowCount <- rowCount + 1
    }
    ###Add the rows of the temporary genetic interaction information to the main data frame
    geneticInteractionsExp <- rbind(geneticInteractionsExp, geneticInteractionsTmp)
  
}
### Add a column with the normalization scheme
geneticInteractionsExp <- geneticInteractionsExp%>%
  mutate(Normalization = "ByExp")

### Combine the data frames
geneticInteractionsAll <- rbind(geneticInteractions, geneticInteractionsExp)

####################################################################################
### Selection Coefficient Calculations
####################################################################################

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
    ### Initialize a row count to be 1
    count <- 1
    ### For each genotype
    for (g in unique(data$Genotype)){
      ### Compute the mean fitness of the mutant
      fm <- mean(filter(data, Genotype == g, Experiment == e, Measurement == m)$NormalizedMeasuredValue)
      ### Compute the mean fitness of wild type
      fwt <- mean(filter(data, Genotype == "Col", Experiment == e, Measurement == m)$NormalizedMeasuredValue)
      ### Calculate the selection coefficient using these values
      TempSelectionCoef <- (fm - fwt) / fwt
      ### Populate the row of the data frame with genotype, selection coefficient, experiment, and measurement
      selectionCoefTmp[count, 1] <- g
      selectionCoefTmp[count, 2] <- TempSelectionCoef
      selectionCoefTmp[count, 3] <- e
      selectionCoefTmp[count, 4] <- m
      ### Increase row count variable by 1
      count <- count + 1
    }
    ### Once each genotype has been iterated through, bind the rows of the temp data frames for these measurements to the main data frame
    selectionCoef <- rbind(selectionCoef, selectionCoefTmp)
  }
}
selectionCoef <- selectionCoef%>%
  arrange(Measurement, Mutant)%>%
  mutate(Normalization = "ByFlat")


####################################################################################
### Selection Coefficient Calculations - Repeat with Normalization by Exp
####################################################################################

normByExp <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/Fitness_NormalizedValuesByExp_OutliersIncluded_07232021.csv", header = TRUE, row.names = 1)

###Initialize an empty data frame to populate with information:
selectionCoefExp <- data.frame(Mutant = rep(NA, 0), 
                           SelectionCoefficient = rep(NA, 0),
                           Measurement = rep(NA, 0))

for (m in c("SN", "TSC", "SPF")){
    ###Create an empty data frame to fill with the information and calcuations:
    selectionCoefTmp <- data.frame(Mutant = rep(NA, 38), 
                                   SelectionCoefficient = rep(NA, 38),
                                   Experiment = rep(NA, 38), 
                                   Measurement = rep(NA, 38))
    ### Initialize a row count to be 1
    count <- 1
    ### For each genotype
    for (g in unique(normByExp$Genotype)){
      ### Compute the mean fitness of the mutant
      fm <- mean(filter(normByExp, Genotype == g, Measurement == m)$NormalizedExp)
      ### Compute the mean fitness of wild type
      fwt <- mean(filter(normByExp, Genotype == "Col", Measurement == m)$NormalizedExp)
      ### Calculate the selection coefficient using these values
      TempSelectionCoef <- (fm - fwt) / fwt
      ### Populate the row of the data frame with genotype, selection coefficient, experiment, and measurement
      selectionCoefTmp[count, 1] <- g
      selectionCoefTmp[count, 2] <- TempSelectionCoef
      selectionCoefTmp[count, 3] <- e
      selectionCoefTmp[count, 4] <- m
      ### Increase row count variable by 1
      count <- count + 1
    }
    ### Once each genotype has been iterated through, bind the rows of the temp data frames for these measurements to the main data frame
    selectionCoefExp <- rbind(selectionCoefExp, selectionCoefTmp)
  
}
selectionCoefExp <- selectionCoefExp%>%
  arrange(Measurement, Mutant)%>%
  mutate(Normalization = "ByExp")

### Finally, combine the calculations based on normalization scheme

selectionCoefAll <- rbind(selectionCoef, selectionCoefExp)

####################################################################################
### Plot Selection Coefficients
####################################################################################

### Filter the selection coefficient to only have the values based on the normalization by flat
selectionCoef <- filter(selectionCoef, Normalization == "ByFlat")
### Include columns with number and number2 to sort genotypes in the plots
selectionCoef <- add_number(selectionCoef)
selectionCoef$Mutant <- reorder(selectionCoef$Mutant, desc(selectionCoef$number))
### Initialize an empty file name list
SelectionFileNameList <- c()
### Initialize an empty list to populate with the plots
SelectionPlots <- list()
### Initialize a count set to 1
count <- 1
### For each measurement 
for (m in c("SN", "TSC", "SPF")){ 
  ### Create a plot title
  tmpTitle <- paste("Selection Coefficient - ", m, sep = "")
  ### Append a file name in the form Selectionm.tiff, where m is either SN, TSC, or SPF to the file name list
  SelectionFileNameList[count] <- paste("Selection", m, ".tiff", sep = "")
  ### Create the plot
  plot <- ggplot(data = filter(selectionCoef, Measurement == m), aes(x = Experiment, y = Mutant, fill = SelectionCoefficient))+
    geom_tile()+
    labs(x = "",
         y = "Genotype",
         title = tmpTitle)+
    theme_tufte(base_family = "Calibri")+
    scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(-0.55, .5), breaks = c(-0.55, 0, 0.5), labels = c("-0.55", "0", "0.5"))
  ### Append the plot to the plot list
  SelectionPlots[[count]] <- plot
  ### Increase the count by 1
  count <- count + 1
  ### Print the plot
  print(plot)
}
### For each of the three plots just created
for (i in 1:3) {
  ### Set the file name to be the ith element in the SelectionFileNameList
  file_name = SelectionFileNameList[i]
  ### Use tiff() and the file_name to set the parameters for saving the plot
  tiff(file_name, units = "in", width = 7, height = 7, res = 500)
  ### Print the ith plot in SelectionPlots list 
  print(SelectionPlots[[i]])
  dev.off()
}

####################################################################################
### Plot Genetic Interactions
####################################################################################


### Filter genetic interactions to only have the values based on the normalization by flat
geneticInteractions <- filter(geneticInteractionsAll, Normalization == "ByFlat")
### Include a column with number and number2 for plotting purposes
geneticInteractions <- add_number2(geneticInteractions)
geneticInteractions$DoubleMutant <- reorder(geneticInteractions$DoubleMutant, desc(geneticInteractions$number))
### Initialize an empty vector to populate with file names
AddEpFileNameList <- c()
### Initialize an empty list to populate with the plots
additiveEpistasisPlots <- list()
### Set a count to be 1
count <- 1
### For each measurement 
for (m in c("SN", "TSC", "SPF")){
  ### Create a title for the plot
  tmpTitle <- paste("Additive Epistasis - ", m, sep = "")
  ### Append this title to the ith element of the file name list
  AddEpFileNameList[count] <- paste("AdditiveEpistasis", m, ".tiff", sep = "")
  ### Create the plot
  plot <- ggplot(data = filter(geneticInteractions, Measurement == m), aes(x = Experiment, y = DoubleMutant, fill = AdditiveEpistasis))+
    geom_tile()+
    labs(x = "",
         y = "Genotype",
         title = tmpTitle)+
    theme_tufte(base_family = "Calibri")+
    scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(-85, 85), breaks = c(-85, 0, 85), labels = c("-85", "0", "85"))
  ### Print the plot
  print(plot)
  ### Append the plot to the plot list
  additiveEpistasisPlots[[count]] <- plot
  ### Increase the count by 1 
  count <- count + 1
}
### For each plot
for (i in 1:3) {
  ### Set file_name to be the ith element of the file name list
  file_name = AddEpFileNameList[i]
  ### Set the parameters of tiff to save the plot as an image
  tiff(file_name, units = "in", width = 7, height = 7, res = 500)
  ### Print the ith plot
  print(additiveEpistasisPlots[[i]])
  dev.off()
}
### Initialize an empty vector for the file names
PropEpFileNameList <- c()
### Initialize an empty list for the plots
proportionalEpistasisPlots <- list()
### Set a count to be 1
count <- 1
### For each measurement 
for (m in c("SN", "TSC", "SPF")){
  ### Create a title for the plot
  tmpTitle <- paste("Proportional Epistasis - ", m, sep = "")
  ### Create a file name for the plot and appent to the file name vector
  PropEpFileNameList[count] <- paste("ProportionalEpistasis", m, ".tiff", sep = "")
  #### Create the plot
  plot <- ggplot(data = filter(geneticInteractions, Measurement == m), aes(x = Experiment, y = DoubleMutant, fill = ProportionalEpistatis))+
    geom_tile()+
    labs(x = "",
         y = "Genotype",
         title = tmpTitle)+
    theme_tufte(base_family = "Calibri")+
    scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(-0.51, 0.74), breaks = c(-0.51, 0, 0.74), labels = c("-0.51", "0", "0.74"))
  ### Print the plot
  print(plot)
  ### Append the plot to the plot list
  proportionalEpistasisPlots[[count]] <- plot
  ### Increase the count by 1
  count <- count + 1
}
### For each plot
for (i in 1:3) {
  ### Extract the ith element from the file name list
  file_name = PropEpFileNameList[i]
  ### Set the parameters for saving the plot as a tiff
  tiff(file_name, units = "in", width = 7, height = 7, res = 500)
  ### Print the ith element of the plot list
  print(proportionalEpistasisPlots[[i]])
  dev.off()
}

####################################################################################
### Save geneticInteractions and selectionCoef as .csv files
####################################################################################

###Remove the number and number2 columns (these were used to sort the plots correctly)
geneticInteractionsFinal <- geneticInteractionsAll%>%
  rename(AdditiveEp_OutliersIncluded = AdditiveEpistasis,
         PropEp_OutliersIncluded = ProportionalEpistatis)%>%
  na.omit()



selectionCoefFinal <- selectionCoefAll%>%
  rename(SC_OutliersIncluded = SelectionCoefficient)%>%
  na.omit()

###Write to .csv files:
write.csv(geneticInteractionsFinal, file = "epistasis_07292021_OutliersIncluded.csv", row.names = FALSE)
write.csv(selectionCoefFinal, file = "selectionCoefficients_07292021_OutliersIncluded.csv", row.names = FALSE)

####################################################################################
### Calculate p-values for fitness measurements (SN, TSC, SPF)
####################################################################################

### Note: p_value and correct_p_value are the same functions used in the DEPI analysis
p_value <- function (data_frame){
  ###Initialize an empty data frame
  out = data.frame()
  ###For each genotype:
  for (i in unique(filter(data_frame, (Genotype != "Col" & Genotype != "mpk5_17"))$Genotype)){
    ###We don't want to make comparisons of WT to itself - this could impact FDR correction
    indiv_data <- data_frame%>% 
      ### Group by measurement and experiment
      group_by(Measurement, Experiment)%>%
      ### Add a column with the p-value
      mutate(p = (wilcox.test(NormalizedMeasuredValue[Genotype == i], NormalizedMeasuredValue[Genotype == "Col"], correct = FALSE, paired = FALSE))$p.value)%>%
      ###Add a column with each genotype
      mutate(Genotype = i)%>%
      select(Experiment, Genotype, Measurement, p)
    ###Add individual information to the main data frame
    out <- rbind(as.data.frame(indiv_data), out)%>%distinct()
  }
  return(out)
}
corrected_p_value <- function(data_frame){
  out <- data_frame%>%
    ###Group by experiment and measurement - we are correcting by the number of genotypes 
    group_by(Experiment, Measurement)%>%
    mutate(p_adj = p.adjust(p, method = "fdr"))
  return(out) 
}
###Convert the normalized data from a tibble to a data frame
data <- as.data.frame(data)
### Apply the p_values function to the data
p_values <- p_value(data)  
### Adjust the p-values using an FDR correction
adjusted_p_values <- p_values%>%
  group_by(Measurement, Experiment)%>%
  mutate(p_adj = round(p.adjust(p, method = 'fdr'), 7))%>%
  rename(p_OutliersIncluded = p,
         p_adj_OutliersIncluded = p_adj)
### Add a column with a boolean indicating if the adjusted p-value is significant or not
adjusted_p_values$significant <- adjusted_p_values$p_adj_OutliersIncluded < 0.05
### Convert this to a .csv file to share with Melissa:
write.csv(adjusted_p_values, file = "FitnessPValues_07222021_OutliersIncluded.csv")



