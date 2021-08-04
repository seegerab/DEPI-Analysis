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
  gather(Measurement, Measured_Value, c("SN", "TSC", "SPF"))
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
### SAVE THESE NORMALIZED VALUES TO USE IN LATER SCRIPTS. REMEMBER, THIS IS WITH REMOVING THE OUTLIERS
#write.csv(NormalizedData, file = "Fitness_NormalizedValues_OutliersRemoved.csv")


####################################################################################
####################################################################################
### Bootstrapping - Col0
####################################################################################
####################################################################################


### Set the graphing parameters to have 3 rows and 3 columns of plots
par(mfrow = c(3, 3))
### Initialize a matrix to populate with the 2.5% quantile and 97.5% quantile for each experiment and measurement
quantile_matrix <- matrix(nrow = 9,
                          ncol = 4,
                          dimnames = list(NULL,
                                           c("Experiment", "Measurement", "Quantile_0.025", "Quantile_0.975"))
                          )
### Initialize a row count variable set to 1
row_count <- 1
### For each experiment
for (j in unique(NormalizedData$Experiment)){
  ### For each measurement
  for (k in unique(NormalizedData$Measurement)){
    ### Filter to extract the normalized measured values of Col0 for that measurement and experiment
    temp_bootstrapping_data <- filter(NormalizedData, Experiment == j, Measurement == k, Genotype == "Col")$NormalizedMeasuredValue
    ### Initialize a vector of 10000 NAs 
    temp_mean_deviation_vector <- rep(NA, 10000)
    ### For 10000 bootstrap iterations
    for (i in 1:10000){
      ### replace the ith element in  temp_mean_deviation vector with the mean of the deviation from the sample to the mean of all Col0
      
      ### (There are a lot of operations strung together; I did this to avoid having to re-assign variables 1000 times to speed up the process)
      
      temp_mean_deviation_vector[i] <-mean((
        ### Sample from the measured values associated with Col0 for the jth experiment and the kth measurement
        sample(temp_bootstrapping_data, 
             ### Sample 1/2 of the number of Col0 plants (Should be about 7, though smaller in some cases)
             size = floor(length(temp_bootstrapping_data) / 2),
             ### Sample with replacement 
             replace = TRUE) 
        ### Take the difference between the sample means and the mean of all Col0 plants for the jth experiment and kth measurement
        - mean(temp_bootstrapping_data)) / mean(temp_bootstrapping_data))
    }
    
    ### Create a histogram of these 1,000 deviations from the mean
    hist(temp_mean_deviation_vector, 
         xlab = " ",
         main = paste(j, k, sep = " "))
    
    ### Add the experiment to the first column
    quantile_matrix[row_count, 1] <- j
    ### Add the measurement to the second column
    quantile_matrix[row_count, 2] <- k
    ### Add the 2.5% percentile of the 10000 bootstrap values to the third column
    quantile_matrix[row_count, 3] <- quantile(temp_mean_deviation_vector, probs = c(0.025))
    ### Add the 97.5% percentile of the 10000 boostrap values to the fourth column
    quantile_matrix[row_count, 4] <- quantile(temp_mean_deviation_vector, probs = c(0.975))
    
    ### Re-assign the row_count variable so the next iteration populates the next row in the matrix
    row_count <- row_count + 1
  }
}



####################################################################################
####################################################################################
### Similar approach - each genotype
####################################################################################
####################################################################################

### Create a vector of the unique genotypes in the experiment

### mpk5-17 is missing from DEPI1. So, remove this genotype for now so we're able to compare genotypes between experiments
genotype_vector <- unique(NormalizedData$Genotype)[unique(NormalizedData$Genotype) != "mpk5_17"]

### Create a column name vector to use in the matrix. This is the unique genotypes, with columns included for the Measurement and Experiment
colname_vector <- append(as.vector(genotype_vector), c("Measurement", "Experiment"))

### Initialize a matrix of NA values to populate with mean(measured_values - mean(Col0 measured value)) for each experiment and measured value
selectionCoef_matrix <- matrix(data = NA, 
                      nrow = 9,
                      ncol = 39,
                      dimnames = list(NULL,
                                      append(as.vector(genotype_vector), c("Experiment", "Measurement"))))
                      

### Initialize matrix_count to be 1; this metric is used as an indicator of the row number to populate in the matrix
matrix_count <- 1

### For each experiment
for (j in unique(NormalizedData$Experiment)){
  ### For each measurement (SN, TSC, SPF)
  for (k in unique(NormalizedData$Measurement)){
    
    ### Filter the data to include the jth experiment and the kth measurement
    temp_data <- filter(NormalizedData, Experiment == j, Measurement == k)
    
    ### Initial an empty vector to populate with the mean(measured_values - mean(Col0 measured value)) for each genotype
    temp_vector <- rep(NA, length(unique(temp_data$Genotype)))
    
    ### For each genotype
    for (i in 1:length(genotype_vector)){
      
      ### Put the calculated value for the ith genotype in the ith element in temp_vector
      temp_vector[i] <- 
        
        ### Round to include 6 decimal places
        round(
        ### Calculate the mean of the difference between each normalized measured value and the mean of Col0
        mean(
        ### Filter the data to the ith genotype, and extract the normalized measured values
        (filter(temp_data, Genotype == genotype_vector[i])$NormalizedMeasuredValue - mean(filter(temp_data, Genotype == "Col")$NormalizedMeasuredValue))
        / mean(filter(temp_data, Genotype == "Col")$NormalizedMeasuredValue)), 6)
    }
    ### Filter the data to Col0, extract the normalized values, and take the mean of these values  
    
    ### Add experiment to the 38th element of temp_vector
    temp_vector[i + 1] <- j
    ### Add measurement to the 39th element of temp_vector
    temp_vector[i + 2] <- k
    ### temp_vector now becomes a row of the matrix
    selectionCoef_matrix[matrix_count, ] <- temp_vector
    
    ### Re-assign the matrix_count variable so the next iteration populates the next row in the matrix
    matrix_count <- matrix_count + 1
      
  }
}

####################################################################################
####################################################################################
### Combine matrices with the quantiles from Col0 and the selection coefficients for each genotype
####################################################################################
####################################################################################


### First, convert the matrices to dataframes in order to use the cbind() command

### Note that the experiment and measurement columns will be duplicated, but this is a good check to make sure the data frames are correctly joined

### The final type is a data frame
exceeds_quantile_df <- cbind(as.data.frame(selectionCoef_matrix), as.data.frame(quantile_matrix))

### The selection coefficents are characters. Convert these columns to numeric. 
exceeds_quantile_df[, 1:37] <- lapply(exceeds_quantile_df[, 1:37], as.numeric)

### The quantiles are characters. Convert these columns to numeric.
exceeds_quantile_df[, 42:43] <- lapply(exceeds_quantile_df[, 42:43], as.numeric)

### Create a duplicate of this data frame to use to plot the data later:
heatmap_data <- exceeds_quantile_df

### For each row in the data frame
for (i in 1:nrow(exceeds_quantile_df)){
  ### For each of the columns that are the selection coefficent calculations for each genotype
  for (j in 1:(ncol(exceeds_quantile_df) - 7)){
    ### If the value in the cell is less than the 2.5% quantile for Col0
    if (exceeds_quantile_df[i,j] < exceeds_quantile_df[i,42]){
      ### Replace the value with the statement -1000
      ### Note: I initially included a character string to replace these values, but am using a numeric to make it easier to plot the data later
      exceeds_quantile_df[i,j] <- -1000
    } 
    ### Else if the value the cell is greater than the 97.5% for Col0
    
    else if (exceeds_quantile_df[i,j] > exceeds_quantile_df[i,43]){
      ### Replace the value with the statement 1000
      ### Note: I initially included a character string to replace these values, but am using a numeric to make it easier to plot the data later
      
      exceeds_quantile_df[i,j] <- 1000
    }
  }
}

### This data frame is a way to see "at a glance" where the genotypes are behaving differently from Col0
### But, a better way to visualize this is with heat maps:

####################################################################################
####################################################################################
### Create heat maps
####################################################################################
####################################################################################

### Convert the data to the correct format to use to create heatmaps in ggplot:
### Remove the duplicate experiment and measurement columns, as well as the columns with the quantiles
heatmap_data <- heatmap_data[, 1:39]

heatmap_data <- heatmap_data%>%
  ### Pivot the data frame so the columns for each genotype become a "key" column
  ### i.e. convert the dataframe from wide to long format
  pivot_longer(cols = c(starts_with("mpk"), "Col"),
                            names_to = "genotype",
                            values_to = "selectionCoef"
                            )%>%
  ### Remove the rows with Col, because the selection coefficient is always 0
  ### There is no biological meaning because Col is compared to itself
  filter(genotype != "Col")
### Upon initially creating the heat map, the genotypes are not ordered in a way that makes sense
### Include a function that will sort the genotypes
add_number <- function(data_frame){
  ###First, if the genotype is Col0 (only genotype with length 4), assign 0 as number
  ###Else, assign number as genotype with "mpk" removed
  ###Example: mpk1 will be 1, mpk1-17 will be 1-17
  data_frame <- data_frame%>%
    mutate(number = ifelse(genotype != "Col0",(stri_sub(genotype, 4, length(genotype))), 0))
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

  return(data_frame)
}

### Apply the add_number function above to the data frame
heatmap_data <- add_number(heatmap_data)
### Reorder genotype by the "number" column. Multiply this by -1, because the order was opposite of what I wanted it to be
heatmap_data$genotype <- reorder(heatmap_data$genotype, heatmap_data$number * -1)

### Filter the data to create 3 dataframes to plot, one for each measurement
heatmap_data_tsc <- filter(heatmap_data, Measurement == "TSC")
heatmap_data_sn <- filter(heatmap_data, Measurement == "SN")
heatmap_data_spf <- filter(heatmap_data, Measurement == "SPF")

### Find the range of selection coefficients so I can set the scale of each plot to be the same, so we can easily compare between them
lowerbound <- floor(min(heatmap_data$selectionCoef) * 100)/100
upperbound <- ceiling(max(heatmap_data$selectionCoef) * 100)/100

### Create TSC heatmap:
ggplot(data = heatmap_data_tsc, aes(x = Experiment, y = genotype, fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for TSC")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.52", "0", "0.46"))

### Create SPF heatmap:
ggplot(data = heatmap_data_spf, aes(x = Experiment, y = genotype, fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SPF")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.52", "0", "0.46"))

### Create SN heatmap:
ggplot(data = heatmap_data_sn, aes(x = Experiment, y = genotype, fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SN")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.52", "0", "0.46"))

####################################################################################
####################################################################################
### Create heat maps - binned data
####################################################################################
####################################################################################

### Here, recreate the heat maps with binned data

### This allows us to highlight the cells that either are less than the 0.025 quantile or exceed the 0.975 quantile

heatmap_data_binned <- exceeds_quantile_df[, 1:39]%>%
  ### Pivot the data frame so the columns for each genotype become a "key" column
  ### i.e. convert the dataframe from wide to long format
  pivot_longer(cols = c(starts_with("mpk"), "Col"),
               names_to = "genotype",
               values_to = "selectionCoef"
  )%>%
  ### Remove the rows with Col, because the selection coefficient is always 0
  ### There is no biological meaning because Col is compared to itself
  filter(genotype != "Col")
### Create the bins 

# set up cut-off values 

### NOTE: -1000 represents the selection coeffients that are less than 0.025 quantile
### 1000 represents the selection coefficients that exceed the 0.975 quantile

breaks <- c(-1200,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2, 0.3, 1200)
# specify interval/bin labels
tags <- c("Lower than 0.025 Quantile","[-0.04, -0.3)", "[-0.3, -0.2)", "[-0.2, -0.1)", "[-0.1, 0)", "[0, 0.1)","[0.1, 0.2)", "[0.2, 0.3)","Exceeds 0.975 Quantile")
# bucketing values into bins
group_tags <- cut(heatmap_data_binned$selectionCoef, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_tags)

### Add a columns "bins" to the data frame with the group_tags
heatmap_data_binned$bins <- group_tags
### Make this column a factor
heatmap_data_binned$bins <- as.factor(heatmap_data_binned$bins)

### Apply the add_number function above to the data frame
heatmap_data_binned <- add_number(heatmap_data_binned)
### Reorder genotype by the "number" column. Multiply this by -1, because the order was opposite of what I wanted it to be
heatmap_data_binned$genotype <- reorder(heatmap_data_binned$genotype, heatmap_data_binned$number * -1)

### Filter to create a data frame for each measurement
heatmap_binned_sn <- filter(heatmap_data_binned, Measurement == "SN")
heatmap_binned_tsc <- filter(heatmap_data_binned, Measurement == "TSC")
heatmap_binned_spf <- filter(heatmap_data_binned, Measurement == "SPF")

### Create binned SN heatmap
ggplot(data = heatmap_binned_sn, aes(x = Experiment, y = genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SN")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),

        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("blue", cm.colors(length(unique(heatmap_binned_sn$bins)) - 2) , "red"))

### Create binned TSC heat map
ggplot(data = heatmap_binned_tsc, aes(x = Experiment, y = genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for TSC")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),

        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("blue", cm.colors(length(unique(heatmap_binned_tsc$bins)) - 1)))

### Create binned SPF heat map
ggplot(data = heatmap_binned_spf, aes(x = Experiment, y = genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SPF")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),

        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("blue", cm.colors(length(unique(heatmap_binned_spf$bins)) - 2), "red"))

####################################################################################
####################################################################################
### Visualize normalizations and begin epistasis calculations
####################################################################################
####################################################################################


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
  # ggsave(file=QuantFileName, arrangeGrob(plot1, plot2), width = 7, height = 7)
  } 
}

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

### Calculate p-values comparing each measurement to wild type

p_value <- function (data_frame){
  ###Initialize an empty data frame
  out = data.frame()
  ###For each genotype:
  for (i in unique(filter(data_frame, (Genotype != "Col" & Genotype != "mpk5_17"))$Genotype)){
    ###We don't want to make comparisons of WT to itself - this could impact FDR correction
    indiv_data <- data_frame%>% 
      group_by(Measurement, Experiment)%>%
      mutate(p = (wilcox.test(NormalizedMeasuredValue[Genotype == i], NormalizedMeasuredValue[Genotype == "Col"], correct = FALSE, paired = FALSE))$p.value)%>%
      
      ###Add a column with each genotype
      mutate(Genotype = i)%>%
      
      select(Experiment, Genotype, Measurement, p)
      ###Add individual information to the main data frame
      out <- rbind(as.data.frame(indiv_data), out)%>%distinct()
    }

    
  return(out)}

corrected_p_value <- function(data_frame){
  out <- data_frame%>%
    ###Group by experiment and measurement - we are correcting by the number of genotypes 
    group_by(Experiment, Measurement)%>%
    mutate(p_adj = p.adjust(p, method = "fdr"))
  
    return(out) 
}

###Convert the normalized data from a tibble to a data frame
NormalizedData <- as.data.frame(NormalizedData)
### Apply the p_values function to the NormalizedData
p_values <- p_value(NormalizedData)  
### Adjust the p-values using an FDR correction
adjusted_p_values <- p_values%>%
  group_by(Measurement, Experiment)%>%
  mutate(p_adj = round(p.adjust(p, method = 'fdr'), 7))

### Add a column with a boolean indicating if the adjusted p-value is significant or not
adjusted_p_values$significant <- adjusted_p_values$p_adj < 0.05

### Convert this to a .csv file to share with Melissa:
write.csv(adjusted_p_values, file = "Fitness_PValues_OutliersRemoved.csv")




###Confirm that this p-value correction is actually doing what I think it is:
test <- filter(p_values, Experiment == "DEPI1",  Measurement == "TSC")%>%
  select(p, Genotype)
test$adjusted <- p.adjust(test$p, method = "fdr")
filter(adjusted_p_values, Experiment == "DEPI1", Measurement == "TSC")$p_adj



library(stats)
P_value <- c(0.0001, 0.001, 0.006, 0.03, 0.095, 0.117, 0.234, 0.552, 0.751, 0.985)

p.adjust(P_value, method="bonferroni") ## [1] 0.001 0.010 0.060 0.300 0.950 1.000 1.000 1.000 1.000 1.000

P_value <- test$p

p.adjust(P_value, method="fdr")


adjusted_p_values%>%
  group_by(Experiment, Measurement)%>%
  summarize(length(unique(p_adj)))


