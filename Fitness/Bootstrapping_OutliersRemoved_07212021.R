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

### Read in the data - this is the data frame with the outliers removed
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/Data/Fitness_NormalizedValues_OutliersRemoved_07212021.csv", header = TRUE)

####################################################################################
### Bootstrapping - Col0
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
for (j in unique(data$Experiment)){
  ### For each measurement
  for (k in unique(data$Measurement)){
    ### Filter to extract the normalized measured values of Col0 for that measurement and experiment
    temp_bootstrapping_data <- filter(data, Experiment == j, Measurement == k, Genotype == "Col")$NormalizedMeasuredValue
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
    # Uncomment the lines of code below to plot (I ran into a problem with populating the matrix when I got an "increase plot window" error)
    # hist(temp_mean_deviation_vector, 
    #      xlab = " ",
    #      main = paste(j, k, sep = " "))
    
    ### Add the experiment to the first column
    quantile_matrix[row_count, 1] <- j
    ### Add the measurement to the second column
    quantile_matrix[row_count, 2] <- k
    ### Add the 2.5% percentile of the 10000 bootstrap values to the third column
    quantile_matrix[row_count, 3] <- quantile(temp_mean_deviation_vector, probs = c(0.025))
    ### Add the 97.5% percentile of the 10000 bootstrap values to the fourth column
    quantile_matrix[row_count, 4] <- quantile(temp_mean_deviation_vector, probs = c(0.975))
    ### Re-assign the row_count variable so the next iteration populates the next row in the matrix
    row_count <- row_count + 1
  }
}
####################################################################################
### Similar approach - each genotype
####################################################################################

### Create a vector of the unique genotypes in the experiment
### mpk5-17 is missing from DEPI1. So, remove this genotype for now so we're able to compare genotypes between experiments
genotype_vector <- unique(data$Genotype)[unique(data$Genotype) != "mpk5_17"]
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
for (j in unique(data$Experiment)){
  ### For each measurement (SN, TSC, SPF)
  for (k in unique(data$Measurement)){
    ### Filter the data to include the jth experiment and the kth measurement
    temp_data <- filter(data, Experiment == j, Measurement == k)
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
### Combine matrices with the quantiles from Col0 and the selection coefficients for each genotype
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

### I ran into problems plotting the heat maps; set dev.off() in case it was set "on" in order to plot
dev.off()
### Create TSC heatmap:
ggplot(data = heatmap_data_tsc, aes(x = Experiment, y = genotype, fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for TSC")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.52", "0", "0.45"))

### Create SPF heatmap:
ggplot(data = heatmap_data_spf, aes(x = Experiment, y = genotype, fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SPF")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.52", "0", "0.45"))

### Create SN heatmap:
ggplot(data = heatmap_data_sn, aes(x = Experiment, y = genotype, fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SN")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.52", "0", "0.45"))


####################################################################################
### Create heat maps - binned data
####################################################################################

### Here, recreate the heat maps with binned data. This allows us to highlight the cells that either are less than the 0.025 quantile or exceed the 0.975 quantile.
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
### Create the bins:
### set up cut-off values 
### NOTE: -1000 represents the selection coeffients that are less than 0.025 quantile and 1000 represents the selection coefficients that exceed the 0.975 quantile
breaks <- c(-1200,-0.4,-0.1,0.1, 0.3, 1200)
# specify interval/bin labels
tags <- c("Lower than 0.025 Quantile","[-0.04, -0.1)", "[-0.1, 0.1)", "[0.1, 0.04)","Exceeds 0.975 Quantile")
# bucketing values into bins
group_tags <- cut(heatmap_data_binned$selectionCoef, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect bins
summary(group_tags)
### Add a columns "bins" to the data frame with the group_tags and convert to a factor
heatmap_data_binned$bins <- as.factor(group_tags)
### Apply the add_number function above to the data frame
heatmap_data_binned <- add_number(heatmap_data_binned)
### Reorder genotype by the "number" column. Multiply this by -1, because the order was opposite of what I wanted it to be
heatmap_data_binned$genotype <- reorder(heatmap_data_binned$genotype, heatmap_data_binned$number * -1)
### Filter to create a data frame for each measurement
heatmap_binned_sn <- filter(heatmap_data_binned, Measurement == "SN")
heatmap_binned_tsc <- filter(heatmap_data_binned, Measurement == "TSC")
heatmap_binned_spf <- filter(heatmap_data_binned, Measurement == "SPF")
### Create a color palette with white in the middle and blue and red at the extremes
### This palettes is specifically for the 5 bins
col.palette <- rev(c("#FF0000", 
                 "#FFCCCB",
                 "#FFFEFE",
                 "#D4EBF2",
                 "#0000FF"))
### Create binned SN heatmap
ggplot(data = heatmap_binned_sn, aes(x = Experiment, y = genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SN")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = col.palette)
### Create binned TSC heat map
ggplot(data = heatmap_binned_tsc, aes(x = Experiment, y = genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for TSC")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = col.palette)
### Create binned SPF heat map
ggplot(data = heatmap_binned_spf, aes(x = Experiment, y = genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Experiment", y = NULL, title = "Selection Coefficients for SPF")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = col.palette)




