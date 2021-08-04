####################################################################################
# 
# The goal of this script is to first compute p-values for the data that is normalized by experiment.
# Then, compare these p-values to the p-values normalized by flat.
# Then, plot the distributions of the epistasis values for normalized by flat compared to normalized by experiment data.
# 
# Finally, use bootstapping to get the 0.025 and 0.975 quantiles of the selection coefficient ()
#
####################################################################################

###Load necessary packages
library(tidyverse)
library(stringi)
library(ggthemes)
library(viridis)
library(extrafont)
library(hexbin)
library(ggplot2)
library(gridExtra)
####################################################################################
### Functions to use to plot and calculate p-values
####################################################################################
###This function adds a column that is used to sort the Genotype in the plots:
add_number <- function(data_frame){
  ###First, if the Genotype is Col0 (only Genotype with length 4), assign 0 as number
  ###Else, assign number as Genotype with "mpk" removed
  ###Example: mpk1 will be 1, mpk1-17 will be 1-17
  data_frame <- data_frame%>%
    mutate(number = ifelse(Genotype != "Col",(stri_sub(Genotype, 4, length(Genotype))), 0))
  ###Next, for all double Genotypes, replace "-" with "0"
  ###Example: 1-17 becomes 1017
  data_frame$number <- as.numeric(gsub("_","0", data_frame$number))
  ###Almost there! There's a problem with two single digit double Genotypes
  ###We need a four digit number to sort correctly 
  ###Example: mpk1-3 -> 1-3 -> 103, but we need it to be 1003 to sort correctly
  data_frame$number[data_frame$number == "103"] <- "1003"
  data_frame$number[data_frame$number == "506"] <- "5006"
  data_frame$number[data_frame$number == "608"] <- "6008"
  data_frame$number[data_frame$number == "609"] <- "6009"
  ###Convert number to a numeric in order to sort
  data_frame$number <- as.numeric(data_frame$number)
  data_frame <- data_frame%>%arrange(number)
  data_frame <- data_frame%>%mutate(number_2 = number)
  data_frame$number_2[nchar(data_frame$number_2) == 4] <- 0
  data_frame$number_2[nchar(data_frame$number_2) == 5] <- 0
  return(data_frame)
}

### Note: p_value and correct_p_value are the same functions used in the DEPI analysis
p_value <- function (data_frame){
  ###Initialize an empty data frame
  out = data.frame()
  ###For each genotype:
  for (i in unique(filter(data_frame, (Genotype != "Col" & Genotype != "mpk5_17"))$Genotype)){
    ###We don't want to make comparisons of WT to itself - this could impact FDR correction
    indiv_data <- data_frame%>% 
      ### Group by measurement 
      group_by(Measurement)%>%
      ### Add a column with the p-value
      mutate(p = (wilcox.test(NormalizedExp[Genotype == i], NormalizedExp[Genotype == "Col"], correct = FALSE, paired = FALSE))$p.value)%>%
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
    group_by(Measurement)%>%
    mutate(p_adj = p.adjust(p, method = "fdr"))
  return(out) 
}

####################################################################################
### P-value calculation
####################################################################################
### Read in the data
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/Fitness_NormalizedValuesByExp_OutliersIncluded_07302021.csv", header = TRUE, row.names = 1)
### Apply the p_values function to the data
p_values <- p_value(data)  
### Adjust the p-values using an FDR correction
adjusted_p_values <- p_values%>%
  group_by(Measurement)%>%
  mutate(p_adj = round(p.adjust(p, method = 'fdr'), 7))%>%
  rename(p_ByExp = p,
         p_adj_ByExp = p_adj)
write.csv(adjusted_p_values, "PValues_NormalizedByFlatThenExp_07232021.csv", row.names = FALSE)
### Add a column with a boolean indicating if the adjusted p-value is significant or not
adjusted_p_values$significant_ByExp <- adjusted_p_values$p_adj_ByExp < 0.05
### Read in the p-value data from the normalization by flat, not exp
pVal_NormalizedFlat <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/FitnessPValues_07222021_OutliersRemoved.csv", header = TRUE, row.names = 1)

####################################################################################
### Compare p-values between normalization schemes
####################################################################################
colnames(adjusted_p_values)
colnames(pVal_NormalizedFlat)
plotData <- left_join(pVal_NormalizedFlat, adjusted_p_values, by = c("Experiment", "Genotype", "Measurement"))%>%
  add_number()

plotData$Genotype <- reorder(plotData$Genotype, desc(plotData$number))

colnames(plotData) <- c("Experiment", "Genotype", "Measurement",
                        "p_ByFlat", "pAdj_ByFlat", "Significant_ByFlat",
                        "p_ByExp", "pAdj_ByExp", "Significant_ByExp", "number", "number2")
### 1.9% of the p-values change from significant to not significant
sum(plotData$Significant_ByFlat ==  TRUE & plotData$Significant_ByExp == FALSE) / nrow(plotData)
### 10.5% of the p-values change from not significant to significant
sum(plotData$Significant_ByFlat ==  FALSE & plotData$Significant_ByExp == TRUE) / nrow(plotData)
### Histogram of the differences between p-values between normalization schemes
hist(plotData$pAdj_ByFlat - plotData$pAdj_ByExp,
     xlab = "Adj. P-Val. by Flat - Adj. P-Val. By Flat then Exp",
     main = "Differences in P-Value by Normalization Scheme")
par(mfrow = c(1,2))
hist(plotData$pAdj_ByExp,
     main = "Adj. P-Val After Normalization \nby Flat then Exp.",
     xlab = "Adj. P-Val")
hist(plotData$pAdj_ByFlat,
     main = "Adj. P-Val. After Normalization \nby Flat",
     xlab = "Adj. P-Val.")

plotData <- pivot_longer(plotData, 
                         cols = c("pAdj_ByFlat", "pAdj_ByExp"),
                         names_to = "Normalization_Method",
                         values_to = "pAdj")
### Add a column for plotting purposes; this is binary based on the significance of the p-value
plotData$sig <- ifelse(plotData$pAdj < 0.05, "p < 0.05", "p > 0.05")
### Plot the p-values from the two normalization schemes for each measurement and experiment
ggplot(data = plotData, aes(x = Experiment, y = Genotype, fill = sig)) +
  geom_tile()+
  facet_grid(~Normalization_Method + Measurement)+
  labs(fill = "Adjusted P-Value")+
  theme_minimal()+
  theme(text=element_text(size=12, family="Calibri"))+
  scale_fill_manual(values = c("red", "black"))

### Correlation of p-values between normalization schemes
scatterplotData <- left_join(pVal_NormalizedFlat, adjusted_p_values, by = c("Experiment", "Genotype", "Measurement"))%>%
  add_number()

ggplot(scatterplotData, aes(x=p_adj_ByExp, y=p_adj_OutliersRemoved)) +
  geom_point() +
  labs(x = "Adj. P, By Flat then Exp.",
       y = "Adj. P, By Flat")+
  theme_minimal()

####################################################################################
### Compare distribution of epistasis values by normalization scheme
####################################################################################

### Read in the data with the outliers included
sc.data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/selectionCoefficients_07292021_OutliersIncluded.csv", header = TRUE)
ep.data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/epistasis_07292021_OutliersIncluded.csv", header = TRUE)

### This gives 3 total plots: a grid of SC for TSC, SPF, and SN
for (i in c("SN", "TSC", "SPF")){
  plot1 <- ggplot(data = sc.data%>%filter(Experiment == "DEPI1" & Normalization == "ByFlat" & Measurement == i),
                  aes(x = SC_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Selection Coef. (Outliers Included)",
         y = NULL,
         title = paste(i, "Selection Coefficient: DEPI1", sep = " "))+
    theme_minimal()
  plot2 <- ggplot(data = sc.data%>%filter(Experiment == "DEPI2" & Normalization == "ByFlat" & Measurement == i), 
                  aes(x = SC_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Selection Coef. (Outliers Included)",
         y = NULL,
         title = paste(i, "Selection Coefficient: DEPI2", sep = " "))+
    theme_minimal()
  plot3 <- ggplot(data = sc.data%>%filter(Experiment == "DEPI3" & Normalization == "ByFlat" & Measurement == i), 
                  aes(x = SC_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Selection Coef. (Outliers Included)",
         y = NULL,
         title = paste(i, "Selection Coefficient: DEPI3", sep = " "))+
    theme_minimal()
  plot4 <- ggplot(data = sc.data%>%filter(Normalization == "ByExp" & Measurement == i),
                  aes(x = SC_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Selection Coef. (Outliers Included)",
         y = NULL,
         title = paste(i, "Selection Coefficient:Normalized by Exp.", sep = " "))+
    theme_minimal()
  
  finalplot <- grid.arrange(arrangeGrob(plot1, plot2, plot3), plot4, ncol = 2)
  print(finalplot)
  tempName <- paste("SC", i, ".png", sep = "")
  ggsave(filename = tempName,
         finalplot,
         path = "C:/Users/Owner/Documents/Research/Shiu_Lab/Shiu_Lab_R/Fitness/EpistasisDist_ComparingNormalizationSchemes")
}


### This gives 3 total plots: a grid of add. ep. for TSC, SPF, and SN
for (i in c("SN", "TSC", "SPF")){
  plot1 <- ggplot(data = ep.data%>%filter(Experiment == "DEPI1" & Normalization == "ByFlat" & Measurement == i),
                  aes(x = AdditiveEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Add. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Additive Epistasis: DEPI1", sep = " "))+
    theme_minimal()
  plot2 <- ggplot(data = ep.data%>%filter(Experiment == "DEPI2" & Normalization == "ByFlat" & Measurement == i), 
                  aes(x = AdditiveEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Add. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Additive Epistasis: DEPI2", sep = " "))+
    theme_minimal()
  plot3 <- ggplot(data = ep.data%>%filter(Experiment == "DEPI3" & Normalization == "ByFlat" & Measurement == i), 
                  aes(x = AdditiveEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Add. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Additive Epistasis: DEPI3", sep = " "))+
    theme_minimal()
  plot4 <- ggplot(data = ep.data%>%filter(Normalization == "ByExp" & Measurement == i),
                  aes(x = AdditiveEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Add. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Additive Epistasis: Normalized by Exp.", sep = " "))+
    theme_minimal()
  
  finalplot <- grid.arrange(arrangeGrob(plot1, plot2, plot3), plot4, ncol = 2)
  print(finalplot)
  tempName <- paste("AddEpistasis", i, ".png", sep = "")
  ggsave(filename = tempName,
         finalplot,
         path = "C:/Users/Owner/Documents/Research/Shiu_Lab/Shiu_Lab_R/Fitness/EpistasisDist_ComparingNormalizationSchemes")  
}


### This gives 3 total plots: a grid of prop. ep. for TSC, SPF, and SN
for (i in c("SN", "TSC", "SPF")){
  plot1 <- ggplot(data = ep.data%>%filter(Experiment == "DEPI1" & Normalization == "ByFlat" & Measurement == i),
                  aes(x = PropEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Prop. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Proportional Epistasis: DEPI1", sep = " "))+
    theme_minimal()
  plot2 <- ggplot(data = ep.data%>%filter(Experiment == "DEPI2" & Normalization == "ByFlat" & Measurement == i), 
                  aes(x = PropEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Prop. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Proportional Epistasis: DEPI2", sep = " "))+
    theme_minimal()
  plot3 <- ggplot(data = ep.data%>%filter(Experiment == "DEPI3" & Normalization == "ByFlat" & Measurement == i), 
                  aes(x = PropEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Prop. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Proportional Epistasis: DEPI3", sep = " "))+
    theme_minimal()
  plot4 <- ggplot(data = ep.data%>%filter(Normalization == "ByExp" & Measurement == i),
                  aes(x = PropEp_OutliersIncluded))+
    geom_histogram(bins = 10)+
    labs(x = "Prop. Ep. (Outliers Included)",
         y = NULL,
         title = paste(i, "Proportional Epistasis: Normalized by Exp.", sep = " "))+
    theme_minimal()
  
  finalplot <- grid.arrange(arrangeGrob(plot1, plot2, plot3), plot4, ncol = 2)
  print(finalplot)
  tempName <- paste("PropEpistasis", i, ".png", sep = "")
  ggsave(filename = tempName,
         finalplot,
         path = "C:/Users/Owner/Documents/Research/Shiu_Lab/Shiu_Lab_R/Fitness/EpistasisDist_ComparingNormalizationSchemes")
}

####################################################################################
### Bootstrapping - Col0
####################################################################################
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/Fitness_NormalizedValuesByExp_OutliersRemoved_07232021.csv", header = TRUE, row.names = 1)
### Set the graphing parameters to have 3 rows and 3 columns of plots
par(mfrow = c(1, 3))
### Initialize a matrix to populate with the 2.5% quantile and 97.5% quantile for each experiment and measurement
quantile_matrix <- matrix(nrow = 3,
                          ncol = 3,
                          dimnames = list(NULL,
                                          c("Measurement", "Quantile_0.025", "Quantile_0.975"))
)
### Initialize a row count variable set to 1
row_count <- 1
### For each measurement
for (k in unique(data$Measurement)){
  ### Filter to extract the normalized measured values of Col0 for that measurement and experiment
  temp_bootstrapping_data <- filter(data, Measurement == k, Genotype == "Col")$NormalizedExp
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
  hist(temp_mean_deviation_vector,
       xlab = " ",
       main = paste(k, sep = " "))
  abline(v = mean(temp_mean_deviation_vector), col = "red", lty = 2)

  ### Add the measurement to the first column
  quantile_matrix[row_count, 1] <- k
  ### Add the 2.5% percentile of the 10000 bootstrap values to the second column
  quantile_matrix[row_count, 2] <- quantile(temp_mean_deviation_vector, probs = c(0.025))
  ### Add the 97.5% percentile of the 10000 bootstrap values to the third column
  quantile_matrix[row_count, 3] <- quantile(temp_mean_deviation_vector, probs = c(0.975))
  ### Re-assign the row_count variable so the next iteration populates the next row in the matrix
  row_count <- row_count + 1
}
####################################################################################
### Similar approach - each genotype
####################################################################################

### Create a vector of the unique genotypes in the experiment
### mpk5-17 is missing from DEPI1. So, remove this genotype for now so we're able to compare genotypes between experiments
genotype_vector <- unique(data$Genotype)
### Initialize a matrix of NA values to populate with mean(measured_values - mean(Col0 measured value)) for each measurement
selectionCoef_matrix <- matrix(data = NA, 
                               nrow = 3,
                               ncol = 38,
                               dimnames = list(NULL,
                                               append(as.vector(genotype_vector), c("Measurement"))))
### Initialize matrix_count to be 1; this metric is used as an indicator of the row number to populate in the matrix
matrix_count <- 1
### For each measurement (SN, TSC, SPF)
for (k in unique(data$Measurement)){
  ### Filter the data to include the jth experiment and the kth measurement
  temp_data <- filter(data, Measurement == k)
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
          (filter(temp_data, Genotype == genotype_vector[i])$NormalizedExp - mean(filter(temp_data, Genotype == "Col")$NormalizedExp))
          / mean(filter(temp_data, Genotype == "Col")$NormalizedExp)), 6)
  }
  ### Add measurement to the 38th element of temp_vector
  temp_vector[i + 1] <- k
  ### temp_vector now becomes a row of the matrix
  selectionCoef_matrix[matrix_count, ] <- temp_vector
  ### Re-assign the matrix_count variable so the next iteration populates the next row in the matrix
  matrix_count <- matrix_count + 1
}

####################################################################################
### Combine matrices with the quantiles from Col0 and the selection coefficients for each genotype
####################################################################################

### First, convert the matrices to dataframes in order to use the cbind() command
### Note that the experiment and measurement columns will be duplicated, but this is a good check to make sure the data frames are correctly joined. (The final data type is a data frame
exceeds_quantile_df <- cbind(as.data.frame(selectionCoef_matrix), as.data.frame(quantile_matrix))
### The selection coefficents and quantiles are characters. Convert these columns to numeric. 
exceeds_quantile_df[, 1:37] <- lapply(exceeds_quantile_df[, 1:37], as.numeric)
heatmapData <- exceeds_quantile_df
for (i in 1:nrow(exceeds_quantile_df)){
  ### For each of the columns that are the selection coefficent calculations for each genotype
  for (j in 1:(ncol(exceeds_quantile_df) - 5)){
    ### If the value in the cell is less than the 2.5% quantile for Col0
    if (exceeds_quantile_df[i,j] < exceeds_quantile_df[i,40]){
      ### Replace the value with the statement -1000
      ### Note: I initially included a character string to replace these values, but am using a numeric to make it easier to plot the data later
      exceeds_quantile_df[i,j] <- -1000
    } 
    ### Else if the value the cell is greater than the 97.5% for Col0
    else if (exceeds_quantile_df[i,j] > exceeds_quantile_df[i,41]){
      ### Replace the value with the statement 1000
      ### Note: I initially included a character string to replace these values, but am using a numeric to make it easier to plot the data later
      exceeds_quantile_df[i,j] <- 1000
    }
  }
}


### Convert the data to the correct format to use to create heatmaps in ggplot:
### Remove the duplicate experiment and measurement columns, as well as the columns with the quantiles
heatmapData <- heatmapData[, 1:38]
heatmapData <- heatmapData%>%
  ### Pivot the data frame so the columns for each genotype become a "key" column
  ### i.e. convert the dataframe from wide to long format
  pivot_longer(cols = c(starts_with("mpk"), "Col"),
               names_to = "genotype",
               values_to = "selectionCoef"
  )%>%
  ### Remove the rows with Col, because the selection coefficient is always 0
  ### There is no biological meaning because Col is compared to itself
  filter(genotype != "Col")
### Rename the genotype column to Genotype in order to use the add_number function
heatmapData <- heatmapData%>%
  rename(Genotype = genotype)
### Apply the add_number function above to the data frame
heatmapData <- add_number(heatmapData)
### Reorder genotype by the "number" column. Multiply this by -1, because the order was opposite of what I wanted it to be
heatmapData$genotype <- reorder(heatmapData$Genotype, heatmapData$number * -1)
### Find the range of selection coefficients so I can set the scale of each plot to be the same, so we can easily compare between them
lowerbound <- floor(min(heatmapData$selectionCoef) * 100)/100
upperbound <- ceiling(max(heatmapData$selectionCoef) * 100)/100
### Create heatmap of all measurements:
ggplot(data = heatmapData, aes(x = as.factor(Measurement), y = as.factor(genotype), fill = selectionCoef)) + 
  labs(fill = "Selection Coefficient", x = NULL, y = NULL, title = "Selection Coefficient - Normalization by Exp.")+
  geom_tile()+
  theme_tufte(base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0, limits = c(lowerbound, upperbound), breaks = c(lowerbound, 0, upperbound), labels = c("-0.39", "0", "0.25"))
####################################################################################
### Create heat maps - binned data
####################################################################################
### Here, recreate the heat maps with binned data
### This allows us to highlight the cells that either are less than the 0.025 quantile or exceed the 0.975 quantile
heatmap_data_binned <- exceeds_quantile_df[, 1:38]%>%
  ### Pivot the data frame so the columns for each genotype become a "key" column
  ### i.e. convert the dataframe from wide to long format
  pivot_longer(cols = c(starts_with("mpk"), "Col"),
               names_to = "genotype",
               values_to = "selectionCoef"
  )%>%
  ### Remove the rows with Col, because the selection coefficient is always 0
  ### There is no biological meaning because Col is compared to itself
  filter(genotype != "Col")### Create the bins:
### set up cut-off values 
### NOTE: -1000 represents the selection coeffients that are less than 0.025 quantile and 1000 represents the selection coefficients that exceed the 0.975 quantile
breaks <- c(-1200,-0.15,-0.05,0.05, 0.15, 1200)
# specify interval/bin labels
tags <- c("Lower than 0.025 Quantile","[-0.15, -0.05)", "[-0.05, 0.05)", "[0.05, 0.15)","Exceeds 0.975 Quantile")
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
heatmap_data_binned <- heatmap_data_binned%>%
  rename(Genotype = genotype)
### Apply the add_number function above to the data frame
heatmap_data_binned <- add_number(heatmap_data_binned)
### Reorder genotype by the "number" column. Multiply this by -1, because the order was opposite of what I wanted it to be
heatmap_data_binned$Genotype <- reorder(heatmap_data_binned$Genotype, heatmap_data_binned$number * -1)
### Create a color palette with white in the middle and blue and red at the extremes
### This palettes is specifically for the 5 bins
col.palette <- rev(c("#FF0000", 
                     "#FFCCCB",
                     "#FFFEFE",
                     "#D4EBF2",
                     "#0000FF"))
### Create binned SN heatmap
ggplot(data = heatmap_data_binned, aes(x = Measurement, y = Genotype, fill = bins)) + 
  labs(fill = "Selection Coefficient", x = "Measurement", y = NULL, title = "Selection Coefficients: Normalized by Exp.")+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = col.palette)


####################################################################################
# Here, do Bootstrapping for each genotype, and compare the distribution to Col0
####################################################################################
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/Fitness_NormalizedValuesByExp_OutliersIncluded_07302021.csv", header = TRUE, row.names = 1)
### For each measurement
for (j in unique(data$Measurement)){
  ### Filter to extract the normalized measured values of Col0 for that measurement and experiment
  temp_bootstrapping_data <- filter(data, Measurement == j, Genotype == "Col")$NormalizedExp
  ### Initialize a vector of 10000 NAs 
  temp_mean_deviation_Col0 <- rep(NA, 10000)
  ### For 10000 bootstrap iterations
  for (i in 1:10000){
    ### replace the ith element in  temp_mean_deviation vector with the mean of the deviation from the sample to the mean of all Col0
    ### (There are a lot of operations strung together; I did this to avoid having to re-assign variables 1000 times to speed up the process)
    temp_mean_deviation_Col0[i] <-mean((
      ### Sample from the measured values associated with Col0 for the jth experiment and the kth measurement
      sample(temp_bootstrapping_data, 
             ### Sample 1/2 of the number of Col0 plants (Should be about 7, though smaller in some cases)
             size = floor(length(temp_bootstrapping_data) / 2),
             ### Sample with replacement 
             replace = TRUE) 
      ### Take the difference between the sample means and the mean of all Col0 plants for the jth experiment and kth measurement
      - mean(temp_bootstrapping_data)) / mean(temp_bootstrapping_data))
  }
  
  for (k in unique(data$Genotype)){
  ### Filter to extract the normalized measured values of Col0 for that measurement and experiment
  temp_bootstrapping_data <- filter(data, Measurement == j, Genotype == k)$NormalizedExp
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
  nf <- layout( matrix(c(1,2), ncol=2) )
  Col0Density <- density(temp_mean_deviation_Col0)
  GenotypeDensity <- density(temp_mean_deviation_vector)
  plot(Col0Density,
       xlab = " ",
       main = paste(j,": Col0 and ", k, " Density", sep = ""),
       col = "red",
       lty = 2)
       # ylim = max(Col0Density$y, GenotypeDensity$y))
  lines(GenotypeDensity,
       xlab = " ",
       main = paste(j, "Col0", sep = " "),
       col = "blue",
       lty = 3)

  }
}





