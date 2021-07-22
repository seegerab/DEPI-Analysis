####################################################################################
# The goal of this script is to evaluate whether outliers need to be removed from the data.
# 
# To do this, I will compare the fitness values for the .csv file with the outliers included and removed.
#
# Next, I will see if the p-values for the fitness measurements differ between the .csv files with outliers included and removed.
#
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
  ###Convert number to a numberic in order to sort
  data_frame$number <- as.numeric(data_frame$number)
  data_frame <- data_frame%>%arrange(number)
  data_frame <- data_frame%>%mutate(number_2 = number)
  data_frame$number_2[nchar(data_frame$number_2) == 4] <- 0
  data_frame$number_2[nchar(data_frame$number_2) == 5] <- 0
  return(data_frame)
}

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

### Read in the 6 .csv files:

pval_OutliersIncluded <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/FitnessPValues_07222021_OutliersIncluded.csv", header = TRUE, row.names = 1)
pval_OutliersRemoved <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/FitnessPValues_07222021_OutliersRemoved.csv", header = TRUE, row.names = 1)
sc_OutliersIncluded <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/selectionCoefficients_07222021_OutliersIncluded.csv", header = TRUE)
sc_OutliersRemoved <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/selectionCoefficients_07222021_OutliersRemoved.csv", header = TRUE)
ep_OutliersIncluded <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/epistasis_07222021_OutliersIncluded.csv", header = TRUE)
ep_OutliersRemoved <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/epistasis_07222021_OutliersRemoved.csv", header = TRUE)


### Merge the data frames with the outliers included and outliers removed. The column names should indicate which column has the outliers removed or included

pval.data <- left_join(pval_OutliersIncluded, pval_OutliersRemoved, by = c("Experiment", "Genotype", "Measurement"))
sc.data <- left_join(sc_OutliersIncluded, sc_OutliersRemoved, by = c("Experiment", "Mutant", "Measurement"))
ep.data <- left_join(ep_OutliersIncluded, ep_OutliersRemoved, by = c("Experiment", "DoubleMutant", "Measurement", "MutantA", "MutantB"))

### Change the df with col name Mutant to have col name Genotype

### Add a column with "number" to sort the Genotypes by 
pval.data <- add_number(pval.data)
sc.data <- sc.data%>%
  rename(Genotype = Mutant)%>%
  add_number()
ep.data <- ep.data%>%
  rename(Genotype = DoubleMutant)%>%
  add_number()

### Re-order the Genotypes by "number
pval.data$Genotype <- reorder(pval.data$Genotype, desc(pval.data$number))
sc.data$Genotype <- reorder(sc.data$Genotype, desc(sc.data$number))
ep.data$Genotype <- reorder(ep.data$Genotype, desc(ep.data$number))

####################################################################################
# 
# P-value comparisons
#
####################################################################################

pval.data<- pval.data%>%
  mutate(pval.diff = p_adj_OutliersIncluded - p_adj_OutliersRemoved)

par(mfrow = c(1,3 ))
hist(filter(pval.data, Measurement == "SN")$pval.diff,
     xlab = "p-val (outliers included) - p-val (outliers removed)",
     main = "SN p-value comparison",
     breaks = 15)
hist(filter(pval.data, Measurement == "SPF")$pval.diff,
     xlab = "p-val (outliers included) - p-val (outliers removed)",
     main = "SPF p-value comparison",
     breaks = 15)
hist(filter(pval.data, Measurement == "TSC")$pval.diff,
     xlab = "p-val (outliers included) - p-val (outliers removed)",
     main = "TSC p-value comparison",
     breaks = 15)
### Plot the differences in p-values by measurement
ggplot(data = pval.data, aes(x = Experiment, y = Genotype, fill = pval.diff))+
  geom_tile()+
  facet_grid(~Measurement)+
  labs(fill = "pval (outliers included) - \npval (outliers removed)")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
  
####################################################################################
# 
# Epistasis and Selection Coefficient Comparisons
#
####################################################################################

ep.data <- ep.data%>%
  mutate(add.ep.diff = AdditiveEp_OutliersIncluded - AdditiveEp_OutliersRemoved)%>%
  mutate(prop.ep.diff = PropEp_OutliersIncluded - PropEp_OutliersRemoved)%>%
  na.omit()
sc.data <- sc.data%>%
  mutate(sc.diff = SC_OutliersIncluded - SC_OutliersRemoved)%>%
  na.omit()
### Plot the differences in additive epistasis by measurement
ggplot(data = ep.data, aes(x = Experiment, y = Genotype, fill = add.ep.diff))+
  geom_tile()+
  facet_grid(~Measurement)+
  labs(fill = "add. ep. (outliers included) - \nadd. ep. (outliers removed)")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
### Plot the differences in proportional epistasis by measurement
ggplot(data = ep.data, aes(x = Experiment, y = Genotype, fill = prop.ep.diff))+
  geom_tile()+
  facet_grid(~Measurement)+
  labs(fill = "prop. ep. (outliers included) - \nprop. ep. (outliers removed)")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
### Plot the differences in proportional epistasis by measurement
ggplot(data = sc.data, aes(x = Experiment, y = Genotype, fill = sc.diff))+
  geom_tile()+
  facet_grid(~Measurement)+
  labs(fill = "selection coef. (outliers included) - \nselection coef. (outliers removed)")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)





  

