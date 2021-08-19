####################################################################################
# 
# Create heat maps with the data normalized by experiment and the data normalized by flat.
# These heat maps will have 4 columns: DEPI1, DEPI2, and DEPI3 (normalized by flat) and the normalized by exp values
#
####################################################################################
library(tidyverse)
library(stringi)
library(ggthemes)
library(extrafont)
### Read in the data
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/csvFiles/epistasis_07292021_OutliersIncluded.csv", head = TRUE)
data$plotCol <- ifelse(data$Normalization == "ByFlat", data$Experiment, "Normalized By Exp.")
data <- data%>%
  rename(Genotype = DoubleMutant)%>%
  select(plotCol, Genotype, AdditiveEp_OutliersIncluded, PropEp_OutliersIncluded, Measurement)
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
### Add a column with a number to sort the genotypes by when creating the plot
data <- add_number(data)
### Sort the genotypes by this number
data$Genotype <- reorder(data$Genotype, desc(data$number))
### Create a plot for TSC
ggplot(data = filter(data, Measurement == "TSC"), aes(x = plotCol, y = Genotype, fill = AdditiveEp_OutliersIncluded))+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "TSC Additive Epistasis - Outliers Included")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
### Create a plot for SN
ggplot(data = filter(data, Measurement == "SN"), aes(x = plotCol, y = Genotype, fill = AdditiveEp_OutliersIncluded))+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "SN Additive Epistasis - Outliers Included")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
### Create a plot for SPF
ggplot(data = filter(data, Measurement == "SPF"), aes(x = plotCol, y = Genotype, fill = AdditiveEp_OutliersIncluded))+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "SPF Additive Epistasis - Outliers Included")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)

### Create a plot for TSC
ggplot(data = filter(data, Measurement == "TSC"), aes(x = plotCol, y = Genotype, fill = PropEp_OutliersIncluded))+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "TSC Proportional Epistasis - Outliers Included")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
### Create a plot for SN
ggplot(data = filter(data, Measurement == "SN"), aes(x = plotCol, y = Genotype, fill = PropEp_OutliersIncluded))+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "SN Proportional Epistasis - Outliers Included")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
### Create a plot for SPF
ggplot(data = filter(data, Measurement == "SPF"), aes(x = plotCol, y = Genotype, fill = PropEp_OutliersIncluded))+
  geom_tile()+
  theme_tufte(
    base_family = "Calibri")+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "SPF Proportional Epistasis - Outliers Included")+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)





























