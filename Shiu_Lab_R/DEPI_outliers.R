# This code is for all collection periods

##### -----Load necessary packages ------
library(dplyr)
library(tidyverse)
library(ggplot2)
###Lemon is used in ggplot2 - facet_rep_grid modification
library(lemon)
library(data.table)
library(ggthemes)
library(extrafont)
###Routliers is used for outliersmad to find outliers
library(Routliers)
library(stringi)

##### ----- Load in the data ------

depi_data <-read.table("DEPI_analysis_Seeger.txt", sep = ",", header = FALSE)
head(depi_data)
###add column names
names(depi_data) <- c("individual_plant_metadata", "genotype", "line", "subline", "border", "flat_number", "measurement_ID", "plant_ID", "measurement", "time_point", "measured_value")
###Some time points have an "X" in front of them
###Proportion of time points that have an "X":
sum(str_detect(depi_data$time_point, "X"))/nrow(depi_data)
###Rows that have an "X" in its time point:
depi_X <- depi_data[which(str_detect(depi_data$time_point, "X")),]
head(depi_X)
###Remove these X values
depi_data$time_point <- as.numeric(gsub("X","", depi_data$time_point))
depi_subset <- depi_data%>%filter(!genotype %in% c("b1", "b3", "b1b3", "ftsz2-1", "ftsz2-2", "ftsz-dbl", "Col0" )|(genotype == "Col0"&subline == "1"))
###The December collection period has "growth", while the February collection period has "size"
###Both of these measurements are recording leaf area, so change both to leafarea
levels(depi_subset$measurement)[levels(depi_subset$measurement)=="size"] <- "leafarea"
levels(depi_subset$measurement)[levels(depi_subset$measurement)=="growth"] <- "leafarea"
###Create two data sets based on collection period and remove border plants 
feb_data <- filter(depi_subset, substr(plant_ID, 1,4) == "0218", border == FALSE)
dec_data <- filter(depi_subset, substr(plant_ID, 1,4) == "1217", border == FALSE)

##### ----- Function - add column with day ------

add_day_col <- function(data_frame){
  unique_time <- sort(unique(data_frame$time_point))
  diff <- c()
  for (i in 1:length(unique_time)){
   if (i == 1){
      diff[1] <- 0
   }else{
      diff[i] <- unique_time[i] - unique_time[i-1]}}
  breaks <-c(0)
  for (i in 1:length(diff)){
    if (diff[i] > 5)
      breaks <- append(breaks, unique_time[i])}
  out <- data.frame()
  for (i in 1:length(breaks)){
    if (i == length(breaks)){
      indiv <- data_frame%>%filter(time_point >= breaks[i])%>%mutate(day = i)
   }else{
    indiv <- data_frame%>%filter(time_point >=breaks[i] & time_point < breaks[i+1])%>%mutate(day = i)}
    
    out <- rbind(as.data.frame(indiv), out)}
  return(out)}
##### ----- Function - shift npq ------
shift_npq <- function(data_frame){
  ###Some npq values are negative
  ###They should all be positive - shift all values by the lowest npq value (Siobhan is looking into this)
  #min_npq<-min((filter(data_frame, measurement == "npq"))$measured_value)
  data_frame$measured_value[data_frame$measurement == "npq"] <- (data_frame$measured_value[data_frame$measurement == "npq"])+abs(min((filter(data_frame, measurement == "npq"))$measured_value))
  return(data_frame)
}

##### ----- Function - replace outliers with NA ------
###Right now, not necessary because we are using nonparametric analysis
remove_outliers<- function(df){
  out <- df%>%
  group_by(genotype, measurement, time_point)%>%
  mutate(measured_value = replace(measured_value, outliers_mad(measured_value, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
  arrange(genotype, measurement, time_point)
  return(out)}
##### ----- Create a pipeline to find p-value ------
###First, verify that this notation will prevent the comparison of WT to itself
for (i in unique(filter(feb_data, genotype != "Col0")$genotype)){
  print(i)
}
###Test works - continue because we know Col0 won't be compared to itself
p_value <- function (data_frame){
  ###Initialize an empty data frame
  out = data.frame()
  ###For each genotype:
  for (i in unique(filter(data_frame, genotype != "Col0")$genotype)){
    ###We don't want to make comparisons of WT to itself - this could impact FDR correction
    indiv_data <- data_frame%>% 
      ###Focus on each time point and measurement
      group_by(time_point, measurement)%>%
      ###Create a column with the number of WT individual plants and the number of plants for each genotype
      ###Use this later to calculate effect size
      mutate(n_genotype = length(measured_value[genotype == i]), n_wt = length(measured_value[genotype =="Col0"]))%>%
      ###Create a column of p-values using a nonparametric Wilcox test
      
      ###Default set to exact = TRUE, because our sample sizes are too small to use a normal approximation
      ###But, when there are ties in the values (i.e. one value appears twice in the ranking process), wilcox.test returns to the normal approximation and spits out a warning message
      ###This may be a problem - include correct = FALSE to stop this from happening 
      
      ###More information:
      ###http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/RProgramming/NonParametricStats.html
      ###https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
      
      ###Paired = FALSE, because the Col0 plants are independet from each genotype
      ###Correct = FALSE turns off the continuity correction
      mutate(p = (wilcox.test(measured_value[genotype == i], measured_value[genotype == "Col0"], correct = FALSE, paired = FALSE))$p.value)%>%
  
      ###Add a column with each genotype
      mutate(genotype = i)%>%
      # mutate(number = unique(filter(feb_data, genotype == i)$number))%>%
      # mutate(number_2 = unique(filter(feb_data, genotype == i)$number_2))%>%
      # 
      select(time_point, genotype, measurement, day, p, n_wt, n_genotype)
  
      
    ###Add individual information to the main data frame
    out <- rbind(as.data.frame(indiv_data), out)}
  return(out)}

##### ----- Create a pipeline to correct p-value and add effect size------
###Adjust the p-values using an FDR correction
corrected_p_value <- function(data_frame){
  out <- data_frame%>%
  ###For some reason, I have multiple copies of each row
  distinct()%>%
  ###Group by time point and measurement - we are correcting by the number of genotypes 
  group_by(time_point, measurement)%>%
  mutate(p_adj = p.adjust(p, method = "fdr"))%>%
  ###Effect size calculated using method from:
  ###https://stats.stackexchange.com/questions/133077/effect-size-to-wilcoxon-signed-rank-test
  
  ###Only report an effect size if the p-value is significant; otherwise, NA
  mutate(effect = ifelse(p < 0.05, (abs(qnorm(p_adj))/sqrt(n_wt+n_genotype)), NA))%>%
  mutate(effect_size = case_when(
    ###Make sure these are the right cut offs for magnitude of effect size
    (effect <0.1)~"small",
    (effect>0.1 & effect < 0.5)~"medium",
    (effect>0.5)~"large"))
  ###Gather the data to make it easier to plot according to whether p was adjusted
  out <- gather(out, type, p, p, p_adj)%>%arrange(genotype, time_point)
  return(out) 
}

##### ----- Function - add number and number_2 column to sort plots later ------
add_number <- function(data_frame){
  data_frame <- data_frame%>%
    mutate(number = ifelse(genotype != "Col0",(stri_sub(genotype, 4, length(genotype))), 0))
  data_frame$number <- as.numeric(gsub("-","0", data_frame$number))
  data_frame$number[data_frame$number == "103"] <- "1003"
  data_frame$number[data_frame$number == "506"] <- "5006"
  data_frame$number[data_frame$number == "608"] <- "6008"
  data_frame$number[data_frame$number == "609"] <- "6009"
  data_frame$number <- as.numeric(data_frame$number)
  data_frame <- data_frame%>%arrange(number)
  data_frame <- data_frame%>%mutate(number_2 = number)
  data_frame$number_2[nchar(data_frame$number_2) == 4] <- 0
  return(data_frame)
}
##### ----- Apply functions to Feb data ------
feb_data <- add_day_col(feb_data)
feb_data <- shift_npq(feb_data)
feb_outliers <- remove_outliers(feb_data)
###Number of outliers
sum(is.na(feb_outliers$measured_value))
###Proportion of measured values that are outliers
sum(is.na(feb_outliers$measured_value))/nrow(feb_outliers)
###This is a relatively high proportion! Proceed with caution
feb_data_p <- p_value(feb_data)
###Tests to ensure that these p-values are correct
###Test 1 
filter(feb_data_p, genotype == "mpk1", time_point == "1", measurement == "npq")$p
a <- filter(feb_data, genotype == "mpk1", time_point == "1", measurement == "npq")$measured_value
b <- filter(feb_data, genotype == "Col0", time_point == "1", measurement == "npq")$measured_value
wilcox.test(a, b, correct = FALSE)$p.value
###Test 2
filter(feb_data_p, genotype == "mpk14-17", time_point == "0", measurement == "leafarea")$p
c <- filter(feb_data, genotype == "mpk14-17", time_point == "0", measurement == "leafarea")$measured_value
d <- filter(feb_data, genotype == "Col0", time_point == "0", measurement == "leafarea")$measured_value
wilcox.test(c,d, correct = FALSE)$p.value
###Test 3
filter(feb_data_p, genotype == "mpk13", time_point == "0", measurement == "phi2")$p
e <- filter(feb_data, genotype == "mpk13", time_point == "0", measurement == "phi2")$measured_value
f <- filter(feb_data, genotype == "Col0", time_point == "0", measurement == "phi2")$measured_value
wilcox.test(e,f, correct = FALSE)$p.value
feb_data_p_corrected <- corrected_p_value(feb_data_p)
###Notice that some effect sizes are infinite; this applies to p-values of 1; qnorm(1) = infinity
###Make the effect sizes factors, so we can create levels and the x-axis will be sorted correctly in ggplot2
feb_data_p_corrected$effect_size <- factor(feb_data_p_corrected$effect_size, levels = c("small", "medium", "large"))
###Distribution of effect sizes:
ggplot(data = subset(feb_data_p_corrected, !is.na(effect_size)), aes(x = effect_size, na.rm = TRUE))+
  geom_bar()
feb_data[(which(feb_data_p_corrected$effect==Inf)),]
###Summary statistics for the p values for each measurement;leaf area is never significant
feb_data_p_corrected%>%filter(type == "p_adj")%>%group_by(measurement)%>%summarize(median = median(p), min = min(p), max = max(p))
###Here, get counts for p-values, grouped by whether they were adjusted or not; notice that adjustment makes many significant p-values insignificant
feb_data_p_corrected%>%group_by(type)%>%summarize(p_lessthan_0.05 = sum(p<0.05), p_0.05_0.1 = sum(0.05<p & p<0.1), p_0.1_0.2 = sum(0.1<p & p<0.2))
###feb_data_plot is the final dataframe to use in subsequent plots
feb_data_plot <- add_number(feb_data_p_corrected)

##### ----- Feb p-value visualizations ------
###Plot the p-values seperated by the adjusted and initial values
ggplot(feb_data_plot, aes(x=p, color=type)) +
  geom_histogram(position="stack", fill = "white", size = 2)+
  facet_grid(~measurement)+
  ###Include a vertical line with the significant p-value alpha = 0.05
  geom_text(aes(x = 0.18, label = "p=0.05", y =  3000), color = "red", size = 8)+
  geom_vline(xintercept=0.05, color = "red", size = 1.5)+
  theme_igray(base_family = "Calibri",
              base_size = 20)+
  labs(title = "February: P-value Distribution")
ggsave("feb_adjusted_unadjusted_p_values.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")


###Create seperate data frames for each measurement to use in heat maps; only use adjusted p-values
feb_npq <- filter(feb_data_plot, measurement == "npq", type == "p")%>%
  ###Create a new column with the p-value "bins"
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
feb_phi2 <- filter(feb_data_plot, measurement == "phi2", type == "p")%>%
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
feb_leafarea <- filter(feb_data_plot, measurement == "leafarea", type == "p")%>%
  mutate(bin = case_when(
    (p <= 0.05)~ 'p<0.05',
    (p > 0.05)~ "p>0.05"))


###In order to use these bins as the fill in a heat map, convert to a factor
feb_npq$bin <- as.factor(feb_npq$bin)
feb_phi2$bin <- as.factor(feb_phi2$bin)
feb_leafarea$bin <- as.factor(feb_leafarea$bin)
###We want p<0.1 to be first in the legend, so refactor with p<0.01 as the first term
feb_npq$bin <- relevel(feb_npq$bin, 'p<0.01')
feb_phi2$bin <- relevel(feb_phi2$bin, 'p<0.01')
feb_leafarea$bin <- relevel(feb_leafarea$bin, 'p<0.05')
###Reorder by number so heat map has WT first, then single, then double mutants
feb_npq$genotype <- reorder(feb_npq$genotype, feb_npq$number)
feb_phi2$genotype <- reorder(feb_phi2$genotype, feb_phi2$number)
feb_leafarea$genotype <- reorder(feb_leafarea$genotype, feb_leafarea$number)

###Create a heat map of p-values for each measurement:
### x is time point
### y is genotype
### fill is p-value with FDR correction

##### ----- Feb p-value heat map - NPQ ------
ggplot(data = feb_npq, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwithout FDR Correction", x = "Hours", y = NULL, title = "February: NPQ P-value, Uncorrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("feb_npq_uncorrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")


##### ----- Feb p-value heat map - leaf area ------
ggplot(data = feb_leafarea, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwithout FDR Correction", x = "Hours", y = NULL, title = "February: Leaf Area P-value, Uncorrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
 # scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "black"))
ggsave("feb_leafarea_uncorrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Feb p-value heat map - phi2 ------
ggplot(data = feb_phi2, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwithout FDR Correction", x = "Hours", y = NULL, title = "February: Phi2 P-value, Uncorrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("feb_phi2_uncorrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")




###Create seperate data frames for each measurement to use in heat maps; only use adjusted p-values
feb_npq_adj <- filter(feb_data_plot, measurement == "npq", type == "p_adj")%>%
  ###Create a new column with the p-value "bins"
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
feb_phi2_adj <- filter(feb_data_plot, measurement == "phi2", type == "p_adj")%>%
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
feb_leafarea_adj <- filter(feb_data_plot, measurement == "leafarea", type == "p_adj")%>%
  mutate(bin = case_when(
    ###The minimum p-value for leaf area is 0.304, so adjust the bins:
    (p <0.4)~ 'p<0.4',
    (p>0.4 & p < 0.6)~ "0.4<p<0.6",
    (p>0.6)~"p>0.6"))


###In order to use these bins as the fill in a heat map, convert to a factor
feb_npq_adj$bin <- as.factor(feb_npq_adj$bin)
feb_phi2_adj$bin <- as.factor(feb_phi2_adj$bin)
feb_leafarea_adj$bin <- as.factor(feb_leafarea_adj$bin)
###We want p<0.1 to be first in the legend, so refactor with p<0.01 as the first term
feb_npq$feb_npq_adj <- relevel(feb_npq_adj$bin, 'p<0.01')
feb_phi2_adj$bin <- relevel(feb_phi2_adj$bin, 'p<0.01')
feb_leafarea_adj$bin <- relevel(feb_leafarea_adj$bin, 'p<0.4')
###Reorder by number so heat map has WT first, then single, then double mutants
feb_npq_adj$genotype <- reorder(feb_npq_adj$genotype, feb_npq_adj$number)
feb_phi2_adj$genotype <- reorder(feb_phi2_adj$genotype, feb_phi2_adj$number)
feb_leafarea_adj$genotype <- reorder(feb_leafarea_adj$genotype, feb_leafarea_adj$number)

###Create a heat map of p-values for each measurement:
### x is time point
### y is genotype
### fill is p-value with FDR correction

##### ----- Feb p-value heat map - NPQ ------
ggplot(data = feb_npq_adj, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwith FDR Correction", x = "Hours", y = NULL, title = "February: NPQ P-value, Corrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("feb_npq_corrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Feb p-value heat map - leaf area ------
###This plot doesn't really tell us anything, because nothing is significantly different from wildtype
ggplot(data = feb_leafarea_adj, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwith FDR Correction", x = "Hours", y = NULL, title = "February: Leaf Area P-value, Corrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("feb_leafarea_corrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Feb p-value heat map - phi2 ------
ggplot(data = feb_phi2_adj, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwith FDR Correction", x = "Hours", y = NULL, title = "February: Phi2 P-value, Corrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("feb_phi2_corrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")



##### ----- Apply functions to Dec data ------
dec_data <- add_day_col(dec_data)
dec_data <- shift_npq(dec_data)
dec_outliers <- remove_outliers(dec_data)
###Number of outliers
sum(is.na(dec_outliers$measured_value))
###Proportion of measured values that are outliers
sum(is.na(dec_outliers$measured_value))/nrow(dec_outliers)
###This is a relatively high proportion! Proceed with caution
dec_data_p <- p_value(dec_data)
###Tests to ensure that these p-values are correct
###Test 1 
filter(dec_data_p, genotype == "mpk1", time_point == "1", measurement == "npq")$p
a <- filter(dec_data, genotype == "mpk1", time_point == "1", measurement == "npq")$measured_value
b <- filter(dec_data, genotype == "Col0", time_point == "1", measurement == "npq")$measured_value
wilcox.test(a, b, correct = FALSE)$p.value
###Test 2
filter(dec_data_p, genotype == "mpk14-17", time_point == "29", measurement == "leafarea")$p
c <- filter(dec_data, genotype == "mpk14-17", time_point == "29", measurement == "leafarea")$measured_value
d <- filter(dec_data, genotype == "Col0", time_point == "29", measurement == "leafarea")$measured_value
wilcox.test(c,d, correct = FALSE)$p.value
###Test 3
filter(dec_data_p, genotype == "mpk13", time_point == "0", measurement == "phi2")$p
e <- filter(dec_data, genotype == "mpk13", time_point == "0", measurement == "phi2")$measured_value
f <- filter(dec_data, genotype == "Col0", time_point == "0", measurement == "phi2")$measured_value
wilcox.test(e,f, correct = FALSE)$p.value
dec_data_p_corrected <- corrected_p_value(dec_data_p)
###Notice that some effect sizes are infinite; this applies to p-values of 1; qnorm(1) = infinity
###Make the effect sizes factors, so we can create levels and the x-axis will be sorted correctly in ggplot2
dec_data_p_corrected$effect_size <- factor(dec_data_p_corrected$effect_size, levels = c("small", "medium", "large"))
###Distribution of effect sizes:
ggplot(data = subset(dec_data_p_corrected, !is.na(effect_size)), aes(x = effect_size, na.rm = TRUE))+
  geom_bar()
dec_data_p_corrected[(which(dec_data_p_corrected$effect==Inf)),]
###Summary statistics for the p values for each measurement;leaf area is never significant
dec_data_p_corrected%>%filter(type == "p_adj")%>%group_by(measurement)%>%summarize(median = median(p), min = min(p), max = max(p))
###Here, get counts for p-values, grouped by whether they were adjusted or not; notice that adjustment makes many significant p-values insignificant
dec_data_p_corrected%>%group_by(type)%>%summarize(p_lessthan_0.05 = sum(p<0.05), p_0.05_0.1 = sum(0.05<p & p<0.1), p_0.1_0.2 = sum(0.1<p & p<0.2))
###feb_data_plot is the final dataframe to use in subsequent plots
dec_data_plot <- add_number(dec_data_p_corrected)

##### ----- Dec p-value visualizations: adjusted ------
###Plot the p-values seperated by the adjusted and initial values
ggplot(dec_data_plot, aes(x=p, color=type)) +
  geom_histogram(position="stack", fill = "white", size = 2)+
  facet_grid(~measurement)+
  ###Include a vertical line with the significant p-value alpha = 0.05
  geom_text(aes(x = 0.18, label = "p=0.05", y =  3000), color = "red", size = 8)+
  geom_vline(xintercept=0.05, color = "red", size = 1.5)+
  theme_igray(base_family = "Calibri",
              base_size = 20)+
  labs(title = "December: P-value Distribution")
ggsave("dec_adjusted_unadjusted_p_values.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

ggplot(filter(dec_data_plot, measurement == "phi2"), aes(x=p, color=type)) +
  geom_histogram(position="stack", fill = "white", size = 2)+
  ###Include a vertical line with the significant p-value alpha = 0.05
  geom_text(aes(x = 0.18, label = "p=0.05", y =  3000), color = "red", size = 8)+
  geom_vline(xintercept=0.05, color = "red", size = 1.5)+
  theme_igray(base_family = "Calibri",
              base_size = 20)+
  labs(title = "December: Phi2 P-value")
ggsave("dec_adjusted_unadjusted_p_values_phi2.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

ggplot(filter(dec_data_plot, measurement == "npq"), aes(x=p, color=type)) +
  geom_histogram(position="stack", fill = "white", size = 2)+
  ###Include a vertical line with the significant p-value alpha = 0.05
  geom_text(aes(x = 0.18, label = "p=0.05", y =  3000), color = "red", size = 8)+
  geom_vline(xintercept=0.05, color = "red", size = 1.5)+
  theme_igray(base_family = "Calibri",
              base_size = 20)+
  labs(title = "December: NPQ P-value")
ggsave("dec_adjusted_unadjusted_p_values_npq.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")


###Create seperate data frames for each measurement to use in heat maps; only use adjusted p-values
dec_npq_adj <- filter(dec_data_plot, measurement == "npq", type == "p_adj")%>%
  ###Create a new column with the p-value "bins"
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
dec_phi2_adj <- filter(dec_data_plot, measurement == "phi2", type == "p_adj")%>%
  mutate(bin = case_when(
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
dec_leafarea_adj <- filter(dec_data_plot, measurement == "leafarea", type == "p_adj")%>%
  mutate(bin = case_when(
    ###The minimum p-value for leaf area is 0.304, so adjust the bins:
    (p <0.4)~ 'p<0.4',
    (p>0.4 & p < 0.6)~ "0.4<p<0.6",
    (p>0.6)~"p>0.6"))
###In order to use these bins as the fill in a heat map, convert to a factor
dec_npq_adj$bin <- as.factor(dec_npq_adj$bin)
dec_phi2_adj$bin <- as.factor(dec_phi2_adj$bin)
dec_leafarea_adj$bin <- as.factor(dec_leafarea_adj$bin)
###We want p<0.1 to be first in the legend, so refactor with p<0.01 as the first term
dec_npq_adj$bin <- relevel(dec_npq_adj$bin, 'p<0.01')
dec_phi2_adj$bin <- relevel(dec_phi2_adj$bin, "0.01<p<0.05")
dec_leafarea_adj$bin <- relevel(dec_leafarea_adj$bin, 'p<0.4')
###Reorder by number so heat map has WT first, then single, then double mutants
dec_npq_adj$genotype <- reorder(dec_npq$genotype, dec_npq_adj$number)
dec_phi2_adj$genotype <- reorder(dec_phi2_adj$genotype, dec_phi2_adj$number)
dec_leafarea_adj$genotype <- reorder(dec_leafarea_adj$genotype, dec_leafarea_adj$number)



###Create a heat map of p-values for each measurement:
### x is time point
### y is genotype
### fill is p-value with FDR correction

##### ----- Dec p-value heat map - NPQ ------
ggplot(data = dec_npq_adj, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwith FDR Correction", x = "Hours", y = NULL, title = "December: NPQ P-value, Corrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("dec_npq_corrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Dec p-value heat map - leaf area ------
###This plot doesn't really tell us anything, because nothing is significantly different from wildtype
ggplot(data = dec_leafarea_adj, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwith FDR Correction", x = "Hours", y = NULL, title = "December: Leaf Area P-value, Corrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("dec_leafarea_corrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Dec p-value heat map - phi2 ------
ggplot(data = dec_phi2_adj, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwith FDR Correction", x = "Hours", y = NULL, title = "December: Phi2 P-value, Corrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "black"))
ggsave("dec_phi2_corrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")



##### ----- Dec p-value visualizations: not adjusted ------
###Plot the p-values seperated by the adjusted and initial values

###Create seperate data frames for each measurement to use in heat maps; don't use adjusted p-values
dec_npq <- filter(dec_data_plot, measurement == "npq", type == "p")%>%
  ###Create a new column with the p-value "bins"
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
dec_phi2 <- filter(dec_data_plot, measurement == "phi2", type == "p")%>%
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
dec_leafarea <- filter(dec_data_plot, measurement == "leafarea", type == "p")%>%
  mutate(bin = case_when(
    (p <0.01)~ 'p<0.01',
    (p>0.01 & p < 0.05)~ "0.01<p<0.05",
    (p>0.05)~"p>0.05"))
###In order to use these bins as the fill in a heat map, convert to a factor
dec_npq$bin <- as.factor(dec_npq$bin)
dec_phi2$bin <- as.factor(dec_phi2$bin)
dec_leafarea$bin <- as.factor(dec_leafarea$bin)
###We want p<0.1 to be first in the legend, so refactor with p<0.01 as the first term
dec_npq$bin <- relevel(dec_npq$bin, 'p<0.01')
dec_phi2$bin <- relevel(dec_phi2$bin, 'p<0.01')
dec_leafarea$bin <- relevel(dec_leafarea$bin, 'p<0.01')
###Reorder by number so heat map has WT first, then single, then double mutants
dec_npq$genotype <- reorder(dec_npq$genotype, dec_npq$number)
dec_phi2$genotype <- reorder(dec_phi2$genotype, dec_phi2$number)
dec_leafarea$genotype <- reorder(dec_leafarea$genotype, dec_leafarea$number)



###Create a heat map of p-values for each measurement:
### x is time point
### y is genotype
### fill is p-value with FDR correction

##### ----- Dec p-value heat map - NPQ ------
ggplot(data = dec_npq, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwithout FDR Correction", x = "Hours", y = NULL, title = "December: NPQ P-value, Uncorrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("dec_npq_uncorrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Dec p-value heat map - leaf area ------
###This plot doesn't really tell us anything, because nothing is significantly different from wildtype
ggplot(data = dec_leafarea, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwithout FDR Correction", x = "Hours", y = NULL, title = "December: Leaf Area P-value, Uncorrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("dec_leafarea_uncorrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")

##### ----- Dec p-value heat map - phi2 ------
ggplot(data = dec_phi2, aes(x = time_point, y = genotype, fill = bin)) + 
  labs(fill = "P-Value, \nwithout FDR Correction", x = "Hours", y = NULL, title = "December: Phi2 P-value, Uncorrected")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  #scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_manual(values = c("red", "orange", "black"))
ggsave("dec_phi2_uncorrected.png", scale = 2, path = "C:/Users/Joan Seeger/Documents/Shiu Lab - R/Current/Visualizations")




##### ----- Normality investigations ------

###Start by conducting a shapiro test for the data
###Group the data by genotype, measurement, and time point - we are only looking at the measured 
###values for each of these groupings to see if they are normally distributed
# 
# ###http://www.sthda.com/english/wiki/normality-test-in-r
# ###How to interpret: we want p > 0.05 - this means that the data doesn't significanlty differ from a normal distribution
# normality_test_feb <- feb_data%>%
#   group_by(genotype, measurement, time_point)%>%
#   mutate(normality_p = shapiro.test(measured_value)$p.value)%>%
#   select(individual_plant_metadata, genotype, measurement, time_point, normality_p)
# 
# ###Look at the distribution of the p-values - we want all p-values to be above 0.05
# hist(normality_test_feb$normality_p)
# abline(v = 0.05, col = "red")
# 
# ###It looks like most are considered normal, but there is a good proportion that aren't
# 
# ###Here is the breakdown of significant/insignificant p-values
# ###Remember - want p > 0.05
# normality_test_feb %>%
#   ungroup()%>%
#   summarize(percent_sig = (sum(normality_p < 0.05)/nrow(normality_test)), percent_insig = (sum(normality_p > 0.05)/nrow(normality_test)))
# 
# ###Let's try and log transform the data to get more samples to achieve normality:
# normality_test_transform_feb <- feb_data%>%
#   group_by(genotype, measurement, time_point)%>%
#   mutate(normality_p = shapiro.test(measured_value)$p.value)%>%
#   mutate(log_normality_p = shapiro.test(log(measured_value))$p.value)
# 
# hist(normality_test_transform_feb$log_normality_p)
# abline(v = 0.05, col = "red")
# 
# normality_test_transform_feb %>%
#   ungroup()%>%
#   summarize(percent_sig = (sum(log_normality_p < 0.05, na.rm =TRUE)/nrow(normality_test_transform)), percent_insig = (sum(log_normality_p > 0.05, na.rm = TRUE)/nrow(normality_test_transform)))
# 

##### ----- Create a for loop to visualize p-values for each combination of single and double mutants ------
# 
# 
# genotype_combinations <- list(c("mpk1-17", "mpk1", "mpk17"), c("mpk1-16", "mpk1", "mpk16"), c("mpk6-9", "mpk6", "mpk9"),
#                               c("mpk17-20","mpk17", "mpk20"), c("mpk14-17", "mpk14", "mpk17"), c("mpk8-17", "mpk8","mpk17"),
#                               c("mpk8-20","mpk8","mpk20"), c("mpk6-18", "mpk6", "mpk18"), c("mpk1-13", "mpk1", "mpk13"), c("mpk17-20", "mpk17", "mpk20"), c("mpk13-20", "mpk13", "mpk20"), c("mpk6-8", "mpk6", "mpk8"), c("mpk9-18", "mpk9", "mpk18"),
#                               c("mpk6-20", "mpk6", "mpk20"), c("mpk14-16", "mpk14", "mpk16"), c("mpk18-20", "mpk18", "mpk20"), c("mpk5-6", "mpk5", "mpk6"), c("mpk14-18", "mpk14", "mpk18"), c("mpk5-6", "mpk5", "mpk6"), c("mpk14-18", "mpk14", "mpk18"), c("mpk5-17", "mpk5", "mpk17"),
#                               c("mpk1-3", "mpk1", "mpk3"), c("mpk1-17", "mpk1", "mpk17"), c("mpk3-16", "mpk3", "mpk16"), c("mpk9-16", "mpk9", "mpk16"), c("mpk14-20", "mpk14", "mpk20"))
# 
# ###Values below the red line are significant
# ###"Unorthodox" way of visualizing p-values, but still interesting to see how they change over time
# for (element in genotype_combinations){
#   indiv <- filter(corrected_p_value_feb, genotype %in% element)
#   plot <- ggplot(data = indiv, aes(x = time_point, y = p))+
#     geom_line(aes(color = genotype), size = 1)+
#     facet_rep_grid(measurement ~ day, scales = "free" , switch = "y", repeat.tick.labels = FALSE)+
#     labs(x = "Hours", y = NULL)+
#     theme_tufte(base_family = "Calibri",
#                 base_size = 18)+
#     theme(strip.background.x = element_blank(),
#           panel.border = element_rect(color = "black", fill = NA, size = 1),
#           axis.line=element_line(),
#           panel.spacing = unit(1, "lines"))+
#     scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
#     scale_color_viridis_d(begin = 0, end = 1, option = 'viridis', aesthetics = c("colour", "fill"))+
#     geom_hline(yintercept=0.05, color = "red", size = 1.5)
#   print(plot)
# }

##### ----- Everything below this needs to be organized! ------


###Investigate differences between experiments

###Create function to find number of plants per genotype
plants_per_genotype <- function(data_frame){
  data_frame%>%
  filter(time_point == "0")%>%
  group_by(genotype)%>%
  summarize(count = n()/3)%>%
  arrange(count)}

dec_plants <- plants_per_genotype(dec_data)
feb_plants <- plants_per_genotype(feb_data)
hist(feb_plants$count, main = "February - Number of plants per genotype", xlab = "Count", breaks = 10)
hist(dec_plants$count, main = "December - Number of plants per genotype", xlab = "Count", breaks = 5, xlim = c(0, 30))

summary(dec_plants$count)


###Better leaf area visualizations:

dec_leafarea_vis <- filter(dec_data_plot, measurement == "leafarea", type == "p_adj")%>%
  mutate(bin = case_when(
    ###The minimum p-value for leaf area is 0.304, so adjust the bins:
    (p <0.4)~ 'p<0.4',
    (p>0.4 & p < 0.6)~ "0.4<p<0.6",
    (p>0.6)~"p>0.6"))%>%
  group_by(day)%>%
  filter(time_point %in% c(min(time_point), max(time_point), ceiling(median(time_point))))



