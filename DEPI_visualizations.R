###################################################
#
# DEPI Visualizations
#
##################################################

#Load necessary packages
library(dplyr)
library(ggplot2)
library(lemon)
library(data.table)
library(ggthemes)
library(extrafont)
library(Routliers)
library(stringi)

##################################
# Read in the data, and make subsets from the data
##################################

depi_data <-read.table("DEPI_analysis_Seeger.txt", sep = ",", header = FALSE)
names(depi_data) <- c("individual_plant_metadata", "genotype", "line", "subline", "border", "flat_number", "measurement_ID", "plant_ID", "measurement", "time_point", "measured_value")
depi_data$time_point <- as.numeric(gsub("X","", depi_data$time_point))
depi_data <- filter(depi_data, border == FALSE)
###Create a seperate data frame with the border plants
depi_border <- filter(depi_data, border == TRUE)
depi_subset <- depi_data%>%filter(!genotype %in% c("b1", "b3", "b1b3", "ftsz2-1", "ftsz2-2", "ftsz-dbl", "Col0" )|genotype == "Col0"&subline == "1")
levels(depi_subset$measurement)[levels(depi_subset$measurement)=="size"] <- "leafarea"
levels(depi_subset$measurement)[levels(depi_subset$measurement)=="growth"] <- "leafarea"
feb_data <- filter(depi_subset, substr(plant_ID, 1,4) == "0218", border == FALSE)
dec_data <- filter(depi_subset, substr(plant_ID, 1,4) == "1217")
feb_data <- feb_data%>%mutate(day = case_when(
  between(time_point,0,20) ~ 1,
  between(time_point,20,42)~2,
  between(time_point,42,68)~3,
  between(time_point,68,90)~4,
  between(time_point,90,115)~5,
  between(time_point,115,140)~6,
  between(time_point,140,162)~7,
  between(time_point,162,185)~8,
  between(time_point,185,212)~9,
  between(time_point,212,235)~10,
  between(time_point,235,260)~11,
  between(time_point,260,300)~12
))

feb_data <- feb_data%>%
  mutate(number = ifelse(genotype != "Col0",(stri_sub(genotype, 4, length(genotype))), 0))
feb_data$number <- as.numeric(gsub("-","0", feb_data$number))
feb_data$number[feb_data$number == "103"] <- "1003"
feb_data$number[feb_data$number == "506"] <- "5006"
feb_data$number[feb_data$number == "608"] <- "6008"
feb_data$number[feb_data$number == "609"] <- "6009"
feb_data$number <- as.numeric(feb_data$number)
feb_data <- feb_data%>%arrange(number)

feb_data <- feb_data%>%mutate(number_2 = number)
feb_data$number_2[nchar(feb_data$number_2) == 4] <- 0

##################################
# Investigate negative npq values
##################################
ggplot(data = filter(feb_data, measurement == "npq"), aes(x=measured_value))+
  geom_histogram(binwidth = 0.05)+
  labs(title = "Distribution of npq values", x = "NPQ")+
  geom_vline(xintercept = 0, color="red")

###What percent of measurements are negative?
nrow(filter(feb_data, measurement == "npq", measured_value <0))/nrow(feb_data)*100

###Shift all of the npq values up by the min value: -.3556
min_npq<-min((filter(feb_data, measurement == "npq"))$measured_value)
feb_data$measured_value[feb_data$measurement == "npq"] <- (feb_data$measured_value[feb_data$measurement == "npq"])+abs(min_npq)

###This is the new distribution of NPQ values:
ggplot(data = filter(feb_data, measurement == "npq"), aes(x=measured_value))+
  geom_histogram(binwidth = 0.01)+
  labs(title = "Distribution of npq values", x = "NPQ")+
  geom_vline(xintercept = 0, color="red")

###Finally, npq values arranged in increasing order, to get an idea of the minimum values
feb_data%>%filter(measurement == "npq")%>%select(measured_value, measurement)%>%arrange(measured_value)
##################################
# Preliminary Visualizations
##################################

#Filtered data to look at leaf area for mpk17, mpk1, mpk1-17 for February
mpk1_17_subset <- feb_data%>%
  filter(genotype %in% c("mpk17", "mpk1", "mpk1-17", "Col0"), measurement == "leafarea")%>%
  arrange(time_point)%>%
  dplyr::select(c("genotype", "measurement", "measured_value", "individual_plant_metadata", "time_point"))%>%
  group_by(genotype, time_point, measurement)%>%
  summarize(leafarea = median(measured_value))
#Scatter plot for leaf area; shows gaps in data (different days, so different light regimes)
ggplot(data = mpk1_17_subset, aes(x = time_point, y = leafarea)) + 
  geom_point(aes(color = genotype))+
  labs(x = "Time (Hours)", y = "Leaf area (Units)") + 
  scale_color_discrete(name = "Genotype")+
  #ylim(0,400)+
  ggtitle("Leaf area: mpk1, mpk17, mpk1-17, and Col-0 (February 2018)")
#Line plot for leaf area 
ggplot(data = mpk1_17_subset, aes(x = time_point, y = leafarea)) + 
  geom_line(aes(color = genotype)) +
  labs(x = "Time (Hours)", y = "Leaf area (Units)") + 
  scale_color_discrete(name = "Genotype")+
  #ylim(0,400)+
  ggtitle("Leaf area: mpk1, mpk17, mpk1-17, and Col-0 (February 2018)")

###############################################
# Create plot mirrored after page 370 in the Cell paper
###############################################

cell_370_data <- feb_data%>%
  filter(genotype ==("Col0"))%>%
  group_by(genotype, time_point, measurement, day, number)%>%
  summarize(median = median(measured_value))
ggplot(data = cell_370_data, aes(x = time_point, y = median))+
  geom_line()+
  #geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper) ,fill="green", alpha=0.2)+
  facet_rep_grid(measurement ~ day, scales = "free" , switch = "y", repeat.tick.labels = FALSE)+
  labs(x = "Hours", y = NULL)+
  theme_tufte(base_family = "Calibri",
              base_size = 18)+
  theme(strip.background.x = element_blank(),
         panel.border = element_rect(color = "black", fill = NA, size = 1),
         axis.line=element_line(),
        panel.spacing = unit(1, "lines"))+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))
  
cell_370_data_comp<- feb_data%>%
  filter(genotype %in% c("Col0", "mpk1", "mpk17", "mpk1-17"))%>%
  group_by(genotype, time_point, measurement, day)%>%
  summarize(median = median(measured_value),CI_upper = (mean(measured_value) + qt((1-0.95)/2, df=length(measured_value)-1)*(sd(measured_value)/sqrt(length(measured_value)-1))), CI_lower = (mean(measured_value) - qt((1-0.95)/2, df=length(measured_value)-1)*(sd(measured_value)/sqrt(length(measured_value)-1))))
###############################################
# Here, the same plot, but comparing 4 genotypes 
###############################################
ggplot(data = cell_370_data_comp, aes(x = time_point, y = median))+
  geom_line(aes(color = genotype), size = 1)+
  facet_rep_grid(measurement ~ day, scales = "free" , switch = "y", repeat.tick.labels = FALSE)+
  #geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper, fill = genotype), alpha=0.2)+
  labs(x = "Hours", y = NULL)+
  theme_tufte(base_family = "Calibri",
              base_size = 18)+
  theme(strip.background.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line=element_line(),
        panel.spacing = unit(1, "lines"))+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  scale_color_viridis_d(begin = 0, end = 1, option = 'viridis', aesthetics = c("colour", "fill"))

###############################################
# Finally, without standard error 
###############################################
ggplot(data = cell_370_data_comp, aes(x = time_point, y = median))+
  geom_line(aes(color = genotype), size = 1)+
  facet_rep_grid(measurement ~ day, scales = "free" , switch = "y", repeat.tick.labels = FALSE)+
  labs(x = "Hours", y = NULL)+
  theme_tufte(base_family = "Calibri",
              base_size = 18)+
  theme(strip.background.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line=element_line(),
        panel.spacing = unit(1, "lines"))+
  
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  scale_color_viridis_d(begin = 0, end = 1, option = 'viridis', aesthetics = c("colour", "fill"))


####################################
# Error is large for leaf area. Why?
####################################

feb_data%>%
  filter(genotype == "Col0", measurement == "leafarea", day %in% c(1:12))%>%
  dplyr::select(c("time_point", "day","plant_ID", "measured_value"))%>%
  group_by(day, plant_ID)%>%
  summarize(min = min(measured_value), max = max(measured_value))%>%
  print(n =108)
###There is one plant that has a way bigger leaf area than the others,
###especially after day 5. Possibly exclude this plant as an outlier. 

###############################################
# Create plot mirrored after page 371 in the Cell paper ----
# This  visualization is correct!
###############################################
cell_371_data<-feb_data%>%
  group_by(time_point, measurement)%>%
  mutate_each(funs(./median(.[genotype == "Col0"])), measured_value)%>%
  group_by(time_point, measurement, genotype)%>%
  mutate(log2_fold = log2(median(measured_value)))

###############################################
# Verify that the log 2 fold values are correct
# (I just picked random genotypes/times/measurements)
###############################################
###Good job!!!! This works!!!!!
a <- filter(feb_data, genotype == "mpk1", time_point == "1", measurement == "npq")$measured_value
b <- filter(feb_data, genotype == "Col0", time_point == "1", measurement == "npq")$measured_value
c <- log2(median(a/median(b)))
filter(cell_371_data, genotype == "mpk1", time_point == "1", measurement == "npq")$log2_fold
d <- filter(feb_data, genotype == "mpk6-18", time_point == "273.9997", measurement == "leafarea")$measured_value
e <- filter(feb_data, genotype == "Col0", time_point == "273.9997", measurement == "leafarea")$measured_value
f <- log2(median(d/median(e)))
filter(maybe, genotype == "mpk6-18", time_point == "273.9997", measurement == "leafarea")$log2

###############################################
# Here, find the start and end points of each day to use for the x-axis labels
###############################################
cell_371_data%>%group_by(day)%>%summarize(min(time_point), max(time_point))
###############################################
# Heat map for all genotypes
###############################################
cell_371_data_leafarea_all <- cell_371_data%>%filter(measurement == "leafarea")%>%arrange(number)
cell_371_data_npq_all <- cell_371_data%>%filter(measurement == "npq")%>%arrange(number)
cell_371_data_phi2_all <- cell_371_data%>%filter(measurement == "phi2")%>%arrange(number)

cell_371_data_leafarea_all$genotype <- reorder(cell_371_data_leafarea_all$genotype, cell_371_data_leafarea_all$number)
cell_371_data_npq_all$genotype <- reorder(cell_371_data_npq_all$genotype, cell_371_data_npq_all$number)
cell_371_data_phi2_all$genotype <- reorder(cell_371_data_phi2_all$genotype, cell_371_data_phi2_all$number)
###############################################
# Phi2 - All genotypes ----
###############################################
ggplot(data = cell_371_data_phi2_all, aes(x = time_point, y = genotype, fill = log2_fold)) + 
  labs(fill = "Log 2 Fold Change", x = "Hours", y = NULL, title = "phi2 Log 2 Fold Change")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
               base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
###############################################
# Leafarea - All genotypes ----
###############################################
ggplot(data = cell_371_data_leafarea_all, aes(x = time_point, y = genotype, fill = log2_fold)) + 
  labs(fill = "Log 2 Fold Change", x = "Hours", y = NULL, title = "Leaf Area Log 2 Fold Change")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
###############################################
# NPQ - All genotypes
###############################################
ggplot(data = cell_371_data_npq_all, aes(x = time_point, y = genotype, fill = log2_fold)) + 
  labs(fill = "Log 2 Fold Change", x = "Hours", y = NULL, title = "NPQ Log 2 Fold Change")+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
###############################################
# Loop to create plots for all combinations of single/double mutants
# Cell 370
###############################################
genotype_combinations <- list(c("mpk1-17", "mpk1", "mpk17"), c("mpk1-16", "mpk1", "mpk16"), c("mpk6-9", "mpk6", "mpk9"),
                              c("mpk17-20","mpk17", "mpk20"), c("mpk14-17", "mpk14", "mpk17"), c("mpk8-17", "mpk8","mpk17"),
                              c("mpk8-20","mpk8","mpk20"), c("mpk6-18", "mpk6", "mpk18"), c("mpk1-13", "mpk1", "mpk13"), c("mpk17-20", "mpk17", "mpk20"), c("mpk13-20", "mpk13", "mpk20"), c("mpk6-8", "mpk6", "mpk8"), c("mpk9-18", "mpk9", "mpk18"),
                              c("mpk6-20", "mpk6", "mpk20"), c("mpk14-16", "mpk14", "mpk16"), c("mpk18-20", "mpk18", "mpk20"), c("mpk5-6", "mpk5", "mpk6"), c("mpk14-18", "mpk14", "mpk18"), c("mpk5-6", "mpk5", "mpk6"), c("mpk14-18", "mpk14", "mpk18"), c("mpk5-17", "mpk5", "mpk17"),
                              c("mpk1-3", "mpk1", "mpk3"), c("mpk1-17", "mpk1", "mpk17"), c("mpk3-16", "mpk3", "mpk16"), c("mpk9-16", "mpk9", "mpk16"), c("mpk14-20", "mpk14", "mpk20"))
for (element in genotype_combinations){
  data <- feb_data%>%
    filter(genotype %in% c(element, "Col0"))%>%
    group_by(genotype, time_point, measurement, day)%>%
    summarize(median = median(measured_value),CI_upper = (mean(measured_value) + qt((1-0.95)/2, df=length(measured_value)-1)*(sd(measured_value)/sqrt(length(measured_value)-1))), CI_lower = (mean(measured_value) - qt((1-0.95)/2, df=length(measured_value)-1)*(sd(measured_value)/sqrt(length(measured_value)-1))))
  plot <- ggplot(data = data, aes(x = time_point, y = median))+
    geom_line(aes(color = genotype), size = 1)+
    facet_rep_grid(measurement ~ day, scales = "free" , switch = "y", repeat.tick.labels = FALSE)+
    labs(x = "Hours", y = NULL)+
    theme_tufte(base_family = "Calibri",
                base_size = 18)+
    theme(strip.background.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.line=element_line(),
          panel.spacing = unit(1, "lines"))+
    scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
    scale_color_viridis_d(begin = 0, end = 1, option = 'viridis', aesthetics = c("colour", "fill"))
  print(plot)}
###############################################
# Loop to create plots for all combinations of single/double mutants
# Cell 371
###############################################
cell_371_data <- data.frame(feb_data%>%
                              group_by(genotype, time_point, measurement, day, number_2)%>%
                              summarize(median = median(measured_value),CI_upper = (mean(measured_value) + qt((1-0.95)/2, df=length(measured_value)-1)*(sd(measured_value)/sqrt(length(measured_value)-1))), CI_lower = (mean(measured_value) - qt((1-0.95)/2, df=length(measured_value)-1)*(sd(measured_value)/sqrt(length(measured_value)-1))))%>%
                              group_by(genotype)%>%
                              group_by(time_point)%>%
                              mutate(log2_fold = (log2((median)/(median[genotype == "Col0"])))))
for (element in genotype_combinations){
  for (m in c("npq", "leafarea", "phi2")){
    data <- filter(cell_371_data, genotype %in% element, measurement == m)
    data$genotype <- reorder(data$genotype, data$number_2)
    plot <- ggplot(data = data, aes(x = time_point, y = genotype, fill = log2_fold)) + 
      labs(fill = "Log 2 Fold Change", x = "Hours", y = NULL, title = paste(m, "Log 2 Fold Change"))+
      geom_tile(width = 10 , height = 20)+
      facet_grid(genotype ~ day, scales = "free", switch = "y")+
      scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
      theme_tufte(base_family = "Calibri",
                  base_size = 20)+
      theme(strip.background.y = element_blank(),
            strip.text.y = element_blank(),
            panel.spacing=unit(0, "lines"))+
      scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
    print(plot)
  }
}


###############################################
# Find the number of plants in each genotype in each flat
###############################################
plants_per_flat <- feb_data%>%
  filter(time_point == "0")%>%
  group_by(genotype, flat_number)%>%
  summarize(count = n())
hist(plants_per_flat$count)

invest <- feb_data%>%
  filter(time_point == "0")%>%
  group_by(genotype)%>%
  summarize(count = n())

time_4 <- feb_data%>%
  filter(time_point == "0")%>%
  group_by(genotype, measurement)%>%
  summarize(count = n())

filter(feb_data, genotype == "mpk1", time_point == "1")

sort(filter(feb_data, genotype == "Col0", measurement == "npq", time_point == "0")$measured_value)
f = sort(filter(feb_data, genotype == "mpk1", measurement == "npq", time_point == "0")$measured_value)/med_test
med_test = median(sort(filter(feb_data, genotype == "Col0", measurement == "npq", time_point == "0")$measured_value)
)

hist(filter(feb_data,measurement == "npq")$measured_value)
abline(v = 0.83)

cell_370_example<- feb_data%>%
  filter(genotype == "Col0", measurement == "npq")%>%
  group_by(genotype, time_point, measurement, day)%>%
  summarize(median = median(measured_value))


ggplot(data = cell_370_example, aes(x = time_point, y = median))+
  geom_line(aes(color = genotype), size = 1)+
  facet_rep_grid(measurement ~ day, scales = "free" , switch = "y", repeat.tick.labels = FALSE)+
  #geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper, fill = genotype), alpha=0.2)+
  labs(x = "Hours", y = NULL)+
  theme_tufte(base_family = "Calibri",
              base_size = 18)+
  theme(strip.background.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line=element_line(),
        panel.spacing = unit(1, "lines"))+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  scale_color_viridis_d(begin = 0, end = 1, option = 'viridis', aesthetics = c("colour", "fill"))



data <- filter(cell_371_data, genotype =="mpk1", measurement == "npq")
data$genotype <- reorder(data$genotype, data$number_2)
plot <- ggplot(data = data, aes(x = time_point, y = genotype, fill = log2_fold)) + 
  labs(fill = "Log 2 Fold Change", x = "Hours", y = NULL, title = paste(m, "Log 2 Fold Change"))+
  geom_tile(width = 10 , height = 20)+
  facet_grid(genotype ~ day, scales = "free", switch = "y")+
  scale_x_continuous(breaks = round(c(0,15,24,39.5,48,63.7,72,87,96,112,120,135,144,159,168,183,192,207,216,231,240,255,264,279),0))+
  theme_tufte(base_family = "Calibri",
              base_size = 20)+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing=unit(0, "lines"))+
  scale_fill_gradient2(low = "blue", high="red", mid = "white", midpoint = 0)
print(plot)

