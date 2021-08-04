library(dplyr)
library(Routliers)

###Read in the function
remove_outliers<- function(df){
  out <- df%>%
    group_by(genotype, measurement, time_point)%>%
    mutate(measured_value = replace(measured_value, outliers_mad(measured_value, b=1.4826, threshold=3.5, na.rm=TRUE)$outliers_pos, NA))%>%
    arrange(genotype, measurement, time_point)
  return(out)}

###Create a temporary data frame with dummy data
column_1 <- rep("mpk1", 6)
column_2 <- rep("npq", 6)
column_3 <- rep(0, 6)
column_4 <- c(1,2,3,2,3,25)
temp_data <- as.data.frame(cbind(column_1, column_2, column_3, column_4))
colnames(temp_data) <- c("genotype", "measurement", "time_point", "measured_value")

temp_data$measured_value <- as.numeric(temp_data$measured_value)
temp_data$time_point <- as.numeric(temp_data$time_point)

###Does the remove_outliers function detect the 25?
remove_outliers(temp_data)

###Yes! 
###Then, why is the plant size still reporting values that I think should be removed?

depi_jan <- read.table("C:/Users/Owner/Documents/Research/Shiu_Lab/Shiu_Lab_R/Data/Correct_January_Created_Jan132021.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

###Filter to see where the leaf area measurements are 0
jan_leafarea <- depi_jan%>%
  filter(measurement == "growth")

depi_0_leafarea <- jan_leafarea%>%
  filter(measured_value == 0)

test_remove <- remove_outliers(jan_leafarea)

nrow(as.data.frame(test_remove%>%filter(measured_value == 0)))
nrow(depi_0_leafarea)

###So, remove_outliers only removes 3 0 measured values for plant size
###Investigate this

test_1 <- depi_jan%>%
  filter(measurement == "growth", genotype == "mpk8-17", time_point == 151)

remove_outliers(test_1)
outliers_mad(test_1$measured_value, b=1.4826, threshold=3.5, na.rm=TRUE)  

###Okay, so the function is correct, but the limits of acceptable range of values includes 0
###This may need to be modified in the future

###I suspect this is the same issue that the other plants with size 0 are encountering
###Potentially modify the threshold value of 3.5?








