####################################################################################
# 
# The goal is this script is to normalize by experiment, after we have already normalized
# by flat.
#
# This code is adapted from the script that normalizes by flat. 
#
# IN: .csv files of the quantile normalized values by flat for outliers removed and included
# OUT: .csv files of measured values quantile normalied by both flat and experiment
#
####################################################################################

### Read in the data that has already been quantile normalized by flat
data <- read.csv("~/Research/Shiu_Lab/Shiu_Lab_R/Data/Fitness_NormalizedValues_OutliersRemoved_07212021.csv", header = TRUE, row.names =  1)

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

####################################################################################
# Quantile normalization - experiment level
####################################################################################

###Create a loop for each measurement 
for (i in c("SN", "TSC", "SPF")){
  ###Filter to include the rows with measurement i
  measurement.data <- data%>%
    filter(Measurement == i)
  ###Initialize an empty data frame
  temp.df <- data.frame()
  ###Loop through each experiment 
  for (k in c("DEPI1", "DEPI2", "DEPI3")){
    ### Create a temporary vector of the measured values for each flat
    exp.measurements <- measurement.data%>%
      ### Filter to experiment k
      filter(Experiment == k)%>%
      select(NormalizedMeasuredValue)
    ### Add the exp.measurements to temp.df
    ### Note as.data.frame is needed because temp.vector is technically a tibble
    temp.df <- addToDF(temp.df, as.data.frame(exp.measurements)[,1])
  }
  ###Normalize across the flat (temp_normalize is a data frame, but temp.df needs to be converted to a matrix to use the normalize.quantiles function)
  temp_normalize <- as.data.frame(normalize.quantiles(as.matrix(temp.df)))%>%
    ###Add columns with the measurement and experiment to verify flat and experiment information
    mutate(measurement = i)
  temp_normalize <- temp_normalize%>%
    rename(DEPI1 = V1,
           DEPI2 = V2,
           DEPI3 = V3)
  ###Create a name to give the normalized data based on the measurement and experiment
  temp_name <- paste(i , "_normalize", sep = "")
  ###Assign the name to the data frame
  assign(temp_name, temp_normalize)
}

### Test that this worked:

####################################################################################
# Verify quantile normalization by experiment is correct!
####################################################################################


###Filter to only include Flat 2 and DEPI1 from the data frame with the outliers removed 
test <- subset(data, Experiment == "DEPI1" & Measurement == "SN")
###Add a column with the normalized value. These normalized values are the second column, with the string of NA values removed from the column
test$NormalizedSN <- SN_normalize[,1][1:470]
###Add a column with the count to the data frame
test$Count <- 1:nrow(test)
###Sort first by quantile normalized value, and then by the un-normalized value. Because the order is maintained when we quantile normalized, this should be the same data frame!
order1 <- test%>%
  arrange(NormalizedSN)
order2 <- test%>%
  arrange(NormalizedMeasuredValue)
### order1 == order2 is a dataframe of boolean values. If the data frames match, all values will be TRUE, so the sum of each column will be the same - the number of rows in the data frame. 
apply(order1 == order2, 2, sum)
### This looks good! Apply the same approach to the entire data frame. First, divide the data by flat and experiment, do the above approach, and re-combine the data.

####################################################################################
# Create final data frame with quantile normalized values by experiment
####################################################################################

###For each experiment, remove NA values from the data and add a column with the normalized values by experiment to the data
### Initialize an empty data frame for the DEPI1 experiment for the SN measurement
SN.temp <- data.frame()
### For each experiment
for (i in 1:3){
  temp.exp <- c("DEPI1", "DEPI2", "DEPI3")[i]
  ### Filter to the ith flat and the DEPI1 experiment
  ###Remove rows that have NA (remember, these are the NA values appended to the end of the rows in order to combine flats that have a different number of plants)
  tempSN <- data%>%
    filter(Measurement == "SN" & Experiment == temp.exp)%>%
    na.omit()
  ### Add a column with the normalized value from the data frame that we created in the loop above
  tempSN$NormalizedExp <- na.omit(SN_normalize[,i])
  SN.temp <- rbind(SN.temp, tempSN)
}

### Once again verify that this is doing what I want it to
### The order of the ID's for the quantile normalized values and the measured values should be the same
NormalizedID <- SN.temp%>%
  group_by(Experiment)%>%
  arrange(Experiment, NormalizedExp)%>%
  select(ID)
MeasuredValue_ID <- SN.temp%>%
  group_by(Experiment)%>%
  arrange(Experiment, Measured_Value)%>%
  select(ID)
### This is true for 100% of the IDs
sum(NormalizedID[,2] == MeasuredValue_ID[,2]) / nrow(SN.temp)


TSC.temp <- data.frame()
### For each experiment
for (i in 1:3){
  temp.exp <- c("DEPI1", "DEPI2", "DEPI3")[i]
  ### Filter to the ith flat and the DEPI1 experiment
  ###Remove rows that have NA (remember, these are the NA values appended to the end of the rows in order to combine flats that have a different number of plants)
  tempTSC <- data%>%
    filter(Measurement == "TSC" & Experiment == temp.exp)%>%
    na.omit()
  ### Add a column with the normalized value from the data frame that we created in the loop above
  tempTSC$NormalizedExp <- na.omit(TSC_normalize[,i])
  TSC.temp <- rbind(TSC.temp, tempTSC)
}

### Plot the data to see the distribution of the measured values before and after normalization by experiment
plot1 <- ggplot(data = TSC.temp, aes(x = NormalizedMeasuredValue , color = Experiment))+
  geom_density()+
  stat_density(geom = "line", position = "identity")+
  theme_linedraw()+
  labs(x = "Quantile Normalized Values by Flat",
       y = "Density",
       title = "TSC Normalized Values by Flat")
plot2 <- ggplot(data = TSC.temp, aes(x = NormalizedExp, color = Experiment))+
  geom_density()+
  stat_density(geom = "line", position = "identity")+
  theme_linedraw()+
  labs(x = "Quantile Normalized Values by Flat then Experiment",
       y = "Density",
       title = "TSC Normalized Values by Flat then Experiment")
grid.arrange(plot1, plot2)


SPF.temp <- data.frame()
### For each experiment
for (i in 1:3){
  temp.exp <- c("DEPI1", "DEPI2", "DEPI3")[i]
  ### Filter to the ith flat and the DEPI1 experiment
  ###Remove rows that have NA (remember, these are the NA values appended to the end of the rows in order to combine flats that have a different number of plants)
  tempSPF <- data%>%
    filter(Measurement == "SPF" & Experiment == temp.exp)%>%
    na.omit()
  ### Add a column with the normalized value from the data frame that we created in the loop above
  tempSPF$NormalizedExp <- na.omit(SPF_normalize[,i])
  SPF.temp <- rbind(SPF.temp, tempSPF)
}

### Plot the data to see the distribution of the measured values before and after normalization by experiment
plot1 <- ggplot(data = SPF.temp, aes(x = NormalizedMeasuredValue , color = Experiment))+
  geom_density()+
  stat_density(geom = "line", position = "identity")+
  theme_linedraw()+
  labs(x = "Quantile Normalized Values by Flat",
       y = "Density",
       title = "SPF Normalized Values by Flat")
plot2 <- ggplot(data = SPF.temp, aes(x = NormalizedExp, color = Experiment))+
  geom_density()+
  stat_density(geom = "line", position = "identity")+
  theme_linedraw()+
  labs(x = "Quantile Normalized Values by Flat then Experiment",
       y = "Density",
       title = "SPF Normalized Values by Flat then Experiment")
grid.arrange(plot1, plot2)

###Combine all of these data frames into one with all normalized values
NormalizedData_byExp <- rbind(SPF.temp,
                        TSC.temp,
                        SN.temp
)


### SAVE THESE NORMALIZED VALUES TO USE IN LATER SCRIPTS. REMEMBER, THIS IS WITH REMOVING THE OUTLIERS
write.csv(NormalizedData_byExp, file = "Fitness_NormalizedValuesByExp_OutliersRemoved_07232021.csv")
