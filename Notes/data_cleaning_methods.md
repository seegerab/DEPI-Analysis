% Data Cleaning Steps for DEPI Experiment

# Outlier Removal

Outliers were removed from the data using the median absolute deviation. For the measured values associated with each genotype, time point, and measurement, measured values that exceeded 3.5 median absolute deviations were excluded as outliers. A constant scale factor of 1.4826 was used because underlying distribution of the measured values is assumed to be normal.

# Quantile Normalization

Once the outliers were removed, quantile normalization was performed. For each measurement the plants on the same flat of each experiment were normalized against each other. To do this, the measured values were ranked for each flat in each experiment. Then, the measured values were reordered in ascending order for each flat, and the mean of the...[will keep working on this January 26] 

