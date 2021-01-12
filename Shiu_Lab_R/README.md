# DEPI

# Visualizations 

This folder contains the figures generated through the December and February analysis contained in the DEPI_Analysis_Dec and DEPI_Analysis_Feb files. These figures are:

* Heat maps of log-two-fold-change of measured values
* Heat maps of both corrected and uncorrected p-values
* Line plots of the measured values


# Fitness

This folder contains:

## fitness_analysis_loop_091319

Analysis of the fitness data, created by Siobhan, using linear models and anova. 
## Fitness_Pipeline

Analysis of the same fitness data. But, this first removes outliers, quantile normalizes the data, and calculates epistasis. Then, the epistasis calculations are plotted in heat maps.

## Fitness_Figures

This folder has heat maps of the additive and proportional epistasis for SN, TSC, and SPF, as well as the selection coefficients.

This folder also has figures of quantile normalization to confirm that the distributions of each flat are the same among each experiment.

# Experiment_Information

## DEPI_Experiment_Information

This file contains a detailed description of the initial data cleaning process. Though, this file does not contain outlier removal or quantile normalization.

This file contains information on the number of genotypes per experiment, the experiment lengths, as well as the number of measured values per measurement per day of each experiment. 

Finally, this file has notes on why the "X"'s were removed from some time points.

## plantsize_p_investigations

This file looks at the net daily growth for each plant, as well as other ways to analyze leaf area. This investigation was because some plants "shrink" throughout the day. 

# DEPI_R_Scripts

## B1B3 Scripts

This is a folder that has 3 scripts - December, January, and February - for the B1B3 gene family. 

These scripts create the log 2 fold change heat maps, line plots, and p-value heat maps, as well as a description of the data cleaning.

This does not include quantile normalization, but does include outlier removal. 

## Ftsz Scripts

This is a folder that has 3 scripts - December, January, and February - for the Ftsz gene family. 

These scripts create the log 2 fold change heat maps, line plots, and p-value heat maps, as well as a description of the data cleaning.

This does not include quantile normalization, but does include outlier removal.

## MPK Scripts

This is a folder that has 3 scripts - December, January, and February - for the map kinase gene family. 

These scripts create the log 2 fold change heat maps, line plots, and p-value heat maps, as well as a description of the data cleaning.

This does not include quantile normalization, but does include outlier removal.

## DEPI_Analysis_All

This is a general comparison of each of the experiments. This script creates violin plots to compare the distribution of measured values between experiments, and also has a comparison of the number of measured values per measurement per day. 

## DEPI_LinePlots_Cleaned_Data

There are three versions of this file: .Rmd, .pdf, and .html.

This file uses the same code from the MPK Scripts folder above, but the data has been quantile normalized by flat, in addition to having the outliers removed. 


# DEPI_Figures_For_Rmd

There are three files in this folder. All are schematics that are included in various .Rmd files to explain a concept. Two are screenshots of Shinhan's notes, and one is a screenshot of a slide I created to present in lab meeting.

These will be embedded in .Rmd files.

# DEPI_Figures

Note that all the figures in this folder are before quantile normalization.

## Dec B1B3 Figures

P-value heat maps, log two fold change heat maps, and line plots for the B1B3 gene family from the December experiment.

## Dec Ftsz Figures

P-value heat maps, log two fold change heat maps, and line plots for the Ftsz gene family.

## December - heat maps, all

The log two fold change heat maps that include all genotypes. There are three plots - one for each measurement. 

## December - individual heat maps

This folder contains three folders, one for each measurement. In each folder are the log two fold change heat maps for each combination of single and double mutants.

## December - line graph

Line plots for each combination of single and double mutants. These line plots are faceted to include all three measurements on each plot. 

## December - p-value

P-value heat maps for each measurement, both corrected and uncorrected.

## Feb B1B3 Figures

P-value heat maps, log two fold change heat maps, and line plots for the B1B3 gene family from the February experiment.

## Feb Ftsz Figures

P-value heat maps, log two fold change heat maps, and line plots for the Ftsz gene family from the February experiment.

## February - heat maps, all

The log two fold change heat maps that include all genotypes. There are three plots - one for each measurement from the February experiment.

## February - individual heat maps

This folder contains three folders, one for each measurement. In each folder are the log two fold change heat maps for each combination of single and double mutants from the February experiment.

## February - line graph

Line plots for each combination of single and double mutants. These line plots are faceted to include all three measurements on each plot from the February experiment. 

## February - p-value

P-value heat maps for each measurement, both corrected and uncorrected from the February experiment.

## Jan B1B3 Figures

P-value heat maps, log two fold change heat maps, and line plots for the B1B3 gene family from the January experiment.

## Other

Density plots of the npq values for December and February which shows the amount of npq measured value that are negative for each experiment.

## VP

Violin plots comparing the December and February experiment distribution of each measured value for a range of genotypes compared to wild type. Note that these plots are created for each day, which could lead to problems for the flucuating days.

Also includes subfolders, and plots not in folders, all of which contain further violin plots.













