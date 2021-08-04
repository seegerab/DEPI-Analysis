#load required libraries 
require(lme4)
require(lmerTest)
require(lsmeans)
require(car)
#if doing plots
require(ggplot2)
library(stringr)

#define file name for text output (summary model results etc.)
output_filename = "DM_output_091319.txt"
#read in data file (includes border rows)
data = read.delim("MAPK_DEPI_data_inside_062619.txt", sep = "\t", header = T)

#make sure flat, row and column are coded as categorical variables
data$Flat<-as.factor(data$Flat)
data$Row<-as.factor(data$Row)
data$Column<-as.factor(data$Column)

#Get a list of the double mutants in this dataframe
all_double_mutants = list()
for (gen in levels(data$Genotype)) {
  if (str_detect(gen, "_") == T) {
    all_double_mutants = c(all_double_mutants, gen)
  }
}

#Loop through the list of double mutants in the dataframe and run stats on those alleles. 


for (dm in unlist(as.character(all_double_mutants))) {
  print(dm)
  #Get the names of the corresponding single mutants for this double mutant
  sm1 = strsplit(dm,'_')[[1]][1]
  sm2 = paste('mpk',strsplit(dm,'_')[[1]][2],sep='')
  sm1_colname = toupper(sm1)
  sm2_colname = toupper(sm2)
  
  #write single mutant names to output file
  sm_names_output = paste(c(as.name(sm1_colname), "_", as.name(sm2_colname)), collapse = '')
  write(sm_names_output, output_filename, append = T)
  
  #Get a list of experiments the double mutant is in
  dm_sub = subset(data, Genotype == dm)
  all_experiments = unlist(as.character(levels(dm_sub$Experiment)))
  
  #Generate a dataframe with the double mutant, corresponding single mutants, and WT ONLY for experiments this 
  #double mutant is in
  sm1_df = subset(data, Genotype == sm1)
  sm2_df = subset(data, Genotype == sm2)
  wt_df = subset(data, Genotype == "Col")
  subset_df = dm_sub
  expt_num = length(all_experiments)
  for (expt in 1:expt_num) {
    experiment = all_experiments[expt]
    subset_df = rbind(subset_df, subset(sm1_df, Experiment == experiment))
    subset_df = rbind(subset_df, subset(sm2_df, Experiment == experiment))
    subset_df = rbind(subset_df, subset(wt_df, Experiment == experiment))
  }
  
  #Rename the relevant gene column names so they are consistent for each loop
  names(subset_df)[names(subset_df) == sm1_colname] = "Gene1"
  names(subset_df)[names(subset_df) == sm2_colname] = "Gene2"
  
  #Model using genotype at each gene  to test for interactions between genes
  #Analyze total seed count (TSC)
  TSC.lmer<-lmer(TSC ~ Gene1 + Gene2 + Gene1:Gene2 + Experiment + Gene1:Experiment + Gene2:Experiment + Gene1:Gene2:Experiment + (1|Experiment:Flat), data=subset_df)
  
  #write lmer summary(Variance and Std.Dev for "Random effects" and no. of observations) to output file
  write("linear mixed-effects model analysis TSC", output_filename, append = T)
  number_observations = as.character(((summary(TSC.lmer))$devcomp[2])$dims[1])
  obs_output = paste(c("Number of observations: ", number_observations), collapse = '')
  write(obs_output, output_filename, append = T)
  rand_effects = as.data.frame(summary(TSC.lmer)$varcor)[,c("grp", "var1", "vcov", "sdcor")]
  colnames(rand_effects) = c("group", "Name", "Variance", "Std.Dev.")
  write("Random effects",output_filename, append = T)
  write.table(rand_effects, output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #Write ANOVA table to output file
  write("\nANOVA results TSC", output_filename, append = T)
  write.table(as.data.frame(anova(TSC.lmer)), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #print out lsmeans and contrasts
  TSC.lsm<- lsmeans(TSC.lmer, ~Gene1:Gene2)
  write("\nlsmeans TSC", output_filename, append = T)
  write.table(as.data.frame(summary(TSC.lsm)), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  write("\nlsmeans contrasts TSC", output_filename, append = T)
  write.table(as.data.frame(summary(contrast(TSC.lsm, "pairwise"))), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #model TSC lsmeans per genotype
  TSCg.lmer<-lmer(TSC ~ Genotype + Experiment + Genotype:Experiment + (1|Genotype:Subline) + (1|Experiment:Flat), data=subset_df)
  
  
  #analyze seeds per fruit (squared to improve normality)
  subset_df$sSPF = as.numeric(as.character(subset_df$sSPF))
  sSPF.lmer<-lmer(sSPF ~ Gene1 + Gene2 + Gene1:Gene2 + Experiment + Gene1:Experiment + Gene2:Experiment + Gene1:Gene2:Experiment + (1|Experiment:Flat), data=subset_df)
  
  #write lmer summary(Variance and Std.Dev for "Random effects" and no. of observations) to output file
  write("\nlinear mixed-effects model analysis sSPF", output_filename, append = T)
  number_observations = as.character(((summary(sSPF.lmer))$devcomp[2])$dims[1])
  obs_output = paste(c("Number of observations: ", number_observations), collapse = '')
  write(obs_output, output_filename, append = T)
  rand_effects = as.data.frame(summary(sSPF.lmer)$varcor)[,c("grp", "var1", "vcov", "sdcor")]
  colnames(rand_effects) = c("group", "Name", "Variance", "Std.Dev.")
  write("Random effects",output_filename, append = T)
  write.table(rand_effects, output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #Write ANOVA table to output file
  write("\nANOVA results sSPF", output_filename, append = T)
  write.table(as.data.frame(anova(sSPF.lmer)), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #print out lsmeans and contrasts
  sSPF.lsm<- lsmeans(sSPF.lmer, ~Gene1:Gene2)
  write("\nlsmeans sSPF", output_filename, append = T)
  write.table(as.data.frame(summary(sSPF.lsm)), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  write("\ncontrasts sSPF", output_filename, append = T)
  write.table(as.data.frame(summary(contrast(sSPF.lsm, "pairwise"))), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #model sqSPF lsmeans per genotype
  sSPFg.lmer<-lmer(sSPF ~ Genotype + Experiment + Genotype:Experiment + (1|Genotype:Subline) + (1|Experiment:Flat), data=subset_df)
  
  #would probably be better to report the backtransformed values
  #square each value: lsmean, lower.CL and upper.CL
  
  #fruit number
  SN.lmer<-lmer(SN ~ Gene1 + Gene2 + Gene1:Gene2 + Experiment + Gene1:Experiment + Gene2:Experiment + Gene1:Gene2:Experiment + (1|Experiment:Flat), data=subset_df)
  
  #from summary, need Variance and Std.Dev for "Random effects" and no. of observations
  write("\nlinear mixed-effects model analysis SN", output_filename, append = T)
  number_observations = as.character(((summary(SN.lmer))$devcomp[2])$dims[1])
  obs_output = paste(c("Number of observations: ", number_observations), collapse = '')
  write(obs_output, output_filename, append = T)
  rand_effects = as.data.frame(summary(SN.lmer)$varcor)[,c("grp", "var1", "vcov", "sdcor")]
  colnames(rand_effects) = c("group", "Name", "Variance", "Std.Dev.")
  write("Random effects",output_filename, append = T)
  write.table(rand_effects, output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #print out ANOVA table
  write("\nANOVA results SN", output_filename, append = T)
  write.table(as.data.frame(anova(SN.lmer)), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  
  #print out lsmeans and contrasts
  SN.lsm<- lsmeans(SN.lmer, ~Gene1:Gene2)
  write("\nlsmeans SN", output_filename, append = T)
  write.table(as.data.frame(summary(SN.lsm)), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  write("\ncontrasts SN", output_filename, append = T)
  write.table(as.data.frame(summary(contrast(SN.lsm, "pairwise"))), output_filename, append = T, sep = "\t", row.names = F, col.names = T, quote = F)
  write("\n\n", output_filename, append = T)
  
  #model SN lsmeans per genotype
  SNg.lmer<-lmer(SN ~ Genotype + Experiment + Genotype:Experiment + (1|Genotype:Subline) + (1|Experiment:Flat), data=subset_df)
  
  #save plots to assess model fit
  pdf(paste(c(sm1_colname, "_", sm2_colname, ".pdf"), collapse = ''),width=14, height=14,pointsize=6)
  #TSC qqplot
  print(qqp(resid(TSC.lmer), main = paste(c("TSC ", sm1_colname, "_", sm2_colname), collapse = '')))
  #TSC residuals, histogram 
  print(hist(resid(TSC.lmer), main = paste(c("TSC ", sm1_colname, "_", sm2_colname), collapse = '')))
  #TSC residuals, dotplot
  print(plot(TSC.lmer, which=1, main = paste(c("TSC ", sm1_colname, "_", sm2_colname), collapse = '')))
  #plot TSC lsmeans for each genotype
  print(plot(lsmeans(TSCg.lmer,"Genotype"), main = paste(c("TSC ", sm1_colname, "_", sm2_colname), collapse = '')))
  #sSPF qqplot
  print(qqp(resid(sSPF.lmer),main = paste(c("sSPF ", sm1_colname, "_", sm2_colname), collapse = '')))
  #sSPF residuals, histogram
  print(hist(resid(sSPF.lmer),main = paste(c("sSPF ", sm1_colname, "_", sm2_colname), collapse = '')))
  #sSPF residuals, dotplot
  print(plot(sSPF.lmer, which=1,main = paste(c("sSPF ", sm1_colname, "_", sm2_colname), collapse = '')))
  #plot sSPF lsmeans for each genotype
  print(plot(lsmeans(sSPFg.lmer,"Genotype"),main = paste(c("sSPF ", sm1_colname, "_", sm2_colname), collapse = '')))
  #SN qqplot
  print(qqp(resid(SN.lmer),main = paste(c("SN ", sm1_colname, "_", sm2_colname), collapse = '')))
  #SN residuals, histogram
  print(hist(resid(SN.lmer),main = paste(c("SN ", sm1_colname, "_", sm2_colname), collapse = '')))
  #SN residuals, dotplot
  print(plot(SN.lmer, which=1,main = paste(c("SN ", sm1_colname, "_", sm2_colname), collapse = '')))
  #plot SN lsmeans per genotype
  print(plot(lsmeans(SNg.lmer,"Genotype"),main = paste(c("SN ", sm1_colname, "_", sm2_colname), collapse = '')))
  dev.off()
  
}
