#include packages
library(tidyverse)
library(readxl)
library(writexl)
library(parallel)
library(PowerTOST) #for OwensT function
library("nleqslv") #for solving non-linear eqs
#Parallisation libraries: https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745
library(doParallel)
library(foreach)
library(parallel)  


#source functions
list.files("functions", full.names = TRUE) %>% map(source)
#define path to data and reference charts
data_path <- "~\\Data\\"
chart_path <- "~\\Charts\\"

#read in outcome data 
df_1218 <- read_excel(paste(data_path, "12-18_Data.xlsx", sep=""))
df_511 <- read_excel(paste(data_path, "5-11_Data.xlsx", sep="")) 
#combine age groups
df_1218$age<-"12-18"
df_511$age<-"5-11"
df <- rbind(df_1218, df_511)
#create StudyID label including age
df$studyID <- paste(df$study, df$age) 


# Set fixed variables
rho <- 0.9 #correlation coefficient for calculating CS
icc <- 0.02 #ICC 

# List of Study IDs for multi-arm studies where we do not use the variance of ref arm as covariance
# and the covariances we will use instead
not_cov <- c("Bonsergent 2013") #cross over trial
alt_cov <-c(0)


################################################################################
# Step 1 - map non-zBMI data onto zBMI
################################################################################

#work out average baseline mean/SD for each outcome/age group 
#(this is used to impute baseline + FU for CS only studies -> NB: only one study uses this)
base_BMI_511 <- av_baseline(df, outcome = "BMI", age_gp = "5-11")
base_BMI_1218 <- av_baseline(df, outcome = "BMI", age_gp = "12-18")
base_perc_511 <- av_baseline(df, outcome = "Percentile", age_gp = "5-11")
base_perc_1218 <- av_baseline(df, outcome = "Percentile", age_gp = "12-18")

#Create dataset to be mapped
df_tomap <- create_ds_map(df, icc, rho, base_BMI_511, base_BMI_1218, base_perc_511, base_perc_1218)

#subset by outcome
df_tomap_b <- subset(df_tomap, measure=="BMI")
df_tomap_p <- subset(df_tomap, measure=="Percentile")

#MAP============================================================================
cl <- parallel::makeCluster(detectCores(), outfile = "MappingLog.txt")
doParallel::registerDoParallel(cl)

## from BMI
# sampling method: BMI - lognormal, age - normal,
b.samp <- foreach::foreach(i = 1:nrow(df_tomap_b), .combine = rbind, .packages = c("readxl")) %dopar% {
  set.seed(i)
  map.bmi.sample(df_tomap_b[i,], bmi_dist = "lognormal", age_dist = "normal", chart_path = chart_path,
                 bmi_age_corr = FALSE, num_samps = 10000)
  }

## from PERCENTILE
# analytic method
p.cdf <- foreach::foreach(i = 1:nrow(df_tomap_p), .combine = rbind, .packages = c("PowerTOST", "nleqslv")) %dopar% {
  map.p.cdf(df_tomap_p[i,])}


parallel::stopCluster(cl)#======================================================

save.image(file = "MappedData.RData") #can read back in to reset data to this point
#load("MappedData.RData")


#Restructure mapped data to match original data----------------
p.cdf <- edit_dfmap(p.cdf)
p.samp <- edit_dfmap(p.samp)
b.samp <- edit_dfmap(b.samp)


################################################################################
# Step 2 - calculate outcome data
################################################################################

#Reported data:----------------------------------------------
#calculate contrast level data (MD and SE)
df <- calc_MD(df, rho, not_cov, alt_cov)
#adjust for clustering
df <- adj_cluster(df, icc)

#define number of arms and number of FU times 
df <- count_arm_FU(df)

#identify any missing data 
miss_studies <- check_missing(df) #all necessary data for covariances

#Mapped data:-----------------------------------------------
#calculate contrast level data (MD and SE)
p.cdf <- calc_MD(p.cdf, rho, not_cov, alt_cov)
p.samp <- calc_MD(p.samp, rho, not_cov, alt_cov)
b.samp <- calc_MD(b.samp, rho, not_cov, alt_cov)

#adjust for clustering
p.cdf <- adj_cluster(p.cdf, icc)
p.samp <- adj_cluster(p.samp, icc)
b.samp <- adj_cluster(b.samp, icc)

#define number of arms and number of FU times
p.cdf <- count_arm_FU(p.cdf)
p.samp <- count_arm_FU(p.samp)
b.samp <- count_arm_FU(b.samp)

#identify any missing data 
miss_studies.p.cdf <- check_missing(p.cdf)#all necessary data for covariances
miss_studies.p.samp <- check_missing(p.samp)#all necessary data for covariances
miss_studies.b.samp <- check_missing(b.samp)#all necessary data for covariances

rm(miss_studies, miss_studies.p.cdf, miss_studies.p.samp, miss_studies.b.samp)

#Subset data by outcome/mapping status
#Reported data:------------------------------------
#zBMI only
df_z <- subset(df, measure=="BMI-z")
#BMI only
df_b <- subset(df, measure=="BMI")
#percentile only
df_p <- subset(df, measure=="Percentile")

#Create full data set (all data on zBMI scale)---------------------
#reported zBMI and zBMI from proportion
df_main <- subset(df, measure=="BMI-z" | measure =="BMI-z from proportion")
#combine with mapped data
df_main_samp_cdf <- do.call("rbind", list(df_main, b.samp, p.cdf)) 

#create data sets for sensitvity analysis------------------------
#zBMI & mapped percenbtile only
df_z_p <- do.call("rbind", list(df_z, p.cdf))

save.image("FinalMappedData.RData")

#END