#include packages
library(tidyverse)
library(readxl)
library(writexl)
library(corpcor)

#Load data for analysis
load("FinalMappedData.RData")
#source functions
list.files("functions", full.names = TRUE) %>% map(source)

################################################################################
#step 3 - intervention coding 
################################################################################

#read in intervention coding
coding <- read_excel(paste(data_path, "InterventionLevelCovariates.xlsx", sep=""))

#create studyID including age
coding$studyID <- paste(coding$study, coding$age)

#Check data
#check all studies appear in both outcome & intervention data sets
mis_match <- check_studies(df,coding) 
#check all arms appear in both outcome & intervention data sets
mis_arms <- check_arms(df, coding) 
#check all intervention covariates are coded correctly
check_code <- check_covs(coding) 

#dichotomize intervention covariates
coding <- dichotomize(coding)


################################################################################
# step 4 - create data for JAGS
################################################################################

#max number of arms
A <- max(df$num_int)#5 -> A-1 = 4

#max number of time points
T <- max(df$num_times)#3

#Outcome data===================================================================

#Define which outcome data to analyse (change as necessary)
data <- df_main_samp_cdf #main analysis

#sensitivity analyses (select as appropriate)
#data <- df_z #zBMI only
#data <- df_b #BMI only
#data <- df_z_p #zBMI and mapped percentile only

#ORDER
data <- data[order(data$Aarm),] 
data <- data[order(data$time, decreasing = TRUE),] #short, med, long
data <- data[order(data$studyID),] #by study

studies <- unique(data$studyID)
N <- length(studies) 

#Covariate data=================================================================

##Read in
#intervention covariates
df_code <- coding[coding$studyID %in% data$studyID, ] 
df_code <- restrict_covs(df_code) #final covariate selection
#put coding in same order as data
df_code <- df_code[order(df_code$Aarm),] #by arm
df_code <- df_code[order(df_code$studyID),] #by study

#study level covariates
fn_511_covs <- "5-11_StudyLevelCovariates.xlsx"
fn_1218_covs <- "12-18_StudyLevelCovariates.xlsx"
df_study <- read_study_covs(studies, data_path, fn_511_covs, fn_1218_covs)

#time level covariates
#risk of bias
fn_511_rob <- "5-11_ROBData.xlsx"
fn_1218_rob <- "12-18_ROBData.xlsx"
df_time <- read_rob(studies, data, data_path, fn_511_rob, fn_1218_rob)
#add FU time covariates (medium/long) to ROB data
df_time <- create_time_covs(data, df_time)

#remove unnecessary columns
data <- subset(data, select = -c(study, age))
df_code <- subset(df_code, select = -c(study, age))

#merge data + df_code + df_time + df_study
df_covs <- merge(df_code, df_time, by = c("studyID"))
df_covs <- merge(df_covs, df_study, by = c("studyID"))
df_covs <- df_covs[order(df_covs$Aarm),]
df_covs <- df_covs[order(df_covs$time, decreasing = TRUE),] #short, med, long
df_covs <- df_covs[order(df_covs$studyID),] #by study

## Define the set of covariates to analyse

#Intervention-level
ind0 <- which(colnames(df_code) == "school")#first cov
ind1 <- which(colnames(df_code) == "environment")#last cov
covs_x <- colnames(df_code)[ind0:ind1] 
#Study-level
covs_z <- c("age", "income_country", "SES")
#FU time level
covs_w <- c("mid","long","rob")

## Interactions=================================================================

#Main analysis: selected interactions
covs_J <- c("electronic age", "diet_activity age", "multi_strat age",
            "integration age", "resonance age", "education age", "income_country age",
            "electronic diet_activity", "multi_strat diet_activity",
            "fun_factor diet_activity", "resonance diet_activity", "rob diet_activity",
            "electronic activity", "tot_duration activity","fun_factor activity", 
            "income_country activity","rob activity")

## SSVS models:
## Step 1: no interactions
# covs_J <- c() 

## Step 2: all interactions with age
# covs_J <- c(covs_x, covs_z, covs_w) #all covs
# covs_J <- covs_J[ !covs_J== "age" ] #remove age
# covs_J <- paste(covs_J, "age") #paste age after every term

## Step 3: all interactions with DA & activity
# # With DA:
# covs_J_DA <- c(covs_x, covs_z, covs_w) #all covs
# covs_J_DA <- covs_J_DA[ !covs_J_DA== "diet_activity" ] #remove DA
# covs_J_DA <- covs_J_DA[ !covs_J_DA== "activity" ] #remove activity
# covs_J_DA <- paste(covs_J_DA, "diet_activity") #paste DA after every term
# # With activity:
# covs_J_act <- c(covs_x, covs_z, covs_w) #all covs
# covs_J_act <- covs_J_act[ !covs_J_act== "activity" ] #remove activity
# covs_J_act <- covs_J_act[ !covs_J_act== "diet_activity" ] #remove DA
# covs_J_act <- paste(covs_J_act, "activity") #paste activity after every term
# # Combine:
# covs_J <- c(covs_J_DA, covs_J_act)

df_covs <- add_ints(df_covs, covs_J)

## Format data for analysis=====================================================

#Centre covariates
non_covs <- c("studyID", "Aarm", "cont_act", "ref", "time")
df_covs <- centre_covs(df_covs, non_covs)

#Calculate difference in covariates between arms for active comparison studies
df_covs <- calc_active_comp(df_covs, covs_J)

#merge covariate and outcome data
data <- merge(data, df_covs, by = c("studyID", "Aarm", "time"))

#remove unnecessary columns
data <- subset(data, select = -c(label, study_time))

#make sure its in the correct order
data<-data[order(data$Aarm),]
data <- data[order(data$time, decreasing = TRUE),] #short, med, long
data <- data[order(data$studyID),] #by study

## Create arrays to read into JAGS==============================================

#Study and FU time indicators---------------------------------------------------
#vectors of length N with number of no. of contrasts and no. of time points per study
numcont <- data %>% select(studyID, num_int) %>% unique(.)
numcont <- numcont$num_int-1

numFU <- data %>% select(studyID, num_times) %>% unique(.)
numFU <- numFU$num_times

#observations-------------------------------------------------------------------
y <- data$MD #univariate array of MDs

#dimension indicators-----------------------------------------------------------
dim <- numFU*numcont
#array containing the index where each study begins in y
yind <- create_yind(N, numFU, numcont)
yend <- yind + dim - 1
tind <- array(dim = c(N, T, A-1))
for(i in 1:N){
  for(t in 1:numFU[i]){
    for(k in 1:numcont[i]){
      #index for theta_{i,t}^{(k)}
      tind[i,t,k] <- yind[i] + (t-1)*numcont[i] + (k-1)
    }
  }
}

#within study covariance matrix-------------------------------------------------
## Correlations:
#define rho_tt between MDs
rhoy_tt <- 0.8 # main analysis
#rhoy_tt <- 0.5 # sensitivity analysis
#rhoy_tt <- 0.95 # sensitivity analysis
rhoy_sm <- rhoy_ml <- rhoy_tt #one degree of separation
rhoy_sl <- rhoy_tt**2 #2 degrees of separation

#define rho_tt between CSs
rhod_tt <- 0.8 # main analysis
#rhod_tt <- 0.5 # sensitivity analysis
#rhod_tt <- 1.0 # sensitivity analysis
rhod_sm <- rhod_ml <- rhod_tt #one degree of separation
rhod_sl <- rhod_tt**2 #2 degrees of separation

#define rho_tt between AFs (assumed to be the same as CS)
rhoAF_tt <- 0.8 # main analysis
#rhoAF_tt <- 0.5 # sensitivity analysis
#rhoAF_tt <- 1.0 # sensitivity analysis
rhoAF_sm <- rhoAF_ml <- rhoAF_tt #one degree of separation
rhoAF_sl <- rhoAF_tt**2 #2 degrees of separation

## Create covariance matrix:
V <- create_cov(data, yind, numcont, numFU, N, 
                rhoy_sm, rhoy_sl, rhoy_ml,
                rhod_sm, rhod_sl, rhod_ml, 
                rhoAF_sm, rhoAF_sl, rhoAF_ml)


##Covariates--------------------------------------------------------------------
#intervention-level
n <- length(covs_x)
x <- create_x(data, yind, numcont, N, A, covs_x)

#study-level  
p <- length(covs_z)
z <- create_z(data, yind, N, covs_z)

#FU time 
q <- length(covs_w) #q is no. of FU time covs (use q not q-1)
w <- create_w(data, yind, numFU, numcont, N, T, covs_w)

#interactions
l <- length(covs_J)
J <- create_J(data, yind, numFU, numcont, N, T, A, covs_J) 

#Sigma (RE cov matrix)----------------------------------------------------------
S <- create_S(numFU, numcont, N)

#control/active indicator-------------------------------------------------------
#indicator to say whether the reference arm is a control or an active intervention
ref_arm <- data %>% select(studyID, cont_act) %>% unique(.)
ref_arm <- ref_arm$cont_act
for(i in 1:N){
  if(ref_arm[i]=="control"){
    ref_arm[i] <- 1
  }else if(ref_arm[i]=="active"){
    ref_arm[i] <- 0
  }
}
ref_arm <- as.numeric(ref_arm)


##check if covariance matrices are positive definite----------------------------
for(i in 1:N){
  if(dim[i]>1){
    
    check <- is.positive.definite(V[i, 1:dim[i], 1:dim[i]])
    if(check==FALSE){
      print(paste(studies[i], "FU = ", numFU[i], "constrasts = ", numcont[i]))
      print(V[i, 1:dim[i], 1:dim[i]])
    }
  }
}



## SAVE workspace===============================================================

#Data for input into JAGS
#Data: y, V, x, z, w, J, S
#indicators: yind, numFU, numcont, ref_arm
#numbers: N, n, p, q, l

#delete unnecessary data
remove(check_code, coding, df_1218, df_511, df, mis_match, mis_arms, 
       alt_cov, not_cov, icc, df_code, df_time, df_study, df_covs)
remove(rho, rhoAF_ml, rhoAF_sl, rhoAF_sm, rhoAF_tt, rhod_ml, rhoy_sl,
       rhoy_sm, rhoy_tt, rhod_sl, rhod_sm, rhod_tt, rhoy_ml)
remove(chart_path, data_path, fn_1218_covs, fn_511_covs,
       i, k, t, check, non_covs, fn_1218_rob, fn_511_rob)
#delete functions
remove(add_ints, adj_cluster, age.sample.norm, age.sample.unif, av_baseline, 
       b.sample.lnorm, b.sample.norm, bmi.from.z, calc_active_comp, calc_MD, 
       centre_covs, check_arms, check_covs, check_missing, check_studies, 
       count_arm_FU, create_cov, create_ds_map, create_J, create_S, 
       create_time_covs, create_w, create_x, create_yind, create_z,
       dichotomize, dichot_orig, edit_dfmap, init.sample,map.bmi.sample, 
       map.bmi.sample.corr, map.bmi.sample.nocorr, map.p.cdf, map.p.sample,
       read_rob,read_study_covs, restrict_covs, sim_eqs, standardise_basecont, 
       z.from.bmi, z.sample.norm)
#extras
remove(b.samp, cl, df_b, df_main, df_z, df_z_p, df_main_samp_cdf, df_p, df_tomap,
       df_tomap_b, df_tomap_p, p.cdf)
remove(base_BMI_1218, base_BMI_511, base_perc_1218, base_perc_511, ind0, ind1)


#Save data for reading in
save.image(file = "DataForAnalysis.RData")

