#dichotomize intervention covariates for main analysis
dichotomize <- function(df_code, cutoff_tot = NA, cutoff_peak = NA){
  #cutoff is the limit at which we will dichotomize duration (either total or peak)
  #if cutoff = NA (default) then use the median of reported values
  
  df_code <- df_code[,c("studyID", "study","age","Aarm","cont_act","ref","school",
                        "home","community","home_activity","delivery", "electronic", 
                        "diet","activity","multi_strat","single_phase","continuous",
                        "tot_duration", "peak_duration", "intensity", "integration", 
                        "flexibility","choice", "fun_factor", "resonance", "behaviour",
                        "education","social","environment","commercial")]
  
  #Already binary---------------------------------------------------------------
  #SETTING
  df_code$school <- ifelse(df_code$school=="Y",1,0)
  df_code$home <- ifelse(df_code$home =="Y",1,0)
  df_code$community <- ifelse(df_code$community=="Y",1,0)
  df_code$home_activity <- ifelse(df_code$home_activity=="Y",1,0)
  
  #REALM TARGETED
  #Code diet/act such that: (0,0)=diet only, (1,0)=act only, (0,1)=diet&act
  diet <- ifelse(df_code$diet=="No" | df_code$diet=="Yes minimally",0,1)
  act <- ifelse(df_code$activity=="No" | df_code$activity=="Yes minimally",0,1)
  
  #for checking collinearity between original covariates, instead use:
  #df_code$diet <- ifelse(df_code$diet=="No" | df_code$diet=="Yes minimally",0,1)
  #df_code$activity <- ifelse(df_code$activity=="No" | df_code$activity=="Yes minimally",0,1)
  
  
  df_code$activity <- ifelse(diet==0 & act==1, 1, 0)
  df_code$diet_activity <- ifelse(diet==1 & act==1, 1, 0)
  
  #MULTI-FACTORNESS
  df_code$multi_strat <- ifelse(df_code$multi_strat=="Y",1,0)
  df_code$single_phase <- ifelse(df_code$single_phase =="Y",1,0)
  df_code$continuous <- ifelse(df_code$continuous=="Y",1,0)
  df_code$intensity <- ifelse(df_code$intensity=="High",1,0)
  
  #FLEX/CHOICE/RES
  df_code$flexibility <- ifelse(df_code$flexibility=="Y",1,0)
  df_code$choice <- ifelse(df_code$choice=="Y",1,0)
  df_code$resonance <- ifelse(df_code$resonance=="Y",1,0)
  
  #MECHANISM
  df_code$behaviour <- ifelse(df_code$behaviour=="Y",1,0)
  df_code$education <- ifelse(df_code$education=="Y",1,0)
  df_code$social <- ifelse(df_code$social=="Y",1,0)
  df_code$environment <- ifelse(df_code$environment=="Y",1,0)
  
  #COMMERCIAL
  df_code$commercial <- ifelse(df_code$commercial=="Y",1,0)
  df_code$commercial <- ifelse(is.na(df_code$commercial),0, df_code$commercial)
  
  #Categorical------------------------------------------------------------------
  #DELIVERY
  df_code$group <- ifelse(df_code$delivery=="Exclusively or mainly as a group" | df_code$delivery=="Both individually and as a group", 1, 0)
  df_code$indiv <- ifelse(df_code$delivery=="Exclusively or mainly individually" | df_code$delivery=="Both individually and as a group", 1, 0)
  
  #ELECTRONIC
  df_code$electronic <- ifelse(df_code$electronic == "No", 0, 1)
  
  #Integration:
  df_code$integration <- ifelse(df_code$integration == "Y",1,0)
  
  #FUN FACTOR
  df_code$fun_factor <- ifelse(df_code$fun_factor == "Fun", 1, 0)
  
  #Continuous-------------------------------------------------------------------
  if(is.na(cutoff_tot)){
    cutoff_tot <- median(df_code$tot_duration)
    print(paste("median of total duration = ", cutoff_tot))
  }
  if(is.na(cutoff_peak)){
    cutoff_peak <- median(df_code$peak_duration)
    print(paste("median of peak duration = ", cutoff_peak))
  }
  df_code$tot_duration <- ifelse(df_code$tot_duration >= cutoff_tot, 1, 0)
  df_code$peak_duration <- ifelse(df_code$peak_duration >= cutoff_peak, 1, 0)
  
  
  #print numbers for the splits - categorical and continuous--------------------
  print(paste("Electronic: Y (signif, major, minor) = ", sum(df_code$electronic), ", N = ", nrow(df_code)-sum(df_code$electronic)))
  print(paste("Integration: Y = ", sum(df_code$integration), ", N (&P) = ", nrow(df_code)-sum(df_code$integration)))
  print(paste("Fun factor: Fun = ", sum(df_code$fun_factor), ", Neutral & Not fun = ", nrow(df_code)-sum(df_code$fun_factor)))
  print(paste("Tot Duration: Long = ", sum(df_code$tot_duration), ", short = ", nrow(df_code)-sum(df_code$tot_duration)))
  print(paste("Peak Duration: Long = ", sum(df_code$peak_duration), ", short = ", nrow(df_code)-sum(df_code$peak_duration)))
  
  diet_act <- sum(df_code$diet_activity)
  act_only <- sum(df_code$activity)
  diet_only <- nrow(df_code)-diet_act-act_only
  print(paste("Diet only: Yes = ", diet_only, ", No = ", nrow(df_code)-diet_only))
  print(paste("Activity only: Yes = ", act_only, ", No= ", nrow(df_code)-act_only))
  print(paste("Diet and Act: Yes = ", diet_act, ", No= ", nrow(df_code)-diet_act))
  
  print(paste("School: Y = ", sum(df_code$school), ", N = ", nrow(df_code)-sum(df_code$school)))
  
  df_code
  
}

#select/redefine intervention covariates for main analysis
restrict_covs <- function(df_code){
  #combine home and home_activity -> call home
  for(i in 1:nrow(df_code)){
    if(df_code$home[i]==1 || df_code$home_activity[i]==1){
      df_code$home[i] <- 1
    }else{
      df_code$home[i] <- 0
    }
  }
  
  #combine choice & flexibility -> call it flex_choice
  df_code$flex_choice <- NA
  for(i in 1:nrow(df_code)){
    if(df_code$flexibility[i]==1 || df_code$choice[i]==1){
      df_code$flex_choice[i] <- 1
    }else{
      df_code$flex_choice[i] <- 0
    }
  }
  
  #redefine group/individual as a binary covariate = indiv (=Y for indiv & both, = N for group only)
  #i.e. remove group covariate
  
  #remove group, flexibility, choice, home_activity, peak_duration, continuous, commerical & single phase
  df_code <- df_code[,c("studyID", "study","age","Aarm","cont_act","ref",
                        "school","home","community","indiv", "electronic","diet_activity","activity",
                        "multi_strat","tot_duration","intensity", "integration", "flex_choice",
                        "fun_factor", "resonance", "behaviour","education","social","environment")]
  df_code
  
}

#Read in and code study level covariates and create a df
read_study_covs <- function(studies, data_path, fn_511_covs, fn_1218_covs){
  #studies contains a list of studyIDs for which we have outcome data
  
  covs511 <- read_excel(paste(data_path, fn_511_covs, sep = ""))
  covs1218 <- read_excel(paste(data_path, fn_1218_covs, sep = ""))
  
  #combine age groups
  covs1218$age<-"12-18"
  covs511$age<-"5-11"
  covs <- rbind(covs1218, covs511)
  covs$studyID <- paste(covs$Study, covs$age)
  
  #remove any studies not in "studies"
  covs <- covs[covs$studyID %in% studies, ]
  
  #sort alphabetically by studyID (to match data)
  covs <- covs[order(covs$studyID),] 
  
  #from covs extract "Cov: Income status of country" & "Cov: Socio-economic status"
  covs <- covs[,c("studyID", "age", "Cov: Income status of country","Cov: Socio-economic status")]
  #rename columns
  cols <- c("studyID", "age", "income_country", "SES")
  colnames(covs) <- cols
  
  #dichotomize study level covariates
  for(i in 1:nrow(covs)){
    if(covs$age[i]=="12-18"){
      covs$age[i] <- 1
    }else{
      covs$age[i] <- 0
    }
    
    if(covs$income_country[i]=="High income country"){
      covs$income_country[i] <- 1
    }else{
      covs$income_country[i] <- 0
    }
    
    if(covs$SES[i]=="Mixed socio-economic status"){
      covs$SES[i] <- 1
    }else{
      covs$SES[i] <- 0
    }
  }
  
  #make covariates numeric
  covs$age <- as.numeric(covs$age)
  covs$income_country <- as.numeric(covs$income_country)
  covs$SES <- as.numeric(covs$SES)
  
  covs
}

#Read in ROB data, code and create a df
read_rob <- function(studies, data, data_path, fn_511_rob, fn_1218_rob){
  #studies contains a list of studyIDs for which we have outcome data
  
  rob511 <- read_excel(paste(data_path, fn_511_rob, sep = ""))
  rob1218 <- read_excel(paste(data_path, fn_1218_rob, sep = ""))
  
  #combine age groups
  rob1218$age<-"12-18"
  rob511$age<-"5-11"
  rob <- rbind(rob1218, rob511)
  rob$studyID <- paste(rob$Study, rob$age)
  
  #remove any studies not in "studies"
  rob <- rob[rob$studyID %in% studies, ]
  
  #only keep overall bias judgement
  rob <- subset(rob, Domain=="Overall bias")
  rob <- rob[,c("studyID", "Outcome", "Domain", "Judgement")]
  
  #split outcome into measure (BMI, zBMI, percentile) and time
  rob$measure <- NA
  rob$time <- NA
  for(i in 1:nrow(rob)){
    out_split <- strsplit(rob$Outcome[i], split = " ")
    rob$measure[i] <- out_split[[1]][1]
    rob$time[i] <- out_split[[1]][2]
  }
  
  #replace zBMI with BMI-z and medium with mid
  for(i in 1:nrow(rob)){
    if(rob$measure[i]=="zBMI"){
      rob$measure[i] <- "BMI-z"
    }
    if(rob$time[i]=="medium"){
      rob$time[i]<- "mid"
    }
  }
  
  meas <- data %>% select(studyID, measure) %>% unique(.)
  meas$basic_measure <- NA
  
  for(i in 1:nrow(meas)){
    if(meas$measure[i]=="BMI" || grepl("BMI-z from BMI", meas$measure[i], fixed=TRUE)){
      meas$basic_measure[i] <- "BMI"
    }else if(meas$measure[i]=="Percentile" || grepl("BMI-z from percentile", meas$measure[i], fixed=TRUE)){
      meas$basic_measure[i] <- "Percentile"
    }else if(meas$measure[i]=="BMI-z"||meas$measure[i]=="BMI-z from proportion"){
      meas$basic_measure[i] <- "BMI-z"
    }else{
      print(meas$studyID[i])
    }
  }
  
  #remove any time points not in data
  rob$uniq <- paste(rob$studyID, rob$time)
  data$uniq <- paste(data$studyID, data$time)
  rob <- rob[rob$uniq %in% data$uniq, ]
  
  #compare basic_measure with measure in rob - remove row if they don't match
  new_rob <- rob[0, ] #empty df
  for(i in 1:nrow(rob)){
    row <- rob[i, ]
    id <- which(meas$studyID==row$studyID)
    if(meas$basic_measure[id]==row$measure){
      #add row to new_rob
      new_rob <- rbind(new_rob, row)
    }
  }
  
  #sort to match data
  new_rob <- new_rob[order(new_rob$time, decreasing = TRUE),] #short, med, long
  new_rob <- new_rob[order(new_rob$studyID),] #by study
  
  #select important columns
  new_rob <- new_rob[, c("studyID", "time", "Judgement")]
  cols <- c("studyID", "time", "rob")
  colnames(new_rob) <- cols
  
  #dichotomize
  for(i in 1:nrow(new_rob)){
    if(new_rob$rob[i]=="High risk"){
      new_rob$rob[i] <- 1
    }else{
      new_rob$rob[i] <- 0
    }
  }
  
  new_rob$rob <- as.numeric(new_rob$rob)
  new_rob
}

#create df with time level covariates
create_time_covs <- function(df_merge, df_time){
  #add columns to df_time for FU time points
  df_time$mid <- NA
  df_time$long <- NA
  
  df_merge$uniq_time <- paste(df_merge$studyID, df_merge$time)
  study_time <- df_merge %>% select(uniq_time, time) %>% unique(.) #length = length(df_time)
  
  for(i in 1:nrow(df_time)){
    if(study_time$time[i]=="short"){
      df_time$mid[i] <- 0
      df_time$long[i] <- 0
    }else if(study_time$time[i]=="mid"){
      df_time$mid[i] <- 1
      df_time$long[i] <- 0
    }else if(study_time$time[i]=="long"){
      df_time$mid[i] <- 0
      df_time$long[i] <- 1
    }else{
      print(paste("time error ", study_time$uniq_time[ind]))
    }
  }
  
  df_time
}

#add interaction terms to the main data df
add_ints <- function(df, covs_J){
  
  #add columns for each interaction in covs_j
  for(j in 1:length(covs_J)){
    df$new <- NA
    names(df)[names(df) == "new"] <- covs_J[j]
  }
  
  for(j in 1:length(covs_J)){
    #find column of interaction
    ind12 <- which(colnames(df)==covs_J[j])
    
    #find column of each covariate in the interaction
    int_split <- strsplit(covs_J[j], split = " ")
    int1 <- int_split[[1]][1]
    int2 <- int_split[[1]][2]
    
    ind1 <- which(colnames(df)==int1)
    ind2 <- which(colnames(df)==int2)
    
    #add data to interaction column
    df[,ind12] <- df[,ind1]*df[,ind2]
  }
  
  df
}

#centre covariates about their mean values
centre_covs <- function(df_covs, non_covs){
  #non_covs are the names of the columns of non-covariates that are not to be centred
  for(i in 1:ncol(df_covs)){
    if(!(colnames(df_covs)[i] %in% non_covs)){
      df_covs[[i]] <- df_covs[[i]]-mean(df_covs[[i]])
    }
  }
  df_covs
}

#for active comparison studies, work out difference in covariates between arms
calc_active_comp <- function(df_covs, covs_J){
  
  #Add a column to label the reference arm
  df_covs$ref_arm <- NA
  
  #Create an indicator for study & time (these will be treated individually)
  #NB: must be in order of time then arm
  df_covs$study_time <- paste(df_covs$studyID, df_covs$time)
  
  # for active intervention studies, identify the reference arm, subtract the 
  # reference arm intervention covariates from the other arms, and delete the 
  # reference arm row
  x<-1
  while(x<=nrow(df_covs)){
    if(df_covs$cont_act[x]=="active"){
      j<-1
      while(df_covs$study_time[x+j]==df_covs$study_time[x]){
        j <- j+1
      }
      n_a <- j #number of arms
      
      #find the reference arm
      for(i in 0:(n_a-1)){
        if(df_covs$ref[x+i]=="Y"){
          ind <- i
        }
      }
      
      for(i in 0:(n_a-1)){
        if(i != ind){
          #set reference arm label = intervention of reference arm
          df_covs$ref_arm[x+i] <- df_covs$Aarm[x+ind]
          
          #take difference between other interventions and df_coding arm--------
          df_covs$school[x+i] <- df_covs$school[x+i] - df_covs$school[x+ind]
          df_covs$home[x+i] <- df_covs$home[x+i] - df_covs$home[x+ind]
          df_covs$community[x+i] <- df_covs$community[x+i] - df_covs$community[x+ind]
          df_covs$diet_activity[x+i] <- df_covs$diet_activity[x+i] - df_covs$diet_activity[x+ind]
          df_covs$activity[x+i] <- df_covs$activity[x+i] - df_covs$activity[x+ind]
          df_covs$multi_strat[x+i] <- df_covs$multi_strat[x+i] - df_covs$multi_strat[x+ind]
          df_covs$intensity[x+i] <- df_covs$intensity[x+i] - df_covs$intensity[x+ind]
          df_covs$flex_choice[x+i] <- df_covs$flex_choice[x+i] - df_covs$flex_choice[x+ind]
          df_covs$resonance[x+i] <- df_covs$resonance[x+i] - df_covs$resonance[x+ind]
          df_covs$behaviour[x+i] <- df_covs$behaviour[x+i] - df_covs$behaviour[x+ind]
          df_covs$education[x+i] <- df_covs$education[x+i] - df_covs$education[x+ind]
          df_covs$social[x+i] <- df_covs$social[x+i] - df_covs$social[x+ind]
          df_covs$environment[x+i] <- df_covs$environment[x+i] - df_covs$environment[x+ind]
          df_covs$indiv[x+i] <- df_covs$indiv[x+i] - df_covs$indiv[x+ind]
          df_covs$electronic[x+i] <- df_covs$electronic[x+i] - df_covs$electronic[x+ind]
          df_covs$integration[x+i] <- df_covs$integration[x+i] - df_covs$integration[x+ind]
          df_covs$fun_factor[x+i] <- df_covs$fun_factor[x+i] - df_covs$fun_factor[x+ind]
          df_covs$tot_duration[x+i] <- df_covs$tot_duration[x+i] - df_covs$tot_duration[x+ind]
          
          #INTERACTIONS-------------------------------------------------
          if(length(covs_J)>0){
            ints <- which(colnames(df_covs) %in% covs_J)
            for(j in 1:length(covs_J)){
              #ints[j] is the index of the column, x+i is the index of the row
              df_covs[[x+i, ints[j]]] <- df_covs[[x+i, ints[j]]] - df_covs[[x+ind, ints[j]]]
            }
          }
        }
      }
      #delete reference arm row
      df_covs <- df_covs[-c(x+ind),]
      #move on to next study/time point
      x <- x + (n_a-1) #because of order of data (time then arm) - if multi-FU, the next loop will be the next time point
    }else{
      df_covs$ref_arm[x] <- "Control"
      x <- x+1
    }
  }
  df_covs
}


