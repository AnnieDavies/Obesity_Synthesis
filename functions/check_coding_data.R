#check outcome data and intervention data contain the same study names 
check_studies <- function(df_out, df_code){
  
  #check all studies in coding are in outcome
  miss_code <- character()
  for(i in 1:nrow(df_code)){
    if(!(df_code$studyID[i] %in% df_out$studyID)){
      print(paste(df_code$studyID[i], " in coding but not in outcome"))
      if(!(df_code$studyID[i] %in% miss_code)){
        miss_code <- append(miss_code, df_code$studyID[i])
      }
    }
  }
  
  #check all studies in outcome are in coding 
  miss_out <- character()
  for(i in 1:nrow(df_out)){
    if(!(df_out$studyID[i] %in% df_code$studyID)){
      print(paste(df_out$studyID[i], " in outcome but not in coding"))
      if(!(df_out$studyID[i] %in% miss_out)){
        miss_out <- append(miss_out, df_out$studyID[i])
      }
    }
  }
  
  output <- list(miss_code, miss_out)
  output
}

#check outcome data and intervention data contain the same arms
check_arms <- function(df_out, df_code){
  
  #check all arms in coding are in outcome
  df_out$study_arm <- paste(df_out$studyID, df_out$Aarm)
  df_out$study_ref <- paste(df_out$studyID, df_out$Barm)
  df_code$study_arm <- paste(df_code$studyID, df_code$Aarm)
  
  miss_code <- character()
  for(i in 1:nrow(df_code)){
    if(!(df_code$study_arm[i] %in% df_out$study_arm)){
      if(!(is.na(df_code$ref[i])) & df_code$ref[i]=="Y"){
        if(!(df_code$study_arm[i] %in% df_out$study_ref)){
          print(paste(df_code$study_arm[i], " ref arm of coding not in Barm of outcome"))
          if(!(df_code$study_arm[i] %in% miss_code)){
            miss_code <- append(miss_code, df_code$study_arm[i])
          }
        }
      }else{
        print(paste(df_code$study_arm[i], " int arm of coding not in Aarm of outcome"))
        if(!(df_code$study_arm[i] %in% miss_code)){
          miss_code <- append(miss_code, df_code$study_arm[i])
        }
      }
    }
  }
  miss_code
}

#check intervention covariates are coded correctly (i.e. take one of their possible values)
check_covs <- function(df_code){
  df_new <- df_code[,c("studyID","age","Aarm")]
  
  #SETTING
  df_new$school <- ifelse(df_code$school!="Y" & df_code$school!="N","ERROR",0)
  df_new$home <- ifelse(df_code$home!="Y" & df_code$home!="N","ERROR",0)
  df_new$community <- ifelse(df_code$community!="Y" & df_code$community!="N","ERROR",0)
  df_new$home_activity <- ifelse(df_code$home_activity!="Y" & df_code$home_activity!="N","ERROR",0)
  
  #REALM TARGETED
  df_new$check_diet <- ifelse(df_code$diet!="No" & df_code$diet!="Yes minimally" & df_code$diet != "Yes exclusively/substantially", "ERROR", 0)
  df_new$check_act <- ifelse(df_code$activity!="No" & df_code$activity!="Yes minimally" & df_code$activity != "Yes exclusively/substantially", "ERROR", 0)
  
  #MULTI-FACTORNESS
  df_new$multi_strat <- ifelse(df_code$multi_strat!="Y" & df_code$multi_strat!="N","ERROR",0)
  df_new$single_phase <- ifelse(df_code$single_phase!="Y" & df_code$single_phase!="N","ERROR",0)
  df_new$continuous <- ifelse(df_code$continuous!="Y" & df_code$continuous!="N","ERROR",0)
  df_new$intensity <- ifelse(df_code$intensity!="High" & df_code$intensity!="Low","ERROR",0)
  
  #FLEX/CHOICE/RES
  df_new$flexibility <- ifelse(df_code$flexibility!="Y" & df_code$flexibility!="N","ERROR",0)
  df_new$choice <- ifelse(df_code$choice!="Y" & df_code$choice!="N","ERROR",0)
  df_new$resonance <- ifelse(df_code$resonance!="Y" & df_code$resonance!="N","ERROR",0)
  
  #MECHANISM
  df_new$behaviour <- ifelse(df_code$behaviour!="Y" & df_code$behaviour!="N","ERROR",0)
  df_new$education <- ifelse(df_code$education!="Y" & df_code$education!="N","ERROR",0)
  df_new$social <- ifelse(df_code$social!="Y" & df_code$social!="N","ERROR",0)
  df_new$environment <- ifelse(df_code$environment!="Y" & df_code$environment!="N","ERROR",0)
  
  #COMMERCIAL
  df_new$commercial <- ifelse(df_code$commercial!="Y" & df_code$commercial!="N" & !(is.na(df_code$commercial)),"ERROR",0)

  #DELIVERY
  df_new$check_delivery <- ifelse(df_code$delivery!="Exclusively or mainly individually" & df_code$delivery!="Both individually and as a group" & df_code$delivery!="Exclusively or mainly as a group", "ERROR", 0)
  
  #ELECTRONIC
  df_new$check_elec <- ifelse(df_code$electronic != "Yes exclusively" & df_code$electronic != "Yes significantly" & df_code$electronic != "Yes as a minor component" & df_code$electronic != "No", "ERROR", 0)
  
  #INTEGRATION
  df_new$check_integrat <- ifelse(df_code$integration != "Y" & df_code$integration != "N" & df_code$integration != "P", "ERROR", 0)
  
  #FUN FACTOR
  df_new$check_ff <- ifelse(df_code$fun_factor != "Better than neutral" & df_code$fun_factor != "Fun" & df_code$fun_factor != "Worse than neutral" & df_code$fun_factor != "Neutral" & df_code$fun_factor != "Boring", "ERROR",0)
  
  df_new
}