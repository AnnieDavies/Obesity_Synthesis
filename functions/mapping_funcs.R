# create a dataframe of data to be mapped onto zBMI
create_ds_map <- function(df, icc, rho, base_BMI_511, base_BMI_1218, base_perc_511, base_perc_1218){
  #df is the original dataframe of outcome data (both age groups)
  #icc is the imputed intra-class correlation coefficient
  #rho is the imputed correlation coefficient between baseline and follow-up
  #base_BMI & base_perc contain the average mean [1] and SD [2] at baseline of 
  #each outcome for the different age groups to use when we have CS only.
  #NB: results in df starting with MAP are data to be mapped
  
  #subset by Map
  df <- subset(df, Map == "Y")
  
  #exclude mapped from proportion and contrast level only (MD)
  df <- subset(df, Map_outcome != "Proportion" & Map_result != "MD")
  
  for(i in 1:nrow(df)){
   
    #unadjust for clustering-----------------------------------------------------
    if(df$Map_ReqUnAdj[i]=="Y"){
      #check if ICC is reported in trial, otherwise use imputed icc
      if(!is.na(df$calcICC[i]) && df$calcICC[i]=="Y"){ 
        del <- df$ICC[i]
      }else{ # Else if Calculated ICC= N: del = icc
        del <- icc
      }
      
      #unadjust SDs
      #check nAc and nBc (no. of clusters per group) are defined
      if(is.na(df$nAc[i]) || is.na(df$nBc[i])){print(paste("n cluster not defined:", df$studyID[i]))}
      
      #mean cluster size-------------------------------------------
      nA1 <- df$nA1[i]
      nB1 <- df$nB1[i]
      nA0 <- df$nA0[i]
      nB0 <- df$nB0[i]
      nAc <- as.numeric(df$nAc[i])
      nBc <- as.numeric(df$nBc[i])
      mB <- nB1/nBc
      mA <- nA1/nAc
      mB0<- nB0/nBc
      mA0 <- nA0/nAc
      
      #un-adjust SDs---------------------------------------------------
      #sqrt(design effect) = SQRT(1+del(m-1))
      desA1 <- sqrt(1+del*(mA-1)) 
      desB1 <- sqrt(1+del*(mB-1))
      desA0 <- sqrt(1+del*(mA0-1)) 
      desB0 <- sqrt(1+del*(mB0-1))
      
      #SD' = SD*sqrt(design effect)
      #SD = SD'/sqrt(design effect)
      
      #unadjust CS
      if(df$Map_result[i]=="CS"||df$Map_result[i]=="Baseline and CS"){
        df$MAP_CSsdA[i] <- df$MAP_CSsdA[i]/desA1
        df$MAP_CSsdB[i] <- df$MAP_CSsdB[i]/desB1
      }
      #unadjust baseline
      if(df$Map_result[i]=="Baseline and FU" || df$Map_result[i]=="Baseline and CS"){
        df$MAP_BFsd0A[i] <- df$MAP_BFsd0A[i]/desA0
        df$MAP_BFsd0B[i] <- df$MAP_BFsd0B[i]/desB0
      }
      #unadjust FU
      if(df$Map_result[i]=="Baseline and FU"){
        df$MAP_BFsd1A[i] <- df$MAP_BFsd1A[i]/desA1
        df$MAP_BFsd1B[i] <- df$MAP_BFsd1B[i]/desB1
      }
      #unadjust adjusted FU
      if(df$Map_result[i]=="Adjusted FU"){
        df$MAP_AFsdA[i] <- df$MAP_AFsdA[i]/desA1
        df$MAP_AFsdB[i] <- df$MAP_AFsdB[i]/desB1
      }
    }#end unadjust for clustering
    
    #calculate baseline & FU from CS--------------------------------------------
    #baseline & CS
    if(df$Map_result[i]=="Baseline and CS"){
      #FU = baseline + CS
      df$MAP_BFmean1A[i] <- df$MAP_BFmean0A[i] + df$MAP_CSmeanA[i]
      df$MAP_BFmean1B[i] <- df$MAP_BFmean0B[i] + df$MAP_CSmeanB[i]
      
      #set FU SD = baseline SD
      df$MAP_BFsd1A[i] <- df$MAP_BFsd0A[i] 
      df$MAP_BFsd1B[i] <- df$MAP_BFsd0B[i] 
    }
    #CS only
    if(df$Map_result[i]=="CS"){
      #use average baseline
      if(df$Map_outcome[i]=="BMI"){
        if(df$age[i]=="5-11"){
          base <- base_BMI_511[1]
          base_SD <- base_BMI_511[2]
        }else if(df$age[i]=="12-18"){
          base <- base_BMI_1218[1]
          base_SD <- base_BMI_1218[2]
        }
      }else if(df$Map_outcome[i]=="Percentile"){
        if(df$age[i]=="5-11"){
          base <- base_perc_511[1]
          base_SD <- base_perc_511[2]
        }else if(df$age[i]=="12-18"){
          base <- base_perc_1218[1]
          base_SD <- base_perc_1218[2]
        }
      }
      #define baselines
      df$MAP_BFmean0A[i] <- base
      df$MAP_BFmean0B[i] <- base
      #FU = baseline + CS
      df$MAP_BFmean1A[i] <- base + df$MAP_CSmeanA[i]
      df$MAP_BFmean1B[i] <- base + df$MAP_CSmeanB[i]
      
      #set FU SD = baseline SD
      df$MAP_BFsd0A[i] <- base_SD
      df$MAP_BFsd0B[i] <- base_SD
      
      df$MAP_BFsd1A[i] <- base_SD
      df$MAP_BFsd1B[i] <- base_SD 
    }
    #move Adjusted FU results to FU position
    if(df$Map_result[i]=="Adjusted FU"){
      df$MAP_BFmean1A[i] <- df$MAP_AFmeanA[i]
      df$MAP_BFsd1A[i] <- df$MAP_AFsdA[i]
      df$MAP_BFmean1B[i] <- df$MAP_AFmeanB[i]
      df$MAP_BFsd1B[i] <- df$MAP_AFsdB[i]
    }
  }
  
  #re-define result based on mapping result
  for(i in 1:nrow(df)){
    if(df$Map_result[i]=="Baseline and FU" || df$Map_result[i]=="Baseline and CS" || df$Map_result[i]=="CS"){
      df$result[i] <- "Baseline and FU"
    }else{
      df$result[i] <- df$Map_result[i]
    }
  }

  df <- df[, c("study", "age", "studyID","Map_outcome","Map_result", "result", "time", "Aarm",
               "Barm", "nA0", "nA1","nAc","nB0", "nB1", "nBc", "Map_ReqUnAdj", "Map_ReqReAdj",
               "MAP_BFmean0A", "MAP_BFsd0A", "MAP_BFmean0B", "MAP_BFsd0B",
               "MAP_BFmean1A", "MAP_BFsd1A", "MAP_BFmean1B", "MAP_BFsd1B",
               "Mean_Age", "SD_Age", "Prop_Male", "fu_months", "country",	"chart_ref", 
               "chart", "calcICC","ICC")]
  
  #rename variables
  cols <- c("study", "age", "studyID","measure","Map_result", "result", "time", "Aarm", "Barm",
            "nA0", "nA1","nAc","nB0", "nB1", "nBc", "Map_ReqUnAdj", "CA",
            "BFmean0A", "BFsd0A", "BFmean0B", "BFsd0B",
            "BFmean1A", "BFsd1A", "BFmean1B", "BFsd1B",
            "Mean_Age", "SD_Age", "Prop_Male", "fu_months", "country",	"chart_ref", 
            "chart", "calcICC","ICC")
  
  colnames(df)<- cols
  df
}


#calculate the average baseline mean & SD across trials for a specified outcome and age group
av_baseline <- function(df, outcome, age_gp){
  #subset by outcome & age & result = baseline + FU 
  df <- subset(df, measure == outcome & age == age_gp & result == "Baseline and FU")
  
  #only include observations that haven't been cluster adjusted
  df <- subset(df, CRCT=="N" | (CRCT=="Y" & CA=="Y"))
  
  #create a df for the reference arm (Barm) at baseline (0)
  df_B <- df[, c("studyID", "time", "Barm", "BFmean0B", "BFsd0B")]
  #remove duplicated baselines (multi-arm/multi-follow-up studies)
  df_B <- df_B %>% distinct(studyID, Barm, .keep_all = TRUE)
  
  #create a df for the intervention arm (Aarm) at baseline (0)
  df_A <- df[, c("studyID", "time", "Aarm", "BFmean0A", "BFsd0A")]
  #remove duplicated baselines (multi-follow-up studies)
  df_A <- df_A %>% distinct(studyID, Aarm, .keep_all = TRUE)
  
  #rename columns to allow merging of dfs
  cols <- c("studyID", "time", "arm", "BFmean0", "BFsd0")
  colnames(df_A) <- cols
  colnames(df_B) <- cols
  
  #merge
  dfAB <- rbind(df_A, df_B)
  
  #return a vector: [1] = mean, [2] = SD
  vec <- c(mean(dfAB$BFmean0), mean(dfAB$BFsd0))
  vec
}


#convert mapped ds into the same format as the original dataset
#Create columns for Adjusted FU & CS and move this data to the correct columns
edit_dfmap <- function(df){
  
  #add columns for CS
  df$CSmeanA <- NA
  df$CSsdA <- NA
  df$CSmeanB <- NA
  df$CSsdB <- NA
  
  #add columns for Adjusted FU
  df$AFmeanA <- NA
  df$AFsdA <- NA
  df$AFmeanB <- NA
  df$AFsdB <- NA
  
  #move Adjusted FU results to AF columns
  for(i in 1:nrow(df)){
    if(df$result[i]=="Adjusted FU"){
      df$AFmeanA[i] <- df$BFmean1A[i]
      df$AFsdA[i] <- df$BFsd1A[i]
      df$AFmeanB[i] <- df$BFmean1B[i]
      df$AFsdB[i] <- df$BFsd1B[i]
      
      df$BFmean1A[i] <- NA
      df$BFsd1A[i] <- NA
      df$BFmean1B[i] <- NA
      df$BFsd1B[i] <- NA
    }
  }
  
  df
}
