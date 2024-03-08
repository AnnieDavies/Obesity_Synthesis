# calculate mean difference and SE for each arm and time point
calc_MD <- function(data, rho, not_cov, alt_cov){
  #data = dataframe, rho = correlation coeff for CS
  
  # Create MD and SE columns
  data$MD <- NA
  data$SE <- NA
  
  #create CSref columns
  data$refCS <- NA #CS in reference arm [or Adj FU]
  
  # Create reference variance columns (for covariances of multi-arm/multi-FU)
  data$refvar <-NA #variance of CS in reference arm var(CS) [or var(m_it(r)) for Adj FU and FU]
  data$varmA <- NA #variance on FU mean in intervention arm var(m_it(k))
  data$varmB <- NA #variance on FU mean in reference arm var(m_it(r))
  data$varmAbase <- NA #variance on baseline mean in intervention arm var(m_ib(k))
  data$varmBbase <- NA #variance on baseline mean in reference arm var(m_ib(r))
  
  #label observation (study, outcome measure and time point)
  data$ref <- paste(data$study, data$measure, data$time, sep = " ") 
  
  for(i in 1:nrow(data)){
    
    # If Result type = MD (contrast level)
    if(data$result[i]=="MD"){
      # Set MD = MD and SE = SE
      data$MD[i] <- data$MDmean[i]
      data$SE[i] <- data$MDse[i]
      
      #NB:check for reference variance in dataset (either CS or Base+FU)
      #This will probably be an imputed value but we will store it in the spreadsheet
      if(!is.na(data$CSsdB[i])){  #First check CS
        sigB <- data$CSsdB[i]
        nB <- data$nB1[i] 
        varB <- sigB**2/nB
        
        #Reference variance
        data$refvar[i] <- varB
      }else if(!is.na(data$BFsd0B[i]) & !is.na(data$BFsd1B[i])){ #Then check baseline & FU
        sigB0 <- data$BFsd0B[i]
        nB0 <- data$nB0[i]
        sigB1 <- data$BFsd1B[i]
        nB1 <- data$nB1[i]
        varB0 <- sigB0**2/nB0
        varB1 <- sigB1**2/nB1
        covB <- (2*rho*sigB0*sigB1)/sqrt(nB0*nB1)
        
        #Reference variance
        data$refvar[i] <- varB0+varB1-covB
        
      }else{ #else print that we are missing this data
        #print(paste("Missing data on variance of CS in reference arm for study ", data$study[i]))
        #NB: this is checked in check_out_data (for multi-arm studies)
      }
      
      #see if there is information about CS in reference arm
      if(!is.na(data$CSmeanB[i])){  #First check CS
        data$refCS[i] <- data$CSmeanB[i]
      }else if(!is.na(data$BFmean0B[i]) & !is.na(data$BFmean1B[i])){ #then baseline and FU
        data$refCS[i] <- data$BFmean1B[i] - data$BFmean0B[i]
      }
      
    }
    # Else if Result type = Adjusted FU, CS or FU, calculate MD in relevant arm level means
    else if(data$result[i]=="Adjusted FU" || data$result[i]=="CS" || data$result[i]=="FU"){ 
      if(data$result[i]=="Adjusted FU"){
        # If Adjusted FU, define mean/SD A/B using Adj FU columns
        muA <- data$AFmeanA[i]
        sigA <- data$AFsdA[i]
        muB <- data$AFmeanB[i]
        sigB <- data$AFsdB[i]
      } else if(data$result[i]=="CS"){
        # Else if CS, define mean/SD A/B using CS columns
        muA <- data$CSmeanA[i]
        sigA <- data$CSsdA[i]
        muB <- data$CSmeanB[i]
        sigB <- data$CSsdB[i]
      } else if(data$result[i]=="FU"){
        # Else if FU, define mean/SD A/B using FU columns
        muA <- data$BFmean1A[i]
        sigA <- data$BFsd1A[i]
        muB <- data$BFmean1B[i]
        sigB <- data$BFsd1B[i]
      }
      else{ print(paste("result AF/CS error:", data$study[i])) }
      #participant numbers
      nA <- data$nA1[i]
      nB <- data$nB1[i] 
      #group variances 
      varA <- sigA**2/nA
      varB <- sigB**2/nB
      #calculate MD and SE
      data$MD[i] <- muA - muB
      data$SE[i] <- sqrt(varA+varB)
      
      #Reference variance
      data$refvar[i] <- varB
      data$refCS[i] <- muB
      
    }
    # Else if Result type = baseline and FU, calculate change score per group and then MD
    else if(data$result[i]=="Baseline and FU"){
      #baseline means/SDs & n
      muA0 <- data$BFmean0A[i]
      sigA0 <- data$BFsd0A[i]
      muB0 <- data$BFmean0B[i]
      sigB0 <- data$BFsd0B[i]
      nA0 <- data$nA0[i]
      nB0 <- data$nB0[i]
      #FU means/SDs & n
      muA1 <- data$BFmean1A[i]
      sigA1 <- data$BFsd1A[i]
      muB1 <- data$BFmean1B[i]
      sigB1 <- data$BFsd1B[i]
      nA1 <- data$nA1[i]
      nB1 <- data$nB1[i]
      
      #calculate variances and covariance between base and FU
      varA0 <- sigA0**2/nA0
      varA1 <- sigA1**2/nA1
      covA <- (2*rho*sigA0*sigA1)/sqrt(nA0*nA1)
      varB0 <- sigB0**2/nB0
      varB1 <- sigB1**2/nB1
      covB <- (2*rho*sigB0*sigB1)/sqrt(nB0*nB1)
      
      #calculate MDs
      data$MD[i] <- (muA1-muA0) - (muB1-muB0)
      data$SE[i] <- sqrt(varA0+varA1-covA+varB0+varB1-covB)
      
      #Reference variance
      data$refvar[i] <- varB0+varB1-covB
      data$refCS[i] <- muB1-muB0
      
    }else{
      # Else: print error
      print(paste("result error:", data$study[i]))
    }
    
    
    #arm level variances
    #I will extract/impute arm level variances where possible for MD/CS studies (and input into data sheets)
    #Differences will be for Adjusted FU & FU only studies (see model development section 8)
    
    #NB: in the data, baseline and reference arm data is repeated over rows where necessary
    #therefore we can use [i] indexing and will get the same covariances where appropriate
    
    if(data$result[i]=="Adjusted FU"){
      sigA <- data$AFsdA[i] #SD on adj FU mean in intervention arm var(m_it(k))
      sigB <- data$AFsdB[i] #SD on adj FU mean in control arm var(m_it(r))
      
      #participant numbers
      nA <- data$nA1[i]
      nB <- data$nB1[i] 
      varA <- sigA**2/nA
      varB <- sigB**2/nB
      
      data$varmA[i] <- varA #variance on FU mean in intervention arm var(m_it(k))
      data$varmB[i] <- varB #variance on FU mean in reference arm var(m_it(r))
      
      #baseline variances are NA for Adj FU (and FU only - see below)
    }else{
      #baseline SDs & n
      sigA0 <- data$BFsd0A[i]
      sigB0 <- data$BFsd0B[i]
      nA0 <- data$nA0[i]
      nB0 <- data$nB0[i]
      #FU SDs & n
      sigA1 <- data$BFsd1A[i]
      sigB1 <- data$BFsd1B[i]
      nA1 <- data$nA1[i]
      nB1 <- data$nB1[i]
      
      #calculate variances 
      varA0 <- sigA0**2/nA0
      varA1 <- sigA1**2/nA1
      varB0 <- sigB0**2/nB0
      varB1 <- sigB1**2/nB1
      
      data$varmA[i] <- varA1 #variance on FU mean in intervention arm var(m_it(k))
      data$varmB[i] <- varB1 #variance on FU mean in reference arm var(m_it(r))
      data$varmAbase[i] <- varA0 #variance on baseline mean in intervention arm var(m_ib(k))
      data$varmBbase[i] <- varB0 #variance on baseline mean in reference arm var(m_ib(r))
    }
    
    #Change reference variance for complex studies (multi-arm)-----------------
    if(data$study[i]%in%not_cov){
      #find the index of the study within the list of not_cov studies
      ind <- which(not_cov %in% data$study[i])
      #covariance is then the correpsonding element of the vector alt_cov (alternative covariances)
      data$refvar[i] <- alt_cov[ind]
    }
  }
  data
}


#adjust SE and all variances for clustering where necessary
adj_cluster <- function(data, icc){
  #add column for adjusted SE/variances
  data$SE2 <- NA
  data$refvar2 <-NA #variance of CS in reference arm var(CS) [or var(m_it(r)) for Adj FU and FU]
  data$varmA2 <- NA #variance on FU mean in intervention arm var(m_it(k))
  data$varmB2 <- NA #variance on FU mean in reference arm var(m_it(r))
  data$varmAbase2 <- NA #variance on baseline mean in intervention arm var(m_ib(k))
  data$varmBbase2 <- NA #variance on baseline mean in reference arm var(m_ib(r))
  
  for(i in 1:nrow(data)){
    # If Requires Cluster Adjusting = Y:
    if(data$CA[i]=="Y"){
      #Define ICC by calculated value or imputed icc------------
      if(!is.na(data$calcICC[i]) && data$calcICC[i]=="Y"){ # If Calculated ICC = Y: del = ICC[]
        del <- data$ICC[i]
      }else if(is.na(data$calcICC[i]) || data$calcICC[i]=="N"){ # Else if Calculated ICC= N: del = icc
        del <- icc
      }else{print(paste("CalcICC error:", data$study[i]))}
      
      #check nAc and nBc (no. of clusters per group) are defined
      if(is.na(data$nAc[i]) || is.na(data$nBc[i])){print(paste("n cluster:", data$study[i]))}
      
      #mean cluster size-------------------------------------------
      nA1 <- data$nA1[i]
      nB1 <- data$nB1[i]
      nA0 <- data$nA0[i]
      nB0 <- data$nB0[i]
      nAc <- as.numeric(data$nAc[i])
      nBc <- as.numeric(data$nBc[i])
      m <- (nA1+nB1)/(nAc+nBc)
      mB <- nB1/nBc
      mA <- nA1/nAc
      mB0<- nB0/nBc
      mA0 <- nA0/nAc
      
      #adjust SE ---------------------------------------------------
      data$SE2[i] <- data$SE[i]*sqrt(1+del*(m-1)) #SE' = SE*SQRT(1+del(m-1))
      
      #adjust variances
      data$refvar2[i] <- data$refvar[i]*(1+del*(mB-1))
      data$varmA2[i] <- data$varmA[i]*(1+del*(mA-1))
      data$varmB2[i] <- data$varmB[i]*(1+del*(mB-1))
      data$varmAbase2[i] <- data$varmAbase[i]*(1+del*(mA0-1))
      data$varmBbase2[i] <- data$varmBbase[i]*(1+del*(mB0-1))
      
    }else if(data$CA[i]=="N"){ # Else if does not require cluster adjusting: SE' = SE
      
      #don't adjust SE/variances---------------------------------------------
      data$SE2[i] <- data$SE[i]
      data$refvar2[i] <- data$refvar[i]
      data$varmA2[i] <- data$varmA[i]
      data$varmB2[i] <- data$varmB[i]
      data$varmAbase2[i] <- data$varmAbase[i]
      data$varmBbase2[i] <- data$varmBbase[i]
      
    }else{print(paste("CA Y/N error:", data$study[i]))}
    
  }
  
  data <- data[, c("studyID","study","age","measure", "time","result","Aarm",
                   "nA1","Barm","nB1","MD","SE","SE2","refCS", "refvar", "refvar2", 
                   "varmA", "varmA2", "varmB","varmB2","varmAbase", "varmAbase2","varmBbase",
                   "varmBbase2")]
  data
  
}