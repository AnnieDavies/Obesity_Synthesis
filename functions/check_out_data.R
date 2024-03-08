#count the number of arms and number of follow-up times per study
count_arm_FU <- function(data){
  data$num_times <- NA #no. of follow up times
  data$num_int <- NA #no. of arms
  
  data$label <- paste(data$study, data$age, data$measure, sep = " ")
  
  x<- 1
  while(x<=nrow(data)){
    intv <- 2 #count no. of arms (minimum 2)
    j <- x
    while(paste(data$label[j], data$time[j]) == paste(data$label[j+1], data$time[j+1])){
      intv <- intv + 1
      j <- j + 1
    }
    
    k <- x
    tim <- 1 #count number of times (minimum 1)
    #check times
    while(paste(data$label[k], data$Aarm[k]) == paste(data$label[k+(intv-1)], data$Aarm[k+(intv-1)])){
      tim <- tim + 1
      k <- k + (intv-1)
    }
    
    #label all rows for that study with number of arms and time points
    for(i in 1:(tim*(intv-1))){
      data$num_int[x+i-1] <- intv
      data$num_times[x+i-1] <- tim
    }
    x <- x + tim*(intv - 1)
  }
  data
}

#check for any missing data on variances/covariances
check_missing <- function(data){
  
  studies <- character() #studies with missing data
  
  for(i in 1:nrow(data)){
    if(data$num_int[i]>2){
      #check for refvar
      if(is.na(data$refvar[i])){
        print(paste("Missing variance on CS in reference arm for study ", data$label[i]))
        if(!(data$study[i]%in% studies)){
          studies <- append(studies, data$study[i])
        }
      }
    }
    if(data$num_times[i]>1){
      #check var(m_it(k))
      if(is.na(data$varmA[i])){
        print(paste("Missing var(m_it(k)) for study ", data$label[i]))
        if(!(data$study[i]%in% studies)){
          studies <- append(studies, data$study[i])
        }
      }
      #check var(m_it(r))
      if(is.na(data$varmB[i])){
        print(paste("Missing var(m_it(r)) for study ", data$label[i]))
        if(!(data$study[i]%in% studies)){
          studies <- append(studies, data$study[i])
        }
      }
      if(data$result[i]!="Adjusted FU"){
        #check var(m_ib(k))
        if(is.na(data$varmAbase[i])){
          print(paste("Missing var(m_ib(k)) for study ", data$label[i]))
          if(!(data$study[i]%in% studies)){
            studies <- append(studies, data$study[i])
          }
        }
        #check var(m_ib(r))
        if(is.na(data$varmBbase[i])){
          print(paste("Missing var(m_ib(r)) for study ", data$label[i]))
          if(!(data$study[i]%in% studies)){
            studies <- append(studies, data$study[i])
          }
        }
      }
    }
  }
  studies
}