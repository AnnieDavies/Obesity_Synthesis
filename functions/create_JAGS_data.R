#Create within study covariance matrix
create_cov <- function(data, yind, numcont, numFU, N, rhoy_sm, rhoy_sl, rhoy_ml,
                       rhod_sm, rhod_sl, rhod_ml, rhoAF_sm, rhoAF_sl, rhoAF_ml){
  #yind is a vector giving the index where each study starts in data
  #numcont is a vector with the number of contrasts (Ai-1) per study
  #numFU is a vector with the number of FU times (Ti) per study
  #N = number of studies
  #rhos give the correlation between s = short, m = medium, l = long term FU
  #rhoy are correlations between MDs, rhod are correlations between CSs and
  #rhoAF are correlations for adjusted FU measurements
  
  #study level dimensions Ti(Ai-1)
  study_dim <- numcont*numFU
  max_dim <- max(study_dim)
  V <- array(dim = c(N, max_dim, max_dim))
  
  for(i in 1:N){
    #create Ti Vit matrices (fixed time point)
    for(t in 1:numFU[i]){
      #Vit has dimensions (Ai-1)x(Ai-1)
      Vit <- matrix(nrow = numcont[i], ncol = numcont[i])
      
      for(k in 1:numcont[i]){
        index <- yind[i] + (t-1)*numcont[i] + (k-1)
        #diagonal elements are equal to SE2^2
        Vit[k,k] <- data$SE2[index]**2
        
        #off diag are equal to refvar
        if((k+1)<=numcont[i]){
          for(kp in (k+1):numcont[i]){
            Vit[k,kp] <- data$refvar2[yind[i]+ (t-1)*numcont[i]] #same for all at the same time point
            Vit[kp,k] <- Vit[k,kp] #symmetric
          }
        }
      }
      
      #add to V
      V[i, (1+(t-1)*numcont[i]):(t*numcont[i]), (1+(t-1)*numcont[i]):(t*numcont[i])] <- Vit
      
      #create Vitt' (off diagonal matrices)
      if((t+1)<=numFU[i]){
        for(tp in (t+1):numFU[i]){
          Vittp <- matrix(nrow = numcont[i], ncol = numcont[i])
          #need to find out whether t and tp are s, m, or l--------------------
          t_ind <- yind[i] + (t-1)*numcont[i] #same for all k
          tp_ind <- yind[i] + (tp-1)*numcont[i] #same for all k
          
          t_time <- data$time[t_ind]
          tp_time <- data$time[tp_ind]
          
          #choose rho based on times
          if(t_time=="short" && tp_time=="mid"){
            rhoy <- rhoy_sm
            rhod <- rhod_sm
            rhoAF <- rhoAF_sm
          }else if (t_time=="short" && tp_time=="long"){
            rhoy <- rhoy_sl
            rhod <- rhod_sl
            rhoAF <- rhoAF_sl
          }else if(t_time=="mid" && tp_time=="long"){
            rhoy <- rhoy_ml
            rhod <- rhod_ml
            rhoAF <- rhoAF_ml
          }else{
            print(paste("Error t_time, study = ", data$studyID[yind[i]], ", yind = ", yind[i], ", i = ", i))
          }
          
          for(k in 1:numcont[i]){
            tk_ind <- yind[i] + (t-1)*numcont[i] + (k-1)
            tpk_ind <- yind[i] + (tp-1)*numcont[i] + (k-1)
            
            #diagonal elements t!=t', k=k' (variance of MDs) - same for CS & AF
            vy_itk <- data$SE2[tk_ind]**2
            vy_itpk <- data$SE2[tpk_ind]**2
            
            Vittp[k,k] <- rhoy*sqrt(vy_itk*vy_itpk)
            
            
            #off diag elements t!=t', k!=k'  
            vd_itk <- data$refvar2[tk_ind]
            vd_itpk <- data$refvar2[tpk_ind]
            
            if((k+1)<=numcont[i]){
              for(kp in (k+1):numcont[i]){
                if(data$result[yind[i]]== "Adjusted FU"){
                  Vittp[k,kp] <- rhoAF*sqrt(vd_itk*vd_itpk) #doesn't depend on k or k'
                }else{
                  Vittp[k,kp] <- rhod*sqrt(vd_itk*vd_itpk) #doesn't depend on k or k'
                }
                Vittp[kp,k] <- Vittp[k,kp] #symmetric
              }
            }
          }
          
          #add to V
          V[i, (1+(t-1)*numcont[i]):(t*numcont[i]), (1+(tp-1)*numcont[i]):(tp*numcont[i])] <- Vittp
          V[i, (1+(tp-1)*numcont[i]):(tp*numcont[i]), (1+(t-1)*numcont[i]):(t*numcont[i])] <- t(Vittp) #Vit't = Vitt'^T
          #NB: transpose doesn't matter bc matrix is symmetric anyway:
          #off diag elements for t!=t' and k!=k' are equal for given t
        }#end of t' loop
      }#end if (t+1)<=Ti
    }#end of t loop
  }#end of i loop
  V
}

#create indicator for where each new study begins in the vector of observations
create_yind <- function(N, numFU, numcont){
  yind <- c()
  sum_ind <- 1
  for(i in 1:N){
    yind <- append(yind, sum_ind)
    sum_ind <- sum_ind + numFU[i]*numcont[i]
  }
  yind
}


#Create matrix with ones on diag and 0.5 on off diag (per study - correct dimensions)
#to create the between study covariance matrix
create_S <- function(numFU, numcont, N){
  study_dim <- numcont*numFU
  max_dim <- max(study_dim)
  S <- array(dim = c(N, max_dim, max_dim)) 
  
  #Sigma = tau^2*S (S=(1,0.5...\\0.5,1,...\\...))
  for(i in 1:N){
    for(j in 1:(numFU[i]*numcont[i])){
      for(k in 1:(numFU[i]*numcont[i])){
        if(j==k){
          S[i,j,k] <- 1
        }else{
          S[i,j,k] <- 0.5
        }
      }
    }
  }
  S
}

#create ragged array of intervention covariates
create_x <- function(df, yind, numcont, N, A, covs){
  #N is number of studies
  #A is maximum number of arms
  #covs is a vector of covariate names (must match colnames in df)
  
  n <- length(covs) #n = number of covariates
  x <- array(dim = c(N, A-1, n))
  
  for(i in 1:N){ #loop over study
    for(k in 1:numcont[i]){#loop over arms
      for(j in 1:n){#loop over covariates
        cind <- which(colnames(df)==covs[j]) #index of column
        #df[[row, column]]
        x[i, k, j] <- df[[yind[i] + (k-1), cind]]
      }
    }
  }
  x
}

#create an array of study level covariates
create_z <- function(df, yind, N, covs){
  #N is number of studies
  #covs is a vector of covariate names (must match colnames in study_covs)
  
  p <- length(covs)
  z <- array(dim = c(N, p)) 
  
  for(i in 1:N){ #loop over study
    for(j in 1:p){#loop over covariates
      cind <- which(colnames(df)==covs[j]) #index of column
      #df[[row, column]]
      z[i, j] <- df[[yind[i], cind]]
    }
  }
  z
}


#create a ragged array of time-dependent covariates
create_w <- function(df, yind, numFU, numcont, N, T, covs){
  #T is maximum number of FU times
  #covs is a vector of covariate names
  
  #q = number of FU time categories - 1 + number of time dep covs
  q <- length(covs)
  w <- array(dim = c(N, T, q)) 
  
  for(i in 1:N){#loop over study
    for(t in 1:numFU[i]){#loop over FU times
      for(j in 1:q){#loop over covariates
        cind <- which(colnames(df)==covs[j]) #index of column 
        w[i,t,j] <- df[[yind[i] + (t-1)*numcont[i], cind]]
      }
    }
  }
  w
}

#create ragged array of interactions
create_J <- function(df, yind, numFU, numcont, N, T, A, covs){
  #N=no. of studies, A=max number of arms, T=max number of FU times
  #covs = vector of interaction terms
  
  l <- length(covs)
  J <- array(dim = c(N,T,A-1,l))
  
  for(i in 1:N){#loop over study
    for(t in 1:numFU[i]){#loop over FU time
      for(k in 1:numcont[i]){#loop over arms
        for(j in 1:l){#loop over interaction term
          cind <- which(colnames(df)==covs[j])
          J[i, t, k, j] <- df[[yind[i] + (t-1)*numcont[i] + (k-1), cind]]
        }
      }
    }
  }
  J
}
