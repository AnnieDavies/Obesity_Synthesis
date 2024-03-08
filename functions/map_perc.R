#map from percentile to zBMI using the analytic method
map.p.cdf <- function(df_row){
  
  #map baselines (labelled 0)============================================
  if(df_row$result == "Baseline and FU"){
    
    #Intervention (A)-------------------
    pA0 <- df_row$BFmean0A/100
    spA0 <- df_row$BFsd0A/100
    
    #solve simultaneous equations
    resA0 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pA0, spA0, method = "Newton")
    
    #re-define 
    df_row$BFmean0A <- resA0$x[1]
    df_row$BFsd0A <- resA0$x[2]
    
    #Reference (B)----------------------
    pB0 <- df_row$BFmean0B/100
    spB0 <- df_row$BFsd0B/100
    
    resB0 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pB0, spB0, method = "Newton")
    
    #re-define 
    df_row$BFmean0B <- resB0$x[1]
    df_row$BFsd0B <- resB0$x[2]
  }
  
  #Map FU (labelled 1)===================================================
  #Intervention (A)------------------------
  pA1 <- df_row$BFmean1A/100
  spA1 <- df_row$BFsd1A/100
  
  #solve simultaneous equations
  resA1 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pA1, spA1, method = "Newton")
  
  #re-define 
  df_row$BFmean1A <- resA1$x[1]
  df_row$BFsd1A <- resA1$x[2]
  
  #Reference (B)----------------------------
  pB1 <- df_row$BFmean1B/100
  spB1 <- df_row$BFsd1B/100
  
  resB1 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pB1, spB1, method = "Newton")
  
  #re-define 
  df_row$BFmean1B <- resB1$x[1]
  df_row$BFsd1B <- resB1$x[2]

  
  df_row$measure <- "BMI-z from percentile"
  
  df_row
  
}

#simultaneous equations to solve
sim_eqs <- function(x, meanp, SDp){
  z = x[1]
  s = x[2]
  eq1 = meanp - pnorm(z/sqrt(1+s^2))
  
  h <- z/sqrt(1 + s^2)
  a <- 1/sqrt(1 + 2*s^2)
  
  eq2 = SDp^2 - (pnorm(z/sqrt(1+s^2))-2*OwensT(h,a)-pnorm(z/sqrt(1+s^2))^2)
  result <- c(eq1, eq2)
  result
}

