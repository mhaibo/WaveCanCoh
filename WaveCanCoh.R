
library(mvLSW)
library(zoo)
library(stats)
library(ggplot2)





# The input of function extract_spectrum, here we extract based on the structure of S = [Sxx,Sxy; Syx, Syy]
extract_spectrum <- function(spectrum,size_XY){
  dim_spec <- dim(spectrum)
  nrows <- dim_spec[1]
  ncols <- dim_spec[2]
  J <- dim_spec[3]
  T <- dim_spec[4]
  Px <- size_XY[1]
  Py <- size_XY[2]
  S_xx <- spectrum[1:Px, 1:Px,,]
  S_xy <- spectrum[1:Px,(Px+1):ncols,,]
  S_yx <- spectrum[(Px+1):nrows,1:Px,,]
  S_yy <- spectrum[(Px+1):nrows,(Px+1):ncols,,]
  S <- list(Sxx=S_xx,Sxy=S_xy,Syx=S_yx,Syy=S_yy)
  return(S)
}


# function for claculating the WaveCanCoh with partitioned sepctrum
# rho[J,T] is the time varying scale-specific canonical coherence between group X and group Y, 
#which is the largest eigenvalue of inv(Sxx) * Sxy * inv(Syy) * Syx, 
# a[J,Px,T] stores the contribution of each channel in X, which is the eigenvector of inv(Sxx) * Sxy * inv(Syy) * Syx, 
#corresponding to the largest eigenvalue
# b[J,Py,T] stores the contribution of each channel in Y,which is the eigenvector of inv(Syy) * Syx * inv(Sxx) * Sxy, 
#corresponding to the largest eigenvalue
WaveCan <- function(Sxx,Sxy,Syx,Syy){
  T <- dim(Sxx)[4]
  J <- dim(Sxx)[3]
  Px <- dim(Sxx)[1] # dim of group X
  Py <- dim(Syy)[1] # dim of group Y
  a <- array(0,dim = c(J,Px,T))
  b <- array(0,dim = c(J,Py,T))
  rho <- array(0,dim = c(J,T))
  
  for (j in 1:J){
    for (t in 1:T){
      
      inv_Sxx <- solve(Sxx[,,j,t])
      inv_Syy <- solve(Syy[,,j,t])
      A <- inv_Sxx %*% Sxy[,,j,t] %*% inv_Syy %*% Syx[,,j,t] 
      B <- inv_Syy %*% Syx[,,j,t] %*% inv_Sxx %*% Sxy[,,j,t] 
      eig_A <- eigen(A)
      eig_B <- eigen(B)
      rho[j,t] <- Re(eig_A$values[1]) # take the real part to avoid numerical issues, the imaginary parts are zero indeed
      a[j,,t] <- Re(eig_A$vectors[,1])
      b[j,,t] <- Re(eig_B$vectors[,1])
    }
  }
  waveCan <- list(rho = rho, a = a, b = b)
  return(waveCan)
}



# calculate the wavelet coherence between two group of multivariate time series X and Y.
# The input includes X and Y, should with formation T x P, filter_number and wavelet_family is used for specifying wavelet functions.
Calculate_WaveCan <- function(X,Y, filter_number = 1, wavelet_family = "DaubExPhase"){
  data_XY <- cbind(X,Y) # combine the X and Y
  Px <- ncol(X) # access the dim of X and Y
  Py <- ncol(Y)
  size_XY <- c(Px,Py)
  
  group_spectrum <- mvEWS(X = data_XY, filter.number = filter_number, family = wavelet_family,  kernel.name = "daniell",
                          optimize = TRUE, bias.correct = TRUE,  tol = 1e-10)  # estimate the group spectrum S= [Sxx,Sxy; Syx, Syy]
  extracted_spectrum <- extract_spectrum(group_spectrum$spectrum,size_XY)  # extract the group spectrum into four corresponding blocks
  
  WaveCan_XY_results <- WaveCan(extracted_spectrum$Sxx,extracted_spectrum$Sxy,extracted_spectrum$Syx,extracted_spectrum$Syy) 
  
  return(WaveCan_XY_results)
}
