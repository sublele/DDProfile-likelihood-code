# These are specific functions for the multivariate
# normal distribution. 

DC.MLE.fn = function(){
  # Observation process
  for (k in 1:ncl){
    for (i in 1:n){
      Y[i,1:D,k] ~ dmnorm(Mu[1:D],Omega[1:D,1:D])
    }}
  # Priors:  
  parms[1:P] ~ dmnorm(MuP[1:P],PrecP[1:P,1:P])
  Mu[1:D] <- parms[1:D]
  # A is lower triangular and positive definite because diagonal values (same as eigenvalues) are strictly positive.
  for (i in 1:D){
    A[index1[i],index2[i]] <- exp(parms[D+i])
  } 
  for (i in (D+1):S){
    A[index1[i],index2[i]] <- parms[D+i]
  }
  for (i in 2:D){
    for (j in 1:(i-1)){
      A[i,j] <- 0
    }
  }
  Sigma <- A %*% t(A)  
  Omega <- inverse(Sigma)
}

DC.profile.fn = function(){
  # Observation process
  for (k in 1:ncl){
    for (i in 1:n){
      Y[i,1:D,k] ~ dmnorm(Mu[1:D],Omega[1:D,1:D])
    }}
  # Priors:  
  parms[1:P] ~ dmnorm(MLE,FI)
  Mu[1:D] <- parms[1:D]
  # A is lower triangular and positive definite because diagonal values (same as eigenvalues) are strictly positive.
  for (i in 1:D){
    A[index1[i],index2[i]] <- exp(parms[D+i])
  } 
  for (i in (D+1):S){
    A[index1[i],index2[i]] <- parms[D+i]
  }
  for (i in 2:D){
    for (j in 1:(i-1)){
      A[i,j] <- 0
    }
  }
  Sigma <- A %*% t(A)  
  Omega <- inverse(Sigma)
}
