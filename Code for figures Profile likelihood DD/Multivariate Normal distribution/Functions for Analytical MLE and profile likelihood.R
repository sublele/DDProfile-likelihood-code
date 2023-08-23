# Functions for analytical profile likelihood

floglike.MN = function(parms,Y,psi.fix){
  # Parms are the mean and the Cholesky decomposition. Need to create the Var matrix from the Cholesky decomposition.
  P = length(parms)
  D = ncol(Y)
  N = nrow(Y)
  mu1 = parms[1:D]
  L = matrix(0,D,D)
  L[lower.tri(L,diag=T)]=parms[(D+1):P]
  for (i in 1:D){
    L[i,i] = exp(L[i,i])}
  Var1 = L %*% t(L)
  out = -sum(dmvnorm(Y,mu1,Var1,log=T))
  return(out)
}

# We can use different psi.fn to get different profile likelihoods. 
fpsi = function(parms,Y,psi.fn, psi.fix,...){
  # Evaluate the constraint function at a given set of parameter value.
  psi = psi.fn(parms,...)
  return(psi)
}
