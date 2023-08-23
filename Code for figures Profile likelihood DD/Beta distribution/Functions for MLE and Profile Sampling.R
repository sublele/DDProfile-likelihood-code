# These are functions to compute the profile likelihood using data cloning and 
# data doubling.

DC.MLE.fn = function(){
  
  # Observations model
  for (k in 1:ncl){
    for (i in 1:n){
      Y[i,k] ~ dbeta(a1,b1)
    }}
  
  # Priors: We use the parameterization that is unconstrained. 
  a1 <- exp(parms[1])
  b1 <- exp(parms[2])
  parms ~ dmnorm(MuP,PrecP)
}

DC.profile.fn = function(){
  # Observations model
  for (k in 1:ncl){
    for (i in 1:n){
      Y[i,k] ~ dbeta(a1,b1)
    }}
  # Prior: This is based on the asymptotic normal distribution
  # of the MLE. These are unconstrained parameters.
  parms ~ dmnorm(MLE,FI)
  a1 <- exp(parms[1])
  b1 <- exp(parms[2])
}


# These are functions to compute the profile likelihood using  
# data doubling when the parameter of interest is a component of the canonical 
# parameter vector. For the Beta distribution, we have to change the
# parameterization to the appropriate vector. Then we only need to write
# the Bayesian model with flat priors. 

Bayes.fn = function(){
  
  # Observations model
  for (k in 1:ncl){
    for (i in 1:n){
      Y[i,k] ~ dbeta(a1,b1)
    }}
  
  # Priors: We use the parameterization so that parameter of interest
  # is the first component of the canonical parameter vector. The flat prior
  # should be on the appropriate scale. We need to monitor only the 
  # parameter of interest.
  
  mu1 ~ dunif(0,1)  # mu1 = a1/(a1+b1)
  c1 ~ dgamma(0.1,0.1)   # a1 + b1 = c1
  a1 <- mu1 * c1
  b1 <- (a1/mu1)-a1
}

