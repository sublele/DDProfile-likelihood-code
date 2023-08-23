DC.MLE.fn= function() {
  for (k in 1:ncl) {
    for(i in 2:(n+1)){
      Y[(i-1), k] ~ dpois(exp(X[i, k])) 
      X[i, k] ~ dnorm(mu[i, k], 1 / sigma^2) 
      mu[i, k] <- X[(i-1), k] + log(lambda) - log(1 + beta * exp(X[(i-1), k]))
    }
    X[1, k] ~ dnorm(mu0, 1 / sigma^2) 
  }
  # Priors on model parameters: They are on the real line.
  
  ln.beta ~ dnorm(0, 0.1) 
  ln.sigma ~ dnorm(0, 0.1)
  ln.tmp ~ dnorm(0, 0.1)
  # Parameters on the natural scale
  beta <- exp(ln.beta)
  sigma <- exp(ln.sigma)
  tmp <- exp(ln.tmp)
  lambda <- tmp + 1
  mu0 <- log(2)  + log(lambda) - log(1 + beta * 2)
}

DC.profile.fn= function() {
  for (k in 1:ncl) {
    for(i in 2:(n+1)){
      Y[(i-1), k] ~ dpois(exp(X[i, k])) 
      X[i, k] ~ dnorm(mu[i, k], 1 / sigma^2) 
      mu[i, k] <- X[(i-1), k] + log(lambda) - log(1 + beta * exp(X[(i-1), k]))
      
    }
    X[1, k] ~ dnorm(mu0, 1 / sigma^2)
  }
  
  # Priors on model parameters (Use the MLE and its variance)
  
  parms ~ dmnorm(MLE,FI)
  beta <- exp(parms[1])
  sigma <- exp(parms[2])
  tmp <- exp(parms[3])
  lambda <- tmp + 1
  mu0 <- log(2)  + log(lambda) - log(1 + beta * 2)
  
}

# Following function is only useful for the case when 
# parameter of interest is a component.

Bayes.fn.lambda= function() {
  for (k in 1:ncl) {
    for(i in 2:(n+1)){
      Y[(i-1), k] ~ dpois(exp(X[i, k])) 
      X[i, k] ~ dnorm(mu[i, k], 1 / sigma^2) 
      mu[i, k] <- X[(i-1), k] + log(lambda) - log(1 + beta * exp(X[(i-1), k]))
    }
    X[1, k] ~ dnorm(mu0, 1 / sigma^2) 
  }
  # Priors on model parameters: 
  
  beta ~ dlnorm(0, 0.1) 
  sigma ~ dlnorm(0, 0.1)
  tmp ~ dlnorm(0, 0.1)
  lambda <- tmp + 1   # This is just shifted flat dist.
  mu0 <- log(2)  + log(lambda) - log(1 + beta * 2)
}

Bayes.fn.K= function() {
  for (k in 1:ncl) {
    for(i in 2:(n+1)){
      Y[(i-1), k] ~ dpois(exp(X[i, k])) 
      X[i, k] ~ dnorm(mu[i, k], 1 / sigma^2) 
      mu[i, k] <- X[(i-1), k] + log(lambda) - log(1 + beta * exp(X[(i-1), k]))
    }
    X[1, k] ~ dnorm(mu0, 1 / sigma^2) 
  }
  # Priors on model parameters: 
  
  tmp ~ dlnorm(0, 0.1) 
  sigma ~ dlnorm(0, 0.1)
  Carry ~ dlnorm(0,0.1)
  lambda <- tmp + 1
  beta <- tmp/Carry
  mu0 <- log(2)  + log(lambda) - log(1 + beta * 2)
}
