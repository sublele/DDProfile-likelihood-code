# Program for computing the log-likelihood for SIR model under
# Normal, Poisson and Negative Binomial errors
# Parameters are unconstrained

mLL_norm = function(beta, gamma, sigma,day, cases,N) {
  beta <- exp(beta) 
  gamma <- exp(gamma)
  sigma <- exp(sigma)
  I0 <- cases[1] # initial number of infectious
  observations <- cases[-1] # the fit is done on the other data points
  predictions <- sir_1(beta = beta, gamma = gamma,
                       S0 = N - I0, I0 = I0, R0 = 0, times = day)
  predictions <- predictions$I[-1] # removing the first point too
  if (any(predictions < 0)) return(NA) # safety
  # returning minus log-likelihood:
  out = -sum(dnorm(x = observations, mean = predictions, sd = sigma, log = TRUE))
  return(out)
}

mLL_pois = function(beta, gamma,day, cases,N) {
    beta <- exp(beta)
    gamma <- exp(gamma)
    I0 <- cases[1] # initial number of infectious
    observations <- cases[-1] # the fit is done on the other data points
    predictions <- sir_1(beta = beta, gamma = gamma,
                         S0 = N - I0, I0 = I0, R0 = 0, times = day)
    predictions <- predictions$I[-1] # removing the first point too
    if (any(predictions < 0)) return(NA) # safety
    # returning minus log-likelihood:
    out = -sum(dpois(x = observations, lambda = predictions, log = TRUE))
    return(out)
  }
mLL_NB = function(beta, gamma,phi,day, cases,N) {
  beta <- exp(beta) 
  gamma <- exp(gamma)
  phi = exp(phi)
  I0 <- cases[1] # initial number of infectious
  observations <- cases[-1] # the fit is done on the other data points
  predictions <- sir_1(beta = beta, gamma = gamma,
                       S0 = N - I0, I0 = I0, R0 = 0, times = day)
  predictions <- predictions$I[-1] # removing the first point too
  if (any(predictions < 0)) return(NA) # safety
  # returning minus log-likelihood:
  out = -sum(dnbinom(x = observations, size = phi, mu = predictions, log = TRUE))
  return(out)
}




