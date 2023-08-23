
#' Generating data doubling MCMC samples for profile likelihood
#'
#' @param Pmodel.fn This is the same model function as is used for the 
#'                  computation of the MLE except the prior distribution 
#'                  is the asymptotic distribution of the MLE. 
#' @param dat.p Change the input data list to include MLE and 
#'              Fisher Information matrix as MuP and SigmaP.
#' @param params  Parameter vector to be monitored
#' @param ...     Optional arguments n.adapt, n.iter, n.chains,n.update
#' @param unchanged Nodes in the model that remain unchanged. 
#'                  For example, n the sample size.
#'
#' @return list with elements matrix of MCMC outpute M1, M2 
#'         and diagnostics Rhat.out
#' @export
#'
DDSample_fn = function(Pmodel.fn,dat.p,params,unchanged,...){
  
  model1.fit = dclone::dc.fit(dat.p,params,Pmodel.fn,n.clones=1,multiply="ncl", unchanged=unchanged,...)
  model2.fit = dclone::dc.fit(dat.p,params,Pmodel.fn,n.clones=2,multiply="ncl", unchanged=unchanged,...)
  
  # Extract the MCMC output
  M1 = as.matrix(model1.fit)
  M2 = as.matrix(model2.fit)
  
  # MCMC convergence diagnostics
  Rhat.out = rbind(dcdiag(model1.fit),dcdiag(model2.fit))
  out = list( M1=M1, M2=M2,
             Rhat.out=Rhat.out)
  return(out)
}



#' Profile likelihood function using the Logistic regression method
#'         
#' @param psi.fn This is a user supplied function to compute the parameter of
#'               of interest from the canonical parameters.              
#' @param M1 This is a matrix of MCMC samples from the posterior with 1 clone
#' @param M2 This is a matrix of MCMC samples from the posterior with 2 clones
#' @param grid.length This is the number of grid points at which PL is computed.
#'                    The range is decided based on the range of psi values. 
#'                    Default value is 100.   
#' @param psi.grid These are the grid points at which PL is computed. Default is 
#'                 NULL. User can supply this in place of the grid length.
#' @param method  Choose between 'q'uartic, 'c'oncave or 's'plines as link function.
#'                Concave fit using 'scam' function in the R package 'mgcv'.
#' @param ...     Optional arguments in the psi.fn if needed.
#'
#' @return        list containing two objects. Profile.like is a matrix with first 
#'                column psi.grid, 2nd column is the relative profile loglikelihood value.
#'                Second object is the profile MLE of the parameter of interest and
#'                associated confidence interval at the prescribed alpha value. This
#'                can be used to make sure the procedure worked properly. Profile MLE
#'                should be the same as the one obtained by transforming the canonical
#'                MLE directly. Default for confidence interval is 90%
#' @export
#'
#' @examples
profile_fn = function(psi.fn,M1,M2,grid.length=100,psi.grid=NULL,method='c',CIalpha=0.90,...){
  
  # Apply the psi function on the MCMC output to get the psi transformation
  psi1 = psi.fn(M1,...)
  psi2 = psi.fn(M2,...)
  M12 = data.frame(Y.m=c(rep(0,length(psi1)),rep(1,length(psi2))),X = c(psi1,psi2))
  
  if(is.null(psi.grid)){
  psi.range = range(M12$X)
  psi.grid = seq(psi.range[1],psi.range[2],length.out=grid.length)}
  
  if (method == 'q'){
  # We fit a fourth degree polynomial (orthogonal) regression
  glm.fit12.2 = glm(Y.m ~ poly(X,2),family=binomial,data=M12)
  glm.fit12.3 = glm(Y.m ~ poly(X,3),family=binomial,data=M12)
  glm.fit12.4 = glm(Y.m ~ poly(X,4),family=binomial,data=M12)
  tmp = AIC(glm.fit12.2,glm.fit12.3,glm.fit12.4)
  
  # Choose the model fit object that has minimum AIC
  if (AIC(glm.fit12.2) == min(tmp[,2])){glm.fit12=glm.fit12.2}
  if (AIC(glm.fit12.3) == min(tmp[,2])){glm.fit12=glm.fit12.3}
  if (AIC(glm.fit12.4) == min(tmp[,2])){glm.fit12=glm.fit12.4}
}
  
  if (method == 'c'){ 
    # Fit a concave function
    glm.fit12= scam::scam(Y.m ~ s(X,bs="cv"),family=binomial,data=M12)}
  
  if (method == 's'){ 
    # We will fit using general splines
    glm.fit12= mgcv::gam(Y.m ~ s(X),family=binomial,data=M12)
  }
  
  # Predict on the assigned grid values (psi.grid)
  newdata = data.frame(X = psi.grid)
  profile.like = predict(glm.fit12,newdata) 
  profile.out = data.frame(psi.grid=psi.grid,profile.like=profile.like)
  
  # Compute the PL estimate along with its LPL value
  tmp = arrange(profile.out,profile.like)
  Profile.MLE = tmp[nrow(tmp),]
  
  # Standardize the logPL by subtracting the value at the maximum
  profile.out$profile.like = profile.out$profile.like - as.numeric(Profile.MLE[2])
  
  # Compute the confidence interval
  cutoff = -qchisq(CIalpha,1)/2
  index = which(profile.out$profile.like >= cutoff)
  L_index = min(index)
  U_index = max(index)
  CI = c(profile.out[L_index,1],profile.out[U_index,1])
  out = list(Profile.out=profile.out, Profile.MLE=Profile.MLE,Profile.CI = c(CIalpha,CI))
  
  return(out)
}


  
  
  