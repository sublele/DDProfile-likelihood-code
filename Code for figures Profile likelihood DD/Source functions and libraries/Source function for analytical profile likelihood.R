# Source functions for analytical profile likelihood when possible

# This is a function that computes the log-likelihood for a model at a parameter
# combination (parms) for given data (Y). The argument 'psi.fix' is unused but needed 
# in order to use package 'alabama' and function 'auglog'. Functions 'floglike' and 'fpsi' 
# are a template. They need to be written for each case separately.  


profile.analytical = function(parms.start,floglike,fpsi,Y,psi.grid,...){
  
  lplike.out = matrix(0,length(psi.grid),2)
  for (i in 1:length(psi.grid)){
    sink("nul")
    tmp = suppressWarnings(alabama::auglag(parms.start,floglike,heq=fpsi,Y=Y,psi.fix=psi.grid[i],...))
    sink()
    lplike.out[i,] = c(psi.grid[i],tmp$val)
  }
  lplike.out = data.frame(psi=lplike.out[,1],profile.like=lplike.out[,2])
  return(lplike.out)
}
