# Functions for analytical profile likelihood 
# (specific to Beta distribution)

floglike.Beta = function(parms,Y,psi.fix){
  out = -sum(dbeta(Y,shape1=exp(parms[1]), shape2= exp(parms[2]),log=T))
  return(out)
}

fpsi = function(parms,Y,psi.fn, psi.fix){
  # Evaluate the constraint function at a given set of parameter value.
  psi = psi.fn(parms)
  return(psi)
}