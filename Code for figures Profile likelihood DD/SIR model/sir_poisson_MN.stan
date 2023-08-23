// Stan model file for SIR model


// Function for SIR

functions {
real[] sir(real t, real[] y, real[] theta,
real[] x_r, int[] x_i) {
real S = y[1];
real I = y[2];
real R = y[3];
real N = x_i[1];
real beta = theta[1];
real gamma = theta[2];
real dS_dt = -beta * I * S / N;
real dI_dt = beta * I * S / N - gamma * I;
//real dS_dt = -beta * I * S ;
//real dI_dt = beta * I * S - gamma * I;
real dR_dt = gamma * I;
return {dS_dt, dI_dt, dR_dt};
}
}

data {
int<lower=1> n_days;
real y0[3];
real t0;
real ts[n_days];
int N;
int nclones;
int cases[n_days, nclones];
vector [2] Mu;
matrix [2,2] Sigma;
}

transformed data {
real x_r[0];
int x_i[1] = { N };
}

parameters {
vector [2] parms;
}

transformed parameters{
real <lower=0> y[n_days, 3];
{
real theta[2];
theta[1] = exp(parms[1]);
theta[2] = exp(parms[2]);
y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
}
}

model {
//priors
parms ~ multi_normal(Mu,Sigma); 

//sampling distribution
//col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people
for (i in 1:nclones){
cases[,i] ~ poisson(col(to_matrix(y), 2));}
}

generated quantities {
real R0 = exp(parms[1]) / exp(parms[2]);
real recovery_time = 1 / exp(parms[2]);
}




