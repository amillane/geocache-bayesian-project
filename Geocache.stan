// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of observations
  int<lower=0> id[N]; // group id 
  int<lower=0> J; // type of terrain 
  int<lower=0> y[N]; // data 
  real<lower=0> a; // alpha prior for alpha
  real<lower=0> b; // beta prior for alpha
  real<lower=0> c; // alpha prior for beta 
  real<lower=0> d;// beta prior for beta 
}
  

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> lambda[J]; //
  real<lower=0> alpha[J];
  real<lower=0> beta[J];
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
   for(j in 1:9){
    lambda[j] ~ gamma(alpha[j], beta[j]);
    alpha[j] ~ gamma(a, b);
    beta[j] ~ gamma(c, d);
  }
  
  for(i in 1:N){
    y[i] ~ poisson(lambda[id[i]]);
    }
}


