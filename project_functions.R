### Functions Used in our project

#log full conditional of lambda

log_full_conditional_lambda <- function(lambda){
  
  if(alpha > 0 & lambda > 0 & beta > 0){
    sum(y*log(lambda)) - (length(y)*lambda) + ((alpha-1)*log(lambda)) - (lambda*beta)
  } else {
    -Inf
  }
}


#log full conditional of alpha
log_full_conditional_alpha <- function(alpha){
  
  if(alpha > 0 & lambda > 0 & beta > 0){
    (alpha*log(beta)) - log(gamma(alpha)) + (alpha*log(lambda)) + ((c-1)*log(alpha)) - (alpha*d)
  } else {
    -Inf
  }
}


#log full conditional of beta
log_full_conditional_beta <- function(alpha,lambda){
  
  if(alpha > 0 & lambda > 0){
    rgamma(1, shape = alpha + e, rate = lambda + f)
  } else {
    -Inf
  }
}


#Log likelihood function
# Define your log-likelihood function
log_likelihood <- function(par, y) {
  value <- (sum(y) * log(par)) - (length(y)*par) - sum(lfactorial(y))
  return(value)  # Minimization, so negate the log-likelihood
}

#log likelihood function for exponential distribution
log_likelihood_exp <- function(par, y) {
  value <- length(y)*log(par) - par*sum(y)
  return(value)  # Minimization, so negate the log-likelihood
}