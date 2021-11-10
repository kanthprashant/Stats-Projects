# calculate mle of binomial distribution -
calculate_binom_mle <- function(sample) {
  # calculate negative log likelihood
  nll <- function(parameters, data) {
    n = parameters[1]
    p = parameters[2]
    -sum(dbinom(x=data, size=n, prob=p, log=TRUE))
  }
  
  mle = optim(par = c(n = 10, p = 0.5), fn = nll, data = sample, 
              control=list(parscale = c(n = 10, p = 0.5)))
  return(mle)
}

# calculate mle of geometric
calculate_geometric_mle <- function(sample) {
  nll <- function(parameters, data) {
    prob = parameters[1]
    -sum(dgeom(x=data, prob=prob, log=TRUE))
  }
  
  mle = optim(par = c(prob=0.4), fn=nll, data=sample, 
              control = list(parscale = c(prob=0.4)))
  return(mle)
  
}

# calculate mle of poisson
calculate_poisson_mle <- function(sample) {
  nll <- function(parameters, data) {
    lambda = parameters[1]
    -sum(dpois(x = data, lambda = lambda, log=TRUE))
  }
  
  mle = optim(par = c(lambda=8), fn = nll, data = sample, 
              control = list(parscale = c(lambda = 8)))
  return(mle)
}

# calculate mle of exponential -
calculate_exp_mle <- function(sample) {
  nll <- function(parameters, data) {
    lambda = parameters[1]
    -sum(dexp(x = data, rate  = lambda, log=TRUE))
  }
  
  mle = optim(par = c(lambda=8), fn = nll, data = sample, 
              control = list(parscale = c(lambda = 8)))
  return(mle)
}

# calculate mle of beta -
calculate_beta_mle <- function(sample) {
  nll <- function(parameters, data) {
    alpha = parameters[1]
    beta = parameters[2]
    -sum(dbeta(x = data, shape1 = alpha, shape2 = beta, log=TRUE))
  }
  
  mle = optim(par = c(alpha = 1, beta = 10), fn = nll, data = sample, 
              control=list(parscale = c(alpha = 1, beta = 10)))
  return(mle)
}

calculate_uniform_mle <- function(sample) {
  nll <- function(parameters, data) {
    min = parameters[1]
    max = parameters[2]
    -sum(dunif(x = data, min = min, max = max, log=TRUE))
  }
  
  mle = optim(par = c(min = -10, max = 10), fn = nll, data = sample, 
              control=list(parscale = c(min = -10, max = 10)))
  return(mle)
}

# calculate mle of normal distribution 
calculate_normal_mle <- function (sample) {
  # calculate negative log likelihood
  nll <- function(parameters, data) {
    mu = parameters[1]
    sigma = parameters[2]
    -sum(dnorm(x = data, mean = mu, sd = sigma, log = TRUE))
  }
  
  mle = optim(par = c(mu=0.2, sigma = 1.5), fn = nll, data = sample, 
              control = list(parscale = c(mu = 0.2, sigma = 1.5)))
  
  return(mle)
}


main <- function (sample_vec,dist_type)
{
    
    
    if(dist_type == "uniform" || dist_type == "Uniform" || dist_type == "UNIFORM")
    {
      uniform_sample = runif(10, -10, 10)
      uniform_mle = calculate_uniform_mle(uniform_sample)
      print(uniform_mle$par)
      
    }
    
    else if(dist_type == "beta" || dist_type == "Beta" || dist_type == "BETA")
    {
      sample = rbeta(10,1,10)
      cat(sample)
      sample_mle = calculate_beta_mle(sample)
      print(sample_mle$par)
      
    }
    else if(dist_type == "normal" || dist_type == "Normal" || dist_type == "NORMAL")
    {
      sample = rnorm(10)
      normal_mle = calculate_normal_mle(sample)
      print(normal_mle$par)
    }
    else if(dist_type=="poisson"||dist_type=="Poisson"||dist_type=="POISSON")
    {
      
      poisson_sample = rpois(10, 8)
      poisson_mle = calculate_poisson_mle(poisson_sample)
      print(poisson_mle$par)
      
    }
    else if(dist_type=="gamma"||dist_type=="Gamma"||dist_type=="GAMMA")
    {
      
      Gamma_dist(sample_vec)
    }
    
    
    
    else if(dist_type=="geometric"||dist_type=="Geometric"||dist_type=="GEOMETRIC")
    {
      
      geom_sample = rgeom(10, 0.4)
      geom_mle = calculate_geometric_mle(geom_sample)
      print(geom_mle$par)
    }
    else if(dist_type=="binomial"||dist_type=="Binomial"||dist_type=="BINOMIAL")
    {
      binom_sample = rbinom(4,10,0.5)
      binom_mle = calculate_binom_mle(binom_sample)
      print(binom_mle$par)
    }
    else if(dist_type=="exponential"||dist_type=="Exponential"||dist_type=="EXPONENTIAL")
    {
      sample = rexp(5)
      cat(sample)
      sample_mle = calculate_exp_mle(sample)
      print(sample_mle$par)
    }
    else if(dist_type=="multinomial"||dist_type=="Multinomial"||dist_type=="MULTINOMIAL")
    {
      Multinomial_Dist(sample_vec)
    }
    
    else if(dist_type =="Multivariate Normal")
    {
      Multivariate_Normal_dist(sample_vec)
    }
    
  }
  
  

