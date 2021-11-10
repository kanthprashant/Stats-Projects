library(EnvStats)
library(Rfast)

# MLE of Binomial Distribution
Binomial_Dist <- function(sample_vec)
{
  # sample_vec = rbinom(n,1,0.5)
  negative_likelihood <- function(p) {
    dbinom(sample_vec, length(sample_vec), p)*-1
  }
  
  nlm(negative_likelihood, 0.5, stepmax = 0.5)
  $estimation
}

# MLE of Uniform Distribution
Unif_Dist <- function(sample_vec)
{
  
  #  sample_vec = runif(n,0,1)
  eunif(sample_vec, method = 'mle')
  $parameters
}

# MLE of Geometric Distribution
Geometric_Dist <- function(sample_vec)
{
  # sample_vec = rgeom(n,0.1)
  egeom(sample_vec, method = "mle")
  $parameters
}

# MLE of Beta Distribution
Beta_Dist <- function(sample_vec)
{
  ebeta(sample_vec)
  $parameters
}

# MLE of Exponential Distribution
Exponential_Dist <- function(sample_vec)
{
  eexp(sample_vec, ci = TRUE, conf = 0.9, method = "mle")
  $parameters
}

# MLE of Multivariate Normal Distribution
Multivariate_Normal_Dist <- function(sample_vec)
{
  mvnorm.mle(sample_vec)
  $loglik
}


# Main Function
main <- function(sample_vec,dist_type)
{
  if(dist_type == "uniform" || dist_type == "Uniform" || dist_type == "UNIFORM")
  {
    Unif_Dist(sample_vec)
  }
  else if(dist_type == "normal" || dist_type == "Normal" || dist_type == "NORMAL")
  {
    Norm_Dist(sample_vec)
  }
  else if(dist_type=="poisson"||dist_type=="Poisson"||dist_type=="POISSON")
  { 
    Poisson_Dist(sample_vec)
  }
  else if(dist_type=="gamma"||dist_type=="Gamma"||dist_type=="GAMMA")
  {
    Gamma_dist(sample_vec)
  }
  else if(dist_type=="beta"||dist_type=="Beta"||dist_type=="BETA")
  {
    Beta_dist(sample_vec)
  }
  else if(dist_type=="geometric"||dist_type=="Geometric"||dist_type=="GEOMETRIC")
  {
    Geometric_Dist(sample_vec)
  }
  else if(dist_type=="binomial"||dist_type=="Binomial"||dist_type=="BINOMIAL")
  {
    Binomial_Dist(sample_vec)
  }
  else if(dist_type=="exponential"||dist_type=="Exponential"||dist_type=="EXPONENTIAL")
  {
    Exponential_Dist(sample_vec)
  }
  else if(dist_type=="multinomial"||dist_type=="Multinomial"||dist_type=="MULTINOMIAL")
  {
    Multinomial_Dist(sample_vec)
  }
  else if(dist_type=="multivariate normal"||dist_type=="Multivariate Normal"||dist_type=="MULTIVARIATE NORMAL")
  {
    Multivariate_Normal_Dist(sample_vec)
  }
}
