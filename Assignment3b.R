Norm_Dist <- function(sample_vec)
  #finds method of moment estimates for normal distribution
{ 
  cat("\nDistribution Type - Normal")
  
  cat("\nsample input data is\n ", sample_vec)
  e1 = cal_moment1(sample_vec)
  e2 = cal_moment2(sample_vec)
  X = e1 #sample mean
  uNorm = X
  
  cat("\n\nfirst maximum likelihood estimate mean estimate in normal distribution ", uNorm)
  sdev = sqrt(e2 - (uNorm*uNorm))
  cat("\n\nsecond method of moments estimate standard deviation estimate in normal distribution ", sdev)
  
}

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
  
  else if(dist_type =="Multivariate Normal"){
    Multivariate_Normal_dist(sample_vec)
  }
  
}

