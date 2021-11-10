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