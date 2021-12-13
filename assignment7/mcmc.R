require(MASS)
library(LaplacesDemon)
library(mvtnorm)
library(CholWishart)

# Below is the example prior density used for input,output generation
# prior <- function(x)
# {
#   m<-c(1,1,1,1,1)
#   c <-0.5*diag(5)
#   pd<-dmvnorm(x,m,c,log = FALSE)
#   return(pd)
# }

likelihood <- function(x,data){
  # m is mean and has 2 dimensions
  # c is of  covariance 2*2 dimensions
  # x is theta and has 5 dimesnions, with first 2 as mean and next 3 as covariance triangle
  m<-c(0,0) 
  m[1] = x[1]
  m[2] = x[2]
  c <- diag(2)
  c[1,1] = x[3]
  c[1,2] = x[4]
  c[2,1] = x[4]
  c[2,2] = x[5]
  return(prod(dmvnorm(data,m,c))) 
  
}


calculate_mcmc <- function(prior,data,theta_init)
  # prior  contains prior density from user , data is observed data with 2 dimesions
  # theta_init is of 5 dimesnions, with furst 2 dimesions as mean vector and rest as covariance triangle
{
  require(MASS)
  x <- theta_init
  th <- x # th is theta
  for (i in 2:10000)
  {
    y <- mvrnorm(1, x, 1.2*diag(5)) # step 1 :proposed distribution,y is proposed value
    temp_y = prior(y)*likelihood(y,data)
    temp_x = prior(x)*likelihood(x,data)
    if (runif(1) < min(1, temp_y/temp_x))  # step 2 calculate(r(x,y) )
      x <- y # step 3 setting the next value with probability r as proposed value, otherwise previous value 
    
    th <- rbind(th,x)
  } 
  
  require(rgl)
  plot3d(th[,1], th[,2], th[,3], xlab='x',ylab='y', zlab='z',col='red', size=1)
  for(i in 1:5)
  hist(th[,i],xlim=c(0,10),probability = TRUE, xlab=paste("dimension of theta is ",i), main="Histogram of values of theta visited by MCMC algorithm")
  cat("final 5  theta points are \n")
  print(th[(9996):(10000),])
}

#calculate_mcmc(prior,mvrnorm(10,c(5,5),diag(2)),c(1,2,3,4,5))
