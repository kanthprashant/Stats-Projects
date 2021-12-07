library(Rlab)

# SPRT Function
SPRT <- function(alpha, beta, h0, h1, d, count) {
  # alpha = P{Deciding for h1 when h0 is True} = alpha
  # beta = P{Deciding for h0 when h1 is True} = beta
  num0 = sum(d == 1)
  num1 = sum(d == 0)
  den0 = sum(d == 1)
  den1 = sum(d == 0)
  # Calculating lambda
  lambda = log(((0.45^num0)*(0.55^num1))/((0.55^den0)*(0.45^den1)))
  s = lambda
  # Calculating a and b values
  a = log(beta/(1-alpha))
  b = log((1-beta)/alpha)
  if(s > b){
    #H1 is true and stop
    message = "H1 is True"
    return_list = list(s = s, message = message, tf = TRUE, count = count)
    return(return_list)
  }
  else if(s < a){
    #H0 is true and stop
    message = "H0 is true"
    count = count+1
    return_list = list(s = s, message = message, tf = TRUE, count = count)
    return(return_list)
  }
  else{
    #Collect another observation
    message = "Collect another observation"
    return_list = list(s = s, message = message, tf = FALSE, count = count)
    return(return_list)
  }
}

# Simulation function
simulation <- function(p, count) {
  svector = c()
  h0 = 0.45
  h1 = 0.55
  alpha = 0.01
  beta = 0.01
  tf = FALSE
  count = count
  d = rbinom(1, 1, p)
  a = log(beta/(1-alpha))
  b = log((1-beta)/alpha)
  # while loop to iterate through sprt
  while(!tf){
    return_list = SPRT(alpha, beta, h0, h1, d, count)
    count = return_list$count
    s = return_list$s
    svector = append(svector, s)
    message = return_list$message
    tf = return_list$tf
    if(tf == FALSE){
      d = append(d, rbinom(1, 1, p))
    }
  }
  # Graph using blue points overlayed by a line 
  plot(svector, type="o", col="blue")
  
  abline(h=a, col="red")
  abline(h=b, col="red")
  
  # Create a title with a red, bold/italic font
  pstr = sprintf("%0.2f", p)
  title(main=pstr, col.main="red", font.main=4)
  ret_list = list(message = message, count = count)
  return(ret_list)
}

count = 0
# Iterating the simulation 100 times
for(i in 1:100){
  #ret_list = simulation(p = 0.3, count)
  #ret_list = simulation(p = 0.56, count)
  ret_list = simulation(p = 0.54, count)
  count = ret_list$count
  cat(ret_list$message, "\n")
  cat("count of H0 being true = ", count, "\n")
}

