library(Rlab)

SPRT <- function(alpha, beta, h0, h1, d) {
  # alpha = P{Deciding for h1 when h0 is True} = alpha
  # beta = P{Deciding for h0 when h1 is True} = beta
  num0 = sum(d == 1)
  num1 = sum(d == 0)
  den0 = sum(d == 1)
  den1 = sum(d == 0)
  lambda = log(((0.45^num0)*(0.55^num1))/((0.55^den0)*(0.45^den1)))
  #cat("lambda: ", lambda, "\n")
  s = lambda
  #cat("s: ", s, "\n")
  a = log(beta/(1-alpha))
  #cat("a: ", a, "\n")
  b = log((1-beta)/alpha)
  #cat("b: ", b, "\n")
  if(s > b){
    #H1 is true and stop
    message = "H1 is True"
    return_list = list(s = s, message = message, tf = TRUE)
    return(return_list)
  }
  else if(s < a){
    #H0 is true and stop
    message = "H0 is true"
    return_list = list(s = s, message = message, tf = TRUE)
    return(return_list)
  }
  else{
    #Collect another observation
    message = "Collect another observation"
    return_list = list(s = s, message = message, tf = FALSE)
    return(return_list)
  }
}

simulation <- function(p) {
  svector = c()
  h0 = 0.45
  h1 = 0.55
  alpha = 0.01
  beta = 0.01
  tf = FALSE
  d = rbinom(1, 1, p)
  a = log(beta/(1-alpha))
  b = log((1-beta)/alpha)
  #cat(a, b)
  # while loop to iterate through sprt
  while(!tf){
    return_list = SPRT(alpha, beta, h0, h1, d)
    s = return_list$s
    svector = append(svector, s)
    message = return_list$message
    tf = return_list$tf
    if(tf == FALSE){
      d = append(d, rbinom(1, 1, p))
    }
  }
  # Graph cars using blue points overlayed by a line 
  plot(svector, type="o", col="blue")
  
  abline(h=a, col="red")
  abline(h=b, col="red")
  
  # Create a title with a red, bold/italic font
  pstr = sprintf("%0.2f", p)
  title(main=pstr, col.main="red", font.main=4)
  return(message)
}

simulation(p = 0.3)
simulation(p = 0.56)
simulation(p = 0.54)

