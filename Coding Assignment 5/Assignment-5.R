library(Rlab)
SPRT <- function(type1.err = 0.01, type2.err = 0.01, h0, h1, rv=dbinom(1, 1, 0.4)) {
  # Bernoulli Random Variable p <= .45 vs P >=.55 using a SPRT
  # At a1 = a0 = 0.01 have the alphas as input variables
  # Have n, the input sequence, and the accepted hypothesis as output variables
  # type 1 err
  # type 2 err
  # h0 = hypothesis0
  # h1 = hypothesis1
  # a = probability of type I error
  # b = probability of type II error
  # k = # successes in n trials

  # LLR 
  a = type2.err / (1 - type1.err)
  b = (1 - type2.err) / type1.err
  
  n = length(rv)
  k = sum(rv, na.rm = TRUE)
  
  k.coeff = (log(h1 / (1-h1))) - (log(h0) / (1-h0))
  n.coeff = ((-1 * log(1-h1)) - (-1 * log(1-h0)) * -1)
  
  llr <- (k * k.coeff) / (n * n.coeff)

  decision <- NULL
  msg <- NULL
  if (llr < a) {
    decision <- TRUE
    msg <- "Accept h0"
  } else if (llr > b) {
    decision <- FALSE
    msg <- "Accept h1"
  } else {
    decision <- NA
    msg <- "Keep testing"
  }
  
  output <- list(n = n, 
                 k = k, 
                 h0 = h0,
                 h1 = h1,
                 wald.a = a,
                 wald.b = b,
                 k.coeff = k.coeff,
                 n.coeff = n.coeff,
                 llr = llr,
                 rv = rv,
                 decision = decision,
                 msg = msg)
  

}


simulation <- function() {
  # embed in a simulation that produces bernoulli random variables one at a time and tests to rejection (rbinom(1, 1, p))
  # run it in a sequence of x's Bernoulli 0.3 and 0.56
  # what would it do for 0.54? why does it give the result you got? 
  trial.3.true = 0
  trial.3.false = 0
  trial.3.na = 0
  
  trial.56.true = 0
  trial.56.false = 0
  trial.56.na = 0
  
  trial.54.true = 0
  trial.54.false = 0
  trial.54.na = 0
  
  # Bernoulli 0.3
  for (x in 1:1000) {
    rv <- rbinom(1, 1, 0.3)
    sprt_result <- SPRT(type1 = 0.01, type2 = 0.01, h0 = 0.45, h1 = 0.55, rv=rv)
    if (is.na(sprt_result$decision)) {
      trial.3.na  = trial.3.na + 1
    } else if (sprt_result == TRUE) {
      trial.3.true  = trial.3.true + 1
    } else {
      trial.3.false  = trial.3.false + 1
    }
  }
  cat("\n")
  cat("For Bernoulli p = 0.3")
  cat("\n")
  cat("Number of trials return FALSE: ", trial.3.false)
  cat("\n")
  cat("Number of trials return TRUE: ", trial.3.true)
  cat("\n")
  cat("Number of trials return NA:", trial.3.na)
  cat("\n")
  
  
  
  # Bernoulli 0.56 
  for (x in 1:1000) {
    rv <- rbinom(1, 1, 0.56)
    sprt_result <- SPRT(type1 = 0.01, type2 = 0.01, h0 = 0.45, h1 = 0.55, rv=rv)
    if (is.na(sprt_result$decision)) {
      trial.56.na  = trial.56.na + 1
    } else if (sprt_result == TRUE) {
      trial.56.true  = trial.56.true + 1
    } else {
      trial.56.false  = trial.56.false + 1
    }
  }
  
  cat("\n")
  cat("For Bernoulli p = 0.56")
  cat("\n")
  cat("Number of trials return FALSE: ", trial.56.false)
  cat("\n")
  cat("Number of trials return TRUE: ", trial.56.true)
  cat("\n")
  cat("Number of trials return NA:", trial.56.na)
  cat("\n")
  
  for (x in 1:1000) {
    rv <- rbinom(1, 1, 0.54)
    sprt_result <- SPRT(type1 = 0.01, type2 = 0.01, h0 = 0.45, h1 = 0.55, rv=rv)
    if (is.na(sprt_result$decision)) {
      trial.54.na  = trial.54.na + 1
    } else if (sprt_result == TRUE) {
      trial.54.true  = trial.54.true + 1
    } else {
      trial.54.false  = trial.54.false + 1
    }
  }
  
  cat("\n")
  cat("For Bernoulli p = 0.54")
  cat("\n")
  cat("Number of trials return FALSE: ", trial.54.false)
  cat("\n")
  cat("Number of trials return TRUE: ", trial.54.true)
  cat("\n")
  cat("Number of trials return NA:", trial.54.na)
  cat("\n")
  
}

simulation()
#for (x in 1:1000) {
#  rv <- rbinom(1, 1, 0.5)
#  sprt_result <- SPRT(type1 = 0.45, type2 = 0.55, h0 = 0.45, h1 = 0.55, rv=rv)
#  print(sprt_result$decision)
#}


