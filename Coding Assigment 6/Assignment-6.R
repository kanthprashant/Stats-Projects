library(ggplot2)
FDR_calculator <- function(pvals, q, independent=FALSE) {
  # order pvals in ascending order
  ordered.pvals <- sort(pvals, decreasing = FALSE)
  bh.vals <- vector()
  m <- length(ordered.pvals)
  
  # for each pvalue, generate bh critical value (i/m)*q 
  for (i in 1:length(ordered.pvals)) {
    bh <- (ordered.pvals[i] / m) * q
    bh.vals <- append(bh.vals, bh)
  }
  
  # find highest p value such that still smaller than critical value 
  highest.p <- 0
  for (i in 1:length(ordered.pvals)) {
    if (ordered.pvals[i] < bh.vals[i]) {
      highest.p = i
    }
  }

  # create significant values, all values lower than highest p-value smaller than critical value are significant
  sig.col <- vector()
  sig.vals <- vector()
  for (i in 1:length(ordered.pvals)) {
    if (i <= highest.p) {
      sig.col <- append(sig.col,1)
      sig.vals <- append(sig.vals, ordered.pvals[i])
    } else {
      sig.col <- append(  sig.col,0)
    }
  }
  
  dataframe <- data.frame(ordered.pvals,   sig.col)
  dataframe['rank'] <- 1:length(ordered.pvals)
  
  fdr <- highest.p/m
  cat("FDR: ", fdr)
  cat("\n")
  cat("Significant Values: ")
  print(sig.vals)
  # blue points are significant values
  ggplot(data=dataframe, aes(x=rank, y=ordered.pvals)) + geom_point(aes(color=sig.col))
}

pvals <- c(.205, .074, .060, .042, .041, .039, .041, .008, .001, .005, .007, .06, .56, .734, 
           .58, .5345, .06, .004, .235, .34)

FDR_calculator(pvals, 20)