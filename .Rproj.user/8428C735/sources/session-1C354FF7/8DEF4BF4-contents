scale_trimmed <- function(X, trim = 0.025){
  X <- as.matrix(X)
  for (col in 1:NCOL(X)) {
    trmean <- mean(X[ ,col], trim = trim )
    deviations <- (X[ ,col] - trmean)^2
    small_qq <- quantile(deviations, probs = trim, na.rm = TRUE)
    upper_qq <- quantile(deviations, probs = (1-trim), na.rm = TRUE)
    trdeviations <- deviations[which((deviations >= small_qq) & (deviations <= upper_qq) ) ]
    trsd <- sqrt(mean(trdeviations))
    if(trsd !=0){
      X[, col] <- (X[, col] - trmean) / trsd
    }
  }
  X
}
