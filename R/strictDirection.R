#' @title Apply directional penalty to network
#' @description Applies lambda penalty to enforce edge directionality
#' @param X Network matrix
#' @param lambda Penalty parameter between 0 and 1
#' @return Modified network matrix
strictDirection <- function(X, lambda = 1){
  S <- as.matrix(X)
  
  # Optimized: only compute upper triangle
  n <- nrow(S)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(abs(S[i,j]) < abs(S[j,i])){
        S[i,j] <- 0
      } else {
        S[j,i] <- 0
      }
    }
  }
  
  O <- (((1-lambda) * X) + (lambda * S))
  O <- Matrix::Matrix(O, sparse = TRUE)
  return(O)
}
