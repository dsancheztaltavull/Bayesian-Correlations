#Bayesian Correlation algorithm


#More information:
#Statistical properties: https://doi.org/10.1371/journal.pone.0163595
#To study scRNA-seq data: https://doi.org/10.1101/714824


BaCo <- function(X){
  
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}



#Example

#create a matrix (or load your own)
X <- matrix(1:1000, ncol=20)

#compute the Bayesian correlation matrix
B <- BaCo(X)

#look at the first 5 columns and 5 rows
B[1:5,1:5]

