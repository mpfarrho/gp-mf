# principal components estimates of the factor model X = F L' + e (with L'L/N = I)
get_pc_factors <- function(X,n_fac){
  N <- ncol(X)

  lambda <- svd(crossprod(X))$u[,1:n_fac,drop=FALSE]*sqrt(N) # loadings: sqrt(N) times the eigenvectors of X'X
  f_hat <- X%*%lambda/N                                      # factors
  mse <- mean((X - f_hat%*%t(lambda))^2)

  return(list(factors = f_hat, lambda = lambda, mse = mse))
}

# factors from data with missing values
get_em_factors <- function(data, n, it_max=50){
  X <- as.matrix(data)
  sl.na <- is.na(X)

  # standardize, extract factors and replace the missings with the common component
  em_step <- function(x){
    x0 <- scale(x)
    x_mu <- attr(x0, "scaled:center")
    x_sd <- attr(x0, "scaled:scale")
    attr(x0, "scaled:center") <- attr(x0, "scaled:scale") <- NULL
    x0[is.na(x0)] <- 0 # unconditional mean of standardized data is zero

    pc <- get_pc_factors(x0, n_fac=n)
    x0[sl.na] <- (pc$factors%*%t(pc$lambda))[sl.na]
    pc$data <- t(t(x0)*x_sd + x_mu)
    return(pc)
  }

  it <- 0
  err <- Inf
  con <- 1e-6 # convergence criterion
  pc <- em_step(X)
  while(it < it_max && err > con){
    f_old <- pc$factors
    pc <- em_step(pc$data)

    err <- abs(mean(apply(f_old^2, 2, mean) - apply(pc$factors^2, 2, mean)))
    it <- it + 1
  }
  pc$data[!sl.na] <- X[!sl.na]

  return(list(data=pc$data, factors=pc$factors, lambda=pc$lambda, iterations=it, mse=pc$mse))
}

# --------------------------------------------------------------------------------------
# number of factors based on the information criteria
get_num_factors <- function(x, kmax, criteria = "IC2"){
  Tn <- nrow(x)
  N <- ncol(x)
  ii <- 1:kmax

  # penalty
  CT <- switch(criteria,
               "IC1" = ii * (N+Tn)/(N*Tn) * log(N*Tn/(N+Tn)),
               "IC2" = ii * (N+Tn)/(N*Tn) * log(min(N,Tn)),
               "IC3" = ii * log(min(N,Tn))/min(N,Tn),
               stop("Choose available criterion ('IC1', 'IC2' or 'IC3')."))

  Fhat0 <- svd(tcrossprod(x))$u
  IC <- rep(NA, kmax+1)
  for(i in ii){
    Fhat <- Fhat0[,1:i,drop=FALSE]
    ehat <- x - Fhat %*% crossprod(Fhat,x)
    IC[i] <- log(mean(ehat^2)) + CT[i]
  }
  IC[kmax+1] <- log(mean(x^2)) # no factors

  ic <- which.min(IC)
  if(ic > kmax) ic <- 0

  return(list(ic = ic, IC = IC))
}
