# ---------------------------------------------------------------------------------------------
# horseshoe prior
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# ---------------------------------------------------------------------------------------------
# weighting function for exponential almon (xalm) with 2 parameters
get_xalm <- function(j,par1,par2){
  jj <- seq(1,j)
  w_tmp <- exp(par1*jj + par2*jj^2)+1e-6
  return(w_tmp/sum(w_tmp))
}

# several useful polynomials
get_poly <- function(PH,L,alpha=NULL,beta=NULL,m=3,type){
  require(pracma)
  
  xseq <- seq(0,1,length.out=PH)
  W <- matrix(1,nrow=PH,ncol=L+1) # holds polynomials
  colnames(W) <- paste0("deg",0:L)
  rownames(W) <- paste0("lag",1:PH)
  
  if(type == "almon"){
    P <- seq(0,1,length.out=PH)
    for(i in 0:L){
      W[,i+1] <- P^i
    }
  } else if(type == "jacobi"){
    P <- matrix(1,nrow=PH,ncol=L+2)
    P[, 2] <-  2*xseq - 1
    
    for (i in 1:L){
      d0 <- (2*i + 2*alpha + 1)*(2*i + 2*alpha + 2)/(2*(i+1)*(i + 2*alpha + 1))
      c0 <- (alpha + i)^2*(2*i + 2*alpha + 2)/( (i+1)*(i + 2*alpha + 1)*(2*i + 2*alpha)) 
      P[, i+2]   <- d0 * P[, 2]*P[, i+1] - c0 * P[, i]
      W[, i+1] <- sqrt((2*i + 1)) %*% P[, i+1]
    }
  } else if(type == "bernstein"){
    for (i in 0:L) {
      for(j in 1:PH){
        W[j,i+1] <- bernsteinb(i, L, xseq[j])
      }
    }
  } else if(type == "fourier"){
    omega <- 2 * pi / (L * m)
    for (i in 2:(L+1)) {
      if (i %% 2 == 0) {
        W[, i] <- sin(i * omega * xseq)
      } else {
        W[, i] <- cos(i * omega * xseq)
      }
    }
  }
  return(W)
}

# ---------------------------------------------------------------------------------------------
# compute various kernels
get_Kernel <- function(r,xi,h,zeta,type="matern"){
  # h = 1/l^2 (inverse length-scale)
  d <- r^2*h
  if(type=="sqexp"){
    KX <- xi * exp(-(d/2))
  }else if(type=="matern"){
    if(zeta == 1) zeta <- 1-1e-8 # rule out special case
    
    if(zeta == Inf || zeta >= 20){
      KX <- xi * exp(-(d/2))
    }else if(zeta==1/2){
      KX <- xi * exp(-sqrt(2*zeta*d))
    }else if(zeta==3/2){
      KX <- xi * (1 + sqrt(3*d)) * exp(-sqrt(3*d))
    }else if(zeta==5/2){
      KX <- xi * (1 + sqrt(5*d) + (5/3)*d) * exp(-sqrt(5*d))
    }else{
      KX <- xi * (2^(1-zeta))/gamma(zeta) * (sqrt(2*zeta*d)^zeta) * besselK(sqrt(2*zeta*d),zeta)
      KX[is.nan(KX)] <- xi
    }
  }
  return(KX)
}

# construct lags
get.lag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- NCOL(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 2 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# ---------------------------------------------------------------------------------------------
# function for quantile weighted CRPS
qwCRPS <-function(true,Qtau,tau,weighting="none"){
  require(pracma)
  
  tau_len <- length(tau)
  true_rep <- rep(true,tau_len)
  QS.vec <- 2*(true_rep-Qtau)*(tau-((true_rep<=Qtau)*1))
  
  weights <- switch(tolower(weighting),
                    "none" = 1,
                    "tails" = (2*tau-1)^2,
                    "right" = tau^2,
                    "left" = (1-tau)^2,
                    "center" = tau*(1-tau))
  wghs <- QS.vec*weights
  return(pracma::trapz(tau,wghs))
}
