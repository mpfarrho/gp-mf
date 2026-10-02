# ---------------------------------------------------------------------------------------------
# nonparametric midas estimation/forecast function
npmf <- function(y,X,Xth,prior_setup=NULL,fcst=TRUE,insamp=FALSE,
                    nburn=1000,nsave=1000,nsave_lim=1000,nthin=1){
  library(forecast)
  library(stochvol)
  library(dbarts)
  library(MASS)
  library(spam)
  library(fields)
  library(abind)
  
  # setup the model
  set.mean <- prior_setup$set.mean
  set.sv <- prior_setup$set.sv
  set.midas <- prior_setup$set.midas
  
  midas.par1 <- prior_setup$midas.par1
  midas.par2 <- prior_setup$midas.par2
  par1_pr <- prior_setup$par1_pr
  par2_pr <- prior_setup$par2_pr
  if(set.midas == "xalm"){ # Gaussian prior
    par1_a <- par1_pr # Prior mean 1
    par2_a <- par2_pr # Prior mean 2
    par1_b <- par2_b <- 0.1^2 # Prior variance
  }
  
  p.y <- prior_setup$p.y
  p.x <- prior_setup$p.x
  cons <- prior_setup$cons
  K_m <- prior_setup$K
  N <- NCOL(y)
  
  if(!(set.mean %in% c("lin","bart","gp"))) stop("Choose available conditional mean ('lin', 'bart' or 'gp').")
  
  # linear initialization
  shrink <- prior_setup$shrink
  
  # BART initialization
  cgm.level <- prior_setup$cgm.level
  cgm.exp <- prior_setup$cgm.exp
  sd.mu <- prior_setup$sd.mu
  num.trees <- prior_setup$num.trees
  
  # GP initialization
  samp.gp.hyper <- prior_setup$samp.gp.hyper
  max.try.solve <- prior_setup$max.try.solve
  xi <- xi_pr <- prior_setup$xi        # marginal variance
  h <- h_pr <- prior_setup$h           # length-scale
  zeta <- zeta_pr <- prior_setup$zeta  # smoothness
  zeta.grid <- prior_setup$zeta.grid   # discrete grid for smoothness
  gp.type   <- prior_setup$gp.type     # Matern versus exp^2
  
  # multiple equations
  midas.par1 <- rep(midas.par1, N)
  midas.par2 <- rep(midas.par2, N)
  
  # -----------------------------------------------------------------------------
  # sampling and prior settings
  ntot <- nburn+nthin*nsave
  thinset <- seq(nburn+1,ntot,by=nthin)
  thincount <- 0
  
  # hard limit for huge objects
  thinset_lim <- thinset[seq(1,length(thinset),by=round(nsave/nsave_lim,0))]
  if(length(thinset_lim)!=nsave_lim) stop("Check thinset limit.")
  thincount_lim <- 0
  
  # standardize data
  ymu <- apply(y,2,mean)
  ysd <- apply(y,2,sd)
  y <- apply(y,2,function(x) (x - mean(x))/sd(x))
  
  if(sd(X)==0){
    Xmu <- 0
    Xsd <- 1
  }else{
    Xmu <- apply(X, 2, mean)
    Xsd <- apply(X, 2, sd)
  }
  
  varlabs <- colnames(X)
  X <- X_og <- (X - t(matrix(Xmu, ncol(X), nrow(X))))/t(matrix(Xsd, ncol(X), nrow(X)))
  Xth <- Xth_og <- (Xth - matrix(Xmu,nrow=1))/matrix(Xsd,nrow=1)
  
  X_ls <- Xth_ls <- list()
  if(!(set.midas %in% c("u","f"))){
    set.midas.sub <- set.midas
    if(set.midas=="xalm"){
      w_midas <- get_xalm(3*(p.x),midas.par1[1],midas.par2[1])
      L <- 0
    } else if(set.midas=="br"){
      w_midas <- rep(1/(3*(p.x)),3*(p.x))
      L <- 0
    } else if(set.midas %in% c("avg","eop")) {
      L <- p.x - 1
      if(set.midas == "avg"){
        w_midas <- kronecker(diag(p.x),matrix(1/3, 3, 1))
      }else if(set.midas == "eop"){
        w_midas <- kronecker(diag(p.x),matrix(c(1,0,0),3,1))
      }
    } else{
      set.midas.sub <- strsplit(set.midas,"_")[[1]][1]
      L <- as.numeric(strsplit(set.midas,"_")[[1]][2])
      if(set.midas.sub == "alm"){
        w_midas <- get_poly(PH=3*p.x,L=L,type="almon")
      } else if(set.midas.sub == "leg"){
        w_midas <- get_poly(PH=3*p.x,L=L,alpha=0,beta=0,type="jacobi")
      } else if(set.midas.sub == "ber"){
        w_midas <- get_poly(PH=3*p.x,L=L,type="bernstein")
      } else if(set.midas.sub == "fou"){
        w_midas <- get_poly(PH=3*p.x,L=L,m=3,type="fourier")
      }
    }
    
    Xrow_id <- matrix(seq((N*p.y+cons+1),ncol(X)),ncol=K_m)
    Xcol_id <- matrix(seq(N*p.y+cons+1,N*p.y+cons+K_m*(L+1)),nrow=(L+1))
    
    WX <- matrix(0,N*p.y+cons+3*K_m*(p.x),N*p.y+cons+K_m*(L+1))
    if((N*p.y + cons) != 0) WX[1:(N*p.y+cons),1:(N*p.y+cons)] <- diag(N*p.y+cons)

    for(i in 1:K_m){
      WX[Xrow_id[,i],Xcol_id[,i]] <- w_midas
    }
    
    WX_p <- WX
    X <- X_og %*% WX
    Xth <- Xth_og %*% WX
    
    # store by equation
    for(nn in 1:N){
      X_ls[[nn]] <- X
      Xth_ls[[nn]] <- Xth
    }
    
    mh_c <- rep(0.01,N) # proposal variance
    midas_acc_count <- rep(0,N) # acceptance counter
  }else{
    set.midas.sub <- set.midas
    w_midas <- rep(0,3*(p.x))
    
    # store by equation
    for(nn in 1:N){
      X_ls[[nn]] <- X
      Xth_ls[[nn]] <- Xth
    }
  }
  
  # dimensions
  T <- NROW(y)
  K <- NCOL(X)
  
  # variance related objects
  a0_hom <- prior_setup$a0_hom
  b0_hom <- prior_setup$b0_hom
  sv_sig <- prior_setup$sv_sig
  
  Ht <- matrix(0,T,N)
  if(set.sv){
    sv_draw <- sv_latent <- sv_priors <- list()
    for(nn in 1:N){
      sv_draw[[nn]] <- list(mu = 0, phi = 0.99, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      sv_latent[[nn]] <- rep(0,T)
      
      sv_priors[[nn]] <- specify_priors(
        mu = sv_normal(mean=0,sd=1),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*sv_sig)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      ) 
    }
  }
  
  # required objects for sampling
  b_draw <- matrix(0,K,N)
  if(set.mean=="lin"){
    b0 <- matrix(0,K,N)
    
    # horseshoe prior
    theta_b <- lambda_b <- nu_b <- matrix(1,K,N)
    tau_b <- rep(1,N)
    zeta_b <- rep(1,N)
    I_T <- diag(T)
    
    if(K>T) fast.samp <- TRUE else fast.samp <- FALSE
  }else if(set.mean=="bart"){
    sig2.bart <- rep(1,N)
    sampler.list <- list()
    
    for(nn in 1:N){
      control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                               keepTrees = FALSE, n.samples = ntot,
                               n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                               n.threads = 1, n.thin = 1L, printEvery = 1,
                               printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                               updateState = FALSE)
      sampler.list[[nn]] <- dbarts(y[,nn]~X, control = control,
                             tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), 
                             n.samples = nsave, weights=rep(1,T), 
                             sigma=1, 
                             resid.prior = chisq(10^50, 0.5)) # fix variances in BART to 1 and use weights
    }
  }else if(set.mean=="gp"){
    # GP setup
    sc_GP <- sqrt(rep(0.01,N)) # proposal scaling variances
    acc_count <- rep(0,N)
    K_offset <- 1e-8 # add offsetting constant in case kernel is not invertible
    
    # compute hyperparameters
    xi <- rep(xi, N)
    h <- rep(h, N)
    zeta <- rep(zeta, N)
    
    xi_a <- rep(0.5,N)
    xi_b <- xi_a/xi_pr
    # h_pr <- h_pr*arima(y,order = c(p.y,0,0))$sigma2 # ballpark scaling
    
    h_prtmp <- rep(NA,N)
    for(nn in 1:N){
      h_prtmp[nn] <- auto.arima(y[,nn])$sigma2 # ballpark scaling
    }
    h_pr <- h_pr*h_prtmp
    h_a <- rep(0.5, N)
    
    h_b <- h_a/h_pr
    zeta_a <- rep(min(zeta.grid),N)
    zeta_b <- rep(max(zeta.grid),N)
  
    # construct kernels and associated objects
    nlog2pi <- (-T/2)*log(2*pi) # normalizing constant for likelihood
    
    # distance measures to compute kernel
    Xdist <- rdist(X)
    Xthdist <- rdist(Xth)
    XXthdist <- rdist(X,Xth)
    
    # compute kernel
    KXX_ls <- KXXsig_ls <- KXXi_ls <- Xdist_ls <- XXthdist_ls <- list()
    KXXdet_ls <- rep(1,N)
    lik <- pri <- rep(NA, N)
    for(nn in 1:N){
      KXX <- get_Kernel(r = Xdist, xi = xi[nn], h = h[nn], zeta = zeta[nn], type = gp.type)
      KXXsig <- KXX + diag(T)*as.numeric(exp(Ht[,nn]))
      KXXi <- solve(KXXsig)
      KXXdet <- determinant(KXXsig)$modulus
      
      KXX_ls[[nn]] <- KXX
      KXXsig_ls[[nn]] <- KXXsig
      KXXi_ls[[nn]] <- KXXi
      KXXdet_ls[nn] <- KXXdet
      
      Xdist_ls[[nn]] <- Xdist
      XXthdist_ls[[nn]] <- XXthdist
      
      lik[nn] <- as.numeric(nlog2pi - KXXdet/2 - (t(y[,nn,drop=FALSE]) %*% KXXi %*% y[,nn,drop=FALSE])/2)
      pri[nn] <- dgamma(h[nn],h_a[nn],h_b[nn],log=TRUE) + dgamma(xi[nn],xi_a[nn],xi_b[nn],log=TRUE) # + dunif(zeta,zeta_a,zeta_b,log=TRUE) 
    }
  }
  
  # conditional mean functions and structural form parameters
  eps <- fx <- matrix(0, T, N)
  if(set.mean=="bart"){
    for(nn in 1:N) fx[, nn] <- sampler.list[[nn]]$predict(X_ls[[nn]]) # initial fit of the sampler
  }
  B0tilde <- matrix(0, N, N)
  sig2.reg <- rep(1,N)
  
  # -----------------------------------------------------------------------------
  # storage
  fx_store        <- array(NA,dim=c(nsave,T,N))
  sig2_store      <- array(NA,dim=c(nsave,T,N))
  
  w_store         <- array(NA,dim=c(nsave,3*p.x,N))
  bpoly_store     <- array(NA,dim=c(nsave,3*p.x,N))
  b_store         <- array(NA,dim=c(nsave_lim,K,N))
  X_store         <- array(NA,dim=c(nsave_lim,T,K,N))
  Xth_store       <- array(NA,dim=c(nsave_lim,K,N))
  
  fcst_store      <- array(NA,dim=c(nsave,N))

  dimnames(fx_store) <- dimnames(sig2_store) <- list(NULL,rownames(y),colnames(y))

  # estimation
  t.start <- Sys.time()
  pb <- txtProgressBar(min = 0, max = ntot, style = 3) # start progress bar
  
  irep <- 1
  for(irep in 1:ntot){
    # ---------------------------------------------------------------------------------------------------------------
    # (0a) sample contemporaneous correlation and uncorrelated endogenous variable
    normalizer <- exp(-Ht/2)
    if(N > 1){
      for(nn in 2:N){
        yy <- (y[,nn] - fx[,nn]) * normalizer[,nn]
        xx <- (-y[,1:(nn-1),drop = FALSE]) * normalizer[,nn]
        
        VB0p <- solve(diag(nn-1) + crossprod(xx))
        bB0p <- VB0p %*% crossprod(xx,yy)
        
        B0tilde[nn,1:(nn-1)] <- bB0p + t(chol(VB0p)) %*% rnorm(nn-1)
      }
      ytilde <- y + y %*% t(B0tilde)
      B0 <- diag(N) + B0tilde
      B0i <- solve(B0)
    }else{
      ytilde <- y
      B0 <- B0i <- diag(1)
    }
    
    # (0b) Sample the "data" in case of applicable MIDAS variant
    if(set.midas.sub == "xalm"){
      for(nn in 1:N){
        par1_prop <- truncnorm::rtruncnorm(1,a=0,b=Inf,mean=midas.par1[nn],sd=mh_c[nn])
        par2_prop <- truncnorm::rtruncnorm(1,a=-Inf,b=0,mean=midas.par2[nn],sd=mh_c[nn])
        # par1_prop <- rnorm(1,mean=midas.par1[nn],sd=mh_c[nn])
        # par2_prop <- rnorm(1,mean=midas.par2[nn],sd=mh_c[nn])

        pri_midas <-   dnorm(midas.par1[nn],par1_a,sqrt(par1_b),log=TRUE) + dnorm(midas.par2[nn],par2_a,sqrt(par2_b),log=TRUE)
        pri_midas_p <- dnorm(par1_prop ,par1_a, sqrt(par1_b),log=TRUE) + dnorm(par2_prop ,par2_a, sqrt(par2_b),log=TRUE)
        
        w_prop <- get_xalm(3*(p.x),par1_prop,par2_prop)
        for(i in 1:K_m){
          WX_p[Xrow_id[,i],Xcol_id[,i]] <- w_prop
        }
        X_p <- X_og %*% WX_p
        
        if(set.mean=="lin"){
          lik_midas <- sum(dnorm(ytilde[, nn], X_ls[[nn]] %*% b_draw[,nn], exp(Ht[,nn] / 2),log=TRUE))
          lik_midas_p <- sum(dnorm(ytilde[, nn], X_p %*% b_draw[,nn], exp(Ht[,nn] / 2),log=TRUE))
        }else if(set.mean=="bart"){
          lik_midas <- sum(dnorm(ytilde[, nn], fx[, nn], exp(Ht[ ,nn]/2), log=TRUE))
          fx_p <- sampler.list[[nn]]$predict(X_p)
          lik_midas_p <- sum(dnorm(ytilde[, nn], fx_p, exp(Ht[, nn]/2), log=TRUE))
        }else if(set.mean=="gp"){
          Xdist_p <- rdist(X_p)
          KXX_midas_p <- get_Kernel(r = Xdist_p, xi = xi[nn], h = h[nn], zeta = zeta[nn], type = gp.type)
          KXXsig_midas_p <- KXX_midas_p + diag(T) * as.numeric(exp(Ht[, nn]))
          KXXi_midas_p <- solve(KXXsig_midas_p)
          KXXdet_midas_p <- determinant(KXXsig_midas_p)$modulus
          
          lik_midas <- as.numeric(nlog2pi - KXXdet_ls[nn]/2 - (t(ytilde[, nn]) %*% KXXi_ls[[nn]] %*% ytilde[, nn])/2)
          lik_midas_p <- as.numeric(nlog2pi - KXXdet_midas_p/2 - (t(ytilde[, nn]) %*% KXXi_midas_p %*% ytilde[, nn])/2)
        }
        
        # correction for asymmetric (truncated normal) proposals
        prop_corr <- pnorm(midas.par1[nn]/mh_c[nn],log.p=TRUE) - pnorm(par1_prop/mh_c[nn],log.p=TRUE) +
          pnorm(-midas.par2[nn]/mh_c[nn],log.p=TRUE) - pnorm(-par2_prop/mh_c[nn],log.p=TRUE)

        likrat_midas <- (lik_midas_p + pri_midas_p) - (lik_midas + pri_midas) + prop_corr
        likrat_midas <- ifelse(is.nan(likrat_midas),-Inf,likrat_midas)

        # accept/reject (BART: reject if the new predictors would leave a tree with an empty leaf)
        acc_midas <- likrat_midas > log(runif(1,0,1))
        if(acc_midas && set.mean=="bart"){
          acc_midas <- isTRUE(sampler.list[[nn]]$setPredictor(X_p, forceUpdate = FALSE))
        }
        if(acc_midas){
          midas_acc_count[nn] <- midas_acc_count[nn] + 1
          midas.par1[nn] <- par1_prop
          midas.par2[nn] <- par2_prop
          
          X <- X_ls[[nn]] <- X_p
          Xth <- Xth_ls[[nn]] <- Xth_og %*% WX_p
          
          # re-configure nonparametric explanatory variables if accepted (BART predictors already set above)
          if(set.mean=="gp"){
            KXX_ls[[nn]] <- KXX_midas_p
            KXXsig_ls[[nn]] <- KXXsig_midas_p
            KXXi_ls[[nn]] <- KXXi_midas_p
            KXXdet_ls[nn] <- KXXdet_midas_p
            
            # update distances for forecasts
            Xdist_ls[[nn]] <- Xdist_p
            XXthdist_ls[[nn]] <- rdist(X,Xth)
          }
        }
        
        acc_prop_midas <- midas_acc_count[nn]/irep
        if(irep < (2/3)*nburn){
          if(acc_prop_midas>0.3){
            mh_c[nn] <- mh_c[nn]*1.01
          }else if(acc_prop_midas<0.2){
            mh_c[nn] <- mh_c[nn]*0.99
          }
        }
      }
    }
    
    # (0c) sample GP hyperparameters if required
    if(set.mean=="gp" & samp.gp.hyper){
      for(nn in 1:N){
        lik <- as.numeric(nlog2pi - KXXdet_ls[nn]/2 - (t(ytilde[,nn]) %*% KXXi_ls[[nn]] %*% ytilde[,nn])/2)
        pri <- dgamma(h[nn],h_a[nn],h_b[nn],log=TRUE) + dgamma(xi[nn],xi_a[nn],xi_b[nn],log=TRUE) # + dunif(zeta,zeta_a,zeta_b,log=TRUE)  
        
        h_p <- exp(rnorm(1,0,sc_GP[nn]))*h[nn]
        xi_p <- exp(rnorm(1,0,sc_GP[nn]))*xi[nn]
        
        # proposal
        if(gp.type == "sqexp") zeta_p <- Inf else if(gp.type == "matern") zeta_p <- sample(zeta.grid,1)
        KXX_p <- get_Kernel(r = Xdist_ls[[nn]], xi = xi_p, h = h_p, zeta = zeta_p, type = gp.type)
        
        # compute required objects for proposal
        KXXsig_p <- KXX_p + diag(T)*as.numeric(exp(Ht[, nn]))
        KXXi_p <- try(solve(KXXsig_p),silent=TRUE)
        
        check.inv <- 0
        while(is(KXXi_p,"try-error") & check.inv < max.try.solve){
          check.inv <- check.inv + 1
          KXXsig_p <- KXXsig_p + K_offset*diag(T)
          KXXi_p <- try(solve(KXXsig_p),silent=TRUE)
          message("Adding constant: ",K_offset*check.inv)
        }
        KXXdet_p <- determinant(KXXsig_p)$modulus
        
        lik_p <- as.numeric(nlog2pi - KXXdet_p/2 - (t(ytilde[,nn]) %*% KXXi_p %*% ytilde[,nn])/2)
        pri_p <- dgamma(h_p,h_a[nn],h_b[nn],log=TRUE) + dgamma(xi_p,xi_a[nn],xi_b[nn],log=TRUE) #+ dunif(zeta_p,zeta_a,zeta_b,log=TRUE)
        
        # compute the acceptance probability
        likrat <- (lik_p + pri_p) - (lik + pri) + 
          log(h_p) - log(h[nn]) + 
          log(xi_p) - log(xi[nn])
        
        likrat <- ifelse(is.nan(likrat),-Inf,likrat)
        if(likrat > log(runif(1,0,1))){
          acc_count[nn] <- acc_count[nn] + 1
          
          h[nn] <- h_p
          xi[nn] <- xi_p
          zeta[nn] <- zeta_p
          
          # update lists (distances stay the same as X is unchanged)
          KXX_ls[[nn]] <- KXX_p
          KXXsig_ls[[nn]] <- KXXsig_p
          KXXi_ls[[nn]] <- KXXi_p
          KXXdet_ls[nn] <- KXXdet_p
        }  
        
        # tune proposal variance for RW-MH
        acc_prob <- acc_count[nn]/irep
        if(irep %% 500 == 0) message(" ",round(100*acc_prob,digits=1),"% (",acc_count[nn],"/",irep,")")
        if(irep < nburn){
          if(acc_prob>0.4){
            sc_GP[nn] <- sc_GP[nn]*1.01
          }else if(acc_prob<0.1){
            sc_GP[nn] <- sc_GP[nn]*0.99
          }
        }
        sc_GP[sc_GP>sqrt(1)] <- sqrt(1)
      }
    }
    
    # ---------------------------------------------------------------------------------------------------------------
    # (1) sample conditional mean
    normalizer <- exp(-Ht/2)
    if(set.mean=="lin"){
      for(nn in 1:N){
        X <- X_ls[[nn]] # set predictors to updated ones in case applicable
        XX <- X * normalizer[,nn]
        yy <- ytilde[, nn] * normalizer[,nn]
        
        if(fast.samp){
          uu_fs <- rnorm(K,0,sqrt(theta_b[ , nn]))
          dd_fs <- rnorm(T)
          vv_fs <- XX %*% uu_fs + dd_fs
          
          XXV0 <- theta_b[, nn] * t(XX)
          XXi <- XX%*%XXV0 + I_T
          
          bb_fs <- (yy - vv_fs)
          ww_fs <- solve(XXi) %*% bb_fs
          b_draw[, nn] <- uu_fs + XXV0 %*% ww_fs
        }else{
          V0i <- diag(K) / theta_b[, nn]
          Vp <- solve(V0i + crossprod(XX))
          bp <- Vp %*% (V0i %*% b0[, nn] + crossprod(XX,yy))
          b_draw[, nn] <- bp + t(chol(Vp)) %*% rnorm(K)
        }
        
        fx[, nn] <- X %*% b_draw[, nn]
        eps[,nn] <- ytilde[, nn] - fx[, nn]
        
        # (1b) shrinkage prior
        if(shrink){
          hs_draw <- get.hs(bdraw=as.numeric(b_draw[,nn]),
                            lambda.hs=as.numeric(lambda_b[,nn]),
                            nu.hs=nu_b[,nn],tau.hs=tau_b[nn],zeta.hs=zeta_b[nn])
          theta_b[, nn] <- hs_draw$psi
          lambda_b[, nn] <- hs_draw$lambda
          nu_b[,nn] <- hs_draw$nu
          tau_b[nn] <- hs_draw$tau
          zeta_b[nn] <- hs_draw$zeta
        } 
      }
      theta_b[theta_b<1e-8] <- 1e-8
    }else if(set.mean=="bart"){
      for(nn in 1:N){
        sampler.list[[nn]]$setWeights(exp(-Ht[, nn]))
        # response stays y[, nn] and ytilde enters via the offset, since setResponse rescales the prior (dbarts issue #80)
        off_nn <- y[, nn] - ytilde[, nn]
        sampler.list[[nn]]$setOffset(off_nn, updateScale = FALSE)
        # predictors already set when updating MIDAS parameters if applicable
        bart.rep <- sampler.list[[nn]]$run(0L, 1L)

        sig2.bart[nn] <- bart.rep$sigma
        fx[, nn] <- bart.rep$train - off_nn
        eps[, nn] <- ytilde[, nn] - fx[, nn]
      }
    }else if(set.mean=="gp"){
      for(nn in 1:N){
        # predictors for kernels already set when updating MIDAS parameters if applicable
        KXXKXXi <- KXX_ls[[nn]] %*% KXXi_ls[[nn]]
        F_mu <- KXXKXXi %*% ytilde[, nn]
        F_sig <- KXX_ls[[nn]] - KXXKXXi %*% KXX_ls[[nn]]
        fx_nn <- try(F_mu + t(chol(F_sig)) %*% rnorm(T), silent=TRUE)
        if(is(fx_nn,"try-error")){
          F_sig <- F_sig + diag(1e-3, T)
          fx_nn <- try(F_mu + t(chol(F_sig)) %*% rnorm(T), silent=TRUE)
          if(is(fx_nn,"try-error")){
            fx_nn <- F_mu
            message("No draw produced. Setting to posterior mean of GP.")
          }
        }
        fx[, nn] <- fx_nn
        eps[, nn] <- ytilde[, nn] - fx[, nn]
      }
    }
    
    # (2) sample conditional variance
    for(nn in 1:N){
      if(set.sv){
        sv_sample_nn <- svsample_fast_cpp(eps[,nn], startpara=sv_draw[[nn]],startlatent=sv_latent[[nn]],priorspec=sv_priors[[nn]])
        sv_draw[[nn]][c("mu", "phi", "sigma")] <- as.list(sv_sample_nn$para[, c("mu", "phi", "sigma")])
        Ht[, nn] <- sv_latent[[nn]] <- sv_sample_nn$latent
      }else{
        sig2.reg[nn] <- -log(rgamma(1, a0_hom + T/2, b0_hom + crossprod(eps[, nn])/2))
        Ht[, nn] <- sig2.reg[nn]
      }  
    }
    
    # update GP objects where the variance shows up
    if(set.mean=="gp"){
      for(nn in 1:N){
        KXXsig_ls[[nn]] <- KXX_ls[[nn]] + as.numeric(exp(Ht[, nn]))*diag(T)
        KXXi_ls[[nn]] <- solve(KXXsig_ls[[nn]]) 
        KXXdet_ls[nn] <- determinant(KXXsig_ls[[nn]])$modulus
      }
    }
    
    # ---------------------------------------------------------------------------------------------------------------
    # (3) Storage and predictions
    if(irep %in% thinset){
      thincount <- thincount+1
      fx_store[thincount,,] <- t(((B0i %*% t(fx)) * ysd) + matrix(ymu, ncol = T, nrow = N))
      sig2_store[thincount,,] <- t((ysd^2) * t(exp(Ht)))

      # reduced number of draws
      if(irep %in% thinset_lim){
        thincount_lim <- thincount_lim+1
        b_store[thincount_lim,,] <- b_draw
        
        for(nn in 1:N){
          X_store[thincount_lim,,,nn] <- X_ls[[nn]]
          Xth_store[thincount_lim,,nn] <- Xth_ls[[nn]] 
        }
      }
      
      # forecasts
      if(fcst){
        hth <- Ht[T,] # predict variances
        if(set.sv){
          for(nn in 1:N){
            hth[nn] <- sv_draw[[nn]]$mu + sv_draw[[nn]]$phi * (hth[nn] - sv_draw[[nn]]$mu) + sv_draw[[nn]]$sigma * rnorm(1)
          }
        }
        eth <- rnorm(N,0,exp(hth/2))
        
        # predict mean
        fxth <- rep(NA, N)
        for(nn in 1:N){
          if(set.mean=="lin"){
            fxth[nn] <- as.numeric(Xth_ls[[nn]] %*% b_draw[,nn])
          }else if(set.mean=="bart"){
            fxth[nn] <- sampler.list[[nn]]$predict(Xth_ls[[nn]])
          }else if(set.mean=="gp"){
            KXthXth <- xi[nn]
            KXXth <- get_Kernel(r = XXthdist_ls[[nn]], xi = xi[nn], h = h[nn], zeta = zeta[nn], type = gp.type)
            KXXthKXXi <- crossprod(KXXth,KXXi_ls[[nn]])
            
            yth_mu <- as.numeric(KXXthKXXi %*% ytilde[, nn])
            yth_sig <- sqrt(as.numeric(KXthXth - KXXthKXXi%*%KXXth))
            fxth[nn] <- rnorm(1, yth_mu, yth_sig)
          }  
        }
        
        # predict variables
        yth <- rep(NA,N)
        yth[1] <- fxth[1] + eth[1]
        if(N > 1){
          for(nn in 2:N){
            yth[nn] <- fxth[nn] + sum(B0tilde[nn,1:(nn-1)] * (-yth[1:(nn-1)])) + eth[nn]  
          }
        }
        # yth <- as.numeric(B0i %*% (fxth + eth)) # (same as using mapping between structural and reduced form when triangular)
        fcst_store[thincount,] <- (yth * ysd) + ymu
      }
    }
    setTxtProgressBar(pb, irep)
  }
  
  # -----------------------------------------------------------------------------------------------------
  # change dimensions with multiple outputs here
  # some other posteriors
  fx.smry <- abind(apply(fx_store, c(2,3), quantile, c(0.05, 0.5, 0.95),na.rm=TRUE), 
                     "y" = t((t(y)*ysd) + matrix(ymu, ncol = T, nrow = N)), along = 1)
  sig2.smry   <- apply(sig2_store, c(2,3), quantile, c(0.05, 0.5, 0.95),na.rm=TRUE)

  # posterior for (approximate) regression coefficients
  b_store <- apply(b_store,c(2,3),quantile,probs=c(0.16,0.5,0.84),na.rm=TRUE)
  X_store <- apply(X_store,c(2,3,4),median,na.rm=TRUE)
  Xth_store <- apply(Xth_store,c(2,3),median,na.rm=TRUE)
  dimnames(b_store) <- list(c("p16","median","p84"), NULL, colnames(y))
  
  # output object
  t.elapsed <- difftime(Sys.time(),t.start,units="secs")
  message("\nEstimation time: ",round(t.elapsed/60,digits=2)," mins.")
  ret_obj <- list("fx"=fx.smry,"sig2"=sig2.smry,
                  "b_post"=b_store,"X_post"=X_store, "Xth_post" = Xth_store,
                  "fcst"=fcst_store,
                  "setup"=list("priors"=prior_setup,"variables"=varlabs)
                  )
  return(ret_obj)
}

