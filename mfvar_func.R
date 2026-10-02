# ---------------------------------------------------------------------------------------------
# mixed-frequency VAR for monthly and quarterly variables
#
# requirements on Yraw (months x n matrix with unique column names, quarterly variables in the
# first n_lf columns); there is no forecast horizon argument, the nowcast/forecast periods are
# determined by the rows of Yraw:
#  - the first row is the first month of a quarter and the number of rows is a multiple of 3
#  - quarterly variables are observed in the third month of a quarter (NA otherwise), from the
#    first quarter of the sample onwards
#  - missing values are only allowed at the end of the sample (ragged edge)
#  - for nowcasts/forecasts, append rows of NAs up to the final month of the target quarter
#  - p >= 4 (the triangular aggregation scheme spans five months)
#
# returns draws of the quarterly variables (triangular aggregates of the latent monthly states,
# on the original scale) for each month after the first p, as nsave x months x n_lf array
mfvar <- function(Yraw, n_lf, p = 5, cons = FALSE, sv = FALSE,
                   nburn = 3000, nsave = 3000, nthin = 3){
  require(Matrix)
  require(stochvol)

  get.hs.min <- function(bdraw,lam.hs,nu.hs,tau.hs,zeta.hs,update.ls = TRUE){
    k <- length(bdraw)
    # Local shrinkage scalings
    if(update.ls){
      lam.hs <- 1/rgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
      nu.hs <-  1/rgamma(k,shape=1,rate=1+1/lam.hs)
    }else{
      lam.hs <- lam.hs
      nu.hs <- nu.hs
    }
    # Global shrinkage parameter
    tau.hs  <- 1/rgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lam.hs)/2) 
    zeta.hs <- 1/rgamma(1,shape=1,rate=1+1/tau.hs)
    
    ret <- list("psi"=(lam.hs*tau.hs),"lam"=lam.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
    return(ret)
  }
  
  n <- ncol(Yraw)
  n_hf <- n - n_lf
  varlabs <- colnames(Yraw)
  itr <- rep("D", n_lf) # triangular aggregation scheme
  
  datelabs <- rownames(Yraw)
  datelabs_p <- datelabs[-c(1:p)]
  
  Ymu <- apply(Yraw, 2, mean, na.rm = TRUE)
  Ysd <- apply(Yraw, 2, sd, na.rm = TRUE)
  Yraw <- scale(Yraw, center = Ymu, scale = Ysd)
  
  Yraw_NA <- Yraw # original data with missings
  Ytmp <- Yraw[seq(3, nrow(Yraw), by = 3), 1:n_lf, drop = FALSE]
  for (nn in 1:n_lf) {
    Yraw[, nn] <- rep(Ytmp[, nn], each = 3)
  }
  
  # define matrix of lags
  if (any(is.na(Yraw))) {
    for (i in 1:ncol(Yraw)) {
      yy <- Yraw[, i]
      if(any(is.na(yy))){
        yylags <- embed(yy[!is.na(yy)],2)
        y1 <- yylags[,1]
        x1 <- cbind(yylags[,-1],1)
        phic <- solve(crossprod(x1)) %*% crossprod(x1,y1)
        if(phic[1] >= 1) phic[1] <- 0.999
        
        yh <- y1[length(y1)]
        yfc_init <- rep(NA,sum(is.na(yy)))
        for(hh in 1:(sum(is.na(yy)))){
          yfc_init[hh] <- yh <- c(yh,1) %*% phic
        }
        
        Yraw[, i] <- c(yy[!is.na(yy)], yfc_init)
      }
    }
  }
  
  # lag structure
  Y_init <- Yraw[1:p, ]
  Ylags <- embed(Yraw, dim = p + 1)
  rownames(Y_init) <- paste0("T", (-p + 1):0)
  rownames(Ylags) <- paste0("T", 1:nrow(Ylags))
  colnames(Ylags) <- paste0(rep(varlabs, p + 1), "_l", rep(0:p, each = n))
  Y <- Ylags[, 1:n]
  
  Y_full <- rbind(Y_init, Y)
  colnames(Y) <- colnames(Yraw)
  
  # match observed low-frequency data to lag structure
  Ytmp <- Yraw_NA[(p + 1):nrow(Yraw_NA), 1:n_lf, drop = FALSE]
  Ytmp_init <- Yraw_NA[1:p, 1:n_lf, drop = FALSE]
  rownames(Ytmp) <- paste0("T", 1:nrow(Ytmp))
  rownames(Ytmp_init) <- paste0("T", (-p + 1):0)
  
  min_lf_obs <- min(which(apply(!is.na(Ytmp), 1, sum) == n_lf))
  min_lf_init <- min(which(apply(!is.na(Ytmp_init), 1, sum) == n_lf))
  
  Z <- Ytmp[seq(min_lf_obs, nrow(Ytmp), by = 3), , drop = FALSE]
  Z_init <- Ytmp_init[seq(min_lf_init, nrow(Ytmp_init), by = 3), , drop = FALSE]
  Z_full <- rbind(Z_init, Z)
  
  ix_t_Y <- as.numeric(gsub("T", "", rownames(Y_full)))
  ix_t_Z <- as.numeric(gsub("T", "", rownames(Z_full)))
  ix_min_Z <- which(ix_t_Z - 5 == min(ix_t_Y))
  
  Z_full <- Z_full[ix_min_Z:nrow(Z_full), , drop = FALSE]
  Z_vec <- c(t(Z_full))
  sl_na_lf <- is.na(Z_vec)
  Z_vec <- Z_vec[!sl_na_lf]
  
  T_obs_lf <- length(Z_vec)
  X <- Ylags[, -c(1:n)]
  if (cons) X <- cbind(X, 1) # intercept in the last column

  # --------------------------------------------------------------------------
  # setup related to data and dimensions
  k <- ncol(X)
  T <- nrow(Y)

  eps_small <- sqrt(.Machine$double.eps)
  
  # some required indexing for missings
  Yraw_NA_full <- Yraw_NA
  Yraw_NA_full[, 1:n_lf] <- NA
  
  id_lat <- is.na(Yraw_NA_full)
  id_lat_vec <- c(t(id_lat))
  
  Tn_lat <- sum(id_lat)

  Y_full_vec <- c(t(Y_full))
  Yo_vec <- Y_full_vec[!id_lat_vec]
  
  # setup for selection matrices
  Slt_ls <- Smt_ls <- list()
  for (tt in 1:(T + p)) {
    n_lat_t <- sum(id_lat[tt, ])
    n_mea_t <- n - n_lat_t
    
    Slt_tmp <- matrix(0, nrow = n, ncol = n_lat_t)
    Smt_tmp <- matrix(0, nrow = n, ncol = n_mea_t)
    lcount <- mcount <- 0
    for (nn in 1:n) {
      if (id_lat[tt, nn] == 1) {
        lcount <- lcount + 1
        Slt_tmp[nn, lcount] <- 1
      } else {
        mcount <- mcount + 1
        Smt_tmp[nn, mcount] <- 1
      }
    }
    Slt_ls[[tt]] <- Slt_tmp
    Smt_ls[[tt]] <- Smt_tmp
  }
  
  Sl <- bdiag(Slt_ls)
  Sm <- bdiag(Smt_ls)

  # labeling of positions
  Ypos <- Y_full
  for (i in 1:n) {
    Ypos[, i] <- paste0(colnames(Y_full)[i], "_", rownames(Y_full))
  }
  Ypos_vec <- c(t(Ypos))

  # variable and period of each low-frequency observation
  Zvar_vec <- rep(1:n_lf, nrow(Z_full))[!sl_na_lf]
  Zt_vec <- rep(as.numeric(gsub("T", "", rownames(Z_full))), each = n_lf)[!sl_na_lf]

  # loadings matrix for intertemporal restrictions
  itr_Dload <- c(1 / 3, 2 / 3, 1, 2 / 3, 1 / 3) / 3
  itr_Lload <- rep(1 / 3, 3)
  
  Lambda <- matrix(0, T_obs_lf, Tn_lat)
  colnames(Lambda) <- Ypos_vec[id_lat_vec]

  for (nn in 1:n_lf) {
    sl_var <- varlabs[nn]
    sl_Lambda_rows <- which(Zvar_vec == nn)
    maxT <- Zt_vec[sl_Lambda_rows]

    for (tt in seq_along(maxT)) {
      if (itr[nn] == "D") {
        Lambda[sl_Lambda_rows[tt], paste0(sl_var, "_T", maxT[tt]:(maxT[tt] - 4))] <- itr_Dload
      } else {
        Lambda[sl_Lambda_rows[tt], paste0(sl_var, "_T", maxT[tt]:(maxT[tt] - 2))] <- itr_Lload
      }
    }
  }
  
  Lambda <- Matrix(Lambda, sparse = TRUE) # intertemporal restriction loadings
  O_mat <- Matrix(1e-12 * diag(T_obs_lf), sparse = TRUE)
  Oi_mat <- solve(O_mat)
  
  # --------------------------------------------------------------------------
  # prior settings
  A_pr <- matrix(0, k, n) 
  s2.ARp <- matrix(1,n,1) # Scaling set to 1 instead of running a set of AR(p) models 
  
  # VAR contemporaneous relationships in state equation
  LQ_pr       <- matrix(0,n,n)
  
  # Minnesota moments 
  own.slct <- matrix(FALSE, k, n)
  other.slct <- matrix(FALSE, k, n)
  for(dd in 1:n){
    own.slct[seq(dd, n*p, by = n),dd]   <- TRUE       # identify own lags
    other.slct[setdiff(1:k, seq(dd, n*p, by = n)),dd] <- TRUE # identify other lags (and intercept)
  }

  tau_A   <- rep(1, n)
  tau.1_A <- rep(1, n)
  tau.2_A <- rep(0.5, n)
  zeta_A  <- rep(1, n)

  # Local scalings according to Minnesota moments
  lam_A <- matrix(1, k, n)
  nu_A  <- matrix(1, k, n)
  
  # Elements of prior variance-covariance matrix
  psi_A <- theta_A  <- matrix(1, k, n)
  
  for(dd in 1:n){
    own.slct.dd <- own.slct[,dd]
    other.slct.dd <- other.slct[,dd]
    lam_ai <- rep(1, k)
    lam_ai[1:(n*p)] <- 1/rep((1:p)^2, each = n)

    lam_ai[1:(n*p)] <- lam_ai[1:(n*p)]*s2.ARp[dd]/rep(s2.ARp, p)
    
    lam_A[,dd] <- nu_A[,dd]   <- lam_ai
    theta_A[own.slct.dd,dd]   <- psi_A[own.slct.dd,dd]   <- tau.1_A[dd]*lam_A[own.slct.dd,dd]
    theta_A[other.slct.dd,dd] <- psi_A[other.slct.dd,dd] <- tau.2_A[dd]*lam_A[other.slct.dd,dd]
  }
  
  # Prior mean and variance of contemporaneous relationships
  psi_LQ   <- theta_LQ <- matrix(1,n,n)
  lam_LQ   <- matrix(1, n, n)
  nu_LQ    <- matrix(1, n, n)
  tau_LQ   <- 1
  zeta_LQ  <- 1
  
  id.ind_LQ <- matrix(1:n^2, n, n)
  id.slct_LQ  <- lower.tri(id.ind_LQ, diag = F)
  
  # Prior of SV processes
  sv_pr <- specify_priors(
    mu  = sv_normal(0,10),     # prior on unconditional mean in the state equation
    phi = sv_beta(shape1 = 25, shape2 = 1.5), #informative prior to push the model towards a random walk in the state equation (for comparability)
    sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.1)), # Gamma prior on the state innovation variance
    nu = sv_infinity(),
    rho = sv_constant(0)
  )
  q_pr <- 3
  Q_pr <- 3
  
  # --------------------------------------------------------------------------
  # initialization of the state equation parameters
  A.ols <- solve(crossprod(X)) %*% crossprod(X,Y)  # Reduced-form OLS coefficients
  eps.x.ols <- Y - X %*% A.ols # Reduced-form shocks
  
  Q.ols <- crossprod(eps.x.ols)/T  # OLS variance-covariance matrix
  LQ.inv.ols <- t(chol(Q.ols))       # Lower Cholesky factor of Q
  LQ.inv.ols <- LQ.inv.ols*t(matrix(1/diag(LQ.inv.ols), n, n))
  
  LQ.ols <- solve(LQ.inv.ols)

  # Initialization of VAR coefficients in state equation
  B_draw      <- A.ols * 0.5
  LQ_draw     <- LQ.ols

  # Initialization of SV processes in state equation
  lnqt_draw  <- t(matrix(log(diag(Q.ols)), n, T))
  sv_draw <- list()
  for (dd in 1:n) sv_draw[[dd]] <- list(mu = 0, phi = 0.95, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  
  # --------------------------------------------------------------------------
  # set up the dynamic coefficients
  Alag_sl <- matrix(seq(1, n * p), n)
  H_ls <- list()
  for (pp in 0:p) {
    if (pp == p) h10 <- 1 else h10 <- (-1)
    Htmp <- matrix(0, T, T + p)
    for (tt in 1:T) {
      Htmp[tt, tt + pp] <- h10
    }
    H_ls[[paste0("Hp", p - pp)]] <- Matrix(Htmp, sparse = TRUE)
  }

  # big precision matrix
  Sigmai_T <- Matrix(0, T*n, T*n)
  sl_Sigma_T <- matrix(1:(T*n), n, T)

  # structural-form shocks
  eps <- matrix(0, T, n)
  
  # --------------------------------------------------------------------------
  # storage etc.
  ntot <- nburn + nsave * nthin
  thin.set <- seq(nburn + 1, ntot, by = nthin)
  
  mcmclabs <- paste0("mcmc",1:nsave)
  savecount <- 0
  
  Yq_store <- array(NA, c(nsave, T, n_lf))
  dimnames(Yq_store) <- list(mcmclabs, datelabs_p, varlabs[1:n_lf])
  
  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  for (irep in 1:ntot) {
    # -------------------------------------------
    # BLOCK 1: update state equation parameters
    # STEP 1.a: sample the VAR equation-by-equation
    for(nn in 1:n){
      normalizer <- exp(-lnqt_draw[,nn]/2)
      YY <- Y[,nn]*normalizer
      if(nn != 1){
        XXo <- cbind(X,-Y[,1:(nn-1)])
        XX <- XXo*normalizer
        VA_inv <- diag(k+(nn-1))/c(theta_A[,nn],theta_LQ[nn,1:(nn-1)])
        a_prior <- c(A_pr[,nn],rep(0,nn-1))
        k_n <- ncol(XX)
      }else{
        XXo <- X
        XX <- XXo*normalizer
        VA_inv <- diag(k)/theta_A[,nn]
        a_prior <- A_pr[,nn]
        k_n <- k
      }
      
      VA_post <- solve(crossprod(XX) + VA_inv)
      A_post <- VA_post %*% (VA_inv %*% a_prior + crossprod(XX,YY))
      b_draw <- A_post + t(chol(VA_post))%*%rnorm(k_n)
      
      B_draw[,nn] <- b_draw[1:k]
      if(nn>1) LQ_draw[nn,1:(nn-1)] <- b_draw[-c(1:k)]
      eps[,nn] <- Y[,nn] - XXo%*%b_draw
      
      # Update equation-specific shrinkage on VAR coefficients
      hs_A <- get.hs.min(bdraw     = (b_draw[1:k] - A_pr[,nn]),
                         lam.hs    = lam_A[,nn],
                         nu.hs     = nu_A[,nn],
                         tau.hs    = tau_A[nn],
                         zeta.hs   = zeta_A[nn],
                         update.ls = TRUE)
      psi_A[,nn] <- hs_A$psi
      lam_A[,nn] <- hs_A$lam; nu_A[,nn]  <- hs_A$nu
      tau_A[nn]  <- hs_A$tau; zeta_A[nn] <- hs_A$zeta
    }
    
    # map back to reduced form
    LQ.inv_draw <- solve(LQ_draw)
    A_draw <- t(LQ.inv_draw%*%t(B_draw))
    
    # Update overall shrinkage on contemp. relationships
    hs_LQ <- get.hs.min(bdraw=as.numeric(LQ_draw[id.slct_LQ]-LQ_pr[id.slct_LQ]),lam.hs = lam_LQ[id.slct_LQ], nu.hs = nu_LQ[id.slct_LQ], tau.hs = tau_LQ, zeta.hs = zeta_LQ, update.ls = TRUE)
    psi_LQ[id.slct_LQ]           <- hs_LQ$psi
    lam_LQ[id.slct_LQ]           <- hs_LQ$lam 
    nu_LQ[id.slct_LQ]            <- hs_LQ$nu
    tau_LQ <- hs_LQ$tau; zeta_LQ <- hs_LQ$zeta
    
    psi_LQ[psi_LQ < 1e-8] <- 1e-8
    psi_LQ[psi_LQ > 1]   <- 1
    
    theta_LQ <- psi_LQ
    
    # shrinkage on VAR coefficients
    psi_A[psi_A < 1e-8] <- 1e-8
    psi_A[psi_A > 1]   <- 1
    theta_A <- psi_A
    
    # STEP 1.c: sample the variances
    if(sv){
      for(dd in 1:n){
        sv_dd <- svsample_general_cpp(eps[,dd], startpara = sv_draw[[dd]], startlatent = lnqt_draw[,dd], priorspec = sv_pr)
        svpara_dd <- sv_dd$para[, c("mu", "phi","sigma")]
        sv_draw[[dd]][c("mu", "phi","sigma")] <- as.list(svpara_dd)

        lnqt_draw[,dd]     <- as.numeric(sv_dd$latent)
      }
    }else{
      for(dd in 1:n){
        qt.dd <- 1/rgamma(1, q_pr + T/2, Q_pr + sum(eps[,dd]^2)/2)
        lnqt_draw[,dd]    <- log(qt.dd)
      }
    }
    lnqt_draw[lnqt_draw < log(1e-3)] <- log(1e-3)
    lnqt_draw[lnqt_draw > log(3)] <- log(3)

    for(tt in 1:T){
      Sigmai_T[sl_Sigma_T[, tt], sl_Sigma_T[, tt]] <- t(LQ_draw) %*% diag(1/exp(lnqt_draw[tt,])) %*% LQ_draw
    }
    
    # -------------------------------------------
    # BLOCK 2: sample all latent states
    if(cons){
      a0_draw <- A_draw[k,]
    }else{
      a0_draw <- rep(0, n)
    }
    
    h <- c(t(matrix(1,T,1) %*% a0_draw)) # set up companion full data matrices
    H <- kronecker(H_ls[[paste0("Hp", 0)]], diag(n))
    for (pp in 1:p) {
      H <- H + kronecker(H_ls[[paste0("Hp", pp)]], t(A_draw)[, Alag_sl[, pp]])
    }
    
    Gm <- H %*% Sm
    Gl <- H %*% Sl
    
    GltSigi <- crossprod(Gl, Sigmai_T)
    Sigbar_i <- GltSigi %*% Gl # precision
    diag(Sigbar_i)[diag(Sigbar_i) < 1] <- 1 # offset in case precision is too small
    
    # unconditional on observed low-frequency info ---
    # explicit computation: mubar <- solve(Sigbar_i) %*% GltSigi %*% (h - Gm %*% Yo_vec)
    Sigbar_i_chol <- chol(Sigbar_i)
    mubar <- solve(Sigbar_i_chol, solve(t(Sigbar_i_chol), GltSigi %*% (h - Gm %*% Yo_vec))) # precision-based
    
    LambdaOi <- crossprod(Lambda, Oi_mat)
    Sigbar_const_i <- LambdaOi %*% Lambda + Sigbar_i # precision
    
    Sigbar_const_i_chol <- try(chol(Sigbar_const_i), silent = TRUE)
    if (is(Sigbar_const_i_chol, "try-error")) {
      diag(Sigbar_const_i) <- diag(Sigbar_const_i) + eps_small # in case of numerical issues, add a small offsetting constant
      Sigbar_const_i_chol <- chol(Sigbar_const_i)
    }
    
    # conditioning on low-frequency info ---
    # explicit computation: mubar_const <- solve(Sigbar_const_i) %*% (LambdaOi %*% Z_vec + Sigbar_i %*% mubar)
    mubar_const <- solve(Sigbar_const_i_chol, solve(t(Sigbar_const_i_chol), LambdaOi %*% Z_vec + Sigbar_i %*% mubar)) # precision-based
    Ym_vec <- mubar_const + solve(Sigbar_const_i_chol, rnorm(Tn_lat))
    
    # reconfigure data for the next iteration
    Y_full_vec <- Sm %*% Yo_vec + Sl %*% Ym_vec
    Y_full_raw <- matrix(Y_full_vec, T + p, n, byrow = TRUE)
    Y_full <- scale(Y_full_raw) # normalize so that each dataset has mean 0 and sd 1
    
    Ylags <- embed(Y_full, dim = p + 1)
    Y <- Ylags[, 1:n]
    X <- Ylags[, -c(1:n)]
    if (cons) X <- cbind(X, 1)

    if (irep %in% thin.set) {
      savecount <- savecount + 1
      Y_tmp <- matrix(NA, T + p, n_lf)
      for(j in 1:n_lf){
        Y_tmp[, j] <- ((stats::filter(Y_full_raw[,j], filter = itr_Dload, method = "convolution", sides = 1)) * Ysd[j]) + Ymu[j]
      }
      
      Yq_store[savecount,,] <- Y_tmp[-c(1:p),]
    }
    setTxtProgressBar(pb, irep)
  }
  return(Yq_store)
}


