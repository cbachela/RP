  
  
  ############################################################################
  ### MINKOWSKI SUM OF TWO ELLIPSOIDS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.01.2020
  # First version:    20.01.2020
  # --------------------------------------------------------------------------
  
  
  # We distinguish two cases:
  # 1) Both ellipsoids have the same shape matrix but different centre
  # 2) Ellipsoids have different shape matrix and different centre
  
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(GPO)
  require(RP)
  require(rgl)
  wd <- "R:/Asset_Management/Research_Projects/Equity/Optimization/Multi_Objective/"
  source( paste0(wd, "Source/mop_functions.R") )
  wd_vpb <- "R:/Asset_Management/Research_Projects/Assignments/VP_Bank/VPB_2020_01/"
  source( paste0(wd_vpb, "Source/vpb_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # Functions
  # --------------------------------------------------------------------------
  
  sigmaFUN <- function(x) { as.numeric( t(x) %*% covmat_sigma %*% x ) }
  trerrFUN <- function(x) { as.numeric( t(x - w_bm) %*% covmat_te %*% (x - w_bm) ) }
  toFUN <- function(x) { sum(abs(x - w_init)) }
  objFUN <- function(x) 
  { 
    a <- as.numeric( t(x) %*% (covmat_sigma + lambda * (covmat_te - covmat_sigma)) %*% x )
    b <- as.numeric( 2 * lambda * (w_bm %*% covmat_te) %*% x )
    g <- as.numeric( t(w_bm) %*% (covmat_te * lambda) %*% w_bm )
    return( a - b + g )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # Case 1)
  # --------------------------------------------------------------------------
  
  d <- 3
  X <- Data[ ,1:d]
  covmat <- covariance(X, covCtrl(method = "avg"))
  covmat_sample <- cov(X)
  
  
  w_init <- setNames( rep(1/ncol(X), ncol(X)), colnames(X) )
  # upper <- setNames( 1:ncol(X) / sum(1:ncol(X)) + min(w_init) / 1 , colnames(X) )
  # upper <- w_init * 0 + 1
  upper <- setNames( c(1, 1, 0.5), colnames(X) )
  w_bm <- setNames(  rev(1:ncol(X) / sum(1:ncol(X))), colnames(X) )
  w_bm <- setNames(  c(1, 3, 9), colnames(X) )
  w_bm <- w_bm / sum(w_bm)
  
  cbind(w_init, upper, w_bm)
  
  
  # Prepare GPS
  GPS <- gps(Data = X)
  addConstraint(GPS) <- boxConstraint(upper = upper)
  setCovariance(GPS) <- covmat
  
  # Minimum variance portfolio
  gpo_mv <- gpo(GPS = GPS)
  
  # Global minimum variance portfolio (long-short)
  GPS_ls <- GPS
  addConstraint(GPS_ls) <- boxConstraint(name = "Unbounded")
  gpo_gmv <- gpo(GPS = GPS_ls)
  
  # Minimum tracking error portfolio
  setSolver(GPS) <- solverCtrl(custom_portfolio = "mintrackingerror",
                               w_bm = w_bm,
                               w_init = w_init)
  setCovariance(GPS) <- covmat_sample
  gpo_mte <- gpo(GPS = GPS)
  
  # Minimum turnover portfolio
  setSolver(GPS) <- solverCtrl(custom_portfolio = "minturnover",
                               progtype = "LP",
                               w_bm = w_bm,
                               w_init = w_init)
  gpo_mto <- gpo(GPS = GPS)
  
  w_mv <- getWeights(gpo_mv)
  w_gmv <- getWeights(gpo_gmv)
  w_mte <- getWeights(gpo_mte)
  w_mto <- getWeights(gpo_mto)
  
  wmat <- cbind(init = w_init, 
                bm = w_bm, 
                # ub = upper, 
                mv = w_mv,
                mte = w_mte, 
                mto = w_mto)
  wmat
  colors <- 1:ncol(wmat)
  barplot( t(wmat), beside = TRUE, col = colors )
  legend("topleft", colnames(wmat), lwd = 2, col = colors, text.col = colors, bty = "n")
  

  
  # # Standardize
  # covmat <- covmat / sigmaFUN(w_mte)
  # covmat_sample <- covmat_sample / trerrFUN(w_mv)
  
  
  
  
  S <- simplex(d = d)
  H <- hyperplane(v = rep(1, d), b = 1)
  E1 <- ellipsoid(q = rep(0, d), Q = covmat, b = sigmaFUN(w_mv) * 1.05)
  E2 <- ellipsoid(q = w_bm, Q = covmat_sample, b = trerrFUN(w_mte) * 1.05)
  
  # E1H <- intersect(x = E1, y = H)
  
  n_sim <- 10^3

  
  # Get boundary points 
  GPS <- gps(Data = X,
             Covariance = covmat,
             Solver = solverCtrl(portfolio = "meanvariance",
                                 progtype = "QCP",
                                 obj_lin = -1 / (w_mv + 0.001),
                                 obj_quad = covmat,
                                 utility = list(riskaversion = 0)))
  addConstraint(GPS) <- varianceConstraint(rhs = sigmaFUN(w_mv) * 1.05,
                                           Qmat = covmat)
  GPO <- gpo(GPS = GPS)
  w_boundary_mv <- getWeights(GPO)
  
  GPS <- gps(Data = X,
             Covariance = covmat_sample,
             Solver = solverCtrl(portfolio = "meanvariance",
                                 progtype = "QCP",
                                 obj_lin = -1 / (w_mte + 0.001),
                                 obj_quad = covmat_sample,
                                 w_bm = w_bm,
                                 utility = list(riskaversion = 0)))
  addConstraint(GPS) <- trackingerrorConstraint(rhs = trerrFUN(w_mte) * 1.05,
                                                Qmat = covmat_sample,
                                                w_bm = w_bm)
  GPO <- gpo(GPS = GPS)
  w_boundary_mte <- getWeights(GPO)
  
  cbind(w_boundary_mv, w_boundary_mte)
  
  
  RPE1 <- as.RP(object = E1, 
                rpCtrl(n_sim = n_sim, 
                       algo = "sphere", 
                       b_pushy = FALSE,
                       x0 = w_mv,
                       x_boundary = w_boundary_mv))
  E1T <- transform(object = RPE1 )
  RPE1T <- rp(E1T, getSpec(E1T))
  samples_E1T <- getSamples(RPE1T)
  RPE1 <- backTransform(RPE1T)
  samples_E1 <- getSamples(RPE1)

  RPE2 <- as.RP(object = E2, 
                rpCtrl(n_sim = n_sim, 
                       algo = "sphere", 
                       b_pushy = FALSE,
                       x0 = w_mte,
                       x_boundary = w_boundary_mte))
  # debugonce(.transformE)
  E2T <- transform(object = RPE2 )
  # debugonce(.rpE)
  RPE2T <- rp(E2T, getSpec(E2T))
  samples_E2T <- getSamples(RPE2T)
  RPE2 <- backTransform(RPE2T)
  samples_E2 <- getSamples(RPE2)
  
  
  
  
  samples_S <- getSamples( rp(S, rpCtrl(n_sim = 10^4)) )
  
  lambda <- 0.3
  samples_E3_emp <- (1 - lambda) * samples_E1 + lambda * samples_E2
  val_S <- apply(samples_S, 1, function(x) { (1 - lambda) * sigmaFUN(x) + lambda * trerrFUN(x) })
  w_inbetween <- getWeights( gpo(GPS = gps(Data = X,
                                           Covariance = (1-lambda) * covmat + lambda * covmat_sample,
                                           Solver = solverCtrl(custom_portfolio = "minVminTE",
                                                               utility = list(lambda = lambda),
                                                               w_bm = w_bm))) )
  # rhs <- ((1 - lambda) * sigmaFUN(w_inbetween) + lambda * trerrFUN(w_inbetween)) * 1.05
  rhs <- objFUN(w_inbetween) * (1.05)
  rhs
  # Identify points which satisfy lqc
  idx <- which(val_S <= rhs)
  idx  
  
  
  
  Qmat <- (1 - lambda) * getShape(E1) + lambda * getShape(E2)
  solve(-2 * Qmat) %*% (-2 * lambda * t(w_bm %*% covmat_sample))
  

  E3 <- ellipsoid(q = (1 - lambda) * getCentre(E1) + lambda * getCentre(E2),
                  Q = (1 - lambda) * getShape(E1) + lambda * getShape(E2),
                  b = (1 - lambda) * getRHS(E1) + lambda * getRHS(E2))
  GPS <- gps(Data = X,
             Covariance = getShape(E3),
             Solver = solverCtrl(portfolio = "meanvariance",
                                 progtype = "QCP",
                                 obj_lin = c(1, rep(0, d-1)),
                                 obj_quad = getShape(E3),
                                 w_bm = w_bm,
                                 utility = list(riskaversion = 0)))
  rhs <- 1.05 * ((1 - lambda) * sigmaFUN(w_inbetween) + lambda * trerrFUN(w_inbetween))
  Qmat <- covmat + lambda * (covmat_sample - covmat)
  qvec <- as.numeric( -2 * lambda * w_bm %*% covmat_sample )
  constant <- as.numeric( t(w_bm) %*% (covmat_sample * lambda) %*% w_bm )
  addConstraint(GPS) <- quadraticConstraint(rhs = rhs - constant,
                                            Qc = Qmat,
                                            q = qvec)
  GPO <- gpo(GPS = GPS)
  w_boundary_inbetween <- getWeights(GPO)
  
  # same same:
  getShape(E3) - Qmat
  rhs - ((1 - lambda) * sigmaFUN(w_boundary_inbetween) + lambda * trerrFUN(w_boundary_inbetween))

  
  
  RPE3 <- as.RP(object = E3,
                rpCtrl(n_sim = n_sim,
                       algo = "sphere",
                       b_pushy = FALSE,
                       x0 = w_inbetween,
                       x_boundary = w_boundary_inbetween))
  E3T <- transform(object = RPE3 )
  RPE3T <- rp(E3T, getSpec(E3T))
  RPE3 <- backTransform(RPE3T)
  samples_E3 <- getSamples(RPE3)
  
  
  
  
  RP1 <- rp( E1, 
             rpCtrl(n_sim = 10^4, 
                    algo = "sphere",
                    x0 = w_mv,
                    x_boundary = w_boundary_mv))
  samples_rp1 <- getSamples(RP1)
  RP2 <- rp( E2, 
             rpCtrl(n_sim = 10^4, 
                    algo = "sphere",
                    x0 = w_mte,
                    x_boundary = w_boundary_mte))
  samples_rp2 <- getSamples(RP2)
  RP3 <- rp( E3, 
             rpCtrl(n_sim = 10^4, 
                    algo = "sphere",
                    x0 = w_inbetween,
                    x_boundary = w_boundary_inbetween))
  samples_rp3 <- getSamples(RP3)
  FUN <- function(x) { t(qvec) %*% x + (1 - lambda) * sigmaFUN(x) + lambda * trerrFUN(x) }
  val_rp3 <- apply(samples_rp3, 1, FUN )
  range(val_rp3)
  rhs
  
  
  
  plot3d(x = samples_S[ ,1], y = samples_S[ ,2], z = samples_S[ ,3], col = "grey")    
  points3d(x = samples_E1[ ,1], y = samples_E1[ ,2], z = samples_E1[ ,3])
  points3d(x = samples_E2[ ,1], y = samples_E2[ ,2], z = samples_E2[ ,3], col = "red")
  points3d(x = samples_E3[ ,1], y = samples_E3[ ,2], z = samples_E3[ ,3], col = "green")
  points3d(x = w_mv[1], y = w_mv[2], z = w_mv[3], size = 10, col = "magenta")
  points3d(x = w_gmv[1], y = w_gmv[2], z = w_gmv[3], size = 10, col = "magenta")
  points3d(x = w_mte[1], y = w_mte[2], z = w_mte[3], size = 10, col = "orange")
  points3d(x = w_bm[1], y = w_bm[2], z = w_bm[3], size = 10, col = "orange")
  points3d(x = w_inbetween[1], y = w_inbetween[2], z = w_inbetween[3], size = 10, col = "blue")
  points3d(x = samples_S[idx, 1], y = samples_S[idx, 2], z = samples_S[idx, 3], size = 5, col = "cyan")   
 
  
  
  
  