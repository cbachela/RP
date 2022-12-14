  
  
  ############################################################################
  ### RP PACKAGE - TESTING RPE CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.01.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  test_har_ellipsoid <- function()
  {
    
    require(rgl)
    require(RP)
    wd <- "H:/R/RP/R/"
    source( paste0(wd, "class_ellipsoid.R") )
    source( paste0(wd, "class_RPE.R") )
    source( paste0(wd, "rp_geometric_walks.R") )
    
    
    d <- 3
    X <- Data[ ,1:d]
    covmat <- cov(X)
    GPO <- gpo( gps(Data = X) )
    w_mv <- getWeights(GPO)
    w_eqw <- rep(1/d, d)
    sigma <- as.numeric( t(w_mv) %*% covmat %*% w_mv )
    trerr <- as.numeric( t(w_mv - w_eqw) %*% covmat %*% (w_mv - w_eqw) )
    gps_tmp <- gps(Data = X, 
                   Solver = solverCtrl(portfolio = "meanvariance",
                                       progtype = "QCP",
                                       utility = list(riskaversion = 0),
                                       obj_lin = -1 / (w_mv + 0.000001)))
    addConstraint(gps_tmp) <- varianceConstraint(rhs = sigma * 1.05, Qmat = covmat)
    gpo_tmp <- gpo(gps_tmp)
    w_boundary <- getWeights(gpo_tmp)
    
    E <- Ellipsoid$new(centre = rep(0, d),
                       shape = covmat,
                       rhs = sigma * 1.05,
                       linear = rep(0, d),
                       constant = 0,
                       interior_point = w_mv)

    n_sim <- 10^3
    samples_sphere <- E$sample( n_sim = n_sim, algo = "sphere" )[[1]]
    # debugonce(har)
    lSamples_har <- E$sample( n_sim = n_sim, algo = "hitnrun" )
    samples_har <- lSamples_har[[2]]
    
    
    # Add linear and constant part (e.g. tracking errror constraint)
    qvec <- -2 * as.numeric( covmat %*% w_eqw )
    constant <- as.numeric( t(w_eqw) %*% covmat %*% w_eqw )
    E2 <- E$copy()
    E2$linear <- qvec
    E2$constant <- constant
    E2$rhs <- as.numeric( sigma + constant + t(qvec) %*% w_mv )
    E2$interior_point <- w_eqw
    # debugonce(har)
    samples_2_har <- E2$sample( n_sim = n_sim, algo = "hitnrun" )[[2]]
    
    
    E$rhs
    range( apply(samples_sphere, 1, function(x) { t(x) %*% covmat %*% x }))
    range( apply(samples_har, 1, function(x) { t(x) %*% covmat %*% x }))
    
    E2$rhs
    t(w_mv - w_eqw) %*% covmat %*% (w_mv - w_eqw)
    range( apply(samples_2_har, 1, function(x) { t(x) %*% covmat %*% x + t(qvec) %*% x + constant }))
    range( apply(samples_2_har, 1, function(x) { t(x - w_eqw) %*% covmat %*% (x - w_eqw) }))
    
    
    
    # Plot
    plot3d(x = samples_sphere[ ,1], y = samples_sphere[ ,2], z = samples_sphere[ ,3])
    points3d(x = samples_har[ ,1], y = samples_har[ ,2], z = samples_har[ ,3], col = 2)
    points3d(x = samples_2_har[ ,1], y = samples_2_har[ ,2], z = samples_2_har[ ,3], col = 3)
    
    
    
    solve(-2 * covmat) %*% qvec - w_eqw   # same same
    
    
    
    
  }
  
  
  
  # --------------------------------------------------------------------------
  test_var_tracking_error <- function()
  {
    
    require(rgl)
    require(RP)
    source("H:/R/RP/R/class_ellipsoid.R")
    source("H:/R/RP/R/class_RPE.R")
    
    d <- 3
    X <- Data[ ,1:d]
    covmat <- cov(X)
    GPO <- gpo( gps(Data = X) )
    w_mv <- getWeights(GPO)
    w_bm <- rep(1/d, d)
    sigmaFUN <- function(x) { as.numeric( t(x) %*% covmat %*% x ) }
    trerrFUN <- function(x) { as.numeric( t(x - w_bm) %*% covmat %*% (x -w_bm) ) }
    sigma <- sigmaFUN(x = w_mv)
    trerr <- trerrFUN(x = (w_mv + w_bm) / 2)
    
    gps_tmp <- gps(Data = X, 
                   Solver = solverCtrl(portfolio = "meanvariance",
                                       progtype = "QCP",
                                       utility = list(riskaversion = 0),
                                       obj_lin = -1 / (w_mv + 0.000001)))
    addConstraint(gps_tmp) <- varianceConstraint(rhs = sigma * 1.05, Qmat = covmat)
    gpo_tmp <- gpo(gps_tmp)
    w_boundary_mv <- getWeights(gpo_tmp)
    
    
    E1 <- Ellipsoid$new(centre = rep(0, d),
                        shape = covmat,
                        rhs = sigma * 1.05,
                        linear = rep(0, d),
                        constant = 0,
                        interior_point = w_mv)
    samples_E1 <- E1$sample( n_sim = 10^5, algo = "sphere" )[[1]]
    
    E2 <- Ellipsoid$new(centre = w_bm,
                        shape = covmat,
                        rhs = trerr * 0.5,
                        linear = rep(0, d),
                        constant = 0,
                        interior_point = w_bm)
    samples_E2 <- E2$sample( n_sim = 10^5, algo = "sphere" )[[1]]
    
    
    S <- simplex(d = d)
    RPS <- rp(S, spec = rpCtrl(n_sim = 10^4))
    samples_S <- getSamples(RPS)
    
    # Plot
    plot3d(x = samples_S[ ,1], y = samples_S[ ,2], z = samples_S[ ,3], col = "grey")
    points3d(x = samples_E1[ ,1], y = samples_E1[ ,2], z = samples_E1[ ,3])
    points3d(x = samples_E2[ ,1], y = samples_E2[ ,2], z = samples_E2[ ,3], col = "green")
    
    
    points3d(x = samples_har[ ,1], y = samples_har[ ,2], z = samples_har[ ,3], col = 2)
    points3d(x = samples_2_har[ ,1], y = samples_2_har[ ,2], z = samples_2_har[ ,3], col = 3)
    
    
    
  }
  
  
  
  # --------------------------------------------------------------------------
  test <- function()
  {
    
    require(RP)
    source("H:/R/RP/R/class_ellipsoid.R")
    source("H:/R/RP/R/class_RPE.R")
    
    d <- 3
    X <- Data[ ,1:d]
    covmat <- cov(X)
    GPO <- gpo( gps(Data = X) )
    w_mv <- getWeights(GPO)
    w_eqw <- rep(1/d, d)
    sigma <- as.numeric( t(w_mv) %*% covmat %*% w_mv )
    gps_tmp <- gps(Data = X, 
                   Solver = solverCtrl(portfolio = "meanvariance",
                                       progtype = "QCP",
                                       utility = list(riskaversion = 0),
                                       obj_lin = -1 / (w_mv + 0.000001)))
    addConstraint(gps_tmp) <- varianceConstraint(rhs = sigma * 1.05, Qmat = covmat)
    gpo_tmp <- gpo(gps_tmp)
    w_boundary <- getWeights(gpo_tmp)
  
    E <- Ellipsoid$new(centre = rep(0, d),
                       shape = covmat,
                       rhs = sigma * 1.05,
                       linear = rep(0, d),
                       interior_point = w_mv)
    
    E  
    E$shape
    samples_E <- E$sample( n_sim = 10^1, algo = "sphere" )
    # debugonce(E$sample)
    # debugonce(har)
    samples_E <- E$sample( n_sim = 10^1, algo = "hitnrun" )
    E_prime <- E$copy()
    E_prime
    
    transformation <- simplexTransform(d = E$dimension() )
    tmat <- transformation$tmat
    z <-transformation$z
    
    E_prime <- E$transform( x_interior = w_mv,
                            x_boundary = w_boundary )
    E_prime
    
    samples_E_prime <- E_prime$sample()
    head(samples_E_prime[[1]])
    
    
    ###
   
    
    K <- ConvexBody( ellipsoid = E,
                     spec = rpCtrl(n_sim = 10^3,
                                   algo = "sphere",
                                   x0 = w_mv,
                                   x_boundary = w_boundary) )
    K  
    # debugonce(K$transform)
    K_prime <- K$transform()
    # K_orig <- K_prime$backtransform()
    
    
    
    
    # debugonce(K$ellipsoid$boundaryPoint)
    K$ellipsoid$boundaryPoint( u = w_eqw - E$centre )
    
    
    
    
   
    
    # debugonce(E$sample)
    E$sample( n_sim = 10, algo = "hitnrun" )
    
    
    
    
  }
  