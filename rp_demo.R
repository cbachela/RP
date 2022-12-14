  
  
  ############################################################################
  ### RP PACKAGE - DEMO
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.07.2019
  # First version:    07.07.2019
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  test.geometricBodies()
  {
    
    require(rgl) # available on CRAN
    require(RP)
    
    # Parameters
    d <- 3   # Use low dimension so that we can plot things later
    n_sim <- 10^4 * 2  # Number of MCMC samples to be returned
    Names <- paste0("x", 1:d) # Give the dimensions names (for plotting)
    x0 <- setNames( rep(1/d, d), Names )  # Interior point used for starting the MCMC
    x_boundary <- c(0.3522251, 0.2398346, 0.4079403)  
    # x_boundary is needed to transform the ellipsoid. Can be obtained by numerical
    # optimization. Will be implemented in the future.
    
    
    # --------------------------------------------------------------------------
    # Step 1: Building geometric bodies
    # --------------------------------------------------------------------------
    
    # Simplex
    S <- simplex(d = d)
    
    # Ellipsoid
    ## Begin by defining a shape matrix as the sample covariance matrix of
    ## some mulitvariate normal data
    n_row <- 10^3
    set.seed(1234)
    X <- matrix( rnorm(n_row * d) / 100, nrow = n_row, ncol = d)
    X[ ,2] <- X[ ,2] * 0.8
    X[ ,3] <- X[ ,3] * 1.3
    colnames(X) <- Names
    covmat <- cov(X)
    ## Set the centre to the be the origin
    centre <- rep(0, d)
    ## Define un upper bound
    b <- as.numeric(t(x0) %*% covmat %*% x0) * (1 + 0.2)
    
    ## Instantiate the ellipsoid
    E <- ellipsoid(q = centre, 
                   Q = covmat, 
                   b = b)
    E
   
    
    # Bounded Simplex
    upper_bounds <- c(0.8, 0.7, 0.6)
    S2 <- boundedSimplex(d = d, v = 1, upper = upper_bounds)
    S2
    
    # Polytope
    A <- matrix( c(1, 1, 0,
                   1, 0, 1),
                 nrow = 2, byrow = TRUE)
    b <- c(0.9, 0.5)
    sense <- c("<=", ">=")
    P <- polytope(A = A,
                  b = b,
                  sense = sense)
    addHyperplane(P) <- hyperplane(v = c(1, 1, 0), b = 0.8)  # renders a previous linear constraint obsolete
    P
    
    
    # Redefine new Polytope as the intersection of a Simplex and a Polytope
    PS <- intersect(x = S, y = P)
    PS
    

    # Convex Body
    B <- convexBody(P = PS, E = E)
    str(B)
    
    
    
    
    # Use get methods to access parts of the objects
    # e.g.:

    str(getPolytope(B))
    str(getEllipsoid(B))
    getAmat(P)
    getRHS(P)
    getShape(E)
    getCentre(E)
    getHyperplane(P)
    getHyperplane(PS)
    
    
    
    
    # --------------------------------------------------------------------------
    # Step 2: Transformation (from variable to sample space)
    # --------------------------------------------------------------------------
    
    S_prime <- transform( as.RP(S, spec = rpCtrl(n_sim = n_sim, x0 = x0)) )
    P_prime <- transform( as.RP(PS, spec = rpCtrl(n_sim = n_sim, x0 = x0)) )
    E_prime <- transform( as.RP(E, spec = rpCtrl(n_sim = n_sim, x0 = x0, x_boundary = x_boundary, algo = "sphere")) )
    # debugonce(.transformB)
    B_prime <- transform( as.RP(B, spec = rpCtrl(n_sim = n_sim, x0 = x0, x_boundary = x_boundary)) )
    
    tmp <- as.RP(B, spec = rpCtrl(n_sim = n_sim, x0 = x0, x_boundary = x_boundary))
    
    
    ### 
    # 25.07.2019, PCA
    
    E_prime1 <- getEllipsoid(B_prime)
    P_prime1 <- getPolytope(B_prime)
    covmat <- getShape(E_prime1)
    evec <- eigen(covmat)$vectors
    t(evec) %*% covmat %*% evec
    
    P_prime2 <- pMult(P = P_prime1, Gmat = t(evec))

    
    spec <- rpCtrl(n_sim = 10^3,
                   algo = "hitnrun",
                   x0 = getSpec(E_prime)$x0,
                   x_boundary = getSpec(E_prime)$x_boundary,
                   transformation = list(tmat = t(evec), 
                                         tmat_inv = evec, 
                                         z = rep(0, dimension(E_prime))))
   
    
    E_prime2 <- eMult(E = E_prime1, Amat = t(evec))
    x_b <- getSpec(E_prime)$x_boundary
    x_b_prime <- t(evec) %*% x_b
    rhs <- t(x_b_prime - getCentre(E_prime2)) %*% getShape(E_prime2) %*% (x_b_prime - getCentre(E_prime2))
    x0_prime <- t(evec) %*% spec$x0
    
    getShape(E_prime1)
    getShape(E_prime2)
    
    
    
    
    
    
    
    
    # debugonce(.rpE)
    # debugonce(dists2Ellipsoid)
    # debugonce(hitnrun)
    spec <- getSpec(E_prime)
    spec$algo <- "hitnrun"
    tmp <- rp( object = E_prime, spec = spec )
    
    
    
    
    
    
    # --------------------------------------------------------------------------
    # Step 3: Sampling (in sampling space)
    # --------------------------------------------------------------------------
  
    # debugonce(.rpP)
    # RPS_prime <- rp( object = S_prime, spec = getSpec(S_prime) )
    RPS_prime <- rp( object = S_prime, spec = S_prime@spec )
    
    # debugonce(.rpP)
    # RPP_prime <- rp( object = P_prime, spec = getSpec(P_prime) )
    RPP_prime <- rp( object = P_prime, spec = P_prime@spec )
    
    # debugonce(.rpE)
    # RPE_prime <- rp( object = E_prime, spec = getSpec(E_prime) )
    RPE_prime <- rp( object = E_prime, spec = E_prime@spec )
    
    # debugonce(.rpB)
    # RPB_prime <- rp( object = B_prime, spec = getSpec(B_prime) )
    RPB_prime <- rp( object = B_prime, spec = B_prime@spec )
    
    
    
    head(getSamples(RPS_prime))
    head(getSamples(RPB_prime))
    
    
   
    
    
    # --------------------------------------------------------------------------
    # Step 4: Back transformation (from sample to variable space)
    # --------------------------------------------------------------------------
    
    # debugonce(.backTransform)
    RPS <- backTransform(RPS_prime)
    RPP <- backTransform(RPP_prime)
    RPE <- backTransform(RPE_prime)
    RPB <- backTransform(RPB_prime)
    
    
    
    
    
    
    # --------------------------------------------------------------------------
    # All in one (i.e. steps 2 - 4) using wrapper rp
    # --------------------------------------------------------------------------
    
    # debugonce(.rpS)
    RPS <- rp( object = S, 
               spec = rpCtrl(n_sim = n_sim) )

    # debugonce(.rpE)
    RPE <- rp( object = E, 
               spec = rpCtrl(algo = "sphere", # could also use 'hitnrun', but is much slower
                             x0 = x0,
                             x_boundary = x_boundary,
                             n_sim = n_sim, 
                             b_pushy = FALSE) )
    
    # debugonce(.rpP)
    RPP <- rp( object = PS, 
               spec = rpCtrl(n_sim = n_sim, 
                             algo = "hitnrun", 
                             x0 = x0) )

    # debugonce(.rpB)
    RPB <- rp( object = B, 
               spec = rpCtrl(n_sim = n_sim, 
                             algo = "hitnrun", 
                             x0 = x0,
                             x_boundary = x_boundary) )
    
    
    
    
   
    
    # --------------------------------------------------------------------------
    # Visalization
    # --------------------------------------------------------------------------
    
    # Extract samples
    samples_S <- getSamples(RPS)
    samples_P <- getSamples(RPP)
    samples_E <- getSamples(RPE) 
    samples_B <- getSamples(RPB)

    
    
    # Plot
    plot3d( x = samples_S[ ,1], y = samples_S[ ,2], z = samples_S[ ,3], 
            xlab = Names[1],  ylab = Names[2], zlab = Names[3], col = "lightgrey" )
    points3d( x = samples_E[ ,1],  y = samples_E[ ,2], z = samples_E[ ,3], col = "orange" )
    points3d( x = samples_P[ ,1],  y = samples_P[ ,2], z = samples_P[ ,3], col = "darkgrey" )
    points3d( x = samples_B[ ,1],  y = samples_B[ ,2], z = samples_B[ ,3], col = "blue", size = 4 )
    # points3d( x = getSpec(RPE)$x_boundary[1],  
    #           y = getSpec(RPE)$x_boundary[2], 
    #           z = getSpec(RPE)$x_boundary[3], 
    #           col = "blue", size = 10 )
    
    
    
    
    # Plot samples obtained in reduced space
    
    samples_S_prime <- getSamples(RPS_prime)
    samples_P_prime <- getSamples(RPP_prime)
    samples_E_prime <- getSamples(RPE_prime) 
    samples_B_prime <- getSamples(RPB_prime)
    
    plot3d( x = samples_S[ ,1], y = samples_S[ ,2], z = samples_S[ ,3], 
            xlab = Names[1],  ylab = Names[2], zlab = Names[3], col = "lightgreen" )
    points3d( x = samples_S_prime[ ,1] * 0,  y = samples_S_prime[ ,1], z = samples_S_prime[ ,2], col = "lightgrey" )
    points3d( x = samples_P_prime[ ,1] * 0,  y = samples_P_prime[ ,1], z = samples_P_prime[ ,2], col = "darkgrey" )
    points3d( x = samples_E_prime[ ,1] * 0,  y = samples_E_prime[ ,1], z = samples_E_prime[ ,2], size = 3, col = "orange" )
    points3d( x = samples_B_prime[ ,1] * 0,  y = samples_B_prime[ ,1], z = samples_B_prime[ ,2], col = "blue" )
    
    
    
  }
  
  
  