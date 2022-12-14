  
  
  ############################################################################
  ### RP PACKAGE - TESTING
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     06.08.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # test.billiard
  # test.high_dim
  
  # --------------------------------------------------------------------------
  test.billiard <- function()
  {
    
    require(RP)
    
    # Parameters
    d <- 3   # Use low dimension so that we can plot things later
    n_sim <- 10^4 * 2  # Number of MCMC samples to be returned
    Names <- paste0("x", 1:d) # Give the dimensions names (for plotting)
    x0 <- setNames( rep(1/d, d), Names )  # Interior point used for starting the MCMC
    x_boundary <- c(0.3522251, 0.2398346, 0.4079403)  
    # x_boundary is needed to transform the ellipsoid. Can be obtained by numerical
    # optimization. Will be implemented in the future.

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
    
    # Transform
    spec <- rpCtrl(n_sim = 100, 
                   thin = 1,
                   x0 = x0, 
                   algo = "billiard",
                   jmp = rep(1, d))
    P_prime <- transform( as.RP(PS, spec) )
    
    
    tmat_inv <- getSpec(P_prime)$transform$tmat_inv
    z <- getSpec(P_prime)$transform$z
    jmp_prime <- spec$jmp - as.numeric( tmat_inv %*% z )
    jmp_prime
    
    
    
      
    # Sample
    # debugonce(.rpP)
    # debugonce(billiard)
    RPP_prime <- rp( object = P_prime, spec = getSpec(P_prime) )
    samples_prime <- getSamples(RPP_prime)
    
    
    plot( unlist(lK) )
    unlist(lK)
    lJmp
    
    
    colors <- fBasics::divPalette(n = nrow(samples_prime), "RdYlGn")
    plot(x = samples_prime[ ,1], y = samples_prime[ ,2], type = "o", pch = 19, col = colors )  
    
    plot( samples_prime[ ,1], type = "p", col = colors )
    
    
    
  }
  
  
  # --------------------------------------------------------------------------
  test.high_dim <- function()
  {
    
    require(RP)
    wd_data <- "H:/Papers Cyril/PhD_Proposal/R/Data/"
    env <- readRDS( file = paste0(wd_data, "europe_ex_ch_survivorshipbias_factor.rds") )
  
    
    
      
  }
  
  
  
  
  
  
  
  
  
  
  
  