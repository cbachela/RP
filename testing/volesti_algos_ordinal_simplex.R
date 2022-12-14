  
  
  ############################################################################
  ### RP PACKAGE - TEST VOLESTI ALGOS ON ORDINALLY CONSTRAINED SIMPLEX
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     17.05.2021
  # First version:    17.05.2021
  # --------------------------------------------------------------------------

  
  require(volesti)
  require(GPO)
  require(RP)
  
  
  ?sample_points
  
  
  
  # --------------------------------------------------------------------------
  # Toy example
  # --------------------------------------------------------------------------
  
  n <- 3
  B <- 10^3
  set.seed(1111)
  x <- rnorm(n)
  z <- (x - mean(x))^2
  # names(z) <- aviationAlphabet()[1:n]
  names(z) <- paste0("X", 1:n)
  
  
  Constraints <- constraints( selection = names(z) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint( lower = setNames(rep(0, n), names(z)),
                                               upper = setNames(rep(1, n), names(z)) )
  addConstraint(Constraints) <- orderingConstraint( Names = names(z),
                                                    ordering = 1:n, 
                                                    eps = 1e-04 )
  P <- as.polytope(Constraints) 
  
  # Interior point used for starting the MCMC
  x0 <- setNames( (1:n) / sum(1:n), names(z) )
  
  # Specifications
  rp_ctrl <- rpCtrl( n_sim = 10^4,
                     thin = 10^3,
                     algo = "volesti",
                     x0 = x0 )
  
  # debugonce( .rpP )
  RPP <- rp( object = P, spec = rp_ctrl )
  
  # # Transformation (from variable to sample space)
  # P_prime <- transform( as.RP(P, spec = rp_ctrl) )
  # 
  # # Sampling (in sampling space)
  # # debugonce(.rpP)
  # RPP_prime <- rp( object = P_prime, spec = getSpec(P_prime) )
  
  samples <- getSamples(RPP)
  lp_samples <- apply( samples, 1, lpNorm, p = 2 )
  
  
  
  # Use Dirichlet
  alpha <- apply(samples, 2, mean)
  Omega <- rdirichlet( n = rp_ctrl$n_sim,
                       alpha = alpha * n )
  lp_omega <- apply( Omega, 1, lpNorm, p = 2 )
  
  cbind( alpha, apply(Omega, 2, mean) )
  
  
  plot( density(lp_samples), xlim = c(0, 1) )  
  lines( density(lp_omega), col = 2 )
  abline( v = lpNorm(alpha, p = 2) )
  abline( v = lpNorm(apply(Omega, 2, mean), p = 2), col = 2 )
  
  
  
  
  # 3D plot
  
  samples_simplex <- getSamples( rp( simplex( d = n ), rpCtrl( n_sim = 10^4 ) ) )
  
  plot3d( x = samples_simplex[ ,1], y = samples_simplex[ ,2], z = samples_simplex[ ,3], 
          xlab = names(z)[1], ylab = names(z)[2], zlab = names(z)[3], col = "lightgrey" )
  points3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "darkgrey" )
  points3d( x = Omega[ ,1], y = Omega[ ,2], z = Omega[ ,3], col = "green" )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Higher dimensions
  # --------------------------------------------------------------------------
  
  n_sim <- 10^3
  n_thin <- 10^1
  tau <- 10
  X <- log( 1 + GPO::Data[ ,1] )
  X <- head(X, 100)
  mu <- mean(X)
  mu_ewma <- mean( weightEWMA( X = X, tau = tau ) )
  n <- nrow(X)
  z <- setNames( as.numeric(X), rownames(X) ) 
 
  lambda <- exp(-log(2) / tau)
  wt <- lambda^(0:(n-1))
  wt <- rev( length(wt) * wt / sum(wt) ) / n
 
  
  # Using Volesti
  Constraints <- constraints( selection = names(z) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint( lower = setNames(rep(0, n), names(z)),
                                               upper = setNames(rep(1, n), names(z)) )
  addConstraint(Constraints) <- orderingConstraint( Names = names(z),
                                                    ordering = 1:n, 
                                                    eps = 0 )
  P <- as.polytope(Constraints) 
  
  # Interior point used for starting the MCMC
  x0 <- wt
  
  
  walks <- c("CDHR", "BiW", "john", "vaidya", "dikin", "billiard")
  env <- new.env()

  for ( i in seq(along = walks) ) {
    
    # Specifications
    rp_ctrl <- rpCtrl( n_sim = n_sim,
                       thin = n_thin,
                       algo = "volesti",
                       volesti = list( density = "uniform",
                                       walk = walks[i] ),
                       x0 = x0 )
    if ( walks[i] == "billiard" ) {
      rp_ctrl$algo = walks[i]
    }
    
    tic <- Sys.time()
    RPP <- rp( object = P, spec = rp_ctrl )
    toc <- Sys.time() - tic
    assign(envir = env, paste0("RPP_", walks[i]), RPP)
    assign(envir = env, paste0("tictoc_", toc), RPP)
    
  }
 
  
  lTheta <- list()
  for ( i in seq(along = walks) ) {
    
    RPP <- get( paste0("RPP_", walks[i]), envir = env )
    samples <- getSamples(RPP)
    lTheta[[i]] <- apply( samples, 1, function(p) { sum( p * z ) } )
    
  }
  names(lTheta) <- walks  
  theta_mat <- do.call( cbind, lTheta )  
  boxplot( theta_mat )  

  
  plot( samples[ ncol(samples)] )
  
  
  
  
  