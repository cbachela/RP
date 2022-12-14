  
  
  ############################################################################
  ### RP PACKAGE - NON-MCMC SAMPLING
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.05.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  # FUNCTIONS:
  #
  # runifS
  # runifE
  # rdirichlet
  # ddirichlet
  
  
  
  # --------------------------------------------------------------------------
  #' @description Sample uniformly from standard simplex
  #' @export
  # --------------------------------------------------------------------------
  runifS <- function(d = NULL, n_sim = 100)
  {
    stopifnot( is.numeric(d) )
    U <- matrix( runif(n_sim * d), nrow = n_sim, ncol = d )
    S <- t( apply( U, 1, function(x) { log(x) / sum(log(x)) } ) )
    return( S ) 
  }
  
  
  
  
  # --------------------------------------------------------------------------
  #' @description Sample uniformly from ellipsoid
  #' @export
  # --------------------------------------------------------------------------
  runifE <- function(S, z_hat, gamma_threshold, n_points, b_pushy = TRUE)
  {
    
    
    
    # Do some checks first
    # ...
    n_dim <- dim(S)[1]
    # Generate random vectors from normal distribution
    U <- matrix( rnorm(n_dim * n_points), 
                 nrow = n_points, 
                 ncol = n_dim )
    #
    # U <- abs(U) # ~~~~~~~~~~~~~~~~
    #
    # Standardise each random vector by its L2-norm to get uniform distributed 
    # points on the surface of a hypersphere
    U <- t(apply(U, 1, function(x) { x / sqrt(sum(x^2)) }))
    # Generate scalar factors to push the points on the surface
    # of the hypersphere to its inside
    r_vec <- runif(n_points, 0, 1)^(1/n_dim)
    if ( isTRUE(b_pushy) ) {
      Y <- t(scale(t(U), FALSE, 1/r_vec))         #// but Y's are not uncorrelated...??
    } else {
      Y <- U
    }
    # Find square root of S using Cholesky factorization
    S <- makePositiveDefinite(m = S)    # ~~~~~~~~~~~~~~~   
    T_mat <- chol(S)
    # T_mat <- chol( solve(S) )
    # Hypersphere to hyperellipsoid mapping
    # X <- Y %*% t(T_mat)
    # X <- Y %*% T_mat
    # X <- Y %*% solve(T_mat)
    X <- Y %*% t(solve(T_mat))
    # Translation around centroid
    ans <- scale( X * sqrt(gamma_threshold), -z_hat, FALSE)     # should it be +z_hat ??
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @description Random generation from the Dirichlet distribution.
  #' @param n [scalar] Number of random vectors to generate.
  #' @param x	[vector] A vector containing a single deviate or matrix containing 
  #'                   one random deviate per row.
  #' @param alpha	[vector] Vector of shape parameters, or matrix of shape parameters 
  #'                       corresponding to the number of draw.
  #' @return Returns a matrix with n rows, each containing a single 
  #'         Dirichlet random deviate.
  #' @export
  # --------------------------------------------------------------------------
  rdirichlet <- function( n, alpha ) 
  {
    # Copy from MCMCpack
    # see ?MCMCpack:rdirichlet
    l <- length(alpha)
    x <- matrix( rgamma(l * n, alpha), ncol = l, byrow = TRUE )
    sm <- x %*% rep(1, l)
    ans <- x / as.vector(sm)
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  ddirichlet <- function( x, alpha ) 
  {
    # Copy from MCMCpack
    # see ?MCMCpack:rdirichlet
    dirichlet1 <- function(x, alpha) 
    {
      logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
      s <- sum((alpha - 1) * log(x))
      exp(sum(s) - logD)
    }
    if (!is.matrix(x)) 
      if (is.data.frame(x)) 
        x <- as.matrix(x)
      else x <- t(x)
      if (!is.matrix(alpha)) 
        alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                        byrow = TRUE)
      if (any(dim(x) != dim(alpha))) 
        stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
      pd <- vector(length = nrow(x))
      for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
                                                             ])
      pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
      pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
      return(pd)
  }
  
