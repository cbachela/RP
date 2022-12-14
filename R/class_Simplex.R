  
  
  ############################################################################
  ### RP PACKAGE - Simplex CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.01.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # Simplex
  # Simplex.initialize
  # Simplex.runif
  
  
  # --------------------------------------------------------------------------
  Simplex <- setRefClass( Class = "Simplex",
                          contains = "Polytope" )
  
  
  
  # --------------------------------------------------------------------------
  Simplex.initialize <- function( d = NULL, ... )
  {
    if ( !is.null(d) ) {
      A <<- rbind( matrix(1, nrow = 1, ncol = d), 
                   diag(d) * (-1), 
                   diag(d) )
      b <<- c(b, rep(0, d), rep(1, d))
      sense <<- c("=", rep("<=", d * 2))
    }
    return( TRUE )
  }
  Simplex$methods( initialize = Simplex.initialize )
  
  
  # --------------------------------------------------------------------------
  Simplex.runif <- function( n_sim = 10^2 )
  {
    d <- .self$dimension()
    U <- matrix( stats::runif(n_sim * d), nrow = n_sim, ncol = d )
    samples <- t( apply( U, 1, function(x) { log(x) / sum(log(x)) } ) )
    return( samples )
  }
  Simplex$methods( runif = Simplex.runif )
  
  
  
  
  
  
  
  