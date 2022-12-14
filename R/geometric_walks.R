  
  
  ############################################################################
  ### RP PACKAGE - GEOMETRIC WALKS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  # FUNCTIONS:
  # hitnrun
  
  
  # --------------------------------------------------------------------------
  hitnrun <- function( K, x = NULL, u = NULL, n_sim = 100 )
  {
    
    if ( is.null(x) ) {
      stop("point x within the body is missing.")
      # x <- K$getInteriorPoint()
    }
    
    # stopifnot( .self$isInterior(x = x) )
    if ( !K$isInterior(x = x) ) {
      x <- x * 0.9                     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    }
    stopifnot( K$isInterior(x = x) )
    
    
    if ( is.null(u) ) {
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
    }
    
    ans <- matrix( 0, nrow = n_sim, ncol = length(x) )
    # U <- ans
    B_lo <- ans
    B_up <- ans
    
    itermax <- n_sim * 100
    i <- iter <- 1
    while ( i <= n_sim && iter < itermax ) {
      
      if ( K$isInterior(x) ) {  # is this necessary? only if we don't trust boundaryDistance
        
        D <- K$boundaryDistance(x = x, u = u)
        
        if ( sum(abs(D)) > 0 ) {
          
          # lower <- D["lower"]
          # upper <- D["upper"]
          lower <- D[1]
          upper <- D[2]
          
          # U[i, ] <- u
          B_lo[i, ] <- x + lower * u
          B_up[i, ] <- x + upper * u
          
          x <- x + (lower + (upper - lower) * runif(1)) * u
          ans[i, ] <- x
          i <- i + 1
        }
      }
      
      # Next random direction
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
      iter <- iter + 1
      
    }
    
    out <- list( S = ans,
                 # U = U,
                 B_lo = B_lo,
                 B_up = B_up )
    
    return( out ) 
  }
