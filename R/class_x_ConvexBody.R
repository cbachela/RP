  
  
  ############################################################################
  ### RP PACKAGE - ConvexBody CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.01.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # ConvexBody
  # ConvexBody.transform
  # ConvexBody.transformPolytope
  # ConvexBody.transformEllipsoid

 
  
  
  
  
  # --------------------------------------------------------------------------
  ConvexBody <- setRefClass( Class = "ConvexBody",
                             fields = list( polytope = "Polytope",
                                            ellipsoid = "Ellipsoid",
                                            spec = "list",
                                            samples = "list" ) )
  
  
  # --------------------------------------------------------------------------
  ConvexBody.transform <- function()
  {
    B_prime <- ConvexBody$new()
    B_prime$polytope <- .self$transformPolytope()
    B_prime$ellipsoid <- .self$transformEllipsoid()
    return( B_prime )
  }
  ConvexBody$methods( transform = ConvexBody.transform )
  
  
  # --------------------------------------------------------------------------
  ConvexBody.transformPolytope <- function( )
  {
    
    d <- polytope$dimension()
    idx_eq <- which(polytope$sense == "=")
    
    if ( length(idx_eq) < 1 ) {
      
      if ( isTRUE(spec$verbose) ) {
        warning("No transformation is performed since the 
                linear system contains no equality.")
        P_prime <- polytope
      }
      
      } else if ( length(idx_eq) == 1 ) {
        
        # Extract or compute transformation
        if ( is.null(spec$transformation$tmat) ) {
          spec$transformation <<- simplexTransform( d = d )
        }
        tmat <- spec$transformation$tmat
        tmat_inv <- spec$transformation$tmat_inv
        z <- spec$transformation$z
        
        # Create new (d-1) dimensional polytope
        A_prime <- (polytope$A %*% tmat_inv)[-idx_eq, -1]
        b_prime <- polytope$b[-idx_eq] - as.numeric( (polytope$A[-idx_eq, ] %*% tmat_inv) %*% z )
        sense_prime <- polytope$sense[-idx_eq]
        
        P_prime <- Polytope$new()
        P_prime$A <- A_prime
        P_prime$b <- b_prime
        P_prime$sense <- sense_prime
      
      } else {
        stop("cannot deal with more than one equality.")
      }
    
    return( P_prime )
  }
  ConvexBody$methods( transformPolytope = ConvexBody.transformPolytope )
  
  
  
  
  # --------------------------------------------------------------------------
  ConvexBody.transformEllipsoid <- function()
  {
    spec_prime <- spec
    
    if ( is.null(spec$transformation$tmat) ) {
      spec$transformation <<- simplexTransform( d = ellipsoid$dimension() )
    }
    tmat <- spec$transformation$tmat 
    tmat_inv <- spec$transformation$tmat_inv
    z <- spec$transformation$z
    x0 <- spec$x0
    x_boundary <- spec$x_boundary
    
    # Checks
    if ( is.null(tmat) ) stop("Transformation matrix 'tmat' is needed.")
    stopifnot( dim(tmat)[2] == ellipsoid$dimension() ) 
    if ( is.null(tmat_inv) ) {
      tmat_inv <- inverse(m = tmat)
    }
    if ( is.null(z) ) stop("Variable 'z' cannot be NULL.")
    # if ( is.null(x_boundary) ) stop("Variable 'x_boundary' cannot be NULL.")
    
    # Transform ellipsoid
    ellipsoid$multiplication( Amat = tmat )
    ellipsoid$addition( bvec = -z )
    
    # Reduce dimension
    q_tmat <- ellipsoid$centre
    M <- ellipsoid$shape
    M <- 0.5 * (M + t(M))
    M_bar <- M[-1, -1]
    m_bar <- M[-1, 1]
    m11 <- M[1, 1]
    h <- q_tmat[1]^2 * (m11 * t(m_bar) %*% M_bar %*% m_bar )
    Q_prime <- as.numeric(1 - h) * M_bar
    q_prime <- as.numeric( q_tmat + q_tmat[1] * t(c(-1, pseudoinverse(M_bar) %*% m_bar)) )[-1]
    
    # Transform threshold 
    if ( is.null(x_boundary) ) {
      # x_boundary <- getBoundaryPoint(E)   # ~~~~~~~~~~~~~~~~ change this later!
      stop("boundary point required.")
    }
    v <- as.numeric( tmat %*% x_boundary - z )[-1]
    b_prime <- as.numeric( t(v - q_prime) %*% Q_prime %*% (v - q_prime) )
    spec_prime$x_boundary <- v
    
    # # Transform starting point
    # if ( is.null(x0) ) {
    #   # x0 <- getInitialPoint(E)            # ~~~~~~~~~~~~~~~~ change this later!
    #   stop("initial point required.")
    # }
    # spec_prime$x0 <- (tmat %*% x0 - z)[-1]
    
    
    K_prime <- ConvexBody$new( ellipsoid = Ellipsoid$new(centre = q_prime,
                                                         shape = Q_prime,
                                                         rhs = b_prime ),
                               spec = spec_prime )
    
    return( K_prime )
  }
  ConvexBody$methods( transformEllipsoid = ConvexBody.transformEllipsoid )
  
  
  
  # # --------------------------------------------------------------------------
  # Body.thin <- function()
  # {
  #   thin <- ifelse( is.null(getSpec(object)$thin), 0, getSpec(object)$thin )
  #   if ( thin > 0 ) {
  #     lSamples <- getSamples(object, "all")
  #     for ( i in 1:length(lSamples) ) {
  #       samples <- lSamples[[i]]
  #       if ( nrow(samples) > thin ) {
  #         idx <- floor(seq(from = 1, to = nrow(samples), length.out = nrow(samples) - thin))
  #         lSamples[[i]] <- samples[idx, ]
  #       }
  #     }
  #     setSamples(object) <- lSamples
  #   }
  #   return( object )
  # }
  # 
  
  
  