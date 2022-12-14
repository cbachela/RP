  
  
  ############################################################################
  ### RP PACKAGE - RPE CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.01.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # RPE
  # transformRPE
  # sampleRPE
  
  
 
  
  # --------------------------------------------------------------------------
  transformRPE <- function()
  {
    spec_prime <- spec
    
    if ( is.null(spec$transformation$tmat) ) {
      spec$transformation <<- simplexTransform(d = ellipsoid$dimension())
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
    
    # Transform starting point
    if ( is.null(x0) ) {
      # x0 <- getInitialPoint(E)            # ~~~~~~~~~~~~~~~~ change this later!
      stop("initial point required.")
    }
    spec_prime$x0 <- (tmat %*% x0 - z)[-1]
    

    K_prime <- ConvexBody$new(ellipsoid = Ellipsoid$new(centre = q_prime,
                                                    shape = Q_prime,
                                                    rhs = b_prime),
                              spec = spec_prime)
    
    return( K_prime )
  }
  
  # --------------------------------------------------------------------------
  sampleRPE <- function( )
  {
    
    if ( spec$algo == "sphere" ) {
      samples <<- runifE( S = ellipsoid$shape,
                          z_hat = ellipsoid$centre,
                          gamma_threshold = ellipsoid$rhs,
                          n_points = spec$n_sim + spec$thin,
                          b_pushy = spec$b_pushy )
      
    } else if ( spec$algo == "hitnrun" ) {
      samples <<- hitnrun( B = ellipsoid, 
                           x = spec$x0, 
                           u = NULL, 
                           n_sim = spec$n_sim + spec$thin )
      
    } else {
      stop("algorithm '", algo, "' not implemented.")
    }
    if ( !is.list(samples) ) {
      samples <<- list(samples)
    }
    
    .self$thin()
    
    return( TRUE )
  }
  

  
  
  
  # --------------------------------------------------------------------------
  RPE <- setRefClass( Class = "RPE",
                      fields = list(ellipsoid = "Ellipsoid",
                                    spec = "list",
                                    samples = "list"),
                      methods = list(transform = transformRPE,
                                     sample = sampleRPE) )
  
  

  # --------------------------------------------------------------------------
  ConvexBody <- setRefClass( Class = "ConvexBody",
                             fields = list(ellipsoid = "Ellipsoid",
                                           spec = "list",
                                           samples = "list"),
                             methods = list(transform = transformRPE,
                                            sample = sampleRPE) )
  
  
  
  
  
