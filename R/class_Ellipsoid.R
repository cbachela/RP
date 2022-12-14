  
  
  ############################################################################
  ### RP PACKAGE - ELLIPSOID CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.01.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # Ellipsoid
  # Ellipsoid.dimension
  # Ellipsoid.isEmpty
  # Ellipsoid.isInterior
  # Ellipsoid.multiplication
  # Ellipsoid.addition
  # Ellipsoid.boundaryPoint
  # Ellipsoid.boundaryDistance
  # Ellipsoid.transform
  # Ellipsoid.sample
 
  # EllipsoidSampler
  # EllipsoidSampler.setCtrl
  # EllipsoidSampler.transform
  # EllipsoidSampler.sample
  
 
  
  
  
  # --------------------------------------------------------------------------
  Ellipsoid <- setRefClass( Class = "Ellipsoid",
                            fields = list( centre = "vector",
                                           shape = "matrix",
                                           # linear = "vector",
                                           # constant = "numeric",
                                           rhs = "numeric" ) )
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.dimension <- function()
  {
    ans = 0
    if ( !.self$isEmpty() ) {
      ans <- ncol( shape )
    }
    return( ans )
  }
  Ellipsoid$methods( dimension = Ellipsoid.dimension )
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.isEmpty <- function()
  {
    ans <- FALSE
    if ( any(is.na(shape)) ) {
      ans <- TRUE
    }
    if ( any(is.na(centre)) ) {
      ans <- TRUE
    }
    return( ans )
  }
  Ellipsoid$methods( isEmpty = Ellipsoid.isEmpty )
  
  # --------------------------------------------------------------------------
  Ellipsoid.isInterior <- function(x, eps = 1e-08)
  {
    if ( !.self$isEmpty() ) {
      # if ( rhs - ( t(x - centre) %*% shape %*% (x - centre) + t(linear) %*% x + constant) < -eps ) {
      if ( rhs - t(x - centre) %*% shape %*% (x - centre) < -eps ) {
        return( FALSE )
      }
    }
    return( TRUE )
  }
  Ellipsoid$methods( isInterior = Ellipsoid.isInterior )
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.multiplication <- function( Amat = NULL )
  {
    # A * E(q, Q) = E(Aq, AQA^T)
    
    if ( !is.null(Amat) ) {
      
      centre_new <- Amat %*% centre
      shape_new <- Amat %*% shape %*% t(Amat)
      shape_new <- 0.5 * (shape_new + t(shape_new))
      
      centre <<- centre_new
      shape <<- shape_new
      return( TRUE )
      
    } else {
      warning("Multiplication cannot be performed 
              since no transformation matrix 'tmat' was provided. ")
    }
  }
  Ellipsoid$methods( multiplication = Ellipsoid.multiplication )
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.addition <- function( bvec = NULL )
  {
    # Operation E + b where E - ellipsoid and b - vector in R^n.
    # E(q, Q) + b = E(q + b, Q)
    
    if ( !is.null(bvec) ) {
      centre <<- centre + bvec
      return( TRUE )
    } else {
      warning("Addition cannot be performed 
              since no transformation vector 'bvec' was provided. ")
    }
  }
  Ellipsoid$methods( addition = Ellipsoid.addition )
  
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.boundaryPoint <- function( u = NULL )
  {
    # Finds a random point on the boundary of the ellipsoid. 
    # Most likely, this point will not be in the intersection with other bodies.
    
    if ( is.null(u) ) {
      u <- rnorm(.self$dimension())
      u <- u / as.numeric(sqrt(sum(u^2)))
    }
    lambda <- .self$boundaryDistance( x = centre, 
                                      u = u )
    p <- centre + lambda[2] * u
    
    # Check
    boundary_check <- t(p) %*% shape %*% p - rhs
    stopifnot( boundary_check < 1e-16 )
    
    return( p )
  }
  Ellipsoid$methods( boundaryPoint = Ellipsoid.boundaryPoint )
  
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.boundaryDistance <- function( x, u )
  {
    if ( !.self$isEmpty() ) {
      
      y <- (x - centre)
      a <- t(u) %*% shape %*% u
      # b <- 2 * t(u) %*% shape %*% y + t(linear) %*% u
      # g <- t(y) %*% shape %*% y + t(linear) %*% y + constant - rhs
      b <- 2 * t(u) %*% shape %*% y
      g <- t(y) %*% shape %*% y - rhs
      d <- sqrt( b^2 - 4 * a * g )
      lambda_minus <- (-b - d) / (2 * a)
      lambda_plus <- (-b + d) / (2 * a)
      lambda <- sort( c(lambda_minus, lambda_plus) )
      
    } else {
      lambda <- NULL
    }
    
    return( lambda )
  }
  Ellipsoid$methods( boundaryDistance = Ellipsoid.boundaryDistance )
  
  
 
  
  # --------------------------------------------------------------------------
  Ellipsoid.transform <- function( tmat = NULL, 
                                   z = NULL, 
                                   boundary_point = NULL )
  {
    
    if ( is.null(tmat) || is.null(z) ) {
      transformation <- simplexTransform(d = .self$dimension() )
      tmat <- transformation$tmat
      z <- transformation$z
    }

    # Transform ellipsoid
    E_prime <- .self$copy()
    E_prime$multiplication( Amat = tmat )
    E_prime$addition( bvec = -z )
    
    # Reduce dimension
    q_tmat <- E_prime$centre
    M <- E_prime$shape
    M <- 0.5 * (M + t(M))
    M_bar <- M[-1, -1]
    m_bar <- M[-1, 1]
    m11 <- M[1, 1]
    h <- q_tmat[1]^2 * (m11 * t(m_bar) %*% M_bar %*% m_bar )
    Q_prime <- as.numeric(1 - h) * M_bar
    q_prime <- as.numeric( q_tmat + q_tmat[1] * t(c(-1, pseudoinverse(M_bar) %*% m_bar)) )[-1]
    E_prime$centre <- q_prime
    E_prime$shape <- Q_prime
    
    # Transform threshold 
    if ( is.null(boundary_point) ) {
      # boundary_point <<- getBoundaryPoint(E)   # ~~~~~~~~~~~~~~~~ change this later!
      stop("boundary point required.")
    }
    v <- as.numeric( tmat %*% boundary_point - z )[-1]
    E_prime$rhs <- as.numeric( t(v - q_prime) %*% Q_prime %*% (v - q_prime) )
   
    return( E_prime )
  }
  Ellipsoid$methods( transform = Ellipsoid.transform )
  
  
  # --------------------------------------------------------------------------
  Ellipsoid.sample <- function( n_sim = 10^1,
                                algo = c("sphere", "hitnrun"),
                                n_thin = 0,
                                x0 = NULL,
                                b_pushy = FALSE )
  {
    algo <- match.arg(algo)
    if ( algo == "sphere" ) {
      samples <- runifE( S = shape,
                         z_hat = centre,
                         gamma_threshold = rhs,
                         n_points = n_sim + n_thin,
                         b_pushy = b_pushy )
    } else if ( algo == "hitnrun" ) {
      if ( is.null(x0) ) {
        x0 <- interior_point
      }
      samples <- har( K = .self,
                      x = x0,
                      u = NULL,
                      n_sim = n_sim + n_thin )

    } else {
      stop("algorithm '", algo, "' not implemented.")
    }
    if ( !is.list(samples) ) {
      samples <- list(samples)
    }

    # Remove n_thin points from the samples (thinning)
    if ( n_thin > 0 ) {
      for ( i in 1:length(samples) ) {
        S <- samples[[i]]
        if ( nrow(S) > n_thin ) {
          idx <- floor(seq(from = 1, to = nrow(S), length.out = nrow(S) - n_thin))
          samples[[i]] <- S[idx, ]
        }
      }
    }

    return( samples )
  }
  
  

  
  
  
  
  # --------------------------------------------------------------------------
  EllipsoidSampler <- setRefClass( Class = "EllipsoidSampler",
                                   fields = list( ellipsoid = "Ellipsoid",
                                                  transformation = "list",
                                                  spec = "list",
                                                  samples = "list" ) )
  
  
  # --------------------------------------------------------------------------
  EllipsoidSampler.setCtrl <- function( n_sim = 10^2,
                                        thin = 1,
                                        # burn = 0,
                                        algo = c("sphere",
                                                 "hitnrun", 
                                                 "billiard", 
                                                 "volesti"),
                                        volesti = list( density = "uniform", 
                                                        walk = "BiW" ),
                                        jmp = NULL,
                                        b_pushy = TRUE,
                                        interior_point = NULL,
                                        boundary_point = NULL,
                                        ellipsis = list(),
                                        verbose = TRUE )
  {
    ans <- list()
    ans$n_sim <- n_sim
    ans$thin <- thin
    ans$algo <- match.arg(algo)
    ans$volesti <- volesti
    ans$jmp <- jmp
    ans$b_pushy <- b_pushy
    ans$interior_point <- interior_point
    ans$boundary_point <- boundary_point
    ans$verbose <- verbose
    spec <<- ans
    
    return( TRUE )
  }
  EllipsoidSampler$methods( setCtrl = EllipsoidSampler.setCtrl )
  
  
  
  
  # --------------------------------------------------------------------------
  EllipsoidSampler.transform <- function()
  {
    d <- ellipsoid$dimension()
    
    # Extract or compute transformation
    if ( is.null(transformation$tmat_inv) ) {
      transformation <<- simplexTransform( d = d )
    }
    tmat <- transformation$tmat
    z <- transformation$z
    
    # Create new (d-1) dimensional ellipsoid
    ellipsoid <<- ellipsoid$transform( tmat = tmat, 
                                       z = z,
                                       boundary_point = spec$boundary_point )
    
    # Transform initial and boundary point if available
    if ( !is.null(spec$interior_point) ) {
      as.numeric( tmat %*% spec$boundary_point - z )[-1]
      spec$interior_point <<- (tmat %*% spec$interior_point - z)[-1]
    }
    if ( !is.null(spec$boundary_point) ) {
      spec$boundary_point <<- (tmat %*% spec$boundary_point - z)[-1]
    }
    
    return( TRUE )
  }
  EllipsoidSampler$methods( transform = EllipsoidSampler.transform )
  
  
  
  # --------------------------------------------------------------------------
  EllipsoidSampler.sample <- function()
  {
    if ( spec$algo == "sphere" ) {
      samples$S <<- runifE( S = ellipsoid$shape,
                            z_hat = ellipsoid$centre,
                            gamma_threshold = ellipsoid$rhs,
                            n_points = spec$n_sim + spec$thin,
                            b_pushy = spec$b_pushy )
    } else if ( spec$algo == "billiard" ) {
      samples <<- billiard( K = ellipsoid,
                            x = spec$interior_point,
                            u = NULL,
                            n_sim = spec$n_sim + spec$thin,
                            jmp = spec$jmp )
      
    } else if ( spec$algo == "hitnrun" ) {
      samples <<- hitnrun( K = ellipsoid, 
                           x = spec$interior_point, 
                           u = NULL, 
                           n_sim = spec$n_sim + spec$thin )
    } else {
      stop("algorithm '", spec$algo, "' not implemented.")
    }
    
    # .self$thin()
    
    return( TRUE )
  }
  EllipsoidSampler$methods( sample = EllipsoidSampler.sample )
  
  
  
  
 