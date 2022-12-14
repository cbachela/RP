  
  
  ############################################################################
  ### RP PACKAGE - Polytope CLASS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.01.2020
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # Polytope
  # Polytope.dimension
  # Polytope.isEmpty
  # Polytope.isInterior
  # Polytope.boundaryDistance
  # Polytope.intersectPolytope
  # Polytope.transform
  # Polytope.transformBack
  # Polytope.sample

  
  # PolytopeSampler
  # PolytopeSampler.setCtrl
  # PolytopeSampler.transform
  # PolytopeSampler.transformBack
  # PolytopeSampler.sample

  
  
  
  
  # --------------------------------------------------------------------------
  Polytope <- setRefClass( Class = "Polytope",
                           fields = list( A = "matrix",
                                          b = "vector",
                                          sense = "character" ) )
  # interior_point = "numeric",
  # boundary_point = "numeric" ) )
  
  
  # --------------------------------------------------------------------------
  Polytope.dimension <- function()
  {
    ans <- 0
    if ( !.self$isEmpty() ) {
      ans <- ncol(A)
    }
    return( ans )
  }
  Polytope$methods( dimension = Polytope.dimension )
  
  
  
  # --------------------------------------------------------------------------
  Polytope.isEmpty <- function()
  {
    ans <- FALSE
    if ( any(is.na(A)) ) {
      ans <- TRUE
    }
    if ( any(is.na(b)) ) {
      ans <- TRUE
    }
    return( ans )
  }
  Polytope$methods( isEmpty = Polytope.isEmpty )
  
  
  
  # --------------------------------------------------------------------------
  Polytope.isInterior <- function( x, eps = 1e-08 )
  {
    if ( !.self$isEmpty() ) {
      if ( any(  b - A %*% x < -eps) ) {
        return( FALSE )
      }
    }
    return( TRUE )
  }
  Polytope$methods( isInterior = Polytope.isInterior )
  
  
  # --------------------------------------------------------------------------
  Polytope.boundaryDistance <- function( x, u )
  {
    if ( !.self$isEmpty() ) {
        
      num <- b - A %*% x
      denom <- A %*% u
      ratio <- num / denom
      lambda <- c(0, 0)
      # idx_neg <- which(denom < 0)
      idx_neg <- which(ratio < 0)
      if ( length(idx_neg) >= 1 ) {
        lambda[1] <- max( c(ratio[idx_neg], -1e+300) )
      } 
      # idx_pos <- which(denom > 0)
      idx_pos <- which(ratio > 0)
      if ( length(idx_pos) >= 1 ) {
        lambda[2] <- min( c(ratio[idx_pos], 1e+300) )
      }  
      
    } else {
      lambda <- NULL
    }
    return( lambda )
  }
  Polytope$methods( boundaryDistance = Polytope.boundaryDistance )
  
  
  
  # --------------------------------------------------------------------------
  Polytope.intersectPolytope <- function( P )
  {
    stopifnot( .self$dimension() == P$dimension() )
    
    A <<- rbind( A, P$A )
    colnames(A) <<- colnames(A)
    b <<- c(b, P$b)
    sense <<- c(sense, P$sense)
    
    return( TRUE )
  }
  Polytope$methods( intersectPolytope = Polytope.intersectPolytope )
  
  
  
  # --------------------------------------------------------------------------
  Polytope.transform <- function( tmat_inv = NULL, z = NULL, idx_eq = NULL )
  {
    A_prime <- (A %*% tmat_inv)[-idx_eq, -1]
    b_prime <- b[-idx_eq] - as.numeric( (A[-idx_eq, ] %*% tmat_inv) %*% z )
    sense_prime <- sense[-idx_eq]
    
    # Instantiate a new object
    P_prime <- Polytope$new()
    P_prime$A <- A_prime
    P_prime$b <- b_prime
    P_prime$sense <- sense_prime
   
    return( P_prime )
  }
  Polytope$methods( transform = Polytope.transform )
  
  
  # # --------------------------------------------------------------------------
  # Polytope.transformBack <- function( tmat_inv = NULL, z = NULL, idx_eq = NULL )
  # {
  #  
  #   A <<- tmat_inv %*% cbind( A, 0) + z
  #   b <<- tmat_inv %*% b + z
  #   sense <<- "xxx"
  #   
  #   return( TRUE )
  # }
  # Polytope$methods( transformBack = Polytope.transformBack )
  
  
  # --------------------------------------------------------------------------
  Polytope.sample <- function( n_sim = 10^1,
                               algo = c("billiard", "hitnrun"),
                               n_thin = 0,
                               x0 = NULL,
                               jmp = NULL,
                               jmp_adaptive = FALSE )
  {
    algo <- match.arg(algo)
    if ( is.null(x0) ) {
      # x0 <- .self$interiorPoint()
    }
    if ( algo == "hitnrun" ) {
      
      samples <- har( K = .self,
                      x = x0,
                      u = NULL,
                      n_sim = n_sim + n_thin )

    } else if ( algo == "billiard" ) {
      
      samples <<- billiard( K = .self,
                            x = x0,
                            u = NULL,
                            jmp = jmp,
                            jmp_adaptive = jmp_adaptive,
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
  Polytope$methods( sample = Polytope.sample )
  
  
  
  
  

  
  
  # --------------------------------------------------------------------------
  PolytopeSampler <- setRefClass( Class = "PolytopeSampler",
                                  fields = list( polytope = "Polytope",
                                                 transformation = "list",
                                                 spec = "list",
                                                 samples = "list" ) )
  
  
  # --------------------------------------------------------------------------
  PolytopeSampler.setCtrl <- function( n_sim = 10^2,
                                       thin = 1,
                                       # burn = 0,
                                       algo = c("hitnrun", 
                                                "billiard", 
                                                "volesti"),
                                       volesti = list( density = "uniform", 
                                                       walk = "BiW" ),
                                       jmp = NULL,
                                       jmp_adaptive = FALSE,
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
    ans$jmp_adaptive <- jmp_adaptive
    ans$interior_point <- interior_point
    ans$boundary_point <- boundary_point
    ans$verbose <- verbose
    spec <<- ans
    
    return( TRUE )
  }
  PolytopeSampler$methods( setCtrl = PolytopeSampler.setCtrl )
  
  
  
  
  # --------------------------------------------------------------------------
  PolytopeSampler.transform <- function()
  {
    d <- polytope$dimension()
    idx_eq <- which(polytope$sense == "=")
    
    if ( length(idx_eq) < 1 ) {
      
      if ( isTRUE(spec$verbose) ) {
        warning("No transformation is performed since the 
                linear system contains no equality.")
      }
      
    } else if ( length(idx_eq) == 1 ) {
      
      # Extract or compute transformation
      if ( is.null(transformation$tmat_inv) ) {
        transformation <<- simplexTransform( d = d )
      }
      tmat_inv <- transformation$tmat_inv
      z <- transformation$z
      
      # Create new (d-1) dimensional polytope
      polytope <<- polytope$transform( tmat_inv = tmat_inv, 
                                       z = z, 
                                       idx_eq = idx_eq )
      
      # Transform initial and boundary point if available
      if ( !is.null(spec$interior_point) ) {
        spec$interior_point <<- (transformation$tmat %*% spec$interior_point - z)[-1]
      }
      if ( !is.null(spec$boundary_point) ) {
        spec$boundary_point <<- (transformation$tmat %*% spec$boundary_point - z)[-1]
      }
      
    } else {
      stop("cannot deal with more than one equality.")
    }
    
    return( TRUE )
  }
  PolytopeSampler$methods( transform = PolytopeSampler.transform )
  
  
  # --------------------------------------------------------------------------
  PolytopeSampler.transformBack <- function()
  {
    
    # Extract transformation variables
    tmat_inv <- transformation$tmat_inv
    z <- transformation$z
    
    # # Back-transform the polytope
    # polytope <<- polytope$transformBack( tmat_inv = tmat_inv, 
    #                                      z = z, 
    #                                      idx_eq = idx_eq )
    
    # Back-transform samples
    FUN <- function(mat) { t(apply( mat, 1, function(y) { tmat_inv %*% (c(0, y) + z) } ) ) }
    samples <<- lapply( samples, FUN = FUN )
    
    # Back-transform spec variables
    if ( !is.null(spec$interior_point) ) {
      spec$interior_point <<-  tmat_inv %*% (c(0, spec$interior_point) + z)
    }
    if ( !is.null(spec$boundary_point) ) {
      spec$boundary_point <<-  tmat_inv %*% (c(0, spec$boundary_point) + z)
    }
    
    return( TRUE )
  }
  PolytopeSampler$methods( transformBack = PolytopeSampler.transformBack )
  
  
  # --------------------------------------------------------------------------
  PolytopeSampler.sample <- function()
  {
    if ( spec$algo %in% c("hitnrun", "billiard") ) {

      samples <<- polytope$sample( n_sim = spec$n_sim,
                                   algo = spec$algo,
                                   n_thin = spec$thin,
                                   x0 = spec$interior_point,
                                   jmp = spec$jmp,
                                   jmp_adaptive = spec$jmp_adaptive )
      
    } else if ( spec$algo == "volesti" ) {
      
      ###  ~~~~~~~~~~~~~~~~~~~~~
      require(volesti)
      P <- Hpolytope$new( A = polytope$A, 
                          b = polytope$b )
      
      lRW <- list( walk = spec$volesti$walk )
      if ( !is.null(spec$interior_point) ) {
        lRW$starting_point <- spec$interior_point
      }
      
      # tmp <- inner_ball( P = P )
      # tmp <- getInteriorPoint( object = RP_prime )
      # debugonce( isInterior) 
      # isInterior(RP_prime, tmp)
      
      samples_tmp <- sample_points( P = P,
                                    n = spec$n_sim + spec$thin,
                                    distribution = list( density = spec$volesti$density ),
                                    random_walk = lRW )
      samples <<- list( S = t(samples_tmp) )
      ### ~~~~~~~~~~~~~~~~~~~~~
      
    } else {
      stop("algorithm '", algo, "' not implemented.")
    }
    if ( !is.list(samples) ) {
      samples <<- list(samples)
    }
    
    # .self$thin()
    
    return( TRUE )
  }
  PolytopeSampler$methods( sample = PolytopeSampler.sample )
  
  
  
  
  
  
  
  
  
 