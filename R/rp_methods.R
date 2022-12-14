
  
  ############################################################################
  ### RP PACKAGE - METHODS 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  #
  # METHODS:
  #
  # rp
  # transform
  # backTransform
  # thin
  # isEmpty
  # dimension
  # dists2Body
  # intersect
  
  
  
  
  # --------------------------------------------------------------------------
  #' rp
  #' @description Random sampling from object
  #' @param object    [CONVEXBODY] Any class that inherits from CONVEXBODY.
  #' @param spec      [list]  A list of specifications passed to lower level functions.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname rp-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "rp", 
             def = function(object, spec) { standardGeneric("rp") },
             package = "RP")
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname rp-methods
  #' @aliases rp, SIMPLEX-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod( f = "rp", 
             signature(object = "SIMPLEX"), 
             definition = function(object, spec) 
             { 
               .rpS( object = object, n_sim = spec$n_sim, thin = spec$thin )
             })
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .rpS <- function( object, n_sim = 10, thin = 0 )
  {
    samples <- runifS( d = dimension(object), 
                       n_sim = n_sim )
    samples <- samples * getRHS(object)[which(getSense(object) == "=")]
    ans <- new( "RPP",
                samples = list(samples),
                spec = list(n_sim = n_sim + thin),
                object )   # inherits class 'VBODY'
    ans <- thin(ans)
    return( ans )
  }
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname rp-methods
  #' @aliases rp, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod( f = "rp", 
             signature(object = "ELLIPSOID"), 
             definition = function(object, spec) 
             { 
               .rpE( object = object, spec = spec )
             })
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .rpE <- function( object, spec )
  {

    if ( spec$algo == "sphere" ) {
      samples <- runifE( S = getShape(object),
                         z_hat = getCentre(object),
                         gamma_threshold = getRHS(object),
                         n_points = spec$n_sim + spec$thin,
                         b_pushy = spec$b_pushy )
      
    } else if ( spec$algo == "hitnrun" ) {
      samples <- hitnrun( B = object, 
                          x = spec$x0, 
                          u = NULL, 
                          n_sim = spec$n_sim + spec$thin )
      
    } else {
      stop("algorithm '", algo, "' not implemented.")
    }
    if ( !is.list(samples) ) {
      samples <- list(samples)
    }
    ans <- new( "RPE",
                samples = samples,
                spec = spec,
                object )
    ans <- thin(ans)
    return( ans )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname rp-methods
  #' @aliases rp, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod( f = "rp",
             signature = "POLYTOPE",
             # definition = .rpP )
             definition = function(object, spec) 
             { 
               .rpP( object = object, spec = spec )
             })
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .rpP <- function( object, spec = rpCtrl() )
  {
    
    # Instantiate RP object
    if ( !inherits(object, "RP") ) {
      object <- new( "RPP",
                     spec = spec,
                     object )
    }
    
    # Tranform RP object to sample space 
    # (If there is no simplex condition then no transformation is done)
    RP_prime <- transform(object)
    spec_prime <- getSpec(RP_prime)
    
    # Run MCMC procedure
    if ( spec$algo == "hitnrun" ) {
      samples_prime <- hitnrun(B = RP_prime,
                               x = spec_prime$x0,
                               u = NULL,
                               n_sim = spec_prime$n_sim + spec_prime$thin)
      
    } else if ( spec$algo == "billiard" ) {
      samples_prime <- billiard(B = RP_prime,
                                x = spec_prime$x0,
                                u = NULL,
                                n_sim = spec_prime$n_sim + spec_prime$thin,
                                jmp = spec_prime$jmp)
    } else if ( spec$algo == "volesti" ) {
      
      ###  ~~~~~~~~~~~~~~~~~~~~~
      require(volesti)
      P <- Hpolytope$new( A = getAmat(RP_prime), 
                          b = getRHS(RP_prime) )
      # samples_tmp <- sample_points(P = P, 
      #                              N = spec_prime$n_sim + spec_prime$thin, 
      #                              distribution = "uniform",
      #                              InnerPoint = spec_prime$x0,
      #                              WalkType = "CDHR")
      lRW <- list( walk = spec_prime$volesti$walk )
      if ( !is.null(spec_prime$x0) ) {
        lRW$starting_point <- spec_prime$x0
      }
      
      tmp <- inner_ball( P = P )
      tmp <- getInteriorPoint( object = RP_prime )
      debugonce( isInterior) 
      isInterior(RP_prime, tmp)
      
      
      samples_tmp <- sample_points( P = P,
                                    n = spec_prime$n_sim + spec_prime$thin,
                                    distribution = list( density = spec_prime$volesti$density ),
                                    random_walk = lRW )
      samples_prime <- list( S = t(samples_tmp) )
      ### ~~~~~~~~~~~~~~~~~~~~~
      
    }
    setSamples(RP_prime) <- samples_prime
    
    # Transform back to variable space
    if ( dimension(RP_prime) < dimension(object) ) {
      RP <- backTransform(RP_prime)
    } else {
      RP <- RP_prime
    }
    RP <- thin(RP)
    
    return( RP )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname rp-methods
  #' @aliases rp, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod( f = "rp",
             signature = "CONVEXBODY",
             definition = function(object, spec) 
             { 
               .rpB( object = object, spec = spec ) 
             })
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .rpB <- function( object, spec )
  {
    
    # Instantiate RP object
    if ( !inherits(object, "RP") ) {
      object <- new( "RPB",
                     spec = spec,
                     object )
    }
    
    # Tranform RP object to sample space 
    # (If there is no simplex condition then no transformation is done)
    RP_prime <- transform(object)
    spec_prime <- getSpec(RP_prime)
    
    # Run MCMC procedure
    if ( spec$algo == "hitnrun" ) {
      
      samples_prime <- hitnrun(B = RP_prime,
                               x = spec_prime$x0,
                               u = NULL,
                               n_sim = spec_prime$n_sim + spec_prime$thin)
      # samples_prime <- list(S = samples_prime$S[-c(1:spec_prime$thin), ],
      #                       B_lo = samples_prime$B_lo[-c(1:spec_prime$thin), ],
      #                       B_up = samples_prime$B_up[-c(1:spec_prime$thin), ])
                            
      
    } else {
      stop("algorithm '", algo, "' not yet implemented.")
    }
  
    setSamples(RP_prime) <- samples_prime
    
    # Transform back to variable space
    if ( dimension(RP_prime) < dimension(object) ) {
      RP <- backTransform(RP_prime)
    } else {
      RP <- RP_prime
    }
    RP <- thin(RP)
    
    return( RP )
  }
  
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' transform
  #' @description Perrforms a (isometric, i.e. distance preserving) transformation
  #'              of body from d-demensional space to (d-1)-demensional space.
  #'              Note that the body needs to be a simplex or a the intersection of a 
  #'              body with a simplex (i.e. a constrained simplex).
  #' @param object [S4 Class] Object of class \code{VBODY}, \code{RP}.
  #' @return Returns an object of the same class as the imput 
  #'         but transformed to (d-1)-dimensional space.
  #' @docType methods
  #' @rdname transform-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "transform", 
             def = function(object) { standardGeneric("transform") },
             package = "RP")
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname transform-methods
  #' @aliases transform, RPP-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "transform", 
            signature = "RPP",
            definition = function(object)
            {
              ans <- .transformP(P = object)
              return( ans )
            })
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .transformP <- function(P)
  {
    
    d <- dimension(P)
    A <- getAmat(P)
    b <- getRHS(P)
    idx_eq <- which(getSense(P) == "=")
    
    if ( length(idx_eq) < 1 ) {
      
      if ( isTRUE(getSpec(P)$verbose) ) {
        warning("No transformation is performed since the 
                linear system contains no equality.")
      }
      
    } else if ( length(idx_eq) == 1 ) {
      
      # Extract or compute transformation
      spec <- getSpec(P)
      if ( is.null(spec$transformation$tmat) ) {
        spec$transformation <- simplexTransform(d = dimension(P))
      }
      tmat <- spec$transformation$tmat
      tmat_inv <- spec$transformation$tmat_inv
      z <- spec$transformation$z
      
      # Create new (d-1) dimensional polytope
      A_new <- (A %*% tmat_inv)[-idx_eq, -1]
      b_new <- b[-idx_eq] - as.numeric( (A[-idx_eq, ] %*% tmat_inv) %*% z )
      P_prime <- polytope(A = A_new,
                          b = b_new,
                          sense = getSense(P)[-idx_eq])
      
      # Transform starting point if available
      if ( !is.null(spec$x0) ) {
        spec$x0 <- (tmat %*% spec$x0 - z)[-1]
      }
      
      setPolytope(P) <- P_prime
      setSpec(P) <- spec
      
    } else {
      stop("cannot deal with more than one equality.")
    }
    
    return( P )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname transform-methods
  #' @aliases transform, RPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "transform", 
            signature = "RPE",
            definition = function(object)
            {
              ans <- .transformE(E = object)
              return( ans )
            })
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .transformE <- function(E)
  {
    
    spec <- getSpec(E)
    if ( is.null(spec$transformation$tmat) ) {
      spec$transformation <- simplexTransform(d = dimension(E))
    }
    tmat <- spec$transformation$tmat 
    tmat_inv <- spec$transformation$tmat_inv
    z <- spec$transformation$z
    x0 <- spec$x0
    x_boundary <- spec$x_boundary
    
    # Checks
    if ( is.null(tmat) ) stop("Transformation matrix 'tmat' is needed.")
    stopifnot( dim(tmat)[2] == dimension(E) ) 
    if ( is.null(tmat_inv) ) {
      tmat_inv <- inverse(m = tmat)
    }
    if ( is.null(z) ) stop("Variable 'z' cannot be NULL.")
    # if ( is.null(x_boundary) ) stop("Variable 'x_boundary' cannot be NULL.")
    
    # Transform ellipsoid
    E_tmat <- eMult(E = E, Amat = tmat)
    E_tmat <- ePlus(E = E_tmat, bvec = -z)
    
    # Reduce dimension
    q_tmat <- getCentre(E_tmat)
    M <- getShape(E_tmat)
    M <- 0.5 * (M + t(M))
    M_bar <- M[-1, -1]
    m_bar <- M[-1, 1]
    m11 <- M[1, 1]
    h <- q_tmat[1]^2 * (m11 * t(m_bar) %*% M_bar %*% m_bar )
    Q_prime <- as.numeric(1 - h) * M_bar
    q_prime <- as.numeric( q_tmat + q_tmat[1] * t(c(-1, pseudoinverse(M_bar) %*% m_bar)) )[-1]
    
    # Transform threshold 
    if ( is.null(x_boundary) ) {
      x_boundary <- getBoundaryPoint(E)   # ~~~~~~~~~~~~~~~~ change this later!
    }
    v <- as.numeric( tmat %*% x_boundary - z )[-1]
    b_prime <- as.numeric( t(v - q_prime) %*% Q_prime %*% (v - q_prime) )
    spec$x_boundary <- v
    
    # Transform starting point
    if ( is.null(spec$x0) ) {
      x0 <- getInitialPoint(E)            # ~~~~~~~~~~~~~~~~ change this later!
    }
    spec$x0 <- (tmat %*% x0 - z)[-1]
    
    E_prime <- ellipsoid(q = q_prime, 
                         Q = Q_prime, 
                         b = b_prime)
    
    setEllipsoid(E) <- E_prime
    setSpec(E) <- spec
    
    return( E )
  }
  
  
  # --------------------------------------------------------------------------
  #' @rdname transform-methods
  #' @aliases transform, RPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "transform", 
            signature = "RPB",
            definition = function(object)
            {
              ans <- .transformB(object = object)
              return( ans )
            })
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .transformB <- function(object)
  {
    P <- getPolytope(object)
    E <- getEllipsoid(object)
    H <- getHyperplane(P)
    
    if ( is.null(H) ) {
      
      if ( isTRUE(getSpec(object)$verbose) ) {
        warning("No transformation is performed since the 
              linear system (polytope) part of the convex body 
              contains no equality.")
      }
      
    } else {
      
      P_prime <- transform(object = as.RP(P, getSpec(object)))
      E_prime <- transform(object = as.RP(E, getSpec(object)))
      B_prime <- convexBody( P = P_prime,
                             E = E_prime )
      setConvexBody(object) <- B_prime
      setSpec(object) <- getSpec(E_prime)
    }
    
    return( object )
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @rdname transform-methods
  #' @aliases transform, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "transform", 
            signature = "POLYTOPE",
            definition = function(object)
            {
              spec <- list(transformation = simplexTransform(d = dimension(object)))
              object <- as.RP(object = object, spec = spec)
              ans <- transform(object = object)
              return( ans )
            })
 
  # --------------------------------------------------------------------------
  #' @rdname transform-methods
  #' @aliases transform, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "transform", 
            signature = "ELLIPSOID",
            definition = function(object)
            {
              spec <- list(transformation = simplexTransform(d = dimension(object)))
              object <- as.RP(object = object, spec = spec)
              ans <- transform(object = object)
              return( ans )
            })
  
  # --------------------------------------------------------------------------
  #' @rdname transform-methods
  #' @aliases transform, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "transform", 
            signature = "CONVEXBODY",
            definition = function(object)
            {
              spec <- list(transformation = simplexTransform(d = dimension(object)))
              object <- as.RP(object = object, spec = spec)
              ans <- transform(object = object)
              return( ans )
            })
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' backTransform
  #' @description Perrforms a (isometric, i.e. distance preserving) transformation
  #'              of body from d-demensional space to (d-1)-demensional space.
  #'              Note that the body needs to be a simplex or a the intersection of a 
  #'              body with a simplex (i.e. a constrained simplex).
  #' @param object [S4 Class] Object of class \code{RP}
  #' @return Returns an object of the same class as the imput 
  #'         but transformed to (d-1)-dimensional space.
  #' @docType methods
  #' @rdname backTransform-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "backTransform", 
             def = function(object) { standardGeneric("backTransform") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname backTransform-methods
  #' @aliases backTransform,RP-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "backTransform", 
            signature = "RP",
            definition = function(object)
            {
              .backTransform(object = object)  
            })
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .backTransform <- function(object)
  {
    spec <- getSpec(object)
    tmat_inv <- spec$transformation$tmat_inv
    z <- spec$transformation$z
    
    # Back-transform samples
    samples_prime <- getSamples(object, name = "all")
    if ( length(samples_prime) > 0 ) { 
      FUN <- function(mat) { t(apply(mat, 1, function(y) { tmat_inv %*% (c(0, y) + z) })) }
      samples <- lapply( samples_prime, FUN = FUN )
    } else {
      samples <- list()
    }
    
    # # Back-transform geometric bodies
    # P_prime <- getPolytope(object)
    # if ( !isEmpty(P_prime) ) {
    #   P <- polytope(A = tmat_inv %*% cbind(getAmat(P_prime), 0) + z,
    #                 b = tmat_inv %*% getRHS(P_prime) + z,
    #                 sense = getSense(P_prime))     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^   
    # }
    # E_prime <- getEllipsoid(object)
    # if ( !isEmpty(E_prime) ) {
    #   E <- ellipsoid()     
    # }
    
    
    # Back-transform spec variables
    if ( !is.null(spec$x0) ) {
      spec$x0 <-  tmat_inv %*% (c(0, spec$x0) + z)
    }
    if ( !is.null(spec$x_boundary) ) {
      spec$x_boundary <-  tmat_inv %*% (c(0, spec$x_boundary) + z)
    }
    
    setSamples(object) <- samples
    setSpec(object) <- spec
    
    return( object )
  }
  
  
  # --------------------------------------------------------------------------
  #' thin
  #' @description Thins the samples of an RP object.
  #' @param object    [S4] Any class that inherits from class \code{RP}.
  #' @return Returns the thinned input object.
  #' @docType methods
  #' @rdname thin-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "thin", 
             def = function(object) { standardGeneric("thin") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname thin-methods
  #' @aliases thin, RP-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "thin", 
            signature = "RP",
            definition = function(object) 
            {
              .thin(object)
            })
  # --------------------------------------------------------------------------
  .thin <- function(object)
  {
    thin <- ifelse( is.null(getSpec(object)$thin), 0, getSpec(object)$thin )
    if ( thin > 0 ) {
      lSamples <- getSamples(object, "all")
      for ( i in 1:length(lSamples) ) {
        samples <- lSamples[[i]]
        if ( nrow(samples) > thin ) {
          idx <- floor(seq(from = 1, to = nrow(samples), length.out = nrow(samples) - thin))
          lSamples[[i]] <- samples[idx, ]
        }
      }
      setSamples(object) <- lSamples
    }
    return( object )
  }

  
  
  
  # --------------------------------------------------------------------------
  #' isEmpty
  #' @description Checks if the input object is empty
  #' @param object    [POLYTOPE] Any class that inherits from POLYTOPE.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname isEmpty-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "isEmpty", 
             def = function(object) { standardGeneric("isEmpty") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname isEmpty-methods
  #' @aliases isEmpty, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "isEmpty", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              
              ans <- FALSE
              if ( any(is.na(getAmat(object))) ) {
                ans <- TRUE
              }
              if ( any(is.na(getRHS(object))) ) {
                ans <- TRUE
              }
              
              return( ans )
            })
  
  # --------------------------------------------------------------------------
  #' @rdname isEmpty-methods
  #' @aliases isEmpty, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "isEmpty", 
            signature = "ELLIPSOID",
            definition = function(object) 
            {
              ans <- FALSE
              if ( any(is.na(getShape(object))) ) {
                ans <- TRUE
              }
              if ( any(is.na(getCentre(object))) ) {
                ans <- TRUE
              }
              return( ans )
            })  
  
  
  
  # --------------------------------------------------------------------------
  #' dimension
  #' @description Returns the dimension of an object
  #' @param object    [POLYTOPE] Any class that inherits from POLYTOPE.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname dimension-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "dimension", 
             def = function(object) { standardGeneric("dimension") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname dimension-methods
  #' @aliases dimension, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "dimension", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              ans = 0
              if ( !isEmpty(object) ) {
                ans <- ncol( getAmat(object) )
              }
              return( ans )
            })
  
  # --------------------------------------------------------------------------
  #' @rdname dimension-methods
  #' @aliases dimension, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "dimension", 
            signature = "ELLIPSOID",
            definition = function(object) 
            {
              ans = 0
              if ( !isEmpty(object) ) {
                ans <- ncol( getShape(object) )
              }
              return( ans )
            })
  
  # --------------------------------------------------------------------------
  #' @rdname dimension-methods
  #' @aliases dimension, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "dimension", 
            signature = "CONVEXBODY",
            definition = function(object) 
            {
              d_p <- dimension( getPolytope(object) )
              d_e <- dimension( getEllipsoid(object) )
              if ( d_p != d_e ) {
                if ( any( c(d_p, d_e) == 0 ) ) {
                  stop("Ellipsoid and Polytope do not have the same dimensions.")
                } else {
                  ans <- max( d_p, d_e )
                }
              } else {
                ans <- d_p
              }
              return( ans )
            })

  
  
  
  
  # --------------------------------------------------------------------------
  #' dists2Body
  #' @description Computes the distances from a point x to the boundary of a body in 
  #'              in directions u and -u.
  #' @param object [S4 Class] Object of class \code{CONVEXBODY}, 
  #'                          \code{POLYTOPE} or \code{ELLIPSOID}.
  #' @param x         [vector] Point within the body.
  #' @param u         [vector] Optional direction vector. 
  #' @return Returns a numeric vector of length two.
  #' @docType methods
  #' @rdname dists2Body-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "dists2Body", 
             def = function(object, x, u = NULL) { standardGeneric("dists2Body") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname dists2Body-methods
  #' @aliases dist2Body, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "dists2Body",
            signature = "POLYTOPE",
            definition = function(object, x, u = NULL)
            {
              d <- dists2Polytope(P = object, x = x, u = u)
              return( d )
            })
  # --------------------------------------------------------------------------
  #' @rdname dists2Body-methods
  #' @aliases dist2Body, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "dists2Body",
            signature = "ELLIPSOID",
            definition = function(object, x, u = NULL)
            {
              d <- dists2Ellipsoid(E = object, x = x, u = u)
              return( d )
            })
  # --------------------------------------------------------------------------
  #' @rdname dists2Body-methods
  #' @aliases dist2Body, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "dists2Body",
            signature = "CONVEXBODY",
            definition = function(object, x, u = NULL)
            {
              d <- dists2ConvexBody(B = object, x = x, u = u)
              return( d )
            })
  
  
  
  
  
  
 
             
  
  
  # --------------------------------------------------------------------------
  #' intersect
  #' @description Creates an obejct 
  #' @param x [S4 Class] Object that inherits from class \code{VBODY}.
  #' @param y [S4 Class] Object that inherits from class \code{VBODY}.
  #' @return Returns an object containing the intersection of the 
  #'         two input objects.
  #' @docType methods
  #' @rdname intersect-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "intersect", 
             def = function(x, y) { standardGeneric("intersect") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname intersect-methods
  #' @aliases intersect, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "intersect", 
            signature = c("POLYTOPE", "POLYTOPE"),
            definition = function(x, y)
            {
              ans <- .intersectPP(p1 = x, p2 = y)
              return( ans )
            })
  # --------------------------------------------------------------------------
  .intersectPP <- function(p1, p2)
  {
    
    stopifnot( dimension(p1) == dimension(p2) )
    
    A <- rbind( getAmat(p1),
                getAmat(p2) )
    colnames(A) <- colnames(getAmat(p1))
    b <- c(getRHS(p1), getRHS(p2))
    sense <- c(getSense(p1), getSense(p2))
    P <- polytope(A = A,
                  b = b,
                  sense = sense)
    return( P )
  }
  
  
  
  
  
  
  
  

  
  
 
  
  