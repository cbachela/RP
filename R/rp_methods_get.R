  
  
  ############################################################################
  ### RP PACKAGE - GET METHODS 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     03.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  #
  # METHODS:
  #
  # getPolytope
  # getEllipsoid
  # getHyperplane
  # getShape
  # getCentre
  # getNormal
  # getRHS
  # getSense
  # getAmat
  # getBoundaryPoint
  # getInteriorPoint
  # getSamples
  # getSpec

  
  
  
  
  
  # --------------------------------------------------------------------------
  #' getPolytope
  #' @description Get polytope object from a \code{VBODY} object.
  #' @param object    [VBODY]
  #' @return Returns an object of class 'VBODY'.
  #' @docType methods
  #' @rdname getPolytope-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getPolytope", 
             def = function(object) { standardGeneric("getPolytope") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getPolytope-methods
  #' @aliases getPolytope, VBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getPolytope", 
            signature = "VBODY",
            definition = function(object) 
            {
              if ( inherits(object, "POLYTOPE") ) {
                ans <- object
              } else if ( inherits(object, "CONVEXBODY") ) {
                ans <- object@P
              } else {
                ans <- polytope()
              }
              return( ans )
            })
  
 
  
  
  
  # --------------------------------------------------------------------------
  #' getEllipsoid
  #' @description Get ellipsoid object from a \code{VBODY} object.
  #' @param object    [VBODY]
  #' @return Returns an object of class 'ELLIPSOID'.
  #' @docType methods
  #' @rdname getEllipsoid-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getEllipsoid", 
             def = function(object) { standardGeneric("getEllipsoid") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getEllipsoid-methods
  #' @aliases getEllipsoid, VBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getEllipsoid", 
            signature = "VBODY",
            definition = function(object) 
            {
              if ( inherits(object, "ELLIPSOID") ) {
                ans <- object
              } else if ( inherits(object, "CONVEXBODY") ) {
                ans <- object@E
              } else {
                ans <- ellipsoid()
              }
              return(ans)
            })
  
 
  

  
  
  # --------------------------------------------------------------------------
  #' getHyperplane
  #' @description Get hyperplane object from a \code{POLYTOPE} object.
  #' @param object    [POLYTOPE]
  #' @return Returns an object of class 'HYPERPLANE'.
  #' @docType methods
  #' @rdname getHyperplane-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getHyperplane", 
             def = function(object) { standardGeneric("getHyperplane") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getHyperplane-methods
  #' @aliases getHyperplane, POLYTOPEY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getHyperplane", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              .getHyperplane(object = object)
            })
  # --------------------------------------------------------------------------
  #' @export
  # -------------------------------------------------------------------------- 
  .getHyperplane <- function(object)
  {
    idx_eq <- which(object@sense == "=")
    H <- NULL
    if ( length(idx_eq) >= 1 ) {
      names(idx_eq) <- rownames(getAmat(object))[idx_eq]
      for ( i in seq(along = idx_eq) ) {
        if ( i == 1 ) {
          H <- hyperplane( v = getAmat(object)[idx_eq[i], ],
                           b = getRHS(object)[idx_eq[i]] )
        } else {
          # H <- list(H, hyperplane( v = getAmat(object)[idx_eq[i], ],
          #                          b = getRHS(object)[idx_eq[i]] ) )
          stop( "object has more than one equality.")
        }
      }
      if ( length(H) > 1 ) {
        names(H) <- names(idx_eq)
      } 
    } 
    return( H )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' getShape
  #' @description Get the shape matrix from an \code{ELLIPSOID} object.
  #' @param object    [ELLIPSOID] Any class that inherits from ELLIPSOID.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getShape-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getShape", 
             def = function(object) { standardGeneric("getShape") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getShape-methods
  #' @aliases getShape, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getShape", 
            signature = "ELLIPSOID",
            definition = function(object)
            {
              return(object@Q)
            })
  
  
  # --------------------------------------------------------------------------
  #' getCentre
  #' @description Get the centroid from an \code{ELLIPSOID} object.
  #' @param object    [ELLIPSOID] Any class that inherits from ELLIPSOID.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getCentre-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getCentre", 
             def = function(object) { standardGeneric("getCentre") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getCentre-methods
  #' @aliases getCentre, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getCentre", 
            signature = "ELLIPSOID",
            definition = function(object) 
            {
              return(object@q)
            })
  
  
  # --------------------------------------------------------------------------
  #' getNormal
  #' @description Get the normal vector from an \code{HYPERPLANE} object.
  #' @param object    [HYPERPLANE] Any class that inherits from HYPERPLANE.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getNormal-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getNormal", 
             def = function(object) { standardGeneric("getNormal") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getNormal-methods
  #' @aliases getNormal, HYPERPLANE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getNormal", 
            signature = "HYPERPLANE",
            definition = function(object) 
            {
              return(object@v)
            })
  
  
  
  # --------------------------------------------------------------------------
  #' getRHS
  #' @description Get the right-hand-side (RHS) variable from a \code{VBODY} object.
  #' @param object    [VBODY] Any class that inherits from VBODY.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getRHS-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getRHS", 
             def = function(object) { standardGeneric("getRHS") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getRHS-methods
  #' @aliases getRHS, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getRHS", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              return(object@b)
            })
  
  # --------------------------------------------------------------------------
  #' @rdname getRHS-methods
  #' @aliases getRHS, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getRHS", 
            signature = "ELLIPSOID",
            definition = function(object) 
            {
              return(object@b)
            })

  # --------------------------------------------------------------------------
  #' @rdname getRHS-methods
  #' @aliases getRHS, HYPERPLANE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getRHS", 
            signature = "HYPERPLANE",
            definition = function(object) 
            {
              return(object@b)
            })

    
  
  # --------------------------------------------------------------------------
  #' getSense
  #' @description Get the right-hand-side (RHS) variable from a \code{VBODY} object.
  #' @param object    [VBODY] Any class that inherits from VBODY.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getSense-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getSense", 
             def = function(object) { standardGeneric("getSense") },
             package = "RP")
  # --------------------------------------------------------------------------
  #' @rdname getSense-methods
  #' @aliases getSense, VBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getSense", 
            signature = "VBODY",
            definition = function(object) 
            {
              return(object@sense)
            })
  
  
  # --------------------------------------------------------------------------
  #' getAmat
  #' @description Get the A-matrix from a \code{POLYTOPE} object.
  #' @param object    [POLYTOPE] Any class that inherits from POLYTOPE.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getAmat-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getAmat", 
             def = function(object) { standardGeneric("getAmat") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getAmat-methods
  #' @aliases getAmat, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getAmat", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              return(object@A)
            })
  
  # --------------------------------------------------------------------------
  #' @rdname getAmat-methods
  #' @aliases getAmat, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getAmat", 
            signature = "CONVEXBODY",
            definition = function(object) 
            {
              return( getAmat(getPolytope(object)) )
            })
  
  
  # --------------------------------------------------------------------------
  #' getBoundaryPoint
  #' @description Get a point on the boundary of a body. 
  #' @param object [S4 Class] Object of class \code{CONVEXBODY}, 
  #'                          \code{POLYTOPE} or \code{ELLIPSOID}.
  #' @return Returns a vector.
  #' @docType methods
  #' @rdname getBoundaryPoint-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getBoundaryPoint", 
             def = function(object) { standardGeneric("getBoundaryPoint") },
             package = "RP")
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .getBoundaryPointP <- function(P)
  {
    x <- getInteriorPoint(P)
    u <- rnorm(length(x))
    u <- u / as.numeric(sqrt(sum(u^2)))
    lambda <- dists2Polytope(P = P, 
                             x = x,
                             u = u)
    p <- x + lambda[2] * u
    
    # Check
    # ...
    
    return( p )
  }
  
  # --------------------------------------------------------------------------
  #' @rdname getBoundaryPoint-methods
  #' @aliases getBoundaryPoint, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getBoundaryPoint", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              x0 <- .getBoundaryPointP(P = object)
              return( x0 )
            })
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .getBoundaryPointE <- function(E)
  {
    x <- getCentre(E)
    u <- rnorm(length(x))
    u <- u / as.numeric(sqrt(sum(u^2)))
    lambda <- dists2Ellipsoid(E = E, 
                              x = x,
                              u = u)
    p <- x + lambda[2] * u
    
    # Check
    boundary_check <- t(p) %*% getShape(E) %*% p - getRHS(E)
    stopifnot( boundary_check < 1e-16 )
    
    return( p )
  }
  
  # --------------------------------------------------------------------------
  #' @rdname getBoundaryPoint-methods
  #' @aliases getBoundaryPoint, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getBoundaryPoint", 
            signature = "ELLIPSOID",
            definition = function(object) 
            {
              x0 <- .getBoundaryPointE(E = object)
              return( x0 )
            })
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  .getBoundaryPointB <- function(B)
  {
    P <- getPolytope(object = B)
    E <- getEllipsoid(object = B)
    
    if ( isEmpty(P) && isEmpty(E) ) {
      stop("Impossible to find boundary point because 
           polytope and ellipsoid are both empty.")
    }
    
    if ( isEmpty(P) ) {
      x0 <- getBoundaryPoint(object = E)
    } else {
    
      if ( isEmpty(E) ) {
        x0 <- getBoundaryPoint(object = P)
        
      } else {
        
        # Compute minimum variance point
        A <- getAmat(P)
        d <- ncol(A)
        selection <- colnames(A)
        if ( is.null(selection) ) {
          selection <- paste0("X", 1:d)
          colnames(A) <- selection
        }
        lincon <- linearConstraint(Amat = A,
                                   rhs = getRHS(P),
                                   sense = P@sense)
        Constraints <- constraints(selection)
        addConstraint(Constraints) <- boxConstraint(name = "Unbounded")
        addConstraint(Constraints) <- lincon
        
        covmat <- getShape(E)
        GPS <- gps(Data = NULL,
                   Covariance = covmat,
                   Constraints = Constraints)
        GPO <- gpo(GPS = GPS)
        w_gmv <- getWeights(GPO)
        w_eqw <- w_gmv * 0 + 1 / length(w_gmv)
        
        sigmaFUN <- function(x) { as.numeric( t(x - centre) %*% covmat %*% (x - centre) ) }
        rhs <- getRHS(E)
        centre <- getCentre(E)

        # obj_fun <- function(alpha)
        # {
        #   sigma_th <- sigmaFUN(x = alpha * w_gmv + (1 - alpha) * w_eqw)
        #   ans <-  abs(sigma_th - rhs)
        #   return( ans )
        # }
        # opt <- optimize(f = obj_fun, interval = c(0, 1))
        # alpha <- opt$minimum
        # x0 <- alpha * w_gmv + (1 - alpha) * w_eqw
        
        # Use closed form solution (NOTE: does not account for centre !! )
        alpha <- getRHS(E) / sigmaFUN(w_gmv) - 1
        tmp <- twoAssetPortfolio(w1 = w_gmv,
                                 w2 = w_eqw,
                                 covmat = covmat,
                                 alpha = alpha)    # ~~~~~~~~~~~~ 
        x0 <- tmp$w
        
        # # Use acceptance-rejection
        # idx <- NULL; iter <- 0
        # while ( length(idx) <= 0 || iter > 100 ) {
        #   PT <- transform(P)
        #   samples_PT <- rp( PT, n_sim = 10^2, x0 = NULL )$S
        #   samples_P <- t(apply(samples_PT, 1, 
        #                        function(y) { PT@info$tmat_inv %*% (c(0, y) + PT@info$z) }))
        #   sigma_values <- apply( samples_P, 1, sigmaFUN )
        #   idx <- which(sigma_values - getRHS(E) > 0)
        #   iter <- iter + 1
        # }
        # w2 <- samples_P[idx[1], ]
        # 
        # alpha <- getRHS(E) / sigmaFUN(w1) - 1
        # tmp <- twoAssetPortfolio(w1 = w1,
        #                          w2 = w2,
        #                          covmat = covmat,
        #                          alpha = alpha)   
        # x0 <- tmp$w
        
      }
    } 
    
    return( x0 )
  }
  
  # --------------------------------------------------------------------------
  #' @rdname getBoundaryPoint-methods
  #' @aliases getBoundaryPoint, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getBoundaryPoint", 
            signature = "CONVEXBODY",
            definition = function(object) 
            {
              #// NOTE: Requires GPO and Gurobi !!!!!!!!!!!!!!!!!!!!!!!!!
              x0 <- .getBoundaryPointB(B = object)
              return( x0 )
            })
            
  
  
  # --------------------------------------------------------------------------
  #' getInteriorPoint
  #' @description Get a point in the interior of a body.
  #' @param object [S4 Class] Object of class \code{CONVEXBODY}, 
  #'                          \code{POLYTOPE} or \code{ELLIPSOID}.
  #' @return Returns a vector.
  #' @docType methods
  #' @rdname getInteriorPoint-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getInteriorPoint", 
             def = function(object) { standardGeneric("getInteriorPoint") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getInteriorPoint-methods
  #' @aliases getInteriorPoint, ELLIPSOID-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getInteriorPoint", 
            signature = "ELLIPSOID",
            definition = function(object) 
            {
              getCentre(object)
            })
  
  # --------------------------------------------------------------------------
  #' @rdname getInteriorPoint-methods
  #' @aliases getInteriorPoint, POLYTOPE-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getInteriorPoint", 
            signature = "POLYTOPE",
            definition = function(object) 
            {
              
              #// NOTE: Requires GPO and Gurobi !!!!!!!!!!!!!!!!!!!!!!!!!
              
              A <- getAmat(object)
              d <- ncol(A)
              selection <- colnames(A)
              if ( is.null(selection) ) {
                selection <- paste0("X", 1:d)
                colnames(A) <- selection
              }
              lincon <- linearConstraint(Amat = A,
                                         rhs = getRHS(object),
                                         sense = object@sense)
              
              covmat <- as.matrix( diag(d), dimnames = list(selection, selection))
              Constraints <- constraints(selection)
              addConstraint(Constraints) <- boxConstraint(name = "Unbounded")
              GPS <- gps(Data = NULL,
                         Covariance = covmat,
                         Constraints = Constraints)
              GPO <- gpo(GPS = GPS)
              x0 <- getWeights(GPO)
              
              return( x0 )
            })
  
  # --------------------------------------------------------------------------
  #' @rdname getInteriorPoint-methods
  #' @aliases getInteriorPoint, CONVEXBODY-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getInteriorPoint", 
            signature = "CONVEXBODY",
            definition = function(object) 
            {
              
              #// NOTE: Requires GPO and Gurobi !!!!!!!!!!!!!!!!!!!!!!!!!
              
              
              P <- getPolytope(object)
              E <- getEllipsoid(object)
              
              if ( !isEmpty(P) ) {
                
                A <- getAmat(P)
                d <- ncol(A)
                selection <- colnames(A)
                if ( is.null(selection) ) {
                  selection <- paste0("X", 1:d)
                  colnames(A) <- selection
                }
                lincon <- linearConstraint(Amat = A,
                                           rhs = getRHS(P),
                                           sense = P@sense)
                
                covmat <- as.matrix( diag(d), dimnames = list(selection, selection))
                Constraints <- constraints(selection)
                addConstraint(Constraints) <- boxConstraint(name = "Unbounded")
                addConstraint(Constraints) <- lincon
                Solver <- solverCtrl()
                
                if ( !isEmpty(E) ) {
                  
                  #// to do: 
                  # variance constraint is not necessary bcs 
                  # minimum variance within polytope always satisfies 
                  # variance bound. But what if the centre of the ellipsoid
                  # is not at the origin (e.g. with a tracking error constraint) ?
                  
                  covmat <- getShape(E)
                  quadcon <- varianceConstraint(rhs = getRHS(E),
                                                Qmat = covmat,
                                                sense = "<=")
                  addConstraint(Constraints) <- quadcon
                  Solver <- solverCtrl(progtype = "QCP")
                }
                
                GPS <- gps(Data = NULL,
                           Solver = Solver,
                           Covariance = covmat,
                           Constraints = Constraints)
                GPO <- gpo(GPS = GPS)
                x0 <- getWeights(GPO)
                
              } else {
                if ( !isEmpty(E) ) {
                  x0 <- getCentre(E)
                } else {
                  stop("Cannot find interior point because 
                       polytope and ellipsoid are both empty.")
                }
              }
              
              return( x0 )
              
            })
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' getSamples
  #' @description Get the samples from an object of class \code{RP}.
  #' @param object [S4 Class] Object of class \code{RP}.
  #' @param name [character] Indicate which samples to return.
  #' @return Returns a matrix.
  #' @docType methods
  #' @rdname getSamples-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getSamples", 
             def = function(object, name = NULL) { standardGeneric("getSamples") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getSamples-methods
  #' @aliases getSamples, RP-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getSamples", 
            signature = "RP",
            definition = function(object, name = NULL)
            {
              name <- ifelse( is.null(name), 1, name )
              samples <- object@samples
              if ( length(samples) > 0 ) {
                if ( name != "all" ) {
                  samples <- samples[[name]]    
                }
              }
              return( samples )
            })
  
  
  # --------------------------------------------------------------------------
  #' getSpec
  #' @description Get the spec slot from an object of class \code{RP}.
  #' @param object [S4 Class] Object of class \code{RP}.
  #' @return Returns a list
  #' @docType methods
  #' @rdname getSpec-methods
  #' @export
  # --------------------------------------------------------------------------
  setGeneric(name = "getSpec", 
             def = function(object) { standardGeneric("getSpec") },
             package = "RP")
  
  # --------------------------------------------------------------------------
  #' @rdname getSpec-methods
  #' @aliases getSpec, RP-method
  #' @export
  # --------------------------------------------------------------------------
  setMethod(f = "getSpec", 
            signature = "RP",
            definition = function(object) 
            {
              return( object@spec )
            })
  
  
  
  
  