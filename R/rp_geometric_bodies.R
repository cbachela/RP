  
  
  ############################################################################
  ### RP PACKAGE - GEOMETRIC BODIES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     03.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # FUNCTIONS:
  # 
  # simplex
  # boundedSimplex
  # polytope
  # ellipsoid
  # convexBody
  # as.polytope
  # as.ellipsoid
  # as.cbody
  # as.RP
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  simplex <- function(d, b = 1)
  {
    
    if ( missing(d) ) stop("Dimension variable 'd' is missing.")
    
    A <- matrix(1, nrow = 1, ncol = d)
    A <- rbind(A, diag(d) * (-1), diag(d))
    b <- c(b, rep(0, d), rep(1, d))
    sense <- c("=", rep("<=", d * 2))
    
    ans <- new("SIMPLEX",
               A = A, 
               b = b, 
               sense = sense)
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  boundedSimplex <- function(d, v = 1, lower = NULL, upper = NULL)
  {
    
    if( missing(d) ) {
      d <- length(lower)
    }
    stopifnot(d > 0)
    
    A <- matrix(1, nrow = 1, ncol = d)
    A <- rbind(A, diag(d) * (-1))
    if( is.null(lower) ) {
      lower <- rep(0, d)
    }
    stopifnot( length(lower) == d )
    b <- c(v, lower)
    sense <- c("=", rep("<=", d))
    
    if( !is.null(upper) ) {
      stopifnot( length(upper) == d )
      A <- rbind(A, diag(d))
      b <- c(b, upper)
      sense <- c(sense, rep("<=", d))
    }
    
    ans <- new("POLYTOPE",
               A = A, 
               b = b, 
               sense = sense)
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  polytope <- function(A = matrix(NA), b = NA, sense = "<=")
  {
    idx_geq <- which(sense == ">=")
    if ( length(idx_geq) > 0 ) {
      A[idx_geq, ] <- A[idx_geq, ] * (-1)
      b[idx_geq] <- b[idx_geq] * (-1)
      sense[idx_geq] <- "<="
    }
    ans <- new("POLYTOPE",
               A = A, 
               b = b, 
               sense = sense)
    return( ans )   
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  hyperplane <- function(v, b)
  {
    ans <- new("HYPERPLANE",
               v = v, 
               b = b)
    return( ans )   
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  ellipsoid <- function(q = NA, Q = matrix(NA), b = 0)
  {
    ans <- new("ELLIPSOID",
               q = q, 
               Q = Q,
               b = b)
    return( ans )   
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  convexBody <- function(P = polytope(), E = ellipsoid())
  {
    ans <- new("CONVEXBODY",
               P = P,
               E = E)
    return( ans )   
  }
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  as.polytope <- function(object)
  {
    if ( class(object) == "SIMPLEX" ) {
      P <-  polytope(A = getAmat(object), 
                     b = getRHS(object),
                     sense = object@sense)
    
    } else if ( class(object) == "CONSTRAINT" ) {
      
      selection <- getConstraints(object, "selection")
      boxcon <- getConstraints(object, "bounds")
      lincon <- getConstraints(object, "linear")
      if ( length(selection) < length(boxcon$upper) ) {
        selection <- names(boxcon$upper)
      }
      
      # P <- simplex(n = length(selection),  
      #              v = 1,                 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      #              lower = boxcon$lower,
      #              upper = boxcon$upper)
      
      n <- length(selection)
      stopifnot( length(boxcon$lower) == n )
      stopifnot( length(boxcon$upper) == n )
      
      # Bounds
      A <- rbind(diag(n) * (-1),
                 diag(n) )
      colnames(A) <- selection
      b <- c(boxcon$lower * (-1), boxcon$upper)
      sense <- rep("<=", n * 2)
      
      P <- polytope(A = A, 
                    b = b, 
                    sense = sense)
      A = b = sense = NULL
      
      # Linear constraints
      if ( !is.null(lincon$Amat) ) {
        
        A <- lincon$Amat
        rhs <- lincon$rhs
        sense <- lincon$sense
        
        idx_geq <- which(sense == ">=")
        if ( length(idx_geq) > 0 ) {
          A[idx_geq, ] <- A[idx_geq, ] * (-1)
          rhs[idx_geq] <- rhs[idx_geq] * (-1)
          sense[idx_geq] <- "<="
        }
        
        P@A <- rbind(P@A, A)
        P@b <- c(P@b, rhs)
        P@sense <- c(P@sense, sense)
      }
      colnames(P@A) <- selection
    }
    
    return( P )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  as.ellipsoid <- function(object)
  {
    quadcon <- getConstraints(object, "quadratic")
    
    if ( !is.null(quadcon$variance) ) {
      E <- ellipsoid(q = rep(0, ncol(quadcon$variance$Qc)),
                     Q = quadcon$variance$Qc,
                     b = quadcon$variance$rhs)
      
    } else {
      E <- ellipsoid()
    }
    
    return( E )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  as.cbody <- function(object)
  {
    P <- as.polytope(object = object)
    E <- as.ellipsoid(object = object)
    CB <- convexBody(P = P,
                     E = E)
    return( CB )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  as.RP <- function(object, spec = NULL)
  {
    if ( is.null(spec) ) {
      spec <- rpCtrl()
    }
    if ( inherits(object, "POLYTOPE") ) {
      ans <- new("RPP", spec = spec, object)
    }
    if ( inherits(object, "ELLIPSOID") ) {
      ans <- new("RPE", spec = spec, object)
    }
    if ( inherits(object, "CONVEXBODY") ) {
      ans <- new("RPB", spec = spec, object)
    }
    return( ans )
  }
  
  
  
  
  
  
  
  