  
  
  ############################################################################
  ### RP PACKAGE - SET METHODS 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     03.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  #
  # METHODS:
  #
  # addHyperplane
  # replaceHyperplane
  # setSamples
  # setSpec
  # setPolytope
  
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Adding one or several inequality equations.
  #' @param object    [POLYTOPE] object of class 'POLYTOPE'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("addHyperplane<-", function(object, value)
  {
    standardGeneric("addHyperplane<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname POLYTOPE-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("addHyperplane", signature = "POLYTOPE",
                   definition = function(object, value)
                   {

                     if ( class(value) == "list" ) {
                       if ( is.null(unlist(value)) ) {
                         return( object )
                       }
                       a <- value$Amat
                       if ( is.null(dim(Atmp)) ) {
                         a <- matrix(Atmp, nrow = 1)
                       }
                       b <- value$rhs
                       sense <- ifelse( is.null(value$sense) == "<=", value$sense)
                     } else if ( class(value) == "HYPERPLANE" ) {
                       a <- getNormal(value)
                       b <- getRHS(value)
                       sense <- "<="
                     } else {
                       stop("Input has to be a list or HPYERPLANE")
                     }
                     
                      # Concatenate existing and new A matrix and rhs
                     object@A <- rbind(object@A, a)
                     object@b <- c(object@b, b)
                     object@sense <- c(object@sense, sense)
                     
                     return( object )
                   })
  
  
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Replacing an existing inequality.
  #' @param object    [POLYTOPE] object of class 'POLYTOPE'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("replaceHyperplane<-", function(P, H, row_name)
  {
    standardGeneric("replaceHyperplane<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname POLYTOPE-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("replaceHyperplane", signature = "POLYTOPE",
                   definition = function(P, H, row_name = NULL)
                   {
                     stopifnot( class(H) == "HYPERPLANE" )
                     if ( isEmptyp(P) ) { stop("Polytope P is empty.") }
                     if ( isNULL(row_name) ) { stop("row_name has to be passed.")}
                     A <- getAmat(P)
                     idx <- which(rownames(A) == row_name)
                     if ( length(idx) != 1 ) {
                       stop("There is no row in the linear system with the corresponding row_name.")
                     } else {
                       P@A[idx, ] <- getNormal(H)
                       P@b[idx] <- getRHS(H)
                     }
                     return( P )
                   })
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Set samples to simplex object.
  #' @param object    [RP] object of class 'RP'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("setSamples<-", function(object, value)
  {
    standardGeneric("setSamples<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname RP-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("setSamples", 
                   signature = "RP",
                   definition = function(object, value)
                   {
                     if ( !is.list(value) ) {
                       value <- list(value)
                     }
                     object@samples <- value
                     return( object )
                   })
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Set specifications to  object.
  #' @param object    [RP] object of class 'RP'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("setSpec<-", function(object, value)
  {
    standardGeneric("setSpec<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname RP-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("setSpec", 
                   signature = "RP",
                   definition = function(object, value)
                   {
                     object@spec <- value
                     return( object )
                   })
  
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Set polytope to object.
  #' @param object    [RP] object of class 'RP'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("setPolytope<-", function(object, value)
  {
    standardGeneric("setPolytope<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname RP-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("setPolytope", 
                   signature = "RP",
                   definition = function(object, value)
                   {
                     stopifnot( inherits(value, "POLYTOPE") )
                     object@A <- getAmat(value)
                     object@b <- getRHS(value)
                     object@sense <- getSense(value)
                     return( object )
                   })
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Set ellipsoid to object.
  #' @param object    [RP] object of class 'RP'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("setEllipsoid<-", function(object, value)
  {
    standardGeneric("setEllipsoid<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname RP-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("setEllipsoid", 
                   signature = "RP",
                   definition = function(object, value)
                   {
                     stopifnot( inherits(value, "ELLIPSOID") )
                     object@q <- getCentre(value)
                     object@Q <- getShape(value)
                     object@b <- getRHS(value)
                     return( object )
                   })
  
  
  # --------------------------------------------------------------------------
  #' Setter Function
  #' @description Set convex body to object.
  #' @param object    [RP] object of class 'RP'.
  #' @return Returns the modified input object.
  #' @export
  # --------------------------------------------------------------------------
  setGeneric("setConvexBody<-", function(object, value)
  {
    standardGeneric("setConvexBody<-")
  })
  # --------------------------------------------------------------------------
  #' @rdname RP-class
  #' @export
  # --------------------------------------------------------------------------
  setReplaceMethod("setConvexBody", 
                   signature = "RP",
                   definition = function(object, value)
                   {
                     stopifnot( inherits(value, "CONVEXBODY") )
                     object@P <- getPolytope(value)
                     object@E <- getEllipsoid(value)
                     return( object )
                   })
