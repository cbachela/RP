  
  
  ############################################################################
  ### RP PACKAGE - CLASSES 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  # CLASSES:
  # 
  # VBODY
  # POLYTOPE
  # SIMPLEX
  # HYPERPLANE
  # ELLIPSOID
  # CONVEXBODY
  # CONSTRAINT
  # RP
  
  
  
  # --------------------------------------------------------------------------
  #' @title Virtual Class
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("VBODY", "VIRTUAL")
 
  
  
  # --------------------------------------------------------------------------
  #' @title S4 class POLYTOPE
  #' @description setClass call for object of class \code{POLYTOPE}.
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("POLYTOPE", 
           slots = c(
               A  = "matrix",
               b  = "vector",
               sense = "character"),
           contains = "VBODY")
  
  
  # --------------------------------------------------------------------------
  #' @title S4 class SIMPLEX
  #' @description setClass call for object of class \code{SIMPLEX}.
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("SIMPLEX", contains = "POLYTOPE" )
  
  
  # --------------------------------------------------------------------------
  #' @title S4 class HYPERPLANE
  #' @description setClass call for object of class \code{HYPERPLANE}.
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("HYPERPLANE", 
           slots = c(
             v  = "vector",
             b  = "numeric") )
  
  
  
  # --------------------------------------------------------------------------
  #' @title S4 class ELLIPSOID
  #' @description setClass call for object of class \code{ELLIPSOID}.
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("ELLIPSOID", 
           slots = c(
               q  = "vector",
               Q  = "matrix",
               b  = "numeric"),
           contains = "VBODY")
  
 
  # --------------------------------------------------------------------------
  #' @title S4 class CONVEXBODY
  #' @description setClass call for object of class \code{CONVEXBODY}.
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("CONVEXBODY",
           slots = c(
               P  = "POLYTOPE",
               E  = "ELLIPSOID"),
           contains = "VBODY")


  
  
  
  
  # --------------------------------------------------------------------------
  #' @title S4 class RP
  #' @description setClass call for object of class \code{RP}.
  #' @import methods
  #' @export
  # --------------------------------------------------------------------------
  setClass("RP", "VIRTUAL")
  
  
  # --------------------------------------------------------------------------
  setClass("RPP", 
           representation(
             samples = "list",
             spec    = "list"),
           contains = c("POLYTOPE", "RP"))
  # --------------------------------------------------------------------------
  setClass("RPE", 
           representation(
             samples = "list",
             spec    = "list"),
           contains = c("ELLIPSOID", "RP"))
  # --------------------------------------------------------------------------
  setClass("RPB", 
           representation(
             samples = "list",
             spec    = "list"),
           contains = c("CONVEXBODY", "RP"))
  
  
  
  
 
  
  