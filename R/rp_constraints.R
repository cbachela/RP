  
  
  #############################################################################
  ### RP PACKAGE - CONSTRAINTS
  #############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.07.2019
  # First version:    01.04.2019
  # --------------------------------------------------------------------------
  #
  # FUNCTIONS
  #
  # constraints
  # simplexConstraint
  # boxConstraint
  # linearConstraint
  # quadraticConstraint
  # groupConstraint
  # bucketorderingConstraint
  # orderingConstraint
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  constraints <- function(selection)
  { 
    if (missing(selection)) {
      stop("Selection vector is required for initialisation.
           Has to be a character vector.")
    } else {
      if (!is.character(selection)) {
        stop("argument 'selection' has to be a character vector.")
      }
      obj <- new("CONSTRAINT", 
                 selection = selection) 
    }
    return(obj)
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  simplexConstraint <- function(rhs = 1,
                                sense = "=",
                                indexOrName = NULL)
  {
    ans <- linearConstraint(name = "simplex",
                            rhs = rhs, 
                            sense = sense, 
                            indexOrName = indexOrName)
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  boxConstraint <- function(name = c("LongOnly", "LongShort", "Unbounded"),
                            lower = NULL,
                            upper = NULL)
  {
    type = "bounds"
    name <- match.arg(name)
    if (name == "Unbounded") {
      lower <- -Inf
      upper <- Inf
    } else {
      if (is.null(lower)) {
        if (is.null(upper)) {
          lower <- ifelse(name == "LongOnly", 0, -1)
          upper <- 1
        } else {
          lower <- upper * 0
        }
      } else {
        if (any(lower < 0) & name == "LongOnly") {
          stop("Inconsistent lower bounds for type 'LongOnly'.
               Change type to LongShort or ensure that lower >= 0.")
        }
        if (is.null(upper)) {
          upper <- lower * 0 + 1
        } 
        }
  }
    ans <- list(name = name, type = type, lower = lower, upper = upper)
    return(ans)
  }  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  linearConstraint <- function(name = NULL,
                               sense = "=", 
                               rhs = Inf, 
                               indexOrName = NULL,
                               Amat = NULL,
                               a_values = NULL,
                               ellipsis = list())
  {
    ans =  list(type = "linear",
                name = name,
                sense = sense,
                rhs = rhs,
                indexOrName = indexOrName,
                Amat = Amat,
                a_values = a_values)
    if (length(ellipsis) > 0) {
      for (i in 1:length(ellipsis)) {
        ans[[names(ellipsis)[i]]] <- ellipsis[[i]]
      }
    }
    return(ans)
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  # Note: Constraint of the form x^TQx + q^T x <= b.
  #       If provided, Qmat needs to positive semi definite.
  #       If Qmat is NULL, the constraint will we be 
  #       x^Tx <= y^2 (second order cone) or 
  #       x^Tx <= y*z, with y and z non-negative .
  quadraticConstraint <- function(name = NULL,
                                  sense = "<=",
                                  rhs = Inf,
                                  indexOrName = NULL,
                                  Qc = NULL,
                                  q = NULL,
                                  ellipsis = list())
  {
    ans =  list(type = "quadratic",
                name = name,
                sense = sense,
                rhs = rhs,
                indexOrName = indexOrName,
                Qc = Qc,
                q = q)
    if (length(ellipsis) > 0) {
      for (i in 1:length(ellipsis)) {
        ans[[names(ellipsis)[i]]] <- ellipsis[[i]]
      }
    }
    return(ans)
  }
  
  
  
  
 
  
  

  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  groupConstraint <- function(groups, rhs = NULL, sense = "<=")
  {
    if (class(groups) != "list") {
      stop("Argument 'groups' has to be a list.")
    }
    if (length(groups) == 0) {
      stop("At least one group has to be provided.")
    }
    Names <- unique(unlist(groups))
    if (sense == "<=") { 
      boundtype <- "_ub"
    } else if (sense == ">=") {
      boundtype <- "_lb"
    } else if (sense == "=") {
      boundtype <- "_eqb"
    } else {
      stop("sense has to be one of '<=', '>=', '='")
    }
    Amat <- matrix(0, nrow = length(groups), ncol = length(Names),
                   dimnames = list(paste0(names(groups), boundtype), 
                                   Names))
    for (i in 1:nrow(Amat)) {
      Amat[i, groups[[i]]] <- 1
    }
    if (is.null(rhs)) {
      if (sense == "<=") rhs <- rep(1, nrow(Amat))
      else rhs = rep(0, nrow(Amat))
    } 
    if (!is.numeric(rhs)) stop("rhs has to be numeric.")
    rhs <- setNames(rhs, rownames(Amat))
    sense <- setNames(rep(sense, nrow(Amat)), rownames(Amat))
    
    ans <- list(name = "group", 
                type = "linear", 
                Amat = Amat, 
                rhs = rhs, 
                sense = sense)
    return(ans)
  }
  
  

  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  bucketorderingConstraint <- function(Names, buckets, eps = 1e-04)
  {
    stopifnot( is.list(buckets) )
    stopifnot( all(unlist(buckets) %in% Names) )
    n <- length(Names)
    m <- 0
    for (i in 1:(length(buckets) - 1)) {
      m <- m + length(buckets[[i]]) * length(buckets[[i+1]])
    }
    Amat <- matrix(0, nrow = m, ncol = n,
                   dimnames = list(NULL, Names))
    l = 1
    for (i in 1:(length(buckets) - 1)) {
      idx_1 <- which(Names %in% buckets[[i]])
      idx_2 <- which(Names %in% unlist(buckets[i + 1]))
      for (j in seq(along = idx_1)) {
        for (k in seq(along = idx_2)) {
          Amat[l, idx_1[j] ] <- 1
          Amat[l, idx_2[k] ] <- -1
          l <- l + 1
        }
      }
    }
    ans <- linearConstraint(name = "bucket_ordering", 
                            rhs = rep(-eps, m), 
                            sense = rep("<=", m), 
                            Amat = Amat)
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  orderingConstraint <- function(Names, ordering, eps = 1e-04)
  {
    
    #//  checkout funciton ordinalConstraint(n, i, j) from hitandrun
    
    
    N <- length(Names)
    Amat <- matrix(0, nrow = (N-1), ncol = N,
                   dimnames = list(NULL, Names))
    for (i in 1:nrow(Amat)) {
      Amat[i, ordering[i]] <- 1
      Amat[i, ordering[i+1]] <- -1
    }
    ans <- linearConstraint(name = "ordering", 
                            rhs = rep(-eps, N-1), 
                            sense = rep("<=", N-1), 
                            Amat = Amat)
    return( ans )
  }
  
  
  