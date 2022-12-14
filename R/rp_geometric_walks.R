  
  
  ############################################################################
  ### RP PACKAGE - GEOMETRIC WALKS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  #
  # FUNCTIONS:
  #
  # rpCtrl
  # hitnrun
  # har
  # billiard
  # dists2Ellipsoid
  # dists2Polytope
  # dists2ConvexBody
  # isInterior
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  rpCtrl <- function( n_sim = 10^2,
                      thin = 1,
                      # burn = 0,
                      algo = c("hitnrun", "billiard", "sphere", "volesti"),
                      volesti = list( density = "uniform", 
                                      walk = "BiW" ),
                      jmp = NULL,
                      x0 = NULL,
                      x_boundary = NULL,
                      b_pushy = FALSE,
                      transformation = list(tmat = NULL, tmat_inv = NULL, z = NULL),
                      ellipsis = list(),
                      verbose = TRUE )
  {
    
    spec <- list()
    spec$n_sim <- n_sim
    spec$thin <- thin
    spec$algo <- match.arg(algo)
    spec$volesti <- volesti
    spec$jmp <- jmp
    spec$x0 <- x0
    spec$x_boundary <- x_boundary
    spec$b_pushy <- b_pushy
    spec$transformation <- transformation
    spec$verbose <- verbose
    
    return( spec )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  hitnrun <- function(B, x, u = NULL, n_sim = 100)
  {
    
    if ( is.null(x) ) {
      # stop("point x within the body is missing.")
      x <- getInteriorPoint(B)
    }
    
    # stopifnot( isInterior(B, x) )
    if ( !isTRUE(isInterior(B, x)) ) {
      x <- x * 0.9                     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    }
    stopifnot( isInterior(B, x) )
    
    
    if ( is.null(u) ) {
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
    }
    
    ans <- matrix( 0, nrow = n_sim, ncol = length(x) )
    # U <- ans
    B_lo <- ans
    B_up <- ans
    
    itermax <- n_sim * 100
    i <- iter <- 1
    while ( i <= n_sim && iter < itermax ) {
      
      if ( isInterior(B, x) ) {
        
        D <- dists2Body(object = B, 
                        x = x, 
                        u = u)
        
        if ( sum(abs(D)) > 0 ) {
          
          # lower <- D["lower"]
          # upper <- D["upper"]
          lower <- D[1]
          upper <- D[2]
          
          # U[i, ] <- u
          B_lo[i, ] <- x + lower * u
          B_up[i, ] <- x + upper * u
          
          x <- x + (lower + (upper - lower) * runif(1)) * u
          ans[i, ] <- x
          i <- i + 1
        }
      }
      
      # Next random direction
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
      iter <- iter + 1
      
    }
    
    out <- list( S = ans,
                 # U = U,
                 B_lo = B_lo,
                 B_up = B_up)
    
    return( out ) 
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  har <- function( K, x = NULL, u = NULL, n_sim = 100 )
  {
    
    if ( is.null(x) ) {
      stop("interrior point is missing.")
      # x <- K$getInteriorPoint()
    }
    
    if ( !K$isInterior(x = x) ) {
      x <- x * 0.9                     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    }
    stopifnot( K$isInterior(x = x) )
    
    
    if ( is.null(u) ) {
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
    }
    
    ans <- matrix( 0, nrow = n_sim, ncol = length(x) )
    # U <- ans
    B_lo <- ans
    B_up <- ans
    
    itermax <- n_sim * 100
    i <- iter <- 1
    while ( i <= n_sim && iter < itermax ) {
      
      if ( K$isInterior(x) ) {  # is this necessary? only if we don't trust boundaryDistance
        
        D <- K$boundaryDistance(x = x, u = u)
        
        if ( sum(abs(D)) > 0 ) {
          
          # lower <- D["lower"]
          # upper <- D["upper"]
          lower <- D[1]
          upper <- D[2]
          
          # U[i, ] <- u
          B_lo[i, ] <- x + lower * u
          B_up[i, ] <- x + upper * u
          
          x <- x + (lower + (upper - lower) * runif(1)) * u
          ans[i, ] <- x
          i <- i + 1
        }
      }
      
      # Next random direction
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
      iter <- iter + 1
      
    }
    
    out <- list( S = ans,
                 # U = U,
                 B_lo = B_lo,
                 B_up = B_up )
    
    return( out ) 
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  billiard <- function( K, x, u = NULL, jmp = NULL, 
                        jmp_adaptive = FALSE, n_sim = 100 )
  {
    
    if ( is.null(x) ) {
      stop("interior is missing.")
      # x <- getInteriorPoint(K)
    }
    
    if ( !K$isInterior(x = x) ) {
      x <- x * 0.9                     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    }
    stopifnot( K$isInterior(x = x) )
    
    if ( is.null(jmp) ) {
      # jmp <- rep(1, length(x))    # ~~~~~~~~~~~~~~~~~~ Write function to get reasonable jmp size
      jmp <- 1
    }
    
    A <- K$A
    b <- K$b 
    x1 <- x
    # x2 <- rnorm(length(x), x, jmp)
    tau <- -jmp * log(runif(1))    # ~~~~~~~~~~~~~~~ 
    u <- rnorm(length(x))
    u <- u / as.numeric(sqrt(sum(u^2)))
    x2 <- x + u * tau 
    
    if ( any(is.na(x2)) )  stop("Direction vector contains NA's")
    x2mat <- matrix(x2, nrow = 1)  # ~~~~~~~~~~~~~~~~~~
    x1mat <- matrix(x1, nrow = 1)  # ~~~~~~~~~~~~~~~~~~
    ans <- matrix( 0, nrow = n_sim, ncol = length(x) )
    boundarypoints <- ans
    itermax <- n_sim * 100
    i <- iter <- 1
    
    ##
    lK <<- list()
    lJmp <<- list()
    ##
    
    
    while ( i <= n_sim && iter < itermax ) {
      
      # if ( isInterior(B, x) ) {
        
        residual <- A %*% x2 - b
        idx <- which(residual > 0)
        
        k <- 1
        while ( length(idx) > 0 ) {
         
          epsilon <- x2 - x1
          alpha <- ( (b[idx] - A[idx, ] %*% x1) /  A[idx, ] %*% epsilon)
          idx_minalpha <- which.min(alpha)
          j <- idx[idx_minalpha]
          d <- -residual[j] / sum(A[j, ]^2)
          x2 <- x2 + 2 * d * A[j, ]
          x2mat <- rbind(x2mat, x2)   # ~~~~~~~~~~~~~~~~  
          residual <- A %*% x2 - b
          idx <- which(residual > 0)
          x1 <- x1 + alpha[idx_minalpha] * epsilon
          x1mat <- rbind(x1mat, x1)    # ~~~~~~~~~~~~~~~~
          k <- k + 1
          
          # if ( k > 10 ) {
          #   # Do hit-and-run instead
          #   D <- dists2Body(object = B, 
          #                   x = x, 
          #                   u = x2)
          #   if ( sum(abs(D)) > 0 ) {
          #     lower <- D[1]
          #     upper <- D[2]          
          #     x2 <- x + (lower + (upper - lower) * runif(1)) * x2
          #     idx <- NULL
          #   }
          # }
          
        }
        
        ans[i, ] <- x2
        boundarypoints[i, ] <- x1
        i <- i + 1
      # }
      
        if ( isTRUE(jmp_adaptive) ) {
          jmp <- 3 * jmp / k      # ~~~~~~~~~~~~~~~~~~~~ ???
          lK[[i-1]] <<- k
          lJmp[[i-1]] <<- jmp
        }
      
        
      # Next random direction
      # x2 <- rnorm(length(x1), x, jmp)
      tau <- -jmp * log(runif(1))    # ~~~~~~~~~~~~~~~ 
      u <- rnorm(length(x))
      u <- u / as.numeric(sqrt(sum(u^2)))
      x2 <- x2 + u * tau 
      iter <- iter + 1
    }
    
    out <- list( S = ans,
                 # x1mat = x1mat )
                 boundarypoints = boundarypoints )
    
    return( out ) 
  }
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  dists2Ellipsoid <- function(E, x, u)
  {
    
    if ( !isEmpty(E) ) {
      
      y <- (x - getCentre(E))
      Sigma <- getShape(E)
      sigma_bar <- getRHS(E)
      a <- t(u) %*% Sigma %*% u
      b <- 2 * t(u) %*% Sigma %*% y
      g <- t(y) %*% Sigma %*% y - sigma_bar
      d <- sqrt( b^2 - 4 * a * g )
      lambda_minus <- (-b - d) / (2 * a)
      lambda_plus <- (-b + d) / (2 * a)
      lambda <- sort( c(lambda_minus, lambda_plus) )
      
    } else {
      lambda <- NULL
    }
    
    return( lambda )
  }
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  dists2Polytope <- function(P, x, u)
  {
    
    if ( !isEmpty(P) ) {
      
      num <- P@b - P@A %*% x
      denom <- P@A %*% u
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
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  dists2ConvexBody <- function(B, x, u)
  {
    lambda_p <- dists2Polytope(P = getPolytope(B), x = x, u = u)
    lambda_e <- dists2Ellipsoid(E = getEllipsoid(B), x = x, u = u)
    lower <- max(lambda_p[1], lambda_e[1])
    upper <- min(lambda_p[2], lambda_e[2])
    # ans <- setNames( c(lower, upper), c("lower", "upper") )
    ans <- c(lower, upper)
    return( ans )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  isInterior <- function(object, x)
  {
    
    eps <- 1e-08
    
    if ( inherits(object, "POLYTOPE") ) {
      P <- getPolytope(object)
      if ( !isEmpty(P) ) {
        if ( any(  getRHS(P) - getAmat(P) %*% x < -eps) ) {
          return( FALSE )
        }
      }
    }
    
    if ( inherits(object, "ELLIPSOID") ) {  
      E <- getEllipsoid(object)
      if ( !isEmpty(E) ) {
        if ( getRHS(E) - t(x - getCentre(E)) %*% getShape(E) %*% (x - getCentre(E)) < -eps ) {
          return( FALSE )
        }
      }
    }
    
    return( TRUE )
  }
  
  
  
  
  