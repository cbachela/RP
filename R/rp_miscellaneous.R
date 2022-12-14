  
  
  ############################################################################
  ### RP PACKAGE - MISCELLANEUOUS FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     07.07.2019
  # First version:    03.06.2019
  # --------------------------------------------------------------------------
  
  
  # FUNCTIONS:
  #
  # inverse
  # pseudoinverse
  # isPositiveDefinite
  # makePositiveDefinite
  # alignVector
  # simplexTransform
  # eMult
  # ePlus
  # twoAssetPortfolio
  # permutations
  # order.distinct
  # angle
  # rad2deg
  # deg2rad
  # n_grid_points
  # fevBias
  # colorCoding
  # varsi
  # lasserre
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  inverse <- function(A)
  {
    stopifnot( nrow(A) == ncol(A) )
    
    if ( isPositiveDefinite(A) ) {
      B <- solve(A)
      ans <- solve(B %*% A) %*% B
    } else {
      B <- pseudoinverse(A)
      ans <- pseudoinverse(B %*% A) %*% B
    }
    ans <- solve(B %*% A) %*% B
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  pseudoinverse <- function(m) 
  {
    msvd <- svd(m)
    if (length(msvd$d) == 0) 
    {
      m_inv <- array(0, dim(m)[2:1])
    } else 
    {
      s <- matrix(0, nrow = length(msvd$d), ncol = length(msvd$d))
      diag(s)[msvd$d != 0] <- 1 / msvd$d[msvd$d != 0]
      m_inv <- msvd$v %*% (s %*% t(msvd$u))
    }
    return( m_inv )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  isPositiveDefinite <- function (m, tol, method = c("eigen", "chol")) 
  {
    method = match.arg(method)
    if (!is.matrix(m)) 
      m = as.matrix(m)
    if (method == "eigen") {
      eval = eigen(m, only.values = TRUE, symmetric = TRUE)$values
      if (is.complex(eval)) {
        warning("Input matrix has complex eigenvalues!")
        return(FALSE)
      }
      if (missing(tol)) 
        tol = max(dim(m)) * max(abs(eval)) * .Machine$double.eps
      if (sum(eval > tol) == length(eval)) 
        return(TRUE)
      else return(FALSE)
    }
    if (method == "chol") {
      val = try(chol(m), silent = TRUE)
      if (class(val) == "try-error") 
        return(FALSE)
      else return(TRUE)
    }
  }
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  makePositiveDefinite <- function (m, tol)
  {
    if (!is.matrix(m)) 
      m = as.matrix(m)
    d = dim(m)[1]
    if (dim(m)[2] != d) 
      stop("Input matrix is not square!")
    es = eigen(m, symmetric = TRUE)
    esv = es$values
    if (missing(tol)) 
      tol = d * max(abs(esv)) * .Machine$double.eps
    delta = 2 * tol
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(m + dm)
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  alignVector <- function(y, x)
  {
    
    # Given vectors y and x \in \mathds{R}, find orthogonal matrix T that 
    # rotates x such that it is parallel to y.
    # i.e. find T (TT^T = T^TT = I) s.t. Tx = ay, 
    # with scalar a = |x|/|y|. |.| denotes the euclidean norm 
    # (i.e. 2-norm, which corresponds to the largest singular value)
    
    
    y = matrix(y, ncol = 1)
    x = matrix(x, ncol = 1)
    
    SVD_y = svd(y, nu = nrow(y), nv = ncol(y))
    SVD_x = svd(x, nu = nrow(x), nv = ncol(x))
    
    ans = SVD_y$u %*% t(SVD_x$u)
    
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  simplexTransform <- function(d)
  {
    v <- rep(1, d)
    b <- 1
    y <- c(1, rep(0, d-1))
    tmat <- alignVector(y = y, x = v)
    tmat_inv <- inverse(tmat)
    z <- as.numeric(b * tmat %*% v) / as.numeric( t(v) %*% v )
    ans <- list(tmat = tmat,
                tmat_inv = tmat_inv,
                z = z)
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  eMult <- function(E, Amat = NULL)  
  {
    
    # A * E(q, Q) = E(Aq, AQA^T)
    
    if ( is.null(Amat) ) {
      return( E )
    }
    
    qvec = getCentre(E)
    Qmat = getShape(E)  
    qvec_new <- Amat %*% qvec
    Qmat_new <- Amat %*% Qmat %*% t(Amat)
    Qmat_new <- 0.5 * (Qmat_new + t(Qmat_new))
    
    ans <- ellipsoid(q = qvec_new,
                     Q = Qmat_new)
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  ePlus <- function(E, bvec)
  {
    # Operation E + b where E - ellipsoid and b - vector in R^n.
    # E(q, Q) + b = E(q + b, Q)
    
    ans <- ellipsoid(q = E@q + bvec,
                     Q = E@Q)
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  pMult <- function(P, Gmat = NULL)  
  {
    
    # G * P(A, b) = P(AG^T, b - AG^T)
    
    if ( is.null(Gmat) ) {
      return( P )
    }
    
    b = getRHS(P)
    Amat = getAmat(P)  
    Amat_new <- Amat %*% t(Gmat)
    # b_new <- b - Amat_new %*% z
    b_new <- b
    
    ans <- polytope(A = Amat_new,
                    b = b_new,
                    sense = getSense(P))
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  twoAssetPortfolio <- function(w1, w2, covmat, alpha)
  {
    w_delta <- w2 - w1
    sigma_1 <- t(w1) %*% covmat %*% w1
    sigma_12 <- t(w_delta) %*% covmat %*% w_delta
    num <- t(-w_delta) %*% covmat %*% w1 + 
      sqrt( (t(w_delta) %*% covmat %*% w1)^2 + alpha * sigma_12 * sigma_1 )
    lambda <- as.numeric( num / sigma_12 )
    w <- w1 + lambda * w_delta
    ans <- list(w = w,
                lambda = lambda)
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  permutations <- function(x) 
  {
    if (length(x) == 1) {
      return(x)
    }
    else {
      res <- matrix(nrow = 0, ncol = length(x))
      for (i in seq_along(x)) {
        res <- rbind(res, cbind(x[i], Recall(x[-i])))
      }
      return(res)
    }
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  order.distinct <- function(x)
  {
    ans <- numeric(length(x))
    dupli <- duplicated(x)
    ordering <- order(x[!dupli])
    ans[!dupli] <- ordering
    dupli_distinct <- unique(x[dupli])
    for (i in seq(along = dupli_distinct)) {
      if ( is.na(dupli_distinct[i]) ) {
        idx <- which(is.na(x))
      } else {
        idx <- which(x == dupli_distinct[i])
      }
      ans[idx] <- rep(ordering[idx[1]], length(idx))
    }
    return( ans )
  }
  
  
  
  # Geometry
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  angle <- function(v1, v2) { acos(sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))) }
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  rad2deg <- function(rad) { (rad * 180) / pi }
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  deg2rad <- function(deg) { (deg * pi) / 180 }
  
  # --------------------------------------------------------------------------
  # N is the number of assets
  # n is the number of points
  # --------------------------------------------------------------------------
  n_grid_points <- function ( N, n ) {
    return( choose( N + n - 2, n - 1) )
  }
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  fevBias <- function(x, q = 4)
  { 
    
    if ( is.null(dim(x)) ) {
      ans <- x^q / sum(x^q) 
    } else {
      ans <- t( apply(x, 1, function(x) { x^q / sum(x^q ) }) )
      colnames(ans) <- colnames(x)
    }
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  colorCoding <- function(Set, 
                          obj_values = NULL, 
                          FUN, 
                          col = c("lightblue", "darkblue"),
                          n_steps = 5)
  {
    txt <- colorRampPalette(colors = col)(n_steps)
    n <- nrow(Set) / length(txt)
    colors <- rep(NA, nrow(Set))
    for ( i in seq(along = txt) ) {
      colors[((n*i)-n+1):(n*i)] <- rep(txt[i], n)
    }
    color <- rep(NA, nrow(Set))
    if ( is.null(obj_values) ) {
      obj_values <- apply(Set, 1, FUN)
    }
    color[order(obj_values)] <- colors
    
    return( color )
  }
  
  
  # --------------------------------------------------------------------------
  #' Varsi's algorithm
  #' @export
  # --------------------------------------------------------------------------
  varsi <- function( mu, b )
  {
    z <- mu - b
    idx_neg <- which(z < 0)
    idx_nonneg <- which(z >= 0)
    J <- length(idx_neg)
    K <- length(idx_nonneg)
    x <- 0
    y <- 0
    if ( J > 0 ) {
      x <- z[idx_neg]
    }
    if ( K > 0 ) {
      y <- z[idx_nonneg]
    }
    a <- c(1, rep(0, K))
    if ( J > 0 ) {
      A <- matrix( NA, nrow = J, ncol = K+1 )
      for ( j in 1:J ) {
        for ( k in 1:K ) {
          a[k+1] <- (y[k] * a[k+1] - x[j] * a[k]) / (y[k] - x[j])
        }
        A[j, ] <- a
      }
    } else {
      A <- a
    }
    
    return( A )
  }

  # varsi <- function( mu, b )
  # {
  #   z <- mu - b
  #   idx_neg <- which(z < 0)
  #   idx_nonneg <- which(z >= 0)
  #   J <- length(idx_neg)
  #   K <- length(idx_nonneg)
  #   if ( J > 0 ) {
  #     x <- z[idx_neg]
  #   }
  #   if ( K > 0 ) {
  #     y <- z[idx_nonneg]
  #   }
  #   a <- c(1, rep(0, K))
  #   if ( J > 0 ) {
  #     for ( j in 1:J ) {
  #       for ( k in 1:K ) {
  #         a[k+1] <- (y[k] * a[k+1] - x[j] * a[k]) / (y[k] - x[j])
  #       }
  #     }
  #   }
  # 
  #   return( tail(a, 1) )
  # }
  
  
  # --------------------------------------------------------------------------
  #' Lasserre's analytic solution
  #' @description Implements Theorem 2.2 in Lasserre (2015)
  #' @references J. B. Lasserre (2015). Volume of slices and section of the simplex in closed form
  #' @export
  # --------------------------------------------------------------------------
  lasserre <- function( a, th )
  {
    
    N <- length(a)
    a <- c(0, a)
    ans <- 0
    for ( i in 1:length(a) ) {
      num <- max(0, th - a[i])^N
      denom <- prod( a[-i] - a[i] )
      ans <- ans + num / denom  
    }
    # ans <- ans  / factorial(N)
    return( ans )    
  }
    
  
  test.lasserre <- function()
  {
    require(RP)
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    b <- quantile(x, 0.6)
    
    tail( as.numeric( varsi( mu = x, b = b ) ), 1 )
    lasserre( a = x, th = b ) 
    lasserre( a = x, th = b ) * factorial(n)
      
    
    
  }
  
  
  
  
  