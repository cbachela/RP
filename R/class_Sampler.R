  # 
  # 
  # ############################################################################
  # ### RP PACKAGE - Sampler CLASS
  # ############################################################################
  # 
  # 
  # # --------------------------------------------------------------------------
  # # Cyril Bachelard
  # # This version:     07.01.2020
  # # First version:    03.06.2019
  # # --------------------------------------------------------------------------
  # 
  # 
  # # PolytopeSampler
  # # PolytopeSampler.setCtrl
  # # PolytopeSampler.transform
  # # PolytopeSampler.sample
  # 
  # # EllipsoidSampler
  # # EllipsoidSampler.setCtrl
  # # EllipsoidSampler.transform
  # # EllipsoidSampler.sample
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # PolytopeSampler <- setRefClass( Class = "PolytopeSampler",
  #                                 fields = list( polytope = "Polytope",
  #                                                transformation = "list",
  #                                                spec = "list",
  #                                                samples = "list" ) )
  # 
  # 
  # # --------------------------------------------------------------------------
  # PolytopeSampler.setCtrl <- function( n_sim = 10^2,
  #                                      thin = 1,
  #                                      # burn = 0,
  #                                      algo = c("hitnrun", 
  #                                               "billiard", 
  #                                               "volesti"),
  #                                      volesti = list( density = "uniform", 
  #                                                      walk = "BiW" ),
  #                                      jmp = NULL,
  #                                      interior_point = NULL,
  #                                      boundary_point = NULL,
  #                                      ellipsis = list(),
  #                                      verbose = TRUE )
  # {
  #   ans <- list()
  #   ans$n_sim <- n_sim
  #   ans$thin <- thin
  #   ans$algo <- match.arg(algo)
  #   ans$volesti <- volesti
  #   ans$jmp <- jmp
  #   ans$interior_point <- interior_point
  #   ans$boundary_point <- boundary_point
  #   ans$verbose <- verbose
  #   spec <<- ans
  #   
  #   return( TRUE )
  # }
  # PolytopeSampler$methods( setCtrl = PolytopeSampler.setCtrl )
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # PolytopeSampler.transform <- function()
  # {
  #   d <- polytope$dimension()
  #   idx_eq <- which(polytope$sense == "=")
  #   
  #   if ( length(idx_eq) < 1 ) {
  #     
  #     if ( isTRUE(spec$verbose) ) {
  #       warning("No transformation is performed since the 
  #                linear system contains no equality.")
  #     }
  #     
  #   } else if ( length(idx_eq) == 1 ) {
  #     
  #     # Extract or compute transformation
  #     if ( is.null(transformation$tmat_inv) ) {
  #       transformation <<- simplexTransform( d = d )
  #     }
  #     tmat_inv <- transformation$tmat_inv
  #     z <- transformation$z
  #     
  #     # Create new (d-1) dimensional polytope
  #     polytope <<- polytope$transform( tmat_inv = tmat_inv, 
  #                                      z = z, 
  #                                      idx_eq = idx_eq )
  #     
  #     # Transform initial and boundary point if available
  #     if ( !is.null(spec$interior_point) ) {
  #       spec$interior_point <<- (transformation$tmat %*% spec$interior_point - z)[-1]
  #     }
  #     if ( !is.null(spec$boundary_point) ) {
  #       spec$boundary_point <<- (transformation$tmat %*% spec$boundary_point - z)[-1]
  #     }
  #     
  #   } else {
  #     stop("cannot deal with more than one equality.")
  #   }
  #   
  #   return( TRUE )
  # }
  # PolytopeSampler$methods( transform = PolytopeSampler.transform )
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # PolytopeSampler.sample <- function()
  # {
  #   if ( spec$algo == "billiard" ) {
  #     samples <<- billiard( K = polytope,
  #                           x = spec$interior_point,
  #                           u = NULL,
  #                           n_sim = spec$n_sim + spec$thin,
  #                           jmp = spec$jmp )
  #     
  #   } else if ( spec$algo == "hitnrun" ) {
  #     samples <<- hitnrun( K = polytope, 
  #                          x = spec$interior_point, 
  #                          u = NULL, 
  #                          n_sim = spec$n_sim + spec$thin )
  #   } else if ( spec$algo == "volest" ) {
  #     
  #     ###  ~~~~~~~~~~~~~~~~~~~~~
  #     require(volesti)
  #     P <- Hpolytope$new( A = polytope$A, 
  #                         b = polytope$b )
  #     
  #     lRW <- list( walk = spec$volesti$walk )
  #     if ( !is.null(spec$interior_point) ) {
  #       lRW$starting_point <- spec$interior_point
  #     }
  #     
  #     # tmp <- inner_ball( P = P )
  #     # tmp <- getInteriorPoint( object = RP_prime )
  #     # debugonce( isInterior) 
  #     # isInterior(RP_prime, tmp)
  #     
  #     samples_tmp <- sample_points( P = P,
  #                                   n = spec$n_sim + spec$thin,
  #                                   distribution = list( density = spec$volesti$density ),
  #                                   random_walk = lRW )
  #     samples <<- list( S = t(samples_tmp) )
  #     ### ~~~~~~~~~~~~~~~~~~~~~
  #     
  #   } else {
  #     stop("algorithm '", algo, "' not implemented.")
  #   }
  #   if ( !is.list(samples) ) {
  #     samples <<- list(samples)
  #   }
  #   
  #   # .self$thin()
  #   
  #   return( TRUE )
  # }
  # PolytopeSampler$methods( sample = PolytopeSampler.sample )
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # # --------------------------------------------------------------------------
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # EllipsoidSampler <- setRefClass( Class = "EllipsoidSampler",
  #                                  fields = list( ellipsoid = "Ellipsoid",
  #                                                 transformation = "list",
  #                                                 spec = "list",
  #                                                 samples = "list" ) )
  # 
  # 
  # # --------------------------------------------------------------------------
  # EllipsoidSampler.setCtrl <- function( n_sim = 10^2,
  #                                       thin = 1,
  #                                       # burn = 0,
  #                                       algo = c("sphere",
  #                                                "hitnrun", 
  #                                                "billiard", 
  #                                                "volesti"),
  #                                       volesti = list( density = "uniform", 
  #                                                       walk = "BiW" ),
  #                                       jmp = NULL,
  #                                       b_pushy = TRUE,
  #                                       interior_point = NULL,
  #                                       boundary_point = NULL,
  #                                       ellipsis = list(),
  #                                       verbose = TRUE )
  # {
  #   ans <- list()
  #   ans$n_sim <- n_sim
  #   ans$thin <- thin
  #   ans$algo <- match.arg(algo)
  #   ans$volesti <- volesti
  #   ans$jmp <- jmp
  #   ans$b_pushy <- b_pushy
  #   ans$interior_point <- interior_point
  #   ans$boundary_point <- boundary_point
  #   ans$verbose <- verbose
  #   spec <<- ans
  #   
  #   return( TRUE )
  # }
  # EllipsoidSampler$methods( setCtrl = EllipsoidSampler.setCtrl )
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # EllipsoidSampler.transform <- function()
  # {
  #   d <- ellipsoid$dimension()
  #   
  #   # Extract or compute transformation
  #   if ( is.null(transformation$tmat_inv) ) {
  #     transformation <<- simplexTransform( d = d )
  #   }
  #   tmat <- transformation$tmat
  #   z <- transformation$z
  #   
  #   # Create new (d-1) dimensional ellipsoid
  #   ellipsoid <<- ellipsoid$transform( tmat = tmat, 
  #                                      z = z,
  #                                      boundary_point = spec$boundary_point )
  #   
  #   # Transform initial and boundary point if available
  #   if ( !is.null(spec$interior_point) ) {
  #     as.numeric( tmat %*% spec$boundary_point - z )[-1]
  #     spec$interior_point <<- (tmat %*% spec$interior_point - z)[-1]
  #   }
  #   if ( !is.null(spec$boundary_point) ) {
  #     spec$boundary_point <<- (tmat %*% spec$boundary_point - z)[-1]
  #   }
  #   
  #   return( TRUE )
  # }
  # EllipsoidSampler$methods( transform = EllipsoidSampler.transform )
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # EllipsoidSampler.sample <- function()
  # {
  #   if ( spec$algo == "sphere" ) {
  #     samples$S <<- runifE( S = ellipsoid$shape,
  #                           z_hat = ellipsoid$centre,
  #                           gamma_threshold = ellipsoid$rhs,
  #                           n_points = spec$n_sim + spec$thin,
  #                           b_pushy = spec$b_pushy )
  #   } else if ( spec$algo == "billiard" ) {
  #     samples <<- billiard( K = ellipsoid,
  #                           x = spec$interior_point,
  #                           u = NULL,
  #                           n_sim = spec$n_sim + spec$thin,
  #                           jmp = spec$jmp )
  #     
  #   } else if ( spec$algo == "hitnrun" ) {
  #     samples <<- hitnrun( K = ellipsoid, 
  #                          x = spec$interior_point, 
  #                          u = NULL, 
  #                          n_sim = spec$n_sim + spec$thin )
  #   } else {
  #     stop("algorithm '", spec$algo, "' not implemented.")
  #   }
  #   
  #   # .self$thin()
  #   
  #   return( TRUE )
  # }
  # EllipsoidSampler$methods( sample = EllipsoidSampler.sample )
  # 
  # 
  # 
  # 
  # 
  #   
  #  
  #   
  #   
  #   
  #   
  #   
  #   
  #   
  #   
  #   
  #   
  #   