  
  
  ############################################################################
  ### CalesChalkisEmiris2019_On the cross-sectional distribution of portfolios returns
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     26.02.2020
  # First version:    26.02.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(boot)
  require(volesti)
  require(rgl)
  require(BCDating)
  require(GPO)
  require(slolz)
  require(momolz)
  require(RP)
  source("H:/R/RP/R/class_ellipsoid.R")
  source("H:/R/RP/R/class_RPE.R")
  
  
  
  wd <- "R:/Asset_Management/Research Projects/Factor_Optimization/"
  r_functions <- paste0(wd, c("1_-_Single_Factor/2_-_Portfolio/Consistent/R_Functions.R",
                              "1_-_Single_Factor/2_-_Portfolio/Consistent/R_Functions_val.R",
                              "1_-_Single_Factor/2_-_Portfolio/Consistent/R_Functions_qlt.R",
                              "1_-_Single_Factor/2_-_Portfolio/Consistent/R_Functions_mom.R",
                              "2_-_Factor_Mix/R_Functions_mix.R"))
  sapply(r_functions, source)
  
  source("H:/Papers Cyril/Seminary_Statistical Learning/kcde_functions.R")
  source("R:/Asset_Management/Research Projects/Factor_Optimization/wip_cyril/multifac_functions.R")
  source("R:/Asset_Management/Research Projects/External Research Projects/Portfolios from Views/Portfolios from Sorts/Sorts.R")
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Replicate example in 
  # CalesChalkisEmiris2019_On the cross-sectional distribution of portfolios returns
  # --------------------------------------------------------------------------
  
  mu <- c(0, 0.01, 0.015)
  
  d <- length(mu)
  n_sim <- 10^4
  Simplex <- simplex(d = d)
  RPS <- rp(Simplex, spec = rpCtrl(n_sim = n_sim))
  samples_S <- getSamples(RPS)
  
  mu_samples <- apply( samples_S, 1, function(x) { t(x) %*% mu } )
  mu_score <- apply(as.matrix(mu_samples, ncol = 1), 1, function(x) { sum(x > mu_samples) / n_sim })
  
  
  plot( x = mu_samples, y = mu_score )
  abline( h = 0.5 )
  abline( v = mu_samples[ which(mu_score == 0.88) ] )
  
  
  

    
  
  plot(density(mu_samples, from = min(mu), to = max(mu)))  
  for ( i in 1:length(mu) ) {
    abline(v = mu[i], col = i)
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  universe <- "usa"
  filename = "MSCI_Model"
  
  
  # Factor simulations - i.e.  OLZ max exposure portfolios
  path <- "R:/Asset_Management/Research_Projects/Equity/Factor_Optimization/Data/max_exposure_series.Rds"
  env_fs <- readRDS(file = path)
  tmp <- env_fs$lSim_nofc[[universe]]
  idx <- grepl("msci_", colnames(tmp))
  X <- tmp[ ,-which(idx == TRUE)]
  X <- X[ ,-which(colnames(X) == "momentum_vol_st")]
  bm <- tmp[ ,"msci_capw"]
  colnames(bm) <- "capw"
  X_fact <- cbind(bm, X)
  X_fact_w <- aggWeekly(X_fact, 
                        compounding = "discrete", 
                        day = "Wed")
  
  
  # Pairwise excess returns (Delta portfolios)
  ## Compute differences in discrete returns between all pairs of factor portfolios
  X_w <- X_fact_w[ ,-which(colnames(X_fact_w) == "capw")]
  K <- (ncol(X_w) * (ncol(X_w) + 1)) / 2 - ncol(X)
  Delta_w <- X_w[ ,rep(1, K)] * NA
  Names <- Names_deltas <- NULL
  k <- 1
  for (i in 1:ncol(X_w)) {
    for (j in 1:ncol(X_w)) {
      if ( j > i ) {
        Names_deltas <- cbind(Names_deltas, paste0(colnames(X)[i], "-", colnames(X)[j]))
        Names <- cbind(Names, paste0(i, "-", j))
        Delta_w[ ,k] <- X_w[ ,i] - X_w[ ,j]
        k <- k + 1
      }
    }
  }
  colnames(Delta_w) <- Names
  
  
  
  
  
  data_ret <- readReturnsRDS(filename = "MSCI_Model", 
                             ccy_est = "LOC", 
                             ccy_sim = "USD")
  
  
  

  
  
  # # Setting
  # lSettings <- universe2lSettings(universe)
  # lSettings$Momentum <- momCtrl(method = "cumretEwma",
  #                               compounding = "continuous",
  #                               scale = NULL,
  #                               width = 52,
  #                               lag = 0,
  #                               zscore_flag = TRUE,
  #                               ellipsis = list(tau = 104))
  # selection_fun <- "db_flag"
  # cons_fun <- c("box", "country", "sector")
  # 
  # # Constraints
  # data_index <- readIndexRDS(universe = lSettings$universe)
  # lSettings$data_settings <- data_index$settings
  # data_cons <- prepareConstraintsData(data_index = data_index$X)
  # data_cons$data_gap <- rodbc.dataGap(filename = lSettings$filename)
  # 
  # # ESG
  # data_esg <- prepareESGData(data_index = data_index$X,
  #                            Industry_Adjusted_Score_min = 5.7,
  #                            Iva_Company_Rating_excl = c(NA, "CCC", "B"),
  #                            UNGC_Compliance_excl = c(NA, "Fail"),
  #                            Overall_Flag_excl = c(NA, "Red"),
  #                            SVVK_excl = c(NA, "Fail"))
  # 
  # # Get rebalancing dates
  # rebdates <- rodbc.datesBacktest(btID = lSettings$btID, filename = filename)
  # rebdates <- intersect(as.character(rebdates$Rebal_Date), rownames(data_ret$X_est))
  # rebdates <- intersect(rebdates, rownames(na.omit(data_esg$esg_LB)))

  
  
  
  # # Regression variance vs esg on single stock level
  # width <- 260
  # for ( today in rebdates ) {
  #   
  #   selection <- colnames(data_esg$esg_mat)[which(!is.na(data_esg$esg_mat[today, ]))]
  #   esg_scores <- data_esg$esg_mat[today, selection]
  #   lridx <- which(rownames(data_ret$X_est) == today)
  #   if ( width == 0 | width > lridx ) {
  #     fridx <- 1
  #   } else {
  #     fridx <- max(1, lridx - width)
  #   }
  #   X_est <- data_ret$X_est[fridx:lridx, ]
  #   
  #   selection <- intersect(colnames(X_est), colnames(esg_scores))
  #   X_est <- X_est[ ,selection]
  #   esg_scores <- setNames( as.numeric(esg_scores[ ,selection]), selection )
  #   sds <- apply(X_est, 2, sd)
  #   
  #   
  #   plot( x = esg_scores, y = sds )
  #   
  #   
  # }
  
  
  
  
  
  
  
  
  
    
  
  # RP
  X <- X_fact[ ,c("value", "quality", "lowsize", "momentum", "lowvola")]
  d <- ncol(X)
  n_sim <- 10^4
  Simplex <- simplex(d = d)
  RPS <- rp(Simplex, spec = rpCtrl(n_sim = n_sim))
  samples_S <- getSamples(RPS)
  samples_S_fev <- fevBias( x = samples_S, q = 4 )
  

  
  mu <- apply(X, 2, mean)
  sds <- apply(X, 2, sd)
  covmat <- cov(X)
  
  mu_samples <- apply( samples_S_fev, 1, function(x) { t(x) %*% mu } )
  mu_score <- apply(as.matrix(mu_samples, ncol = 1), 1, function(x) { sum(x > mu_samples) / n_sim })
  sds_samples <- apply( samples_S_fev, 1, function(x) { t(x) %*% covmat %*% x } )
  sds_score <- apply(as.matrix(sds_samples, ncol = 1), 1, function(x) { sum(x > sds_samples) / n_sim })
  
  
  PDF <- density(mu_samples)
  CDF <- cumsum( PDF$y / sum(PDF$y) )
  
  
  plot( x = mu_samples, y = mu_score )
  abline( h = 0.5 )
  abline( v = mu_samples[ which(mu_score == 0.5) ] )
  
  
  plot( x = sds_samples, y = sds_score )
  abline( h = 0.5 )
  abline( v = sds_samples[ which(sds_score == 0.5) ] )
  
  
  mu_samples[ which(mu_score == 0.5) ]
  median(mu_samples)  
  
  
  
  q_vec <- seq(from = 0.01, to = 1, by = 0.01)
  a_vec <- q_vec * NA
  for ( i in seq(along = q_vec) ) {
    A <- varsi( mu = mu, b = quantile(mu, q_vec[i]) )
    a_vec[i] <- tail(as.numeric(A), 1)  
  }
  
    
  
  
  plot( x = mu_samples, y = mu_score )
  points( x = quantile(mu, q_vec), y = a_vec, col = 2 )
  
  
  
  plot( x =  sds_samples, y = mu_samples )
  
  
  
  
  # --------------------------------------------------------------------------
  # Varsi vs. (Bayesian) bootstrap
  # --------------------------------------------------------------------------

  require(volesti)
  require(RP)
  
  d <- 10^3
  n_sim <- 10^4
  
  set.seed(1234)
  # mu <- rgamma(n = d, shape = 2, rate = 0.1)
  # mu1 <- rgamma(n = d, shape = 20, rate = 0.1)
  # mu2 <- rgamma(n = d, shape = 0.1, rate = 1)
  # mu <- sample( c(mu1, mu2) )[1:d]
  
  mean_true <- 0
  sd_true <- 1
  x <- rnorm(n = 10^8, mean = mean_true, sd = sd_true) / 100
  mu <- rnorm(n = d, mean = mean_true, sd = sd_true) / 100
  
  plot(density(mu))
  
  
  S <- RP::simplex(d = d)
  RPS <- rp(S, spec = rpCtrl(n_sim = n_sim))
  samples_S <- getSamples(RPS)
  samples_S_fev <- fevBias( x = samples_S, q = 40 )
  samples_P <- rdirichlet( n = n_sim, alpha = rep(1 / d, d) )
  samples_P2 <- rdirichlet( n = n_sim, alpha = rep(1 / d, d) * d )
  b <- quantile( mu, 0.5 )
  H <- hyperplane( v = mu, b = b )
  
  # debugonce(varsi)
  A <- varsi( mu = mu, b = H@b )
  A
  B <- frustum_of_simplex( a = mu, z0 = H@b )
  B
  tail(as.numeric(A), 2)
  
  
  
  # Loop over quantiles to obtain cdf
  # q_vec <- seq(from = 0.01, to = 1, by = 0.01)
  q_vec <- seq( from = 0.45, to = 0.65, length.out = 10^3 )
  a_vec <- q_vec * NA
  for ( i in seq(along = q_vec) ) {
    A <- varsi( mu = mu, 
                b = quantile(mu, q_vec[i]) )
    a_vec[i] <- tail(as.numeric(A), 1)  
  }
  FUN <- function(z0) { frustum_of_simplex( a = mu, z0 = z0) }
  b_vec <- unlist( lapply( quantile(mu, q_vec), FUN ) )
  
  plot( x = quantile(mu, q_vec), y = b_vec )
  points( x = quantile(mu, q_vec), y = a_vec, col = 2 )
  
  
  
  
  # sampling distribution of the variance:
  # q_vec <- seq(from = 0, to = 1, length.out = 10^3)
  # tmp <- unlist( lapply( quantile((mu - mean_true)^2, q_vec), FUN ) )
  # plot( x = quantile((mu - mean_true)^2, q_vec), y = tmp, type = "o" )
  
  
  
  # Bayesian bootstrap
  mu_samples <- apply( samples_S, 1, function(x) { t(x) %*% mu } )
  mu_score <- apply(as.matrix(mu_samples, ncol = 1), 1, function(x) { sum(x > mu_samples) / n_sim })
  mu_samples_fev <- apply( samples_S_fev, 1, function(x) { t(x) %*% mu } )
  mu_score_fev <- apply(as.matrix(mu_samples_fev, ncol = 1), 1, function(x) { sum(x > mu_samples_fev) / n_sim })
  mu_samples_P <- apply( samples_P, 1, function(x) { t(x) %*% mu } )
  mu_score_P <- apply(as.matrix(mu_samples_P, ncol = 1), 1, function(x) { sum(x > mu_samples_fev) / n_sim })
  mu_samples_P2 <- apply( samples_P2, 1, function(x) { t(x) %*% mu } )
  mu_score_P2 <- apply(as.matrix(mu_samples_P2, ncol = 1), 1, function(x) { sum(x > mu_samples_fev) / n_sim })
  
  sd(mu_samples)
  sd_true / d
  
  
  # Ordinary bootstrap
  Boot <- boot::boot( data = mu,
                      statistic = meanBoot,
                      R = n_sim )
  
  PDF_boot <- density(Boot$t)
  CDF_boot <- cumsum( PDF_boot$y / sum(PDF_boot$y) )
  
  PDF <- density(mu_samples)
  CDF <- cumsum( PDF$y / sum(PDF$y) )
  
  PDF_mu <- density(mu)
  CDF_mu <- cumsum( PDF_mu$y / sum(PDF_mu$y) )
  
  PDF_P <- density(mu_samples_P)
  CDF_P <- cumsum( PDF_P$y / sum(PDF_P$y) )
  PDF_P2 <- density(mu_samples_P2)
  CDF_P2 <- cumsum( PDF_P2$y / sum(PDF_P2$y) )
  
  
  # Compare results
  plot( x = quantile(mu, q_vec), y = a_vec, type = "o" )
  points( x = quantile(mu, q_vec), y = b_vec, col = "yellow" )
  points( x = mu_samples, y = mu_score, col = 2 )
  points( x = mu_samples_fev, y = mu_score_fev, col = 3 )
  points( x = mu_samples_P, y = mu_score_P, col = "orange" )
  points( x = mu_samples_P2, y = mu_score_P2, col = "tomato" )
  points( x = PDF_P$x, y = CDF_P, col = "brown" )
  points( x = PDF_P2$x, y = CDF_P2, col = "brown" )
  # abline( v = mean(mu_samples) )
  # abline( v = median(mu_samples) )
  # points( x = PDF$x, y = CDF, col = 4 )
  points( x = PDF_mu$x, y = CDF_mu, col = 5 )
  points( x = PDF_boot$x, y = CDF_boot, col = 6 )
  points( x = quantile(mu, q_vec), y = a_vec, type = "o" )  
  

  plot( density(mu), col = 4 )  
  lines( density(mu_samples), col = 1 )  
  lines( density(mu_samples_fev), col = 3 )
  
  
  # Numerical derivative
  d_vec <- b_vec * NA
  x_vec <- quantile(mu, q_vec)
  for ( i in 3:length(b_vec) ) {
    d_vec[i-1] <- (b_vec[i] - b_vec[i-2]) / abs(x_vec[i] - x_vec[i-2])
  }
  plot( x = x_vec, y = d_vec )
  
  # Faster:
  idx <- 3:length(b_vec)
  d_vec <- (b_vec[idx] - b_vec[idx-2]) / abs(x_vec[idx] - x_vec[idx-2])
 
  plot( x = x_vec[idx-1], y = d_vec )
  points( x = PDF_boot$x, y = PDF_boot$y, col = 4 )
  
  
  
  
  # Note: 
  # Dirichlet distribution is log-concave. Can we apply Varsi on dirichlet dist?
  # Probably not, but can use MCMC (see Chalkis Emeris 2020).
  # --> conditional density estimation, posterior mean and stuff inbetween
  
  
  
  
  # --------------------------------------------------------------------------
  # Sampling distribution of variance of normally distributed r.v.
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(slolz)
  require(RP)
  
  d <- 5
  n_sim <- 10^4
  
  mean_true <- 0
  sd_true <- 1
  x_pop <- rnorm(n = 10^8, mean = mean_true, sd = sd_true)
  x <- rnorm(n = d, mean = mean_true, sd = sd_true)
  z <- (x - mean_true)^2
  
  
  # Varsi
  q_vec <- seq( from = 0.01, to = 0.99, length.out = 10^3 )
  FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  p_vec <- unlist( lapply( quantile(z, q_vec), FUN ) )
  
  # Derive density
  x_vec <- quantile(z, q_vec)
  idx <- 3:length(p_vec)
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(x_vec[idx] - x_vec[idx-2])
  
  
  plot( x = quantile(z, q_vec), y = p_vec )
  plot( x = x_vec[idx-1], y = d_vec )
  
  
  # Bootstrap
  Boot <- boot::boot( data = z,
                      statistic = meanBoot,
                      R = n_sim )
  Boot <- boot::boot( data = x, 
                      # data = x_pop,
                      statistic = function(x, idx)
                      {
                        if ( NCOL(x) > 1 ) {
                          ans <- apply(x[idx, ], 2, var)
                        } else {
                          ans <- var(x[idx])
                        }
                        return( ans )
                      },
                      R = n_sim )
  
  pdf_boot <- density(Boot$t)
  cdf_boot <- cumsum( pdf_boot$y / sum(pdf_boot$y) )
  
  
  Boot$t0; sd_true^2
  
  
  
  # Plots
  
  plot( x = quantile(z, q_vec), y = p_vec )
  plot( x = x_vec[idx-1], y = d_vec )
  
  
  plot( x = pdf_boot$x, y = cdf_boot )
  points(  x = quantile(z, q_vec), y = p_vec, col = 2 )

  plot( pdf_boot )
  points( x = x_vec[idx-1], y = d_vec, col = 2 )
  abline( v = Boot$t0 )
  abline( v = sd_true^2 )  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Sampling distribution of variance of heavy tailed distributed r.v.
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(slolz)
  require(RP)
  
  d <- 5000
  n_sim <- 10^4
  
  mean_true <- 0
  sd_true <- 1
  x_pop <- rcauchy( n = 10^8, location = mean_true, scale = 100 )
  x <- rcauchy( n = d, location = 0, scale = 100 )
  z <- (x - mean_true)^2
  
  
  # Varsi
  q_vec <- seq( from = 0.00001, to = 0.999999, length.out = 10^3 )
  FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  p_vec <- unlist( lapply( quantile(z, q_vec), FUN ) )
  
  # Derive density
  x_vec <- quantile(z, q_vec)
  idx <- 3:length(p_vec)
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(x_vec[idx] - x_vec[idx-2])
  
  
  plot( x = quantile(z, q_vec), y = p_vec )
  plot( x = x_vec[idx-1], y = d_vec )
  
  
  # Bootstrap
  Boot <- boot::boot( data = z,
                      statistic = meanBoot,
                      R = n_sim )
  Boot <- boot::boot( data = x,
                      statistic = function(x, idx)
                      {
                        if ( NCOL(x) > 1 ) {
                          ans <- apply(x[idx, ], 2, var)
                        } else {
                          ans <- var(x[idx])
                        }
                        return( ans )
                      },
                      R = n_sim )
  
  pdf_boot <- density(Boot$t)
  cdf_boot <- cumsum( pdf_boot$y / sum(pdf_boot$y) )
  
  
  Boot$t0; sd_true^2
  
  
  
  # Plots
  
  plot( x = quantile(z, q_vec), y = p_vec )
  plot( x = x_vec[idx-1], y = d_vec )
  
  
  plot( x = pdf_boot$x, y = cdf_boot )
  points(  x = quantile(z, q_vec), y = p_vec, col = 2 )
  
  plot( pdf_boot )
  points( x = x_vec[idx-1], y = d_vec, col = 2 )
  abline( v = Boot$t0 )
  abline( v = sd_true^2 ) 
  
  
  
  
  # --------------------------------------------------------------------------
  # Replicate example in Bertsimas and Sturt (2020)
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(slolz)
  require(RP)
  
  
  d <- 81
  B <- 10^6
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, rep(1, d - length(z)))
  
  
  tic <- Sys.time()
  Boot <- boot::boot( data = z,
                      statistic = meanBoot,
                      R = B )
  (toc_boot <- Sys.time() - tic)
  
  
  plot( density(Boot$t) )
  quantile( Boot$t, 0.025 )  
  
  
  # Varsi
  q_vec <- seq( from = 1e-12, to = 0.915, length.out = 10^5 )
  FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  tic <- Sys.time()
  p_vec <- unlist( lapply( quantile(z, q_vec), FUN ) )
  (toc_varsi <- Sys.time() - tic)
  
  plot( x = quantile(z, q_vec), y = p_vec )

  idx <- which(p_vec >= 0.025)[1]
  ( quantile(z, q_vec)[idx] + quantile(z, q_vec)[idx-1] ) / 2
  
  
  # Find the quantile H(0.025) directly
  p_th <- 0.025
  err <- 1
  tol <- 1e-12
  z0 <- quantile(z, p_th)
  z_max <- max(z)
  i <- 1
  tic <- Sys.time()
  while( error > tol ) {
    p <- frustum_of_simplex( a = z, z0 = z0 )
    error <- abs(p - p_th)
    if ( p > p_th ) {
      z0 <- z0 - 1/(2^i) * z_max
    } else {
      z0 <- z0 + 1/(2^i) * z_max
    }
    i <- i + 1
  }
  (toc_varsi_direct <- Sys.time() - tic)
  i
  
  z0
  p
  
  
  
  
  