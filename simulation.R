#### Simulation feasibility study --------

#### Packages and settings -----------------------------------------------

library(dplyr)

setwd("E:/MSc_dissertation/data")
set.seed(2026)

number_of_members <- 180
number_of_roll_calls <- 60
number_of_repetitions <- 5

lsm_number_of_members <- 80
lsm_number_of_repetitions <- 3
lsm_maximum_training_dyads <- 2500

results_directory <- "simulation_results"

if (!dir.exists(results_directory)) {
  dir.create(results_directory)
}

old_files <- c(
  "simulation_results.csv",
  "simulation_summary.csv",
  "simulation_lsm_results.csv",
  "simulation_lsm_summary.csv",
  "simulation_cyclic_policy.csv",
  "simulation_multi_bloc.csv"
)

for (old_file in old_files) {
  old_path <- file.path(results_directory, old_file)
  if (file.exists(old_path)) {
    file.remove(old_path)
  }
}

#### MDS simulation -------------------------------------------------------

simulation_results <- data.frame()

for (scenario in c("Multi-bloc", "Cyclic-policy")) {
  
  for (repetition in 1:number_of_repetitions) {
    
    #Generate simulated voting data
    
    n <- number_of_members
    p <- number_of_roll_calls
    
    if (scenario == "Multi-bloc") {
      
      bloc <- rep(1:3, length.out = n)
      
      bloc_angles <- c( 0, 2 * pi / 3, 4 * pi / 3 )
      
      member_angles <- bloc_angles[bloc] + rnorm(n, 0, 0.20)
      member_radius <- rnorm(n, 1, 0.15)
      
      latent_members <- cbind( member_radius * cos(member_angles), member_radius * sin(member_angles) )
      
      item_angles <- runif(p, 0, 2 * pi)
      
      item_directions <- cbind( cos(item_angles), sin(item_angles) )
      
      eta <- outer( latent_members[, 1], item_directions[, 1] ) +
        outer( latent_members[, 2], item_directions[, 2] )
      
    } else {
      
      member_angles <- runif(n, 0, 2 * pi)
      
      latent_members <- cbind( cos(member_angles), sin(member_angles) )
      
      item_angles <- runif(p, 0, 2 * pi)
      
      item_directions <- cbind( cos(item_angles), sin(item_angles) )
      
      eta <- 1.8 * ( outer(latent_members[, 1], item_directions[, 1]) + outer(latent_members[, 2], item_directions[, 2]) )
    }
    
    eta <- eta + matrix( rnorm(n * p, 0, 0.25), nrow = n, ncol = p )
    
    probability <- plogis(eta)
    
    votes <- matrix( rbinom(n * p, 1, as.vector(probability)), nrow = n, ncol = p )
    
    #### 2.2 Construct pairwise disagreement -------------------------------
    
    disagreement <- matrix(0, nrow = n, ncol = n)
    
    for (i in 1:(n - 1)) {
      for (k in (i + 1):n) {
        
        agreement_ik <- mean(votes[i, ] == votes[k, ])
        
        disagreement[i, k] <- 1 - agreement_ik
        disagreement[k, i] <- disagreement[i, k]
      }
    }
    
    #### 2.3 Fit 1D and 2D Euclidean MDS -----------------------------------
    
    classical_2d <- cmdscale( as.dist(disagreement), k = 2 )
    
    classical_1d <- matrix( as.numeric( cmdscale( as.dist(disagreement), k = 1 ) ), ncol = 1 )
    
    #### 2.4 Fit circular MDS -----------------------------------------------
    
    circular_initial <- cmdscale( as.dist(disagreement), k = 2 )
    
    circular_initial_angles <- atan2( circular_initial[, 2], circular_initial[, 1] ) %% (2 * pi)
    
    circular_target <- pi * disagreement / max(disagreement)
    
    # optim() requires an objective function.
    circular_loss <- function(free_angles) {
      
      theta <- c( 0, free_angles %% (2 * pi) )
      
      angle_difference <- abs( outer(theta, theta, "-") )
      
      circular_fitted_distance <- pmin( angle_difference, 2 * pi - angle_difference )
      
      mean( ( circular_fitted_distance[ upper.tri(circular_fitted_distance) ] - circular_target[ upper.tri(circular_target) ] )^2 )
    }
    
    circular_fit <- optim( circular_initial_angles[-1], circular_loss, method = "BFGS", control = list(maxit = 250) )
    
    circular_angles <- c( 0, circular_fit$par %% (2 * pi) )
    
    circular_coordinates <- cbind( cos(circular_angles), sin(circular_angles) )
    
    circular_converged <- circular_fit$convergence == 0
    
    #### 2.5 Fit soft circular MDS ------------------------------------------
    
    soft_target <- pi * disagreement / max(disagreement)
    
    soft_angles <- atan2( circular_coordinates[, 2], circular_coordinates[, 1] )
    
    soft_start_radii <- rep(1, n)
    
    # optim() requires an objective function.
    soft_circular_loss <- function(log_radii) {
      
      radii <- exp(log_radii)
      
      coordinates <- cbind( radii * cos(soft_angles), radii * sin(soft_angles) )
      
      fitted_distance <- as.matrix( dist(coordinates) )
      
      distance_loss <- mean( ( fitted_distance[ upper.tri(fitted_distance) ] - soft_target[ upper.tri(soft_target) ] )^2 )
      
      radius_penalty <- 5 * mean( (radii - mean(radii))^2 )
      
      distance_loss + radius_penalty
    }
    
    soft_circular_fit <- optim( log(soft_start_radii), soft_circular_loss, method = "BFGS", control = list(maxit = 250) )
    
    soft_radii <- exp( soft_circular_fit$par )
    
    soft_circular_coordinates <- cbind( soft_radii * cos(soft_angles), soft_radii * sin(soft_angles) )
    
    soft_circular_converged <-
      soft_circular_fit$convergence == 0
    
    #### 2.6 Fit spherical MDS ----------------------------------------------
    
    spherical_coordinates <- NULL
    spherical_converged <- FALSE
    
    if (requireNamespace("smacof", quietly = TRUE)) {
      
      spherical_fit <- tryCatch(
        smacof::smacofSphere(
          as.dist(disagreement),
          ndim = 3,
          type = "interval",
          algorithm = "dual",
          init = "torgerson",
          itmax = 1000,
          eps = 1e-3,
          verbose = FALSE
        ),
        error = function(error) NULL
      )
      
      if (!is.null(spherical_fit)) {
        
        if (spherical_fit$niter < 1000) {
          
          spherical_raw <- spherical_fit$conf[
            , 1:3,
            drop = FALSE
          ]
          
          spherical_radii <- sqrt( rowSums(spherical_raw^2) )
          
          if (
            all(is.finite(spherical_radii)) &&
            all(spherical_radii > 1e-12)
          ) {
            
            spherical_coordinates <- sweep( spherical_raw, 1, spherical_radii, FUN = "/" )
            
            spherical_converged <- TRUE
          }
        }
      }
    }
    
    #### 2.7 Put the MDS coordinates together -------------------------------
    
    mds_coordinates <- list(
      "1D Euclidean MDS" = classical_1d,
      "2D Euclidean MDS" = classical_2d,
      "Circular MDS" = circular_coordinates,
      "Soft Circular MDS" = soft_circular_coordinates,
      "Spherical MDS" = spherical_coordinates
    )
    
    mds_convergence <- c(
      "1D Euclidean MDS" = TRUE,
      "2D Euclidean MDS" = TRUE,
      "Circular MDS" = circular_converged,
      "Soft Circular MDS" = soft_circular_converged,
      "Spherical MDS" = spherical_converged
    )
    
    #### 2.8 Evaluate each MDS method ---------------------------------------
    
    for (method in names(mds_coordinates)) {
      
      coordinates <- mds_coordinates[[method]]
      
      if (is.null(coordinates)) {
        
        normalized_stress_value <- NA_real_
        held_out_log_loss_value <- NA_real_
        convergence_value <- "no"
        
      } else {
        
        #### Fitted geometry distances
        
        if (method == "Circular MDS") {
          
          fitted_angles <- atan2( coordinates[, 2], coordinates[, 1] ) %% (2 * pi)
          
          fitted_angle_difference <- abs( outer( fitted_angles, fitted_angles, "-" ) )
          
          fitted_distance <- pmin( fitted_angle_difference, 2 * pi - fitted_angle_difference )
          
        } else if (method == "Spherical MDS") {
          
          fitted_inner_product <- tcrossprod( coordinates )
          
          fitted_inner_product <- pmax( pmin(fitted_inner_product, 1), -1 )
          
          fitted_distance <- acos( fitted_inner_product )
          
        } else {
          
          fitted_distance <- as.matrix( dist(coordinates) )
        }
        
        #### Normalized stress
        
        observed_values <- disagreement[
          upper.tri(disagreement)
        ]
        
        fitted_values <- fitted_distance[
          upper.tri(fitted_distance)
        ]
        
        scale_factor <- sum( observed_values * fitted_values ) / sum( fitted_values^2 )
        
        normalized_stress_value <- sqrt( sum( ( observed_values - scale_factor * fitted_values )^2 ) / sum(observed_values^2) )
        
        #### Held-out vote-entry log-loss
        
        decoder_coordinates <- as.matrix( coordinates )[
          ,
          1:min(2, ncol(as.matrix(coordinates))),
          drop = FALSE
        ]
        
        coordinate_data <- as.data.frame( decoder_coordinates )
        
        names(coordinate_data) <- paste0( "x", 1:ncol(coordinate_data) )
        
        roll_call_losses <- numeric(0)
        
        for (roll_call in 1:ncol(votes)) {
          
          test_members <- sample( 1:nrow(votes), size = max( 1, floor(0.20 * nrow(votes)) ) )
          
          train_members <- setdiff( 1:nrow(votes), test_members )
          
          training_data <- coordinate_data[
            train_members,
            ,
            drop = FALSE
          ]
          
          training_data$vote <- votes[
            train_members,
            roll_call
          ]
          
          decoder_fit <- tryCatch( glm( vote ~ ., data = training_data, family = binomial() ), error = function(error) NULL )
          
          if (!is.null(decoder_fit)) {
            
            predicted_probability <- predict( decoder_fit, newdata = coordinate_data[ test_members, , drop = FALSE ], type = "response" )
            
            predicted_probability <- pmin( 1 - 1e-6, pmax( 1e-6, predicted_probability ) )
            
            vote_test <- votes[
              test_members,
              roll_call
            ]
            
            current_loss <- -mean( vote_test * log(predicted_probability) + (1 - vote_test) * log(1 - predicted_probability) )
            
            roll_call_losses <- c( roll_call_losses, current_loss )
          }
        }
        
        if (length(roll_call_losses) == 0) {
          held_out_log_loss_value <- NA_real_
        } else {
          held_out_log_loss_value <- mean( roll_call_losses )
        }
        
        if (isTRUE(mds_convergence[method])) {
          convergence_value <- "yes"
        } else {
          convergence_value <- "no"
        }
      }
      
      one_result <- data.frame(
        scenario = scenario,
        repetition = repetition,
        method = method,
        normalized_stress = normalized_stress_value,
        held_out_log_loss = held_out_log_loss_value,
        converged = convergence_value,
        stringsAsFactors = FALSE
      )
      
      simulation_results <- rbind( simulation_results, one_result )
    }
    
    message( scenario, " repetition ", repetition, " of ", number_of_repetitions, " complete." )
  }
}

#### Summarise MDS simulation results ------------------------------------

simulation_summary <- simulation_results |>
  group_by( scenario, method ) |>
  summarise(
    normalized_stress = mean( normalized_stress, na.rm = TRUE ),
    held_out_log_loss = mean( held_out_log_loss, na.rm = TRUE ),
    .groups = "drop"
  )

#### LSM simulation -------------------------------------------------------

lsm_results <- data.frame()

for (scenario in c("Multi-bloc", "Cyclic-policy")) {
  
  for (repetition in 1:lsm_number_of_repetitions) {
    
    #### 4.1 Generate voting data again -------------------------------------
    
    n <- lsm_number_of_members
    p <- number_of_roll_calls
    
    if (scenario == "Multi-bloc") {
      
      bloc <- rep(1:3, length.out = n)
      
      bloc_angles <- c( 0, 2 * pi / 3, 4 * pi / 3 )
      
      member_angles <- bloc_angles[bloc] + rnorm(n, 0, 0.20)
      member_radius <- rnorm(n, 1, 0.15)
      
      latent_members <- cbind( member_radius * cos(member_angles), member_radius * sin(member_angles) )
      
      item_angles <- runif( p, 0, 2 * pi )
      
      item_directions <- cbind( cos(item_angles), sin(item_angles) )
      
      eta <- outer( latent_members[, 1], item_directions[, 1] ) +
        outer( latent_members[, 2], item_directions[, 2] )
      
    } else {
      
      member_angles <- runif( n, 0, 2 * pi )
      
      latent_members <- cbind( cos(member_angles), sin(member_angles) )
      
      item_angles <- runif( p, 0, 2 * pi )
      
      item_directions <- cbind( cos(item_angles), sin(item_angles) )
      
      eta <- 1.8 * ( outer( latent_members[, 1], item_directions[, 1] ) + outer( latent_members[, 2], item_directions[, 2] ) )
    }
    
    eta <- eta + matrix( rnorm(n * p, 0, 0.25), nrow = n, ncol = p )
    
    probability <- plogis(eta)
    
    votes <- matrix( rbinom( n * p, 1, as.vector(probability) ), nrow = n, ncol = p )
    
    #### 4.2 Construct co-voting adjacency matrix ---------------------------
    
    agreement <- matrix( 0, nrow = n, ncol = n )
    
    for (i in 1:(n - 1)) {
      for (k in (i + 1):n) {
        
        agreement[i, k] <- mean( votes[i, ] == votes[k, ] )
        
        agreement[k, i] <- agreement[i, k]
      }
    }
    
    agreement_threshold <- median( agreement[ upper.tri(agreement) ] )
    
    adjacency <- matrix( 0, nrow = n, ncol = n )
    
    adjacency[
      agreement >= agreement_threshold
    ] <- 1
    
    diag(adjacency) <- 0
    
    #### 4.3 Training and held-out dyads ------------------------------------
    
    all_dyads <- which( upper.tri(adjacency), arr.ind = TRUE )
    
    test_rows <- sample( 1:nrow(all_dyads), floor(0.20 * nrow(all_dyads)) )
    
    test_dyads <- all_dyads[
      test_rows,
      ,
      drop = FALSE
    ]
    
    training_dyads <- all_dyads[
      -test_rows,
      ,
      drop = FALSE
    ]
    
    if (
      nrow(training_dyads) >
      lsm_maximum_training_dyads
    ) {
      
      training_dyads <- training_dyads[
        sample( 1:nrow(training_dyads), lsm_maximum_training_dyads ),
        ,
        drop = FALSE
      ]
    }
    
    #### 4.4 Fit 1D Euclidean LSM ------------------------------------------
    
    use <- training_dyads
    ndim <- 1
    
    initial <- cmdscale( as.dist(1 - adjacency), k = ndim )
    
    initial <- matrix( initial, ncol = ndim )
    
    initial <- scale( initial, center = TRUE, scale = FALSE )
    
    network_density <- mean( adjacency[ upper.tri(adjacency) ] )
    
    euclidean_1d_start <- c( qlogis( pmin( 0.999, pmax(0.001, network_density) ) ), log(1), as.vector(initial) )
    
    euclidean_1d_nll <- function(par) {
      
      alpha <- par[1]
      lambda <- exp(par[2])
      
      coordinates <- matrix( par[-c(1, 2)], nrow = n, ncol = ndim )
      
      distances <- sqrt( rowSums( ( coordinates[ use[, 1], , drop = FALSE ] - coordinates[ use[, 2], , drop = FALSE ] )^2 ) )
      
      eta_fit <- alpha - lambda * distances
      
      y <- adjacency[
        cbind( use[, 1], use[, 2] )
      ]
      
      fitted_probability <- plogis( pmin( 30, pmax( -30, eta_fit ) ) )
      
      fitted_probability <- pmin( 1 - 1e-12, pmax( 1e-12, fitted_probability ) )
      
      -sum( y * log(fitted_probability) + (1 - y) * log1p(-fitted_probability) )
    }
    
    euclidean_1d_fit <- tryCatch(
      optim(
        euclidean_1d_start,
        euclidean_1d_nll,
        method = "L-BFGS-B",
        lower = c( -8, -5, rep(-5, n * ndim) ),
        upper = c( 8, 5, rep(5, n * ndim) ),
        control = list( maxit = 600, factr = 1e8, pgtol = 1e-4 )
      ),
      error = function(error) NULL
    )
    
    euclidean_1d_coordinates <- NULL
    euclidean_1d_alpha <- NA_real_
    euclidean_1d_lambda <- NA_real_
    euclidean_1d_logLik <- NA_real_
    euclidean_1d_converged <- FALSE
    
    if (!is.null(euclidean_1d_fit)) {
      
      if (all(is.finite(euclidean_1d_fit$par))) {
        
        euclidean_1d_coordinates <- matrix( euclidean_1d_fit$par[ -c(1, 2) ], nrow = n, ncol = ndim )
        
        euclidean_1d_alpha <-
          euclidean_1d_fit$par[1]
        
        euclidean_1d_lambda <-
          exp(euclidean_1d_fit$par[2])
        
        euclidean_1d_logLik <-
          -euclidean_1d_fit$value
        
        euclidean_1d_converged <-
          euclidean_1d_fit$convergence == 0
      }
    }
    
    #### 4.5 Fit 2D Euclidean LSM ------------------------------------------
    
    use <- training_dyads
    ndim <- 2
    
    initial <- cmdscale( as.dist(1 - adjacency), k = ndim )
    initial <- scale( initial, center = TRUE, scale = FALSE )
    
    euclidean_2d_start <- c( qlogis( pmin( 0.999, pmax(0.001, network_density) ) ), log(1), as.vector(initial) )
    euclidean_2d_nll <- function(par) {
      
      alpha <- par[1]
      lambda <- exp(par[2])
      
      coordinates <- matrix( par[-c(1, 2)], nrow = n, ncol = ndim )
      distances <- sqrt( rowSums( ( coordinates[ use[, 1], , drop = FALSE ] - coordinates[ use[, 2], , drop = FALSE ] )^2 ) )
      eta_fit <- alpha - lambda * distances
      
      y <- adjacency[
        cbind( use[, 1], use[, 2] )
      ]
      
      fitted_probability <- plogis( pmin( 30, pmax( -30, eta_fit ) ) )
      fitted_probability <- pmin( 1 - 1e-12, pmax( 1e-12, fitted_probability ) )
      
      -sum( y * log(fitted_probability) + (1 - y) * log1p(-fitted_probability) )
    }
    
    euclidean_2d_fit <- tryCatch(
      optim(
        euclidean_2d_start,
        euclidean_2d_nll,
        method = "L-BFGS-B",
        lower = c( -8, -5, rep(-5, n * ndim) ),
        upper = c( 8, 5, rep(5, n * ndim) ),
        control = list( maxit = 600, factr = 1e8, pgtol = 1e-4 )
      ),
      error = function(error) NULL
    )
    
    euclidean_2d_coordinates <- NULL
    euclidean_2d_alpha <- NA_real_
    euclidean_2d_lambda <- NA_real_
    euclidean_2d_logLik <- NA_real_
    euclidean_2d_converged <- FALSE
    
    if (!is.null(euclidean_2d_fit)) {
      
      if (all(is.finite(euclidean_2d_fit$par))) {
        
        euclidean_2d_coordinates <- matrix( euclidean_2d_fit$par[ -c(1, 2) ], nrow = n, ncol = ndim )
        euclidean_2d_alpha <- euclidean_2d_fit$par[1]
        euclidean_2d_lambda <- exp(euclidean_2d_fit$par[2])
        euclidean_2d_logLik <- -euclidean_2d_fit$value
        euclidean_2d_converged <- euclidean_2d_fit$convergence == 0
      }
    }
    
    #### 4.6 Fit circular LSM -----------------------------------------------
    
    use <- training_dyads
    
    circular_initial <- cmdscale( as.dist(1 - adjacency), k = 2 )
    circular_start_angles <- atan2( circular_initial[, 2], circular_initial[, 1] )
    circular_lsm_start <- c( qlogis( pmin( 0.999, pmax(0.001, network_density) ) ), 
                             log(1), circular_start_angles[-1] )
    circular_lsm_nll <- function(par) {
      
      theta <- c( 0, par[-c(1, 2)] )
      alpha <- par[1]
      lambda <- exp(par[2])
      
      eta_fit <- alpha + lambda * cos( theta[use[, 1]] - theta[use[, 2]] )
      
      y <- adjacency[ cbind( use[, 1], use[, 2] )]
      
      fitted_probability <- plogis( pmin( 30, pmax( -30, eta_fit ) ) )
      fitted_probability <- pmin( 1 - 1e-12, pmax( 1e-12, fitted_probability ) )
      
      -sum( y * log(fitted_probability) + (1 - y) * log1p(-fitted_probability) )
    }
    
    circular_lsm_fit <- tryCatch(
      optim(
        circular_lsm_start,
        circular_lsm_nll,
        method = "L-BFGS-B",
        lower = c( -8, -5, rep(-pi, n - 1) ),
        upper = c( 8, 5, rep(pi, n - 1) ),
        control = list( maxit = 600, factr = 1e8, pgtol = 1e-4 )
      ),
      error = function(error) NULL
    )
    
    circular_lsm_coordinates <- NULL
    circular_lsm_alpha <- NA_real_
    circular_lsm_lambda <- NA_real_
    circular_lsm_logLik <- NA_real_
    circular_lsm_converged <- FALSE
    
    if (!is.null(circular_lsm_fit)) {
      
      if (all(is.finite(circular_lsm_fit$par))) {
        
        circular_lsm_theta <- c( 0, circular_lsm_fit$par[ -c(1, 2) ] )
        circular_lsm_coordinates <- cbind( cos(circular_lsm_theta), sin(circular_lsm_theta) )
        circular_lsm_alpha <- circular_lsm_fit$par[1]
        circular_lsm_lambda <- exp(circular_lsm_fit$par[2])
        circular_lsm_logLik <- -circular_lsm_fit$value
        circular_lsm_converged <- circular_lsm_fit$convergence == 0
      }
    }
    
    #### 4.7 Fit soft circular LSM ------------------------------------------
    
    use <- training_dyads
    
    soft_lsm_initial <- cmdscale( as.dist(1 - adjacency), k = 2 )
    soft_lsm_initial <- scale( soft_lsm_initial, center = TRUE, scale = FALSE )
    soft_lsm_start <- c( qlogis( pmin( 0.999, pmax(0.001, network_density) ) ), log(1), as.vector(soft_lsm_initial) )
    
    soft_lsm_nll <- function(par) {
      
      alpha <- par[1]
      lambda <- exp(par[2])
      
      raw_coordinates <- matrix( par[-c(1, 2)], nrow = n, ncol = 2 )
      raw_norm <- sqrt( sum(raw_coordinates^2) )
      
      if (
        !is.finite(raw_norm) ||
        raw_norm <= 1e-10
      ) {
        return(1e12)
      }
      
      coordinates <- sqrt(n) * raw_coordinates / raw_norm
      
      radii <- sqrt( rowSums(coordinates^2) )
      distances <- sqrt( rowSums( ( coordinates[ use[, 1], , drop = FALSE ] - coordinates[ use[, 2], , drop = FALSE ] )^2 ) )
      
      y <- adjacency[cbind( use[, 1], use[, 2] )]
      
      eta_fit <- alpha - lambda * distances
      
      fitted_probability <- plogis( pmin( 30, pmax( -30, eta_fit ) ) )
      fitted_probability <- pmin( 1 - 1e-12, pmax( 1e-12, fitted_probability ) )
      
      negative_log_likelihood <- -sum( y * log(fitted_probability) + (1 - y) * log1p(-fitted_probability) )
      radius_penalty <- 5 * mean( ( radii - mean(radii) )^2 )
      negative_log_likelihood + radius_penalty
    }
    
    soft_lsm_fit <- tryCatch(
      optim(
        soft_lsm_start,
        soft_lsm_nll,
        method = "L-BFGS-B",
        lower = c( -8, -5, rep(-5, 2 * n) ),
        upper = c( 8, 5, rep(5, 2 * n) ),
        control = list( maxit = 600, factr = 1e8, pgtol = 1e-4 )
      ),
      error = function(error) NULL
    )
    
    soft_lsm_coordinates <- NULL
    soft_lsm_alpha <- NA_real_
    soft_lsm_lambda <- NA_real_
    soft_lsm_logLik <- NA_real_
    soft_lsm_converged <- FALSE
    
    if (!is.null(soft_lsm_fit)) {
      
      if (all(is.finite(soft_lsm_fit$par))) {
        
        raw_coordinates <- matrix( soft_lsm_fit$par[ -c(1, 2) ], nrow = n, ncol = 2 )
        soft_lsm_coordinates <- sqrt(n) * raw_coordinates /sqrt( sum(raw_coordinates^2) )
        soft_lsm_distances <- sqrt(
          rowSums( ( soft_lsm_coordinates[ use[, 1], , drop = FALSE ] - soft_lsm_coordinates[ use[, 2], , drop = FALSE ] )^2 )
        )
        
        soft_lsm_y <- adjacency[cbind( use[, 1], use[, 2] )]
        soft_lsm_eta <- soft_lsm_fit$par[1] -
          exp(soft_lsm_fit$par[2]) *
          soft_lsm_distances
        
        soft_lsm_probability <- plogis( pmin( 30, pmax( -30, soft_lsm_eta ) ) )
        soft_lsm_probability <- pmin( 1 - 1e-12, pmax( 1e-12, soft_lsm_probability ) )
        soft_lsm_logLik <- sum( soft_lsm_y * log(soft_lsm_probability) + (1 - soft_lsm_y) * log1p(-soft_lsm_probability) )
        soft_lsm_alpha <- soft_lsm_fit$par[1]
        soft_lsm_lambda <- exp(soft_lsm_fit$par[2])
        soft_lsm_converged <- soft_lsm_fit$convergence == 0
      }
    }
    
    #### 4.8 Fit spherical LSM ----------------------------------------------
    
    use <- training_dyads
    
    spherical_lsm_initial <- cmdscale( as.dist(1 - adjacency), k = 3 )
    spherical_lsm_initial <- spherical_lsm_initial / pmax( sqrt( rowSums( spherical_lsm_initial^2 ) ), 1e-8 )
    
    spherical_start_angles <- cbind(
      atan2( spherical_lsm_initial[, 2], spherical_lsm_initial[, 1] ),
      acos( pmin( 1, pmax( -1, spherical_lsm_initial[, 3] ) ) )
    )
    
    spherical_lsm_start <- c(
      qlogis( pmin( 0.999, pmax(0.001, network_density) ) ),
      log(1),
      as.vector( t( spherical_start_angles[ -1, , drop = FALSE ] ) )
    )
    
    spherical_lsm_nll <- function(par) {
      
      alpha <- par[1]
      lambda <- exp(par[2])
      
      angles <- rbind( spherical_start_angles[1, ], matrix( par[-c(1, 2)], ncol = 2, byrow = TRUE ) )
      coordinates <- cbind( sin(angles[, 2]) * cos(angles[, 1]), sin(angles[, 2]) * sin(angles[, 1]), cos(angles[, 2]) )
      similarity <- rowSums( coordinates[ use[, 1], , drop = FALSE ] * coordinates[ use[, 2], , drop = FALSE ] )
      
      eta_fit <- alpha + lambda * similarity
      
      y <- adjacency[ cbind( use[, 1], use[, 2] )]
      
      fitted_probability <- plogis( pmin( 30, pmax( -30, eta_fit ) ) )
      fitted_probability <- pmin( 1 - 1e-12, pmax( 1e-12, fitted_probability ) )
      
      -sum( y * log(fitted_probability) + (1 - y) * log1p(-fitted_probability) )
    }
    
    spherical_lower <- c( -8, -5, rep( c(-pi, 0.001), n - 1 ) )
    spherical_upper <- c( 8, 5, rep( c(pi, pi - 0.001), n - 1 ) )
    spherical_lsm_fit <- tryCatch(
      optim(
        spherical_lsm_start,
        spherical_lsm_nll,
        method = "L-BFGS-B",
        lower = spherical_lower,
        upper = spherical_upper,
        control = list( maxit = 600, factr = 1e8, pgtol = 1e-4 )
      ),
      error = function(error) NULL
    )
    
    spherical_lsm_coordinates <- NULL
    spherical_lsm_alpha <- NA_real_
    spherical_lsm_lambda <- NA_real_
    spherical_lsm_logLik <- NA_real_
    spherical_lsm_converged <- FALSE
    
    if (!is.null(spherical_lsm_fit)) {
      
      if (all(is.finite(spherical_lsm_fit$par))) {
        
        spherical_final_angles <- rbind(
          spherical_start_angles[1, ],
          matrix( spherical_lsm_fit$par[ -c(1, 2) ], ncol = 2, byrow = TRUE )
        )
        
        spherical_lsm_coordinates <- cbind(
          sin(spherical_final_angles[, 2]) *
            cos(spherical_final_angles[, 1]),
          sin(spherical_final_angles[, 2]) *
            sin(spherical_final_angles[, 1]),
          cos(spherical_final_angles[, 2])
        )
        
        spherical_lsm_alpha <- spherical_lsm_fit$par[1]
        
        spherical_lsm_lambda <- exp(spherical_lsm_fit$par[2])
        
        spherical_lsm_logLik <- -spherical_lsm_fit$value
        
        spherical_lsm_converged <- spherical_lsm_fit$convergence == 0
      }
    }
    
    #### 4.9 Put LSM fits into one list -------------------------------------
    
    lsm_fit_information <- list(
      "1D Euclidean LSM" = list(
        coordinates = euclidean_1d_coordinates,
        alpha = euclidean_1d_alpha,
        lambda = euclidean_1d_lambda,
        logLik = euclidean_1d_logLik,
        converged = euclidean_1d_converged,
        geometry = "Euclidean"
      ),
      "2D Euclidean LSM" = list(
        coordinates = euclidean_2d_coordinates,
        alpha = euclidean_2d_alpha,
        lambda = euclidean_2d_lambda,
        logLik = euclidean_2d_logLik,
        converged = euclidean_2d_converged,
        geometry = "Euclidean"
      ),
      "Circular LSM" = list(
        coordinates = circular_lsm_coordinates,
        alpha = circular_lsm_alpha,
        lambda = circular_lsm_lambda,
        logLik = circular_lsm_logLik,
        converged = circular_lsm_converged,
        geometry = "Circular"
      ),
      "Soft Circular LSM" = list(
        coordinates = soft_lsm_coordinates,
        alpha = soft_lsm_alpha,
        lambda = soft_lsm_lambda,
        logLik = soft_lsm_logLik,
        converged = soft_lsm_converged,
        geometry = "Soft circular"
      ),
      "Spherical LSM" = list(
        coordinates = spherical_lsm_coordinates,
        alpha = spherical_lsm_alpha,
        lambda = spherical_lsm_lambda,
        logLik = spherical_lsm_logLik,
        converged = spherical_lsm_converged,
        geometry = "Spherical"
      )
    )
    
    #### 4.10 Evaluate each LSM ---------------------------------------------
    
    for (method in names(lsm_fit_information)) {
      
      fit_info <- lsm_fit_information[[method]]
      
      coordinates <- fit_info$coordinates
      
      model_converged <-
        isTRUE(fit_info$converged) &&
        !is.null(coordinates)
      
      if (method == "1D Euclidean LSM") {
        parameter_count <- lsm_number_of_members + 2
      } else if (method == "2D Euclidean LSM") {
        parameter_count <- 2 * lsm_number_of_members + 2
      } else if (method == "Circular LSM") {
        parameter_count <- lsm_number_of_members + 1
      } else if (method == "Soft Circular LSM") {
        parameter_count <- 2 * lsm_number_of_members + 2
      } else {
        parameter_count <- 2 * lsm_number_of_members
      }
      
      #### BIC
      
      if (
        model_converged &&
        method != "Soft Circular LSM"
      ) {
        
        bic_value <- -2 *
          fit_info$logLik +
          parameter_count *
          log( nrow(training_dyads) )
        
      } else {
        
        bic_value <- NA_real_
      }
      
      #### Held-out tie log-loss
      
      if (model_converged) {
        
        if (
          fit_info$geometry == "Circular"
        ) {
          
          theta_eval <- atan2( coordinates[, 2], coordinates[, 1] )
          similarity_matrix <- cos( outer( theta_eval, theta_eval, "-" ) )
          
          eta_test <- fit_info$alpha +
            fit_info$lambda *
            similarity_matrix[
              cbind( test_dyads[, 1], test_dyads[, 2] )
            ]
          
        } else if ( fit_info$geometry == "Spherical" ) {
          
          similarity_matrix <- tcrossprod( coordinates )
          
          eta_test <- fit_info$alpha +
            fit_info$lambda *
            similarity_matrix[
              cbind( test_dyads[, 1], test_dyads[, 2] )
            ]
          
        } else {
          
          euclidean_distance_matrix <- as.matrix( dist(coordinates) )
          
          eta_test <- fit_info$alpha -
            fit_info$lambda *
            euclidean_distance_matrix[
              cbind( test_dyads[, 1], test_dyads[, 2] )
            ]
        }
        
        y_test <- adjacency[
          cbind( test_dyads[, 1], test_dyads[, 2] )
        ]
        
        probability_test <- plogis( eta_test )
        probability_test <- pmin( 1 - 1e-12, pmax( 1e-12, probability_test ) )
        held_out_tie_log_loss_value <- -mean( y_test * log(probability_test) + (1 - y_test) * log1p(-probability_test) )
        
      } else {
        
        held_out_tie_log_loss_value <- NA_real_
      }
      
      #### Normalized stress for the LSM coordinates
      
      if (model_converged) {
        
        observed_lsm_distance <- 1 - adjacency
        
        if (method == "Circular LSM") {
          
          theta_stress <- atan2( coordinates[, 2], coordinates[, 1] ) %% (2 * pi)
          angle_difference <- abs( outer( theta_stress, theta_stress, "-" ) )
          fitted_lsm_distance <- pmin( angle_difference, 2 * pi - angle_difference )
          
        } else if (method == "Spherical LSM") {
          
          inner_product <- tcrossprod( coordinates )
          inner_product <- pmax( pmin(inner_product, 1), -1 )
          fitted_lsm_distance <- acos( inner_product )
          
        } else {
          
          fitted_lsm_distance <- as.matrix( dist(coordinates) )
        }
        
        observed_values <- observed_lsm_distance[upper.tri(observed_lsm_distance)]
        fitted_values <- fitted_lsm_distance[upper.tri(fitted_lsm_distance)]
        stress_scale <- sum( observed_values * fitted_values ) / sum( fitted_values^2 )
        
        normalized_stress_value <- sqrt( sum( ( observed_values - stress_scale * fitted_values )^2 ) / sum( observed_values^2 ) )
        
      } else {
        
        normalized_stress_value <- NA_real_
      }
      
      one_lsm_result <- data.frame(
        scenario = scenario,
        repetition = repetition,
        method = method,
        geometry = fit_info$geometry,
        members = lsm_number_of_members,
        training_dyads = nrow(training_dyads),
        held_out_dyads = nrow(test_dyads),
        bic = bic_value,
        held_out_tie_log_loss =
          held_out_tie_log_loss_value,
        normalized_stress =
          normalized_stress_value,
        converged = ifelse( model_converged, "yes", "no" ),
        stringsAsFactors = FALSE
      )
      
      lsm_results <- rbind( lsm_results, one_lsm_result )
    }
    
    message( scenario, " LSM repetition ", repetition, " of ", lsm_number_of_repetitions, " complete." )
  }
}

#### Summarise LSM results -----------------------------------------------
lsm_summary <- lsm_results |>
  group_by( scenario, method, geometry ) |>
  summarise(
    repetitions = n(),
    successful_repetitions =
      sum(converged == "yes"),
    bic = if ( any(is.finite(bic)) ) {
      mean( bic[is.finite(bic)] )
    } else {
      NA_real_
    },
    held_out_tie_log_loss = if ( any(is.finite(held_out_tie_log_loss)) ) {
      mean( held_out_tie_log_loss[ is.finite(held_out_tie_log_loss) ] )
    } else {
      NA_real_
    },
    normalized_stress = if ( any(is.finite(normalized_stress)) ) {
      mean( normalized_stress[ is.finite(normalized_stress) ] )
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Prepare the paper tables
simulation_mds_table <- simulation_results |>
  mutate(
    Family = "MDS",
    Dimension = ifelse( method %in% c( "1D Euclidean MDS", "Circular MDS" ), "1D", "2D" ),
    Method = method,
    `Held-out log-loss` =
      held_out_log_loss,
    `Normalized stress` =
      normalized_stress,
    BIC = NA_real_
  ) |>
  group_by( scenario, Family, Method, Dimension ) |>
  summarise(
    `Held-out log-loss` =
      mean( `Held-out log-loss`, na.rm = TRUE ),
    `Normalized stress` =
      mean( `Normalized stress`, na.rm = TRUE ),
    BIC = NA_real_,
    .groups = "drop"
  )

simulation_lsm_table <- lsm_results |>
  mutate(
    Family = "LSM",
    Dimension = ifelse( method %in% c( "1D Euclidean LSM", "Circular LSM" ), "1D", "2D" ),
    Method = method,
    `Held-out log-loss` =
      held_out_tie_log_loss,
    `Normalized stress` =
      normalized_stress,
    BIC = bic
  ) |>
  group_by( scenario, Family, Method, Dimension ) |>
  summarise(
    `Held-out log-loss` = if ( any( is.finite( `Held-out log-loss` ) ) ) {
      mean( `Held-out log-loss`[ is.finite( `Held-out log-loss` ) ] )
    } else {
      NA_real_
    },
    `Normalized stress` = if ( any( is.finite( `Normalized stress` ) ) ) {
      mean( `Normalized stress`[ is.finite( `Normalized stress` ) ] )
    } else {
      NA_real_
    },
    BIC = if ( any(is.finite(BIC)) ) {
      mean( BIC[is.finite(BIC)] )
    } else {
      NA_real_
    },
    .groups = "drop"
  )

scenario_comparison_table <- bind_rows( simulation_mds_table, simulation_lsm_table ) |>
  select( scenario, Family, Method, Dimension, `Held-out log-loss`, `Normalized stress`, BIC ) |>
  mutate(
    Method = factor(
      Method,
      levels = c(
        "1D Euclidean MDS",
        "2D Euclidean MDS",
        "Circular MDS",
        "Soft Circular MDS",
        "Spherical MDS",
        "1D Euclidean LSM",
        "2D Euclidean LSM",
        "Circular LSM",
        "Soft Circular LSM",
        "Spherical LSM"
      )
    )
  ) |>
  arrange( scenario, match( Family, c("MDS", "LSM") ), Method ) |>
  mutate( Method = as.character(Method) )

# Keep reported results to four decimal places
simulation_results <- simulation_results |>
  mutate(across(c(normalized_stress, held_out_log_loss), round, digits = 4))
simulation_summary <- simulation_summary |>
  mutate(across(c(normalized_stress, held_out_log_loss), round, digits = 4))
lsm_results <- lsm_results |>
  mutate(across(c(bic, held_out_tie_log_loss, normalized_stress), round, digits = 4))
lsm_summary <- lsm_summary |>
  mutate(across(c(bic, held_out_tie_log_loss, normalized_stress), round, digits = 4))
scenario_comparison_table <- scenario_comparison_table |>
  mutate(across(c(`Held-out log-loss`, `Normalized stress`, BIC), round, digits = 4))

# Save outputs
write.csv( simulation_results, file.path( results_directory, "simulation_results.csv" ), row.names = FALSE )
write.csv( simulation_summary, file.path( results_directory, "simulation_summary.csv" ), row.names = FALSE )
write.csv( lsm_results, file.path( results_directory, "simulation_lsm_results.csv" ), row.names = FALSE )
write.csv( lsm_summary, file.path( results_directory, "simulation_lsm_summary.csv" ), row.names = FALSE )

cyclic_policy_table <- scenario_comparison_table |>
  filter( scenario == "Cyclic-policy" ) |>
  select(-scenario)

write.csv(
  cyclic_policy_table,
  file.path( results_directory, "simulation_cyclic_policy.csv" ),
  row.names = FALSE,
  na = "NA"
)

multi_bloc_table <- scenario_comparison_table |>
  filter( scenario == "Multi-bloc" ) |>
  select(-scenario)

write.csv( multi_bloc_table, file.path( results_directory, "simulation_multi_bloc.csv" ), row.names = FALSE, na = "NA" )

#### Print summaries ------------------------------------------------------

print(simulation_summary)
print(lsm_summary)