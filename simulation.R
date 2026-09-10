#### Simulation feasibility study for latent geometry -----------------------
# This script evaluates whether circular or spherical geometry can be useful
# when the data contain multi-bloc or genuinely periodic latent variation.

# Packages and settings ---------------------------------------------------
setwd("E:/MSc_dissertation/data")
set.seed(2026)

number_of_members <- 180
number_of_roll_calls <- 60
number_of_repetitions <- 5
lsm_number_of_members <- 80
lsm_number_of_repetitions <- 3
lsm_maximum_training_dyads <- 2500

results_directory <- "simulation_results"
if (!dir.exists(results_directory)) dir.create(results_directory)
for (old_file in c("simulation_results.csv", "simulation_summary.csv", "simulation_lsm_results.csv", "simulation_lsm_summary.csv")) {
  old_path <- file.path(results_directory, old_file)
  if (file.exists(old_path)) file.remove(old_path)
}

# Generate binary voting data
generate_votes <- function(scenario, n, p) {
  if (scenario == "Multi-bloc") {
    bloc <- rep(1:3, length.out = n)
    bloc_angles <- c(0, 2 * pi / 3, 4 * pi / 3)
    member_angles <- bloc_angles[bloc] + rnorm(n, 0, 0.20)
    member_radius <- rnorm(n, 1, 0.15)
    latent_members <- cbind(member_radius * cos(member_angles), member_radius * sin(member_angles))
    item_angles <- runif(p, 0, 2 * pi)
    item_directions <- cbind(cos(item_angles), sin(item_angles))
    eta <- outer(latent_members[, 1], item_directions[, 1]) + outer(latent_members[, 2], item_directions[, 2])
  } else {
    member_angles <- runif(n, 0, 2 * pi)
    latent_members <- cbind(cos(member_angles), sin(member_angles))
    item_angles <- runif(p, 0, 2 * pi)
    item_directions <- cbind(cos(item_angles), sin(item_angles))
    eta <- 1.8 * (outer(latent_members[, 1], item_directions[, 1]) + outer(latent_members[, 2], item_directions[, 2]))
  }
  eta <- eta + matrix(rnorm(n * p, 0, 0.25), nrow = n, ncol = p)
  probability <- plogis(eta)
  votes <- matrix(rbinom(n * p, 1, as.vector(probability)), nrow = n, ncol = p)
  list(votes = votes, latent_members = latent_members)
}

# Construct pairwise disagreement
make_disagreement <- function(votes) {
  n <- nrow(votes)
  disagreement <- matrix(0, n, n)
  for (i in seq_len(n - 1)) for (k in (i + 1):n) {
    agreement <- mean(votes[i, ] == votes[k, ])
    disagreement[i, k] <- 1 - agreement
    disagreement[k, i] <- disagreement[i, k]
  }
  disagreement
}

# Geometry fits
normalized_stress <- function(observed, fitted) {
  observed_values <- observed[upper.tri(observed)]
  fitted_values <- fitted[upper.tri(fitted)]
  scale_factor <- sum(observed_values * fitted_values) / sum(fitted_values ^ 2)
  sqrt(sum((observed_values - scale_factor * fitted_values)^2) / sum(observed_values^2))
}

fit_circular_mds <- function(disagreement) {
  n <- nrow(disagreement)
  initial <- cmdscale(as.dist(disagreement), k = 2)
  initial_angles <- atan2(initial[, 2], initial[, 1]) %% (2 * pi)
  target <- pi * disagreement / max(disagreement)
  loss <- function(free_angles) {
    theta <- c(0, free_angles %% (2 * pi))
    difference <- abs(outer(theta, theta, "-"))
    fitted <- pmin(difference, 2 * pi - difference)
    mean((fitted[upper.tri(fitted)] - target[upper.tri(target)])^2)
  }
  fit <- optim(initial_angles[-1], loss, method = "BFGS", control = list(maxit = 250))
  angles <- c(0, fit$par %% (2 * pi))
  list(coordinates = cbind(cos(angles), sin(angles)), converged = fit$convergence == 0)
}

fit_soft_circular_mds <- function(disagreement, circular_coordinates) {
  target <- pi * disagreement / max(disagreement)
  angles <- atan2(circular_coordinates[, 2], circular_coordinates[, 1])
  start <- rep(1, nrow(disagreement))
  loss <- function(log_radii) {
    radii <- exp(log_radii)
    coordinates <- cbind(radii * cos(angles), radii * sin(angles))
    fitted <- as.matrix(dist(coordinates))
    mean((fitted[upper.tri(fitted)] - target[upper.tri(target)])^2) + 5 * mean((radii - mean(radii))^2)
  }
  fit <- optim(log(start), loss, method = "BFGS", control = list(maxit = 250))
  radii <- exp(fit$par)
  list(coordinates = cbind(radii * cos(angles), radii * sin(angles)), converged = fit$convergence == 0)
}

fit_spherical_mds <- function(disagreement) {
  if (!requireNamespace("smacof", quietly = TRUE)) return(list(coordinates = NULL, converged = FALSE))
  fit <- tryCatch(smacof::smacofSphere(as.dist(disagreement), ndim = 3, type = "interval", algorithm = "dual", init = "torgerson", itmax = 1000, eps = 1e-3, verbose = FALSE), error = function(error) NULL)
  if (is.null(fit) || fit$niter >= 1000) return(list(coordinates = NULL, converged = FALSE))
  raw_coordinates <- fit$conf[, 1:3, drop = FALSE]
  radii <- sqrt(rowSums(raw_coordinates ^ 2))
  if (any(!is.finite(radii)) || any(radii <= 1e-12)) return(list(coordinates = NULL, converged = FALSE))
  list(coordinates = sweep(raw_coordinates, 1, radii, FUN = "/"), converged = TRUE)
}

fitted_geometry_distance <- function(coordinates, method) {
  if (method == "Circular MDS") {
    angles <- atan2(coordinates[, 2], coordinates[, 1]) %% (2 * pi)
    differences <- abs(outer(angles, angles, "-"))
    return(pmin(differences, 2 * pi - differences))
  }
  if (method == "Spherical MDS") return(acos(pmax(pmin(tcrossprod(coordinates), 1), -1)))
  as.matrix(dist(coordinates))
}

#### Simplified latent-space models ---------------------------------------
bernoulli_negative_log_likelihood <- function(y, eta) {
  probability <- pmin(1 - 1e-12, pmax(1e-12, plogis(pmin(30, pmax(-30, eta)))))
  -sum(y * log(probability) + (1 - y) * log1p(-probability))
}

make_covoting_adjacency <- function(votes) {
  agreement <- matrix(0, nrow(votes), nrow(votes))
  for (i in seq_len(nrow(votes) - 1)) for (k in (i + 1):nrow(votes)) {
    agreement[i, k] <- agreement[k, i] <- mean(votes[i, ] == votes[k, ])
  }
  threshold <- stats::median(agreement[upper.tri(agreement)])
  adjacency <- matrix(0, nrow(votes), nrow(votes)); adjacency[agreement >= threshold] <- 1
  diag(adjacency) <- 0
  adjacency
}

fit_euclidean_lsm <- function(adjacency, ndim = 2, dyads = NULL, max_dyads = 6000) {
  n <- nrow(adjacency); use <- if (is.null(dyads)) which(upper.tri(adjacency), arr.ind = TRUE) else dyads
  if (nrow(use) > max_dyads) use <- use[sample(seq_len(nrow(use)), max_dyads), , drop = FALSE]
  initial <- cmdscale(as.dist(1 - adjacency), k = ndim)
  initial <- scale(initial, center = TRUE, scale = FALSE)
  start <- c(qlogis(pmin(0.999, pmax(0.001, mean(adjacency[upper.tri(adjacency)])))), log(1), as.vector(initial))
  nll <- function(par) {
    alpha <- par[1]; lambda <- exp(par[2]); coordinates <- matrix(par[-c(1, 2)], nrow = n, ncol = ndim)
    distances <- sqrt(rowSums((coordinates[use[, 1], , drop = FALSE] - coordinates[use[, 2], , drop = FALSE])^2))
    eta <- alpha - lambda * distances
    bernoulli_negative_log_likelihood(adjacency[cbind(use[, 1], use[, 2])], eta)
  }
  fit <- tryCatch(optim(start, nll, method = "L-BFGS-B", lower = c(-8, -5, rep(-5, n * ndim)), upper = c(8, 5, rep(5, n * ndim)), control = list(maxit = 600, factr = 1e8, pgtol = 1e-4)), error = function(error) NULL)
  if (is.null(fit) || !all(is.finite(fit$par))) return(list(coordinates = NULL, alpha = NA_real_, logLik = NA_real_, converged = FALSE))
  coordinates <- matrix(fit$par[-c(1, 2)], nrow = n, ncol = ndim)
  list(coordinates = coordinates, alpha = fit$par[1], lambda = exp(fit$par[2]), logLik = -fit$value, converged = fit$convergence == 0)
}

fit_circular_lsm <- function(adjacency, dyads = NULL, max_dyads = 6000) {
  n <- nrow(adjacency); use <- if (is.null(dyads)) which(upper.tri(adjacency), arr.ind = TRUE) else dyads
  if (nrow(use) > max_dyads) use <- use[sample(seq_len(nrow(use)), max_dyads), , drop = FALSE]
  initial <- cmdscale(as.dist(1 - adjacency), k = 2)
  start_angles <- atan2(initial[, 2], initial[, 1]); start <- c(qlogis(pmin(0.999, pmax(0.001, mean(adjacency[upper.tri(adjacency)])))), log(1), start_angles[-1])
  nll <- function(par) {
    theta <- c(0, par[-c(1, 2)]); alpha <- par[1]; lambda <- exp(par[2]); eta <- alpha + lambda * cos(theta[use[, 1]] - theta[use[, 2]])
    y <- adjacency[cbind(use[, 1], use[, 2])]; bernoulli_negative_log_likelihood(y, eta)
  }
  fit <- tryCatch(optim(start, nll, method = "L-BFGS-B", lower = c(-8, -5, rep(-pi, n - 1)), upper = c(8, 5, rep(pi, n - 1)), control = list(maxit = 600, factr = 1e8, pgtol = 1e-4)), error = function(error) NULL)
  if (is.null(fit) || !all(is.finite(fit$par))) return(list(coordinates = NULL, alpha = NA_real_, logLik = NA_real_, converged = FALSE))
  theta <- c(0, fit$par[-c(1, 2)]); list(coordinates = cbind(cos(theta), sin(theta)), alpha = fit$par[1], lambda = exp(fit$par[2]), logLik = -fit$value, converged = fit$convergence == 0)
}

fit_soft_circular_lsm <- function(adjacency, dyads = NULL, max_dyads = 6000, penalty = 5) {
  n <- nrow(adjacency); use <- if (is.null(dyads)) which(upper.tri(adjacency), arr.ind = TRUE) else dyads
  if (nrow(use) > max_dyads) use <- use[sample(seq_len(nrow(use)), max_dyads), , drop = FALSE]
  initial <- cmdscale(as.dist(1 - adjacency), k = 2); initial <- scale(initial, center = TRUE, scale = FALSE)
  start <- c(qlogis(pmin(0.999, pmax(0.001, mean(adjacency[upper.tri(adjacency)])))), log(1), as.vector(initial))
  nll <- function(par) {
    alpha <- par[1]; lambda <- exp(par[2]); raw_coordinates <- matrix(par[-c(1, 2)], nrow = n, ncol = 2); raw_norm <- sqrt(sum(raw_coordinates ^ 2)); if (!is.finite(raw_norm) || raw_norm <= 1e-10) return(1e12); coordinates <- sqrt(n) * raw_coordinates / raw_norm; radii <- sqrt(rowSums(coordinates ^ 2)); distances <- sqrt(rowSums((coordinates[use[, 1], ] - coordinates[use[, 2], ])^2)); y <- adjacency[cbind(use[, 1], use[, 2])]
    bernoulli_negative_log_likelihood(y, alpha - lambda * distances) + penalty * mean((radii - mean(radii)) ^ 2)
  }
  fit <- tryCatch(optim(start, nll, method = "L-BFGS-B", lower = c(-8, -5, rep(-5, 2 * n)), upper = c(8, 5, rep(5, 2 * n)), control = list(maxit = 600, factr = 1e8, pgtol = 1e-4)), error = function(error) NULL)
  if (is.null(fit) || !all(is.finite(fit$par))) return(list(coordinates = NULL, alpha = NA_real_, logLik = NA_real_, converged = FALSE))
  raw_coordinates <- matrix(fit$par[-c(1, 2)], nrow = n, ncol = 2); coordinates <- sqrt(n) * raw_coordinates / sqrt(sum(raw_coordinates ^ 2)); distances <- sqrt(rowSums((coordinates[use[, 1], ] - coordinates[use[, 2], ])^2)); y <- adjacency[cbind(use[, 1], use[, 2])]
  log_likelihood <- -bernoulli_negative_log_likelihood(y, fit$par[1] - exp(fit$par[2]) * distances)
  list(coordinates = coordinates, alpha = fit$par[1], lambda = exp(fit$par[2]), logLik = log_likelihood, converged = fit$convergence == 0)
}

fit_spherical_lsm <- function(adjacency, dyads = NULL, max_dyads = 6000) {
  n <- nrow(adjacency); use <- if (is.null(dyads)) which(upper.tri(adjacency), arr.ind = TRUE) else dyads
  if (nrow(use) > max_dyads) use <- use[sample(seq_len(nrow(use)), max_dyads), , drop = FALSE]
  initial <- cmdscale(as.dist(1 - adjacency), k = 3); initial <- initial / pmax(sqrt(rowSums(initial^2)), 1e-8)
  start_angles <- cbind(atan2(initial[, 2], initial[, 1]), acos(pmin(1, pmax(-1, initial[, 3])))); start <- c(qlogis(pmin(0.999, pmax(0.001, mean(adjacency[upper.tri(adjacency)])))), log(1), as.vector(t(start_angles[-1, , drop = FALSE])))
  nll <- function(par) {
    alpha <- par[1]; lambda <- exp(par[2]); angles <- rbind(start_angles[1, ], matrix(par[-c(1, 2)], ncol = 2, byrow = TRUE)); coordinates <- cbind(sin(angles[, 2]) * cos(angles[, 1]), sin(angles[, 2]) * sin(angles[, 1]), cos(angles[, 2])); eta <- alpha + lambda * rowSums(coordinates[use[, 1], ] * coordinates[use[, 2], ]); y <- adjacency[cbind(use[, 1], use[, 2])]; bernoulli_negative_log_likelihood(y, eta)
  }
  spherical_lower <- c(-8, -5, rep(c(-pi, 0.001), n - 1)); spherical_upper <- c(8, 5, rep(c(pi, pi - 0.001), n - 1))
  fit <- tryCatch(optim(start, nll, method = "L-BFGS-B", lower = spherical_lower, upper = spherical_upper, control = list(maxit = 600, factr = 1e8, pgtol = 1e-4)), error = function(error) NULL)
  if (is.null(fit) || !all(is.finite(fit$par))) return(list(coordinates = NULL, alpha = NA_real_, logLik = NA_real_, converged = FALSE))
  angles <- rbind(start_angles[1, ], matrix(fit$par[-c(1, 2)], ncol = 2, byrow = TRUE)); coordinates <- cbind(sin(angles[, 2]) * cos(angles[, 1]), sin(angles[, 2]) * sin(angles[, 1]), cos(angles[, 2])); list(coordinates = coordinates, alpha = fit$par[1], lambda = exp(fit$par[2]), logLik = -fit$value, converged = fit$convergence == 0)
}

evaluate_lsm_ties <- function(fit, adjacency, geometry, test_dyads) {
  if (is.null(fit$coordinates)) return(NA_real_)
  use <- test_dyads
  coordinates <- fit$coordinates; if (geometry == "Spherical") similarity <- tcrossprod(coordinates) else if (geometry == "Circular") similarity <- cos(outer(atan2(coordinates[, 2], coordinates[, 1]), atan2(coordinates[, 2], coordinates[, 1]), "-")) else similarity <- -as.matrix(dist(coordinates))
  eta <- if (geometry %in% c("Circular", "Spherical")) fit$alpha + fit$lambda * similarity[cbind(use[, 1], use[, 2])] else fit$alpha + fit$lambda * similarity[cbind(use[, 1], use[, 2])]
  y <- adjacency[cbind(use[, 1], use[, 2])]; -mean(y * log(plogis(eta)) + (1 - y) * log1p(-plogis(eta)))
}

#### Conditional held-out vote decoder -------------------------------------
decoder_log_loss <- function(coordinates, votes) {
  if (is.null(coordinates)) return(NA_real_)
  coordinates <- as.matrix(coordinates)[, 1:min(2, ncol(coordinates)), drop = FALSE]
  coordinate_data <- as.data.frame(coordinates)
  names(coordinate_data) <- paste0("x", seq_len(ncol(coordinate_data)))
  losses <- numeric(0)
  for (roll_call in seq_len(ncol(votes))) {
    test <- sample(seq_len(nrow(votes)), size = max(1, floor(0.20 * nrow(votes))))
    train <- setdiff(seq_len(nrow(votes)), test)
    training_data <- coordinate_data[train, , drop = FALSE]
    training_data$vote <- votes[train, roll_call]
    fit <- tryCatch(glm(vote ~ ., data = training_data, family = binomial()), error = function(error) NULL)
    if (is.null(fit)) next
    probability <- pmin(1 - 1e-6, pmax(1e-6, predict(fit, newdata = coordinate_data[test, , drop = FALSE], type = "response")))
    losses <- c(losses, -mean(votes[test, roll_call] * log(probability) + (1 - votes[test, roll_call]) * log(1 - probability)))
  }
  if (length(losses) == 0) NA_real_ else mean(losses)
}

#### Run both simulation scenarios -----------------------------------------
simulation_results <- list()
result_index <- 1

for (scenario in c("Multi-bloc", "Cyclic-policy")) {
  for (repetition in seq_len(number_of_repetitions)) {
    simulated <- generate_votes(scenario, number_of_members, number_of_roll_calls)
    disagreement <- make_disagreement(simulated$votes)
    classical_2d <- cmdscale(as.dist(disagreement), k = 2)
    classical_1d <- matrix(as.numeric(cmdscale(as.dist(disagreement), k = 1)), ncol = 1)
    circular <- fit_circular_mds(disagreement)
    soft_circular <- fit_soft_circular_mds(disagreement, circular$coordinates)
    spherical <- fit_spherical_mds(disagreement)
    fitted_coordinates <- list("1D Euclidean MDS" = classical_1d, "2D Euclidean MDS" = classical_2d, "Circular MDS" = circular$coordinates, "Soft Circular MDS" = soft_circular$coordinates, "Spherical MDS" = spherical$coordinates)
    convergence <- c("1D Euclidean MDS" = TRUE, "2D Euclidean MDS" = TRUE, "Circular MDS" = circular$converged, "Soft Circular MDS" = soft_circular$converged, "Spherical MDS" = spherical$converged)
    for (method in names(fitted_coordinates)) {
      coordinates <- fitted_coordinates[[method]]
      simulation_results[[result_index]] <- data.frame(scenario = scenario, repetition = repetition, method = method, normalized_stress = if (is.null(coordinates)) NA_real_ else normalized_stress(disagreement, fitted_geometry_distance(coordinates, method)), held_out_log_loss = decoder_log_loss(coordinates, simulated$votes), converged = if (is.null(coordinates) || !isTRUE(convergence[method])) "no" else "yes", stringsAsFactors = FALSE)
      result_index <- result_index + 1
    }
    message(scenario, " repetition ", repetition, " of ", number_of_repetitions, " complete.")
  }
}

simulation_results <- do.call(rbind, simulation_results)
simulation_summary <- aggregate(cbind(normalized_stress, held_out_log_loss) ~ scenario + method, data = simulation_results, FUN = function(x) mean(x, na.rm = TRUE))

#### LSM simulation comparison ---------------------------------------------
lsm_results <- list(); lsm_index <- 1
for (scenario in c("Multi-bloc", "Cyclic-policy")) for (repetition in seq_len(lsm_number_of_repetitions)) {
  simulated <- generate_votes(scenario, lsm_number_of_members, number_of_roll_calls); adjacency <- make_covoting_adjacency(simulated$votes); all_dyads <- which(upper.tri(adjacency), arr.ind = TRUE)
  test_rows <- sample(seq_len(nrow(all_dyads)), floor(0.20 * nrow(all_dyads))); test_dyads <- all_dyads[test_rows, , drop = FALSE]; training_dyads <- all_dyads[-test_rows, , drop = FALSE]
  if (nrow(training_dyads) > lsm_maximum_training_dyads) training_dyads <- training_dyads[sample(seq_len(nrow(training_dyads)), lsm_maximum_training_dyads), , drop = FALSE]
  lsm_fits <- list("1D Euclidean LSM" = fit_euclidean_lsm(adjacency, ndim = 1, dyads = training_dyads), "2D Euclidean LSM" = fit_euclidean_lsm(adjacency, ndim = 2, dyads = training_dyads), "Circular LSM" = fit_circular_lsm(adjacency, dyads = training_dyads), "Soft Circular LSM" = fit_soft_circular_lsm(adjacency, dyads = training_dyads), "Spherical LSM" = fit_spherical_lsm(adjacency, dyads = training_dyads))
  lsm_geometry <- c("1D Euclidean LSM" = "Euclidean", "2D Euclidean LSM" = "Euclidean", "Circular LSM" = "Circular", "Soft Circular LSM" = "Soft circular", "Spherical LSM" = "Spherical"); lsm_dimensions <- c("1D Euclidean LSM" = 1, "2D Euclidean LSM" = 2, "Circular LSM" = 1, "Soft Circular LSM" = 2, "Spherical LSM" = 2)
  for (method in names(lsm_fits)) {
    fit <- lsm_fits[[method]]; coordinates <- fit$coordinates; converged <- isTRUE(fit$converged) && !is.null(coordinates); geometry <- unname(lsm_geometry[method]); k <- if (method == "1D Euclidean LSM") lsm_number_of_members + 2 else if (method == "2D Euclidean LSM") 2 * lsm_number_of_members + 2 else if (method == "Circular LSM") lsm_number_of_members + 1 else if (method == "Soft Circular LSM") 2 * lsm_number_of_members + 2 else 2 * lsm_number_of_members
    stress_method <- if (method == "Circular LSM") "Circular MDS" else if (method == "Spherical LSM") "Spherical MDS" else "2D Euclidean MDS"
    lsm_results[[lsm_index]] <- data.frame(scenario = scenario, repetition = repetition, method = method, geometry = geometry, members = lsm_number_of_members, training_dyads = nrow(training_dyads), held_out_dyads = nrow(test_dyads), bic = if (converged && method != "Soft Circular LSM") -2 * fit$logLik + k * log(nrow(training_dyads)) else NA_real_, held_out_tie_log_loss = if (converged) evaluate_lsm_ties(fit, adjacency, if (geometry == "Soft circular") "Euclidean" else geometry, test_dyads) else NA_real_, normalized_stress = if (converged) normalized_stress(1 - adjacency, fitted_geometry_distance(coordinates, stress_method)) else NA_real_, converged = if (converged) "yes" else "no", stringsAsFactors = FALSE); lsm_index <- lsm_index + 1
  }
  message(scenario, " LSM repetition ", repetition, " of ", lsm_number_of_repetitions, " complete.")
}
lsm_results <- do.call(rbind, lsm_results)
lsm_groups <- split(seq_len(nrow(lsm_results)), interaction(lsm_results$scenario, lsm_results$method, lsm_results$geometry, drop = TRUE))
lsm_summary <- do.call(rbind, lapply(lsm_groups, function(index) {
  group_results <- lsm_results[index, , drop = FALSE]
  finite_mean <- function(values) if (any(is.finite(values))) mean(values[is.finite(values)]) else NA_real_
  data.frame(scenario = group_results$scenario[1], method = group_results$method[1], geometry = group_results$geometry[1], repetitions = nrow(group_results), successful_repetitions = sum(group_results$converged == "yes"), bic = finite_mean(group_results$bic), held_out_tie_log_loss = finite_mean(group_results$held_out_tie_log_loss), normalized_stress = finite_mean(group_results$normalized_stress), stringsAsFactors = FALSE)
}))
rownames(lsm_summary) <- NULL

#### Scenario-specific paper tables ----------------------------------------
simulation_mds_table <- simulation_results |>
  mutate(
    Family = "MDS",
    Dimension = if_else(method == "1D Euclidean MDS" | method == "Circular MDS", "1D", "2D"),
    Method = method,
    `Held-out log-loss` = held_out_log_loss,
    `Normalized stress` = normalized_stress,
    BIC = NA_real_
  ) |>
  group_by(scenario, Family, Method, Dimension) |>
  summarise(`Held-out log-loss` = if (any(is.finite(`Held-out log-loss`))) mean(`Held-out log-loss`[is.finite(`Held-out log-loss`)]) else NA_real_, `Normalized stress` = if (any(is.finite(`Normalized stress`))) mean(`Normalized stress`[is.finite(`Normalized stress`)]) else NA_real_, BIC = NA_real_, .groups = "drop")

simulation_lsm_table <- lsm_results |>
  mutate(
    Family = "LSM",
    Dimension = if_else(method %in% c("1D Euclidean LSM", "Circular LSM"), "1D", "2D"),
    Method = method,
    `Held-out log-loss` = held_out_tie_log_loss,
    `Normalized stress` = normalized_stress,
    BIC = bic
  ) |>
  group_by(scenario, Family, Method, Dimension) |>
  summarise(`Held-out log-loss` = if (any(is.finite(`Held-out log-loss`))) mean(`Held-out log-loss`[is.finite(`Held-out log-loss`)]) else NA_real_, `Normalized stress` = if (any(is.finite(`Normalized stress`))) mean(`Normalized stress`[is.finite(`Normalized stress`)]) else NA_real_, BIC = if (any(is.finite(BIC))) mean(BIC[is.finite(BIC)]) else NA_real_, .groups = "drop")

scenario_comparison_table <- bind_rows(simulation_mds_table, simulation_lsm_table) |>
  select(scenario, Family, Method, Dimension, `Held-out log-loss`, `Normalized stress`, BIC) |>
  mutate(
    Method = factor(
      Method,
      levels = c(
        "1D Euclidean MDS", "2D Euclidean MDS", "Circular MDS",
        "Soft Circular MDS", "Spherical MDS",
        "1D Euclidean LSM", "2D Euclidean LSM", "Circular LSM",
        "Soft Circular LSM", "Spherical LSM"
      )
    )
  ) |>
  arrange(scenario, match(Family, c("MDS", "LSM")), Method) |>
  mutate(Method = as.character(Method))

for (scenario_name in c("Cyclic-policy", "Multi-bloc")) {
  scenario_file_name <- if (scenario_name == "Cyclic-policy") "simulation_cyclic_policy.csv" else "simulation_multi_bloc.csv"
  write.csv(scenario_comparison_table |> filter(scenario == scenario_name) |> select(-scenario), file.path(results_directory, scenario_file_name), row.names = FALSE, na = "NA")
}

print(simulation_summary)
print(lsm_summary)
message("Simulation feasibility study complete. MDS models use held-out vote-entry log-loss; LSMs use held-out member-pair tie log-loss and BIC. Spherical MDS results are blank when smacof is unavailable or fails its iteration check.")
