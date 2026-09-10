#### Voteview CRS policy-profile geometry comparison ------------------------
# This extension uses 119th House Members' Votes and CRS policy metadata.
# It compares dimension-reduction representations of member-by-policy-area
# profiles. It does not fit LSMs: a policy profile is not an observed network.

# Packages and working directory
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(smacof)
library(FactoMineR)
library(Gifi)

setwd("E:/MSc_dissertation/data")
set.seed(123)

results_directory <- "voteview_policy_profile_results"
if (!dir.exists(results_directory)) dir.create(results_directory)

#### Read Members' Votes and CRS metadata -----------------------------------

votes_119 <- read_csv("H119_votes.csv", show_col_types = FALSE) |>
  filter(congress == 119, chamber == "House") |>
  mutate(vote_binary = case_when(
    cast_code %in% c(1, 2, 3) ~ 1,
    cast_code %in% c(4, 5, 6) ~ 0,
    TRUE ~ NA_real_
  ))

members_119 <- read_csv("H119_members.csv", show_col_types = FALSE) |>
  filter(congress == 119, chamber == "House") |>
  mutate(party = case_when(
    party_code == 100 ~ "Democrat",
    party_code == 200 ~ "Republican",
    TRUE ~ "Other"
  ))

rollcalls_119 <- jsonlite::fromJSON("H119_rollcalls.json", simplifyDataFrame = TRUE) |>
  as_tibble() |>
  filter(congress == 119, chamber == "House") |>
  transmute(congress, chamber, rollnumber, crs_policy_area)

if (!"crs_policy_area" %in% names(rollcalls_119)) {
  stop("H119_rollcalls.json does not contain crs_policy_area.")
}

#### Keep policy areas with enough labelled roll calls -----------------------

# A minimum of seven roll calls avoids profiles defined by only one or two votes.
minimum_roll_calls_per_area <- 7

area_coverage <- rollcalls_119 |>
  filter(!is.na(crs_policy_area), crs_policy_area != "") |>
  count(crs_policy_area, name = "number_of_roll_calls") |>
  arrange(desc(number_of_roll_calls))

selected_areas <- area_coverage |>
  filter(number_of_roll_calls >= minimum_roll_calls_per_area) |>
  pull(crs_policy_area)

if (length(selected_areas) != 18) {
  warning("Expected 18 CRS areas at the seven-roll-call threshold; retained ",
          length(selected_areas), ". Inspect the printed area coverage before reporting results.")
}

print(area_coverage)

selected_rollcalls <- rollcalls_119 |> filter(crs_policy_area %in% selected_areas)

# Create member-by-policy-area residual voting profiles
labelled_votes <- votes_119 |>
  inner_join(selected_rollcalls, by = c("congress", "chamber", "rollnumber")) |>
  filter(icpsr %in% members_119$icpsr) |>
  group_by(rollnumber) |>
  mutate(roll_call_yea_share = mean(vote_binary, na.rm = TRUE)) |>
  ungroup() |>
  mutate(vote_residual = vote_binary - roll_call_yea_share)

informative_rollcalls <- labelled_votes |>
  group_by(rollnumber, crs_policy_area) |>
  summarise(
    number_observed = sum(!is.na(vote_binary)),
    yea_share = mean(vote_binary, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(number_observed >= 0.80 * max(number_observed), yea_share >= 0.05, yea_share <= 0.95)

labelled_votes <- labelled_votes |>
  semi_join(informative_rollcalls, by = c("rollnumber", "crs_policy_area"))

profile_long <- labelled_votes |>
  group_by(icpsr, crs_policy_area) |>
  summarise(
    number_of_votes = sum(!is.na(vote_residual)),
    policy_profile = mean(vote_residual, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(policy_profile = ifelse(is.nan(policy_profile), NA_real_, policy_profile))

profile_matrix <- profile_long |>
  select(icpsr, crs_policy_area, policy_profile) |>
  pivot_wider(names_from = crs_policy_area, values_from = policy_profile) |>
  arrange(match(icpsr, members_119$icpsr)) |>
  column_to_rownames("icpsr") |>
  as.matrix()

# Apply the same completeness rule to policy profiles
repeat {
  previous_dimensions <- dim(profile_matrix)
  member_completeness <- rowMeans(!is.na(profile_matrix))
  profile_matrix <- profile_matrix[member_completeness >= 0.80, , drop = FALSE]
  area_completeness <- colMeans(!is.na(profile_matrix))
  profile_matrix <- profile_matrix[, area_completeness >= 0.80, drop = FALSE]
  if (identical(dim(profile_matrix), previous_dimensions)) break
}

if (nrow(profile_matrix) < 10 || ncol(profile_matrix) < 3) {
  stop("Too few members or policy areas remain after completeness filtering.")
}

member_ids <- rownames(profile_matrix)
members_profile <- members_119 |>
  filter(icpsr %in% as.numeric(member_ids)) |>
  arrange(match(icpsr, as.numeric(member_ids)))

if (!identical(as.character(members_profile$icpsr), member_ids)) {
  stop("Member metadata could not be aligned to the policy profile matrix.")
}

# Matrix methods require a complete numeric matrix. Column means are used only
# after the observed profile and completeness checks have been recorded.
profile_imputed <- profile_matrix
for (area_index in seq_len(ncol(profile_imputed))) {
  missing_members <- is.na(profile_imputed[, area_index])
  profile_imputed[missing_members, area_index] <- mean(profile_imputed[, area_index], na.rm = TRUE)
}

# Each policy area receives equal weight in Euclidean profile distances.
profile_scaled <- scale(profile_imputed, center = TRUE, scale = TRUE)
profile_scaled <- as.matrix(profile_scaled)

if (any(!is.finite(profile_scaled))) {
  stop("The scaled policy-profile matrix contains non-finite values.")
}

profile_summary <- data.frame(
  number_of_members = nrow(profile_scaled),
  number_of_policy_areas = ncol(profile_scaled),
  number_of_labelled_roll_calls = nrow(informative_rollcalls),
  missing_profile_proportion_before_imputation = mean(is.na(profile_matrix)),
  stringsAsFactors = FALSE
)
print(profile_summary)

# Prepare ordered profiles for categorical methods
profile_ordinal <- as.data.frame(profile_scaled)
ordinal_labels <- c("Lower", "Middle", "Upper")

for (area_name in names(profile_ordinal)) {
  
  values <- profile_ordinal[[area_name]]
  
  # Divide non-missing observations into approximately equal-sized groups.
  group_index <- dplyr::ntile(values, n = 3)
  
  profile_ordinal[[area_name]] <- ordered(
    ordinal_labels[group_index],
    levels = ordinal_labels
  )
}
rownames(profile_ordinal) <- member_ids

# Construct the common profile-distance matrix
profile_distance <- as.matrix(dist(profile_scaled))
profile_distance_values <- profile_distance[upper.tri(profile_distance)]

if (!all(is.finite(profile_distance_values))) {
  stop("Profile distances contain non-finite values.")
}

if (length(profile_distance_values) == 0 ||
    !any(profile_distance_values > 0)) {
  stop("Profile distances are completely degenerate.")
}

normalized_stress <- function(observed_distance, fitted_distance) {
  if (is.null(observed_distance) || is.null(fitted_distance) ||
      !identical(dim(observed_distance), dim(fitted_distance))) {
    return(NA_real_)
  }
  observed_values <- observed_distance[upper.tri(observed_distance)]
  fitted_values <- fitted_distance[upper.tri(fitted_distance)]
  use <- is.finite(observed_values) & is.finite(fitted_values)
  if (sum(use) < 2 || sum(observed_values[use] ^ 2) <= 1e-12 ||
      sum(fitted_values[use] ^ 2) <= 1e-12) return(NA_real_)
  scale_factor <- sum(observed_values[use] * fitted_values[use]) / sum(fitted_values[use] ^ 2)
  sqrt(sum((observed_values[use] - scale_factor * fitted_values[use]) ^ 2) /
         sum(observed_values[use] ^ 2))
}

#### Fit Euclidean, circular, and spherical MDS -----------------------------
message("Fitting Euclidean MDS representations ...")

euclidean_1d_mds_fit <- cmdscale(as.dist(profile_distance), k = 1, eig = TRUE)
euclidean_1d_mds_coordinates <- matrix(
  euclidean_1d_mds_fit$points[, 1],
  ncol = 1,
  dimnames = list(member_ids, "Dimension_1")
)

classical_mds_fit <- cmdscale(as.dist(profile_distance), k = 2, eig = TRUE)
classical_mds_coordinates <- classical_mds_fit$points[, 1:2, drop = FALSE]
rownames(classical_mds_coordinates) <- member_ids

metric_mds_fit <- tryCatch(
  smacof::smacofSym(as.dist(profile_distance), ndim = 2, type = "ratio", itmax = 1000),
  error = function(error) error
)
metric_mds_coordinates <- if (inherits(metric_mds_fit, "error")) NULL else metric_mds_fit$conf[, 1:2, drop = FALSE]
if (!is.null(metric_mds_coordinates)) rownames(metric_mds_coordinates) <- member_ids

ordinal_mds_fit <- tryCatch(
  smacof::smacofSym(as.dist(profile_distance), ndim = 2, type = "ordinal", itmax = 1000),
  error = function(error) error
)
ordinal_mds_coordinates <- if (inherits(ordinal_mds_fit, "error")) NULL else ordinal_mds_fit$conf[, 1:2, drop = FALSE]
if (!is.null(ordinal_mds_coordinates)) rownames(ordinal_mds_coordinates) <- member_ids

message("Fitting strict circular MDS ...")

# Circular geodesic distances lie in [0, pi], so the observed target is scaled
# to the same interval before fitting angles.
angular_target <- pi * profile_distance / max(profile_distance)
initial_angles <- atan2(classical_mds_coordinates[, 2], classical_mds_coordinates[, 1])
initial_angles <- initial_angles - initial_angles[1]
initial_angles <- atan2(sin(initial_angles), cos(initial_angles))

circular_pairs <- which(upper.tri(profile_distance), arr.ind = TRUE)
circular_target_values <- angular_target[circular_pairs]

# Use only the upper-triangular dyads.
circular_objective <- function(free_angles) {
  theta <- c(0, free_angles)
  angle_difference <- theta[circular_pairs[, 1]] - theta[circular_pairs[, 2]]
  absolute_difference <- abs(angle_difference)
  fitted_distance <- pmin(absolute_difference, 2 * pi - absolute_difference)
  sum((fitted_distance - circular_target_values) ^ 2)
}

# Analytic gradient for the upper-triangular objective.
circular_gradient <- function(free_angles) {
  theta <- c(0, free_angles)
  angle_difference <- theta[circular_pairs[, 1]] - theta[circular_pairs[, 2]]
  absolute_difference <- abs(angle_difference)
  fitted_distance <- pmin(absolute_difference, 2 * pi - absolute_difference)

  distance_derivative <- sign(angle_difference)
  distance_derivative[absolute_difference > pi] <-
    -distance_derivative[absolute_difference > pi]
  distance_derivative[absolute_difference == 0 |
                        absolute_difference == pi] <- 0

  pair_gradient <- 2 * (fitted_distance - circular_target_values) *
    distance_derivative
  gradient <- numeric(length(theta))
  first_endpoint <- tapply(pair_gradient, circular_pairs[, 1], sum)
  second_endpoint <- tapply(pair_gradient, circular_pairs[, 2], sum)
  gradient[as.integer(names(first_endpoint))] <- first_endpoint
  gradient[as.integer(names(second_endpoint))] <-
    gradient[as.integer(names(second_endpoint))] - second_endpoint
  gradient[-1]
}

circular_mds_fit <- tryCatch(
  optim(
    par = initial_angles[-1],
    fn = circular_objective,
    gr = circular_gradient,
    method = "L-BFGS-B",
    lower = rep(-pi, length(initial_angles) - 1),
    upper = rep(pi, length(initial_angles) - 1),
    control = list(maxit = 500, factr = 1e7, pgtol = 1e-6)
  ),
  error = function(error) error
)

circular_mds_coordinates <- if (inherits(circular_mds_fit, "error")) NULL else {
  circular_angles <- c(0, circular_mds_fit$par %% (2 * pi))
  coordinates <- cbind(Dimension_1 = cos(circular_angles), Dimension_2 = sin(circular_angles))
  rownames(coordinates) <- member_ids
  coordinates
}
circular_mds_converged <- !inherits(circular_mds_fit, "error") && circular_mds_fit$convergence == 0

#### Fit soft circular MDS --------------------------------------------------

# Soft circular MDS retains the fitted circular angles but estimates a separate
# radius for each member. A coefficient-of-variation penalty prevents extreme
# radial dispersion while allowing departures from a strict circle.
soft_circular_lambda <- 5
soft_circular_angles <- if (circular_mds_converged) {
  c(0, circular_mds_fit$par %% (2 * pi))
} else {
  rep(NA_real_, length(member_ids))
}
soft_circular_cosines <- cos(
  soft_circular_angles[circular_pairs[, 1]] -
    soft_circular_angles[circular_pairs[, 2]]
)

soft_circular_objective <- function(log_radii) {
  radii <- exp(log_radii)
  radius_i <- radii[circular_pairs[, 1]]
  radius_k <- radii[circular_pairs[, 2]]
  fitted_distance <- sqrt(pmax(
    radius_i ^ 2 + radius_k ^ 2 - 2 * radius_i * radius_k * soft_circular_cosines,
    1e-12
  ))
  radius_cv_squared <- (sd(radii) / mean(radii)) ^ 2
  mean((fitted_distance - circular_target_values) ^ 2) +
    soft_circular_lambda * radius_cv_squared
}

soft_circular_gradient <- function(log_radii) {
  radii <- exp(log_radii)
  radius_i <- radii[circular_pairs[, 1]]
  radius_k <- radii[circular_pairs[, 2]]
  fitted_distance <- sqrt(pmax(
    radius_i ^ 2 + radius_k ^ 2 - 2 * radius_i * radius_k * soft_circular_cosines,
    1e-12
  ))
  residual <- fitted_distance - circular_target_values

  derivative_i <- radius_i *
    (radius_i - radius_k * soft_circular_cosines) / fitted_distance
  derivative_k <- radius_k *
    (radius_k - radius_i * soft_circular_cosines) / fitted_distance
  pair_gradient_i <- 2 * residual * derivative_i / nrow(circular_pairs)
  pair_gradient_k <- 2 * residual * derivative_k / nrow(circular_pairs)

  gradient <- numeric(length(radii))
  first_endpoint <- tapply(pair_gradient_i, circular_pairs[, 1], sum)
  second_endpoint <- tapply(pair_gradient_k, circular_pairs[, 2], sum)
  gradient[as.integer(names(first_endpoint))] <- first_endpoint
  gradient[as.integer(names(second_endpoint))] <-
    gradient[as.integer(names(second_endpoint))] + second_endpoint

  radius_mean <- mean(radii)
  radius_sum_squares <- sum((radii - radius_mean) ^ 2)
  penalty_derivative <- soft_circular_lambda * radii * (
    2 * (radii - radius_mean) / ((length(radii) - 1) * radius_mean ^ 2) -
      2 * radius_sum_squares / ((length(radii) - 1) * length(radii) * radius_mean ^ 3)
  )

  gradient + penalty_derivative
}

soft_circular_mds_fit <- if (!circular_mds_converged) {
  structure(list(message = "Circular MDS did not converge."), class = "error")
} else {
  tryCatch(
    optim(
      par = rep(0, length(member_ids)),
      fn = soft_circular_objective,
      gr = soft_circular_gradient,
      method = "L-BFGS-B",
      lower = rep(log(0.2), length(member_ids)),
      upper = rep(log(3), length(member_ids)),
      control = list(maxit = 500, factr = 1e7, pgtol = 1e-6)
    ),
    error = function(error) error
  )
}

soft_circular_mds_coordinates <- if (inherits(soft_circular_mds_fit, "error") ||
                                      soft_circular_mds_fit$convergence != 0) {
  NULL
} else {
  soft_circular_radii <- exp(soft_circular_mds_fit$par)
  coordinates <- cbind(
    Dimension_1 = soft_circular_radii * cos(soft_circular_angles),
    Dimension_2 = soft_circular_radii * sin(soft_circular_angles)
  )
  rownames(coordinates) <- member_ids
  coordinates
}

message("Fitting strict spherical MDS ...")

spherical_mds_itmax <- 1000
spherical_mds_fit <- tryCatch(
  smacof::smacofSphere(
    # Use three Cartesian coordinates for the spherical representation.
    as.dist(profile_distance), ndim = 3, type = "interval", algorithm = "dual",
    init = "torgerson", penalty = 100, itmax = spherical_mds_itmax, eps = 1e-4, verbose = TRUE
  ),
  error = function(error) error
)

spherical_mds_coordinates <- if (inherits(spherical_mds_fit, "error") ||
                                  is.null(spherical_mds_fit$conf) ||
                                  is.null(dim(spherical_mds_fit$conf)) ||
                                  ncol(spherical_mds_fit$conf) < 3) {
  NULL
} else {
  spherical_mds_fit$conf[, 1:3, drop = FALSE]
}
if (!is.null(spherical_mds_coordinates)) {
  spherical_radii <- sqrt(rowSums(spherical_mds_coordinates ^ 2))
  if (any(!is.finite(spherical_radii)) || any(spherical_radii <= 1e-12)) {
    spherical_mds_coordinates <- NULL
  } else {
    spherical_mds_coordinates <- sweep(
      spherical_mds_coordinates, 1, spherical_radii, "/"
    )
    rownames(spherical_mds_coordinates) <- member_ids
  }
}

spherical_pairwise_distances <- if (is.null(spherical_mds_coordinates)) numeric(0) else {
  spherical_cosines <- tcrossprod(spherical_mds_coordinates)
  spherical_cosines[spherical_cosines > 1] <- 1
  spherical_cosines[spherical_cosines < -1] <- -1
  acos(spherical_cosines)[upper.tri(profile_distance)]
}
spherical_distance_spread <- if (length(spherical_pairwise_distances) < 2 || any(!is.finite(spherical_pairwise_distances))) {
  NA_real_
} else {
  spherical_distance_mean <- mean(spherical_pairwise_distances)
  sqrt(sum((spherical_pairwise_distances - spherical_distance_mean) ^ 2) /
         (length(spherical_pairwise_distances) - 1))
}
spherical_mds_converged <- !inherits(spherical_mds_fit, "error") &&
  is.finite(spherical_mds_fit$stress) && spherical_mds_fit$niter < spherical_mds_itmax &&
  !is.null(spherical_mds_coordinates) &&
  is.finite(spherical_distance_spread) && spherical_distance_spread > 1e-10

#### Fit PCA, correspondence, optimal-scaling, and PGA methods --------------
message("Fitting profile dimension-reduction methods ...")

extract_object_scores <- function(fit) {
  scores <- fit$objectscores
  if (is.null(scores)) scores <- fit$objscores
  if (is.null(scores) || is.null(dim(scores)) || ncol(scores) < 2) return(NULL)
  as.matrix(scores[, 1:2, drop = FALSE])
}

pca_fit <- prcomp(profile_scaled, center = TRUE, scale. = FALSE)
pca_coordinates <- pca_fit$x[, 1:2, drop = FALSE]

kernel_pca_fit <- if (requireNamespace("kernlab", quietly = TRUE)) {
  tryCatch(kernlab::kpca(profile_scaled, kernel = "rbfdot", kpar = list(sigma = 0.05), features = 2),
           error = function(error) error)
} else NULL
kernel_pca_coordinates <- if (is.null(kernel_pca_fit) || inherits(kernel_pca_fit, "error")) NULL else {
  as.matrix(kernlab::rotated(kernel_pca_fit))[, 1:2, drop = FALSE]
}

profile_ca_input <- sweep(profile_scaled, 2, apply(profile_scaled, 2, min), FUN = "-") + 1e-6
ca_fit <- tryCatch(FactoMineR::CA(profile_ca_input, graph = FALSE), error = function(error) error)
ca_coordinates <- if (inherits(ca_fit, "error")) NULL else ca_fit$row$coord[, 1:2, drop = FALSE]

mca_fit <- tryCatch(FactoMineR::MCA(profile_ordinal, ncp = 2, graph = FALSE), error = function(error) error)
mca_coordinates <- if (inherits(mca_fit, "error")) NULL else mca_fit$ind$coord[, 1:2, drop = FALSE]

nlpca_fit <- tryCatch(Gifi::princals(profile_ordinal, ndim = 2, ordinal = TRUE), error = function(error) error)
nlpca_coordinates <- if (inherits(nlpca_fit, "error")) NULL else extract_object_scores(nlpca_fit)

homals_fit <- tryCatch(Gifi::homals(profile_ordinal, ndim = 2), error = function(error) error)
homals_coordinates <- if (inherits(homals_fit, "error")) NULL else extract_object_scores(homals_fit)

odd_profiles <- profile_ordinal[, seq(1, ncol(profile_ordinal), by = 2), drop = FALSE]
even_profiles <- profile_ordinal[, seq(2, ncol(profile_ordinal), by = 2), drop = FALSE]
nlcca_fit <- tryCatch({
  odd_fit <- Gifi::princals(odd_profiles, ndim = 2, ordinal = TRUE)
  even_fit <- Gifi::princals(even_profiles, ndim = 2, ordinal = TRUE)
  odd_scores <- extract_object_scores(odd_fit)
  even_scores <- extract_object_scores(even_fit)
  if (is.null(odd_scores) || is.null(even_scores)) stop("Gifi did not return two-dimensional object scores.")
  canonical <- cancor(odd_scores, even_scores)
  list(coordinates = (scale(odd_scores %*% canonical$xcoef[, 1:2]) +
    scale(even_scores %*% canonical$ycoef[, 1:2])) / 2)
}, error = function(error) error)
nlcca_coordinates <- if (inherits(nlcca_fit, "error")) NULL else nlcca_fit$coordinates

# PGA compares directions of centred policy profiles on the unit hypersphere.
pga_centred <- sweep(profile_scaled, 2, colMeans(profile_scaled), FUN = "-")
pga_norms <- sqrt(rowSums(pga_centred ^ 2))
pga_unit_profiles <- sweep(pga_centred, 1, pmax(pga_norms, 1e-12), FUN = "/")
pga_mean <- colMeans(pga_unit_profiles)
pga_mean <- pga_mean / sqrt(sum(pga_mean ^ 2))
pga_coefficients <- as.vector(pga_unit_profiles %*% pga_mean)
pga_tangent <- pga_unit_profiles - pga_coefficients * matrix(pga_mean, nrow(pga_unit_profiles), ncol(pga_unit_profiles), byrow = TRUE)
pga_fit <- prcomp(pga_tangent, center = TRUE, scale. = FALSE)
pga_coordinates <- pga_fit$x[, 1:2, drop = FALSE]

#### Calculate a common distance-preservation comparison --------------------
method_coordinates <- list(
  "PCA" = pca_coordinates,
  "Kernel PCA" = kernel_pca_coordinates,
  "NLPCA" = nlpca_coordinates,
  "1D Euclidean MDS" = euclidean_1d_mds_coordinates,
  "Classical MDS" = classical_mds_coordinates,
  "Metric MDS" = metric_mds_coordinates,
  "Ordinal MDS" = ordinal_mds_coordinates,
  "Circular MDS" = circular_mds_coordinates,
  "Soft Circular MDS" = soft_circular_mds_coordinates,
  "Spherical MDS" = spherical_mds_coordinates,
  "CA" = ca_coordinates,
  "MCA" = mca_coordinates,
  "HOMALS" = homals_coordinates,
  "NLCCA" = nlcca_coordinates,
  "PGA" = pga_coordinates
)

method_metadata <- data.frame(
  method = names(method_coordinates),
  family = c(rep("PCA family", 3), rep("MDS family", 7), rep("Correspondence / optimal scaling", 4), "PGA"),
  geometry = c("Euclidean", "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Circular", "Soft circular", "Spherical", "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Tangent-space"),
  stringsAsFactors = FALSE
)

method_converged <- c(
  "PCA" = TRUE,
  "Kernel PCA" = !is.null(kernel_pca_coordinates),
  "NLPCA" = !is.null(nlpca_coordinates),
  "1D Euclidean MDS" = TRUE,
  "Classical MDS" = TRUE,
  "Metric MDS" = !is.null(metric_mds_coordinates),
  "Ordinal MDS" = !is.null(ordinal_mds_coordinates),
  "Circular MDS" = circular_mds_converged,
  "Soft Circular MDS" = !is.null(soft_circular_mds_coordinates),
  "Spherical MDS" = spherical_mds_converged,
  "CA" = !is.null(ca_coordinates),
  "MCA" = !is.null(mca_coordinates),
  "HOMALS" = !is.null(homals_coordinates),
  "NLCCA" = !is.null(nlcca_coordinates),
  "PGA" = TRUE
)

coordinate_distance <- function(coordinates, geometry) {
  if (is.null(coordinates) || any(!is.finite(coordinates))) return(NULL)
  if (geometry == "Circular") {
    angle <- atan2(coordinates[, 2], coordinates[, 1])
    difference <- abs(outer(angle, angle, "-"))
    return(pmin(difference, 2 * pi - difference))
  }
  if (geometry == "Soft circular") {
    return(as.matrix(dist(coordinates)))
  }
  if (geometry == "Spherical") {
    coordinate_radii <- sqrt(rowSums(coordinates ^ 2))
    if (any(!is.finite(coordinate_radii)) || any(coordinate_radii <= 1e-12)) {
      return(NULL)
    }
    unit_coordinates <- sweep(coordinates, 1, coordinate_radii, "/")
    cosine_similarity <- tcrossprod(unit_coordinates)
    cosine_similarity[cosine_similarity > 1] <- 1
    cosine_similarity[cosine_similarity < -1] <- -1
    return(acos(cosine_similarity))
  }
  as.matrix(dist(coordinates))
}

method_metrics <- do.call(rbind, lapply(seq_len(nrow(method_metadata)), function(method_index) {
  method_name <- method_metadata$method[method_index]
  geometry <- method_metadata$geometry[method_index]
  coordinates <- method_coordinates[[method_name]]
  fitted_distance <- coordinate_distance(coordinates, geometry)
  stress <- normalized_stress(profile_distance, fitted_distance)
  converged <- isTRUE(method_converged[method_name])

  data.frame(
    method = method_name,
    family = method_metadata$family[method_index],
    geometry = geometry,
    normalized_stress = stress,
    converged = converged,
    status = ifelse(converged && is.finite(stress), "ok", "failed"),
    stringsAsFactors = FALSE
  )
})) |>
  arrange(normalized_stress)

geometry_comparison <- method_metrics |>
  filter(method %in% c("Metric MDS", "Circular MDS", "Soft Circular MDS", "Spherical MDS")) |>
  select(method, geometry, normalized_stress, status) |>
  arrange(normalized_stress)

#### Save the dimension-reduction comparison ------------------
policy_profile_comparison_table <- method_metrics |>
  transmute(
    Family = case_when(
      family == "PCA family" ~ "PCA",
      family == "MDS family" ~ "MDS",
      family == "Correspondence / optimal scaling" ~ "CA",
      family == "PGA" ~ "PCA"
    ),
    Method = recode(
      method,
      "PCA" = "Linear PCA",
      "NLPCA" = "Nonlinear PCA",
      "MCA" = "Multiple CA",
      "NLCCA" = "Nonlinear CCA"
    ),
    Dimension = if_else(method %in% c("1D Euclidean MDS", "Circular MDS"), "1D", "2D"),
    `Normalized stress` = normalized_stress,
    method_order = match(
      method,
      c(
        "PCA", "Kernel PCA", "NLPCA", "PGA",
        "1D Euclidean MDS", "Classical MDS", "Metric MDS", "Ordinal MDS",
        "Circular MDS", "Soft Circular MDS", "Spherical MDS",
        "CA", "MCA", "HOMALS", "NLCCA"
      )
    )
  ) |>
  arrange(match(Family, c("PCA", "MDS", "CA")), method_order) |>
  select(Family, Method, Dimension, `Normalized stress`)

write_csv(
  policy_profile_comparison_table,
  file.path(results_directory, "policy_profile_dimension_reduction_comparison.csv")
)

print(policy_profile_comparison_table)
