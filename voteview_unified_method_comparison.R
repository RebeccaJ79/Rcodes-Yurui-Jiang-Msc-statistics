#### Unified Members' Votes method comparison --------------------------------
# This script compares methods from one cleaned 119th House Members' Votes sample.
# Each method uses the representation required by its statistical formulation.
# Run each section from top to bottom in RStudio.

#### Packages and working directory ------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(smacof)
library(FactoMineR)
library(Gifi)

setwd("E:/MSc_dissertation/data")

set.seed(123)

results_directory <- "voteview_unified_results"

if (!dir.exists(results_directory)) {
  dir.create(results_directory)
}

# Keep one paper-ready comparison table. Earlier auxiliary comparison CSVs are
# removed so that the results directory does not contain competing table versions.
obsolete_comparison_csvs <- c("geometry_comparison.csv", "within_family_comparison.csv", "dimension_reduction_vs_lsm.csv", "one_vs_two_dimensional.csv")
for (obsolete_file in obsolete_comparison_csvs) {
  obsolete_path <- file.path(results_directory, obsolete_file)
  if (file.exists(obsolete_path)) file.remove(obsolete_path)
}

#### Read and filter the full 119th House Voteview dataset ------------------

votes_119 <- read_csv(
  "H119_votes.csv",
  show_col_types = FALSE
) |>
  filter(
    chamber == "House"
  ) |>
  mutate(
    vote_binary = case_when(
      cast_code %in% c(1, 2, 3) ~ 1,
      cast_code %in% c(4, 5, 6) ~ 0,
      TRUE ~ NA_real_
    )
  )

members_119 <- read_csv(
  "H119_members.csv",
  show_col_types = FALSE
) |>
  filter(
    chamber == "House"
  ) |>
  mutate(
    party = case_when(
      party_code == 100 ~ "Democrat",
      party_code == 200 ~ "Republican",
      TRUE ~ "Other"
    )
  )

#### Keep informative roll calls ---------------------------------------------

roll_call_summary <- votes_119 |>
  group_by(
    rollnumber
  ) |>
  summarise(
    number_observed = sum(!is.na(vote_binary)),
    yea_share = mean(vote_binary, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(
    number_observed >= 0.80 * max(number_observed),
    yea_share >= 0.05,
    yea_share <= 0.95
  )

#### Construct the full legislator-by-roll-call vote matrix -----------------

vote_matrix_full <- votes_119 |>
  filter(
    icpsr %in% members_119$icpsr,
    rollnumber %in% roll_call_summary$rollnumber
  ) |>
  select(
    icpsr,
    rollnumber,
    vote_binary
  ) |>
  distinct() |>
  pivot_wider(
    names_from = rollnumber,
    values_from = vote_binary,
    names_prefix = "vote_"
  )

#### Read roll-call metadata for audit only ----------------------------------

# Metadata is retained for reproducibility and coverage reporting. NOMINATE
# coordinates are deliberately excluded: they must not enter any model input.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required to audit H119_rollcalls.json.")
}

rollcalls_119_raw <- jsonlite::fromJSON(
  "H119_rollcalls.json",
  simplifyDataFrame = TRUE
) |>
  as_tibble() |>
  filter(congress == 119, chamber == "House")

required_rollcall_keys <- c("congress", "chamber", "rollnumber")

if (!all(required_rollcall_keys %in% names(rollcalls_119_raw))) {
  stop("H119_rollcalls.json is missing one or more roll-call key columns.")
}

if (anyDuplicated(rollcalls_119_raw[required_rollcall_keys]) > 0) {
  stop("H119_rollcalls.json contains duplicate congress--chamber--rollnumber keys.")
}

metadata_columns <- intersect(
  c("congress", "chamber", "rollnumber", "date", "session", "clerk_rollnumber",
    "yea_count", "nay_count", "bill_number", "vote_result", "vote_desc",
    "vote_question", "dtl_desc", "issue_codes", "peltzman_codes", "clausen_codes",
    "crs_policy_area", "crs_subjects", "congress_url", "source_documents"),
  names(rollcalls_119_raw)
)

rollcalls_119 <- rollcalls_119_raw |>
  select(all_of(metadata_columns))

conflicting_vote_records <- votes_119 |>
  filter(
    icpsr %in% members_119$icpsr,
    rollnumber %in% roll_call_summary$rollnumber,
    !is.na(vote_binary)
  ) |>
  group_by(icpsr, rollnumber) |>
  summarise(number_of_distinct_votes = n_distinct(vote_binary), .groups = "drop") |>
  filter(number_of_distinct_votes > 1)

if (nrow(conflicting_vote_records) > 0) {
  stop("Some legislator--roll-call pairs have conflicting binary vote records.")
}

vote_mat <- vote_matrix_full |>
  column_to_rownames("icpsr") |>
  as.matrix()

#### Apply completeness filtering --------------------------------------------

repeat {
  old_dimensions <- dim(vote_mat)
  member_completeness <- rowMeans(!is.na(vote_mat))
  vote_mat <- vote_mat[member_completeness >= 0.80, , drop = FALSE]
  roll_call_completeness <- colMeans(!is.na(vote_mat))
  vote_mat <- vote_mat[, roll_call_completeness >= 0.80, drop = FALSE]

  if (identical(dim(vote_mat), old_dimensions)) {
    break
  }
}

if (nrow(vote_mat) < 2 || ncol(vote_mat) < 2) {
  stop("Completeness filtering retained too few legislators or roll calls.")
}

if (any(!is.na(vote_mat) & !vote_mat %in% c(0, 1))) {
  stop("The final vote matrix must contain only 0, 1, or NA.")
}

if (any(rowMeans(!is.na(vote_mat)) < 0.80) ||
    any(colMeans(!is.na(vote_mat)) < 0.80)) {
  stop("The final vote matrix does not satisfy the 80% completeness rule.")
}

member_ids <- rownames(vote_mat)

members_full <- members_119 |>
  filter(
    icpsr %in% as.numeric(member_ids)
  ) |>
  arrange(
    match(icpsr, as.numeric(member_ids))
  )

if (!identical(as.character(members_full$icpsr), member_ids)) {
  stop("Member metadata could not be aligned to the final vote matrix.")
}

retained_rollnumbers <- as.numeric(sub("^vote_", "", colnames(vote_mat)))

if (any(!is.finite(retained_rollnumbers))) {
  stop("The final vote matrix has roll-call names that cannot be parsed.")
}

rollcalls_119 <- rollcalls_119 |>
  filter(rollnumber %in% retained_rollnumbers) |>
  arrange(match(rollnumber, retained_rollnumbers))

if (nrow(rollcalls_119) != ncol(vote_mat) ||
    !identical(as.numeric(rollcalls_119$rollnumber), retained_rollnumbers)) {
  stop("Roll-call metadata could not be aligned to the final vote matrix.")
}

rollcall_metadata_audit <- data.frame(
  number_of_retained_roll_calls = nrow(rollcalls_119),
  number_with_crs_policy_area = sum(!is.na(rollcalls_119$crs_policy_area) &
                                      rollcalls_119$crs_policy_area != ""),
  proportion_with_crs_policy_area = mean(!is.na(rollcalls_119$crs_policy_area) &
                                           rollcalls_119$crs_policy_area != ""),
  number_with_crs_subjects = sum(!is.na(rollcalls_119$crs_subjects) &
                                   rollcalls_119$crs_subjects != ""),
  proportion_with_crs_subjects = mean(!is.na(rollcalls_119$crs_subjects) &
                                        rollcalls_119$crs_subjects != ""),
  stringsAsFactors = FALSE
)

data_representation_map <- data.frame(
  method_family = c("Dimension reduction", "MDS", "LSM"),
  model_input = c("member x roll-call vote profile", "pairwise voting disagreement",
                  "derived member-member co-voting network"),
  source = rep("common cleaned Members' Votes sample", 3),
  stringsAsFactors = FALSE
)

dim(vote_mat)
mean(is.na(vote_mat))

message(
  "Retained ", nrow(vote_mat),
  " legislators and ", ncol(vote_mat),
  " roll calls."
)

#### Impute missing votes for matrix-based methods --------------------------

vote_mat_imputed <- vote_mat

for (roll_call_index in seq_len(ncol(vote_mat_imputed))) {
  roll_call_mean <- mean(
    vote_mat_imputed[, roll_call_index],
    na.rm = TRUE
  )

  vote_mat_imputed[
    is.na(vote_mat_imputed[, roll_call_index]),
    roll_call_index
  ] <- roll_call_mean
}

vote_factor_data <- as.data.frame(
  lapply(
    as.data.frame(vote_mat),
    function(vote_column) {
      factor(
        ifelse(
          is.na(vote_column),
          NA,
          ifelse(vote_column == 1, "Yea", "Nay")
        ),
        levels = c("Nay", "Yea")
      )
    }
  )
)

rownames(vote_factor_data) <- member_ids

#### Construct the voting-disagreement matrix --------------------------------

number_of_members <- nrow(vote_mat)

disagreement_matrix <- matrix(
  0,
  nrow = number_of_members,
  ncol = number_of_members,
  dimnames = list(member_ids, member_ids)
)

for (member_i in seq_len(number_of_members - 1)) {
  for (member_k in (member_i + 1):number_of_members) {
    jointly_observed <-
      !is.na(vote_mat[member_i, ]) &
      !is.na(vote_mat[member_k, ])

    if (sum(jointly_observed) >= 5) {
      agreement_rate <- mean(
        vote_mat[member_i, jointly_observed] ==
          vote_mat[member_k, jointly_observed]
      )

      disagreement_matrix[member_i, member_k] <- 1 - agreement_rate
      disagreement_matrix[member_k, member_i] <- 1 - agreement_rate
    } else {
      disagreement_matrix[member_i, member_k] <- NA_real_
      disagreement_matrix[member_k, member_i] <- NA_real_
    }
  }
}

median_disagreement <- median(
  disagreement_matrix[upper.tri(disagreement_matrix)],
  na.rm = TRUE
)

disagreement_matrix[
  is.na(disagreement_matrix)
] <- median_disagreement

diag(disagreement_matrix) <- 0

#### Euclidean, circular, and spherical MDS ---------------------------------

message("Fitting Euclidean MDS ...")

# Rescale the observed disagreement matrix to the angular interval [0, pi].
angular_target <- pi * disagreement_matrix /
  max(disagreement_matrix)

#### Euclidean MDS -----------------------------------------------------------

euclidean_mds <- cmdscale(
  as.dist(disagreement_matrix),
  k = 2,
  eig = TRUE
)

euclidean_coordinates <- euclidean_mds$points

euclidean_fitted_distance <- as.matrix(
  dist(euclidean_coordinates)
)

#### One-dimensional Euclidean MDS -----------------------------------------

euclidean_mds_1d <- cmdscale(as.dist(disagreement_matrix), k = 1, eig = TRUE)
euclidean_1d_coordinates <- matrix(euclidean_mds_1d$points[, 1], ncol = 1, dimnames = list(member_ids, "Dimension_1"))

#### Metric and ordinal MDS -------------------------------------------------

metric_mds_fit <- tryCatch(
  smacof::smacofSym(as.dist(disagreement_matrix), ndim = 2, type = "ratio"),
  error = function(error) error
)

metric_mds_coordinates <- if (inherits(metric_mds_fit, "error")) NULL else {
  metric_mds_fit$conf[, 1:2, drop = FALSE]
}

ordinal_mds_fit <- tryCatch(
  smacof::smacofSym(as.dist(disagreement_matrix), ndim = 2, type = "ordinal"),
  error = function(error) error
)

ordinal_mds_coordinates <- if (inherits(ordinal_mds_fit, "error")) NULL else {
  ordinal_mds_fit$conf[, 1:2, drop = FALSE]
}

#### Circular MDS ------------------------------------------------------------

message("Fitting circular MDS ...")

# Initialise circular angles from the Euclidean MDS configuration.
initial_angles <- atan2(
  euclidean_coordinates[, 2],
  euclidean_coordinates[, 1]
) %% (2 * pi)

# Fix the first angle at zero to remove arbitrary rotational freedom.
circular_stress_function <- function(free_angles) {
  theta <- c(0, free_angles %% (2 * pi))

  angular_difference <- abs(
    outer(theta, theta, "-")
  )

  circular_distance <- pmin(
    angular_difference,
    2 * pi - angular_difference
  )

  sum(
    (
      circular_distance[upper.tri(circular_distance)] -
        angular_target[upper.tri(angular_target)]
    ) ^ 2
  )
}

# Analytic gradient avoids one full distance-matrix calculation per angle.
circular_stress_gradient <- function(free_angles) {
  theta <- c(0, free_angles %% (2 * pi))
  angular_difference <- outer(theta, theta, "-")
  absolute_difference <- abs(angular_difference)
  circular_distance <- pmin(
    absolute_difference,
    2 * pi - absolute_difference
  )

  distance_derivative <- sign(angular_difference)
  distance_derivative[absolute_difference > pi] <-
    -distance_derivative[absolute_difference > pi]
  distance_derivative[absolute_difference == 0] <- 0
  distance_derivative[lower.tri(distance_derivative)] <-
    -t(distance_derivative)[lower.tri(distance_derivative)]
  diag(distance_derivative) <- 0

  residual <- circular_distance - angular_target
  gradient_all <- rowSums(2 * residual * distance_derivative)

  gradient_all[-1]
}

circular_mds_fit <- optim(
  par = initial_angles[-1],
  fn = circular_stress_function,
  gr = circular_stress_gradient,
  method = "BFGS",
  control = list(maxit = 800, reltol = 1e-7)
)

if (circular_mds_fit$convergence != 0) {
  warning("Circular MDS did not report convergence; inspect circular_mds_fit.")
}

circular_angles <- c(
  0,
  circular_mds_fit$par %% (2 * pi)
)

circular_coordinates <- cbind(
  Dimension_1 = cos(circular_angles),
  Dimension_2 = sin(circular_angles)
)

rownames(circular_coordinates) <- member_ids

circular_angular_difference <- abs(
  outer(circular_angles, circular_angles, "-")
)

circular_fitted_distance <- pmin(
  circular_angular_difference,
  2 * pi - circular_angular_difference
)

#### Soft circular MDS -------------------------------------------------------

soft_circular_lambda <- 5
soft_circular_start <- pmin(3, pmax(0.2, abs(euclidean_coordinates[, 1]) + abs(euclidean_coordinates[, 2])))
soft_circular_dyads <- which(upper.tri(angular_target), arr.ind = TRUE)
soft_circular_observed <- angular_target[cbind(soft_circular_dyads[, 1], soft_circular_dyads[, 2])]
soft_circular_cosine <- cos(circular_angles[soft_circular_dyads[, 1]] - circular_angles[soft_circular_dyads[, 2]])
soft_circular_loss <- function(log_radii) {
  radii <- exp(log_radii)
  radius_i <- radii[soft_circular_dyads[, 1]]
  radius_k <- radii[soft_circular_dyads[, 2]]
  fitted_distances <- sqrt(pmax(radius_i ^ 2 + radius_k ^ 2 - 2 * radius_i * radius_k * soft_circular_cosine, 1e-12))
  mean((fitted_distances - soft_circular_observed) ^ 2) + soft_circular_lambda * (sd(radii) / mean(radii)) ^ 2
}
soft_circular_fit <- tryCatch(optim(log(soft_circular_start), soft_circular_loss, method = "L-BFGS-B", lower = rep(log(0.2), length(soft_circular_start)), upper = rep(log(3), length(soft_circular_start)), control = list(maxit = 150, factr = 1e8, pgtol = 1e-4)), error = function(error) error)
soft_circular_coordinates <- if (inherits(soft_circular_fit, "error") || soft_circular_fit$convergence != 0) NULL else {
  soft_circular_radii <- exp(soft_circular_fit$par)
  cbind(Dimension_1 = soft_circular_radii * cos(circular_angles), Dimension_2 = soft_circular_radii * sin(circular_angles))
}
if (!is.null(soft_circular_coordinates)) rownames(soft_circular_coordinates) <- member_ids

#### Spherical MDS -----------------------------------------------------------

message("Fitting spherical MDS with smacofSphere ...")

spherical_mds_itmax <- 3000

spherical_mds_fit <- tryCatch(
  smacof::smacofSphere(
    as.dist(disagreement_matrix),
    ndim = 3,
    type = "interval",
    algorithm = "dual",
    init = "torgerson",
    penalty = 100,
    itmax = spherical_mds_itmax,
    eps = 1e-3,
    verbose = FALSE
  ),
  error = function(error) error
)

spherical_coordinates_raw <- if (inherits(spherical_mds_fit, "error")) {
  NULL
} else {
  spherical_mds_fit$conf[, 1:3, drop = FALSE]
}

spherical_radii <- if (is.null(spherical_coordinates_raw)) numeric(0) else sqrt(rowSums(spherical_coordinates_raw ^ 2))
spherical_radius_coefficient_of_variation <- if (length(spherical_radii) < 2 || any(!is.finite(spherical_radii)) || mean(spherical_radii) <= 1e-12) NA_real_ else sd(spherical_radii) / mean(spherical_radii)
spherical_coordinates <- if (is.null(spherical_coordinates_raw) || any(!is.finite(spherical_radii)) || any(spherical_radii <= 1e-12)) NULL else sweep(spherical_coordinates_raw, 1, spherical_radii, FUN = "/")

if (!is.null(spherical_coordinates)) {
  rownames(spherical_coordinates) <- member_ids
}

spherical_fitted_distance <- if (is.null(spherical_coordinates)) {
  NULL
} else {
  acos(pmax(pmin(tcrossprod(spherical_coordinates), 1), -1))
}

spherical_distance_values <- if (is.null(spherical_fitted_distance)) {
  numeric(0)
} else {
  spherical_fitted_distance[upper.tri(spherical_fitted_distance)]
}

spherical_distance_spread <- if (length(spherical_distance_values) < 2 ||
                                  any(!is.finite(spherical_distance_values))) {
  NA_real_
} else {
  spherical_distance_mean <- mean(spherical_distance_values)
  sqrt(sum((spherical_distance_values - spherical_distance_mean) ^ 2) /
         (length(spherical_distance_values) - 1))
}

spherical_mds_degenerate <- !is.finite(spherical_distance_spread) ||
  spherical_distance_spread <= 1e-10 ||
  !is.finite(spherical_radius_coefficient_of_variation) ||
  spherical_radius_coefficient_of_variation > 0.01
spherical_mds_converged <- !inherits(spherical_mds_fit, "error") &&
  is.finite(spherical_mds_fit$stress) &&
  spherical_mds_fit$niter < spherical_mds_itmax &&
  !spherical_mds_degenerate

if (!spherical_mds_converged) {
  warning("Spherical SMACOF did not return a valid non-degenerate configuration; it is excluded from comparison.")
}

#### Additional full-data dimension-reduction methods ------------------------

message("Fitting PCA, CA, MCA, NLPCA, HOMALS, NLCCA, and PGA ...")

pca_full_fit <- prcomp(
  vote_mat_imputed,
  center = TRUE,
  scale. = FALSE
)

pca_full_coordinates <- pca_full_fit$x[, 1:2, drop = FALSE]

kernel_pca_full_fit <- if (requireNamespace("kernlab", quietly = TRUE)) {
  tryCatch(
    kernlab::kpca(as.matrix(vote_mat_imputed), kernel = "rbfdot",
                  kpar = list(sigma = 0.05), features = 2),
    error = function(error) error
  )
} else {
  NULL
}

kernel_pca_full_coordinates <- if (is.null(kernel_pca_full_fit) ||
                                   inherits(kernel_pca_full_fit, "error")) {
  NULL
} else {
  as.matrix(kernlab::rotated(kernel_pca_full_fit))[, 1:2, drop = FALSE]
}

ca_full_fit <- FactoMineR::CA(
  vote_mat_imputed,
  graph = FALSE
)

ca_full_coordinates <- ca_full_fit$row$coord[, 1:2, drop = FALSE]

mca_full_fit <- tryCatch(
  FactoMineR::MCA(
    vote_factor_data,
    ncp = 2,
    graph = FALSE
  ),
  error = function(error) error
)

mca_full_coordinates <- if (inherits(mca_full_fit, "error")) {
  NULL
} else {
  mca_full_fit$ind$coord[, 1:2, drop = FALSE]
}

nlpca_full_fit <- tryCatch(
  Gifi::princals(
    vote_factor_data,
    ndim = 2,
    ordinal = TRUE
  ),
  error = function(error) error
)

nlpca_full_coordinates <- if (inherits(nlpca_full_fit, "error")) {
  NULL
} else {
  nlpca_full_fit$objectscores[, 1:2, drop = FALSE]
}

homals_full_fit <- tryCatch(
  Gifi::homals(
    vote_factor_data,
    ndim = 2
  ),
  error = function(error) error
)

homals_full_coordinates <- if (inherits(homals_full_fit, "error")) {
  NULL
} else {
  homals_full_fit$objectscores[, 1:2, drop = FALSE]
}

odd_full <- vote_factor_data[, seq(1, ncol(vote_factor_data), by = 2), drop = FALSE]
even_full <- vote_factor_data[, seq(2, ncol(vote_factor_data), by = 2), drop = FALSE]

nlcca_full_fit <- tryCatch({
  odd_fit <- Gifi::princals(odd_full, ndim = 2, ordinal = TRUE)
  even_fit <- Gifi::princals(even_full, ndim = 2, ordinal = TRUE)
  odd_scores <- odd_fit$objectscores[, 1:2, drop = FALSE]
  even_scores <- even_fit$objectscores[, 1:2, drop = FALSE]
  canonical <- cancor(odd_scores, even_scores)
  list(
    coordinates = (scale(odd_scores %*% canonical$xcoef[, 1:2]) +
      scale(even_scores %*% canonical$ycoef[, 1:2])) / 2
  )
}, error = function(error) error)

nlcca_full_coordinates <- if (inherits(nlcca_full_fit, "error")) {
  NULL
} else {
  nlcca_full_fit$coordinates
}

pga_center <- colMeans(vote_mat_imputed)
pga_centered <- sweep(vote_mat_imputed, 2, pga_center, FUN = "-")
pga_norms <- sqrt(rowSums(pga_centered ^ 2))
pga_unit_profiles <- pga_centered / pmax(pga_norms, 1e-12)
pga_mean <- colMeans(pga_unit_profiles)
pga_mean <- pga_mean / sqrt(sum(pga_mean ^ 2))
pga_projection_coefficients <- as.vector(
  pga_unit_profiles %*% pga_mean
)

pga_tangent <- pga_unit_profiles -
  pga_projection_coefficients *
  matrix(
    pga_mean,
    nrow = nrow(pga_unit_profiles),
    ncol = ncol(pga_unit_profiles),
    byrow = TRUE
  )
pga_full_fit <- prcomp(pga_tangent, center = TRUE, scale. = FALSE)
pga_full_coordinates <- pga_full_fit$x[, 1:2, drop = FALSE]

#### Normalized stress -------------------------------------------------------

# This function aligns each fitted distance matrix to the scale of the
# observed disagreement matrix before calculating stress.
normalized_stress <- function(observed_distance, fitted_distance) {
  observed_values <- observed_distance[upper.tri(observed_distance)]
  fitted_values <- fitted_distance[upper.tri(fitted_distance)]

  scale_factor <- sum(observed_values * fitted_values) /
    sum(fitted_values ^ 2)

  sqrt(
    sum((observed_values - scale_factor * fitted_values) ^ 2) /
      sum(observed_values ^ 2)
  )
}

euclidean_stress <- normalized_stress(
  disagreement_matrix,
  euclidean_fitted_distance
)

euclidean_1d_stress <- normalized_stress(disagreement_matrix, as.matrix(dist(euclidean_1d_coordinates)))
soft_circular_fitted_distance <- if (is.null(soft_circular_coordinates)) NULL else as.matrix(dist(soft_circular_coordinates))
soft_circular_stress <- if (is.null(soft_circular_fitted_distance)) NA_real_ else normalized_stress(disagreement_matrix, soft_circular_fitted_distance)

circular_stress <- normalized_stress(
  disagreement_matrix,
  circular_fitted_distance
)

spherical_stress <- if (is.null(spherical_fitted_distance) || spherical_mds_degenerate) {
  NA_real_
} else {
  normalized_stress(disagreement_matrix, spherical_fitted_distance)
}

geometry_stress_results <- data.frame(
  method = c("1D Euclidean MDS", "Classical MDS", "Soft Circular MDS", "Circular MDS", "Spherical MDS"),
  geometry = c("Euclidean", "Euclidean", "Soft circular", "Circular", "Spherical"),
  normalized_stress = c(
    euclidean_1d_stress,
    euclidean_stress,
    soft_circular_stress,
    circular_stress,
    spherical_stress
  )
)

print(geometry_stress_results)

#### Circularity index and projection loss ----------------------------------

euclidean_centre <- colMeans(euclidean_coordinates)

euclidean_centered <- sweep(
  euclidean_coordinates,
  2,
  euclidean_centre,
  FUN = "-"
)

euclidean_radius <- sqrt(
  rowSums(euclidean_centered ^ 2)
)

circularity_index <- sd(euclidean_radius) /
  mean(euclidean_radius)

euclidean_projection <- euclidean_centered /
  euclidean_radius

projection_fitted_distance <- as.matrix(
  dist(euclidean_projection)
)

projection_stress <- normalized_stress(
  disagreement_matrix,
  projection_fitted_distance
)

circularity_results <- data.frame(
  circularity_index = circularity_index,
  euclidean_stress = euclidean_stress,
  projected_circle_stress = projection_stress,
  projection_stress_increase = projection_stress - euclidean_stress
)

print(circularity_results)

#### Geometry comparison plot ------------------------------------------------

geometry_plot <- ggplot(
  geometry_stress_results,
  aes(
    x = normalized_stress,
    y = reorder(method, normalized_stress),
    colour = geometry
  )
) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "Euclidean, circular, and spherical MDS comparison",
    subtitle = "Lower normalized stress indicates better distance fit",
    x = "Normalized stress",
    y = NULL,
    colour = "Geometry"
  ) +
  theme_minimal()

ggsave(
  file.path(results_directory, "geometry_comparison.png"),
  geometry_plot,
  width = 9,
  height = 4,
  dpi = 220
)

#### Latent-space models on the full data -----------------------------------

message("Fitting full-data latent-space models ...")

agreement_full <- 1 - disagreement_matrix
agreement_threshold_full <- quantile(
  agreement_full[upper.tri(agreement_full)],
  probs = 0.75,
  na.rm = TRUE
)

covoting_adjacency_full <- ifelse(
  agreement_full >= agreement_threshold_full,
  1,
  0
)

diag(covoting_adjacency_full) <- 0

if (!identical(covoting_adjacency_full, t(covoting_adjacency_full)) ||
    any(diag(covoting_adjacency_full) != 0)) {
  stop("The derived co-voting network must be symmetric and have no self-ties.")
}

# This object records the common cleaned source data and the two derived
# representations. It is saved before model fitting so it can be reused.
voteview_analysis_data <- list(
  vote_matrix = vote_mat,
  observed_mask = !is.na(vote_mat),
  member_metadata = members_full,
  rollcall_metadata = rollcalls_119,
  member_distance = disagreement_matrix,
  covoting_network = covoting_adjacency_full
)

full_lsm_results <- list()

if (requireNamespace("network", quietly = TRUE) &&
    requireNamespace("latentnet", quietly = TRUE)) {
  message("Constructing the full co-voting network for latentnet ...")

  covoting_network_full <- network::network(
    covoting_adjacency_full,
    directed = FALSE,
    matrix.type = "adjacency"
  )

  message("Network constructed. Fitting latentnet extensions with screening MCMC ...")

  latentnet_screening_control <- latentnet::control.ergmm(
    sample.size = 500,
    burnin = 1000,
    interval = 10,
    pilot.runs = 1,
    mle.maxit = 50,
    refine.user.start = FALSE
  )

  # Start the MCMC from the already available Euclidean MDS positions.  This
  # avoids an expensive conditional-posterior-mode search on the full network.
  latentnet_start <- list(Z = euclidean_coordinates)

  # The latentnet model terms are evaluated inside a formula.  Give that
  # formula an environment whose parent is the latentnet namespace, matching
  # the unqualified euclidean() and bilinear() usage in voteview_graph.R.
  latentnet_formula_environment <- new.env(parent = asNamespace("latentnet"))
  latentnet_formula_environment$covoting_network_full <- covoting_network_full
  euclidean_lsm_formula <- stats::as.formula(
    "covoting_network_full ~ euclidean(d = 2)",
    env = latentnet_formula_environment
  )
  bilinear_lsm_formula <- stats::as.formula(
    "covoting_network_full ~ bilinear(d = 2)",
    env = latentnet_formula_environment
  )

  message("Fitting latentnet Euclidean LSM ...")

  euclidean_lsm_full <- tryCatch(
    latentnet::ergmm(
      euclidean_lsm_formula,
      family = "Bernoulli",
      control = latentnet_screening_control,
      tofit = c("mcmc", "procrustes"),
      user.start = latentnet_start,
      seed = 123,
      verbose = TRUE
    ),
    error = function(error) error
  )

  message("Fitting latentnet Bilinear LSM ...")

  bilinear_lsm_full <- tryCatch(
    latentnet::ergmm(
      bilinear_lsm_formula,
      family = "Bernoulli",
      control = latentnet_screening_control,
      tofit = c("mcmc", "procrustes"),
      user.start = latentnet_start,
      seed = 123,
      verbose = TRUE
    ),
    error = function(error) error
  )

  full_lsm_results$euclidean_lsm <- euclidean_lsm_full
  full_lsm_results$bilinear_lsm <- bilinear_lsm_full

  if (inherits(euclidean_lsm_full, "error")) {
    message("Latentnet Euclidean LSM failed: ", conditionMessage(euclidean_lsm_full))
  }

  if (inherits(bilinear_lsm_full, "error")) {
    message("Latentnet Bilinear LSM failed: ", conditionMessage(bilinear_lsm_full))
  }
} else {
  message("Packages 'network' and/or 'latentnet' are unavailable; LSM fits are marked failed.")
  full_lsm_results$euclidean_lsm <- NULL
  full_lsm_results$bilinear_lsm <- NULL
}

extract_latentnet_coordinates <- function(fit, member_ids) {
  if (is.null(fit) || inherits(fit, "error")) {
    return(NULL)
  }

  coordinates <- tryCatch(
    apply(fit$sample$Z, c(2, 3), mean),
    error = function(error) NULL
  )

  if (is.null(coordinates) || nrow(coordinates) != length(member_ids)) {
    return(NULL)
  }

  rownames(coordinates) <- member_ids
  coordinates
}

euclidean_lsm_full_coordinates <- extract_latentnet_coordinates(
  full_lsm_results$euclidean_lsm,
  member_ids
)

bilinear_lsm_full_coordinates <- extract_latentnet_coordinates(
  full_lsm_results$bilinear_lsm,
  member_ids
)

#### Direct Euclidean, circular, and spherical LSMs -------------------------

# The first position is fixed at the origin and the remaining positions have
# fixed Frobenius norm. These constraints remove translation and scale
# non-identifiability while retaining an estimable distance coefficient.
full_euclidean_lsm_coordinates_from_parameters <- function(parameters) {
  raw_coordinates <- matrix(parameters[3:length(parameters)], ncol = 2)
  raw_norm <- sqrt(sum(raw_coordinates ^ 2))

  if (!is.finite(raw_norm) || raw_norm < 1e-10) {
    return(NULL)
  }

  scaled_coordinates <- sqrt(number_of_members) * raw_coordinates / raw_norm
  rbind(c(Dimension_1 = 0, Dimension_2 = 0), scaled_coordinates)
}

full_euclidean_lsm_nll <- function(parameters) {
  coordinates <- full_euclidean_lsm_coordinates_from_parameters(parameters)

  if (is.null(coordinates)) {
    return(.Machine$double.xmax / 100)
  }

  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  distance_matrix <- as.matrix(dist(coordinates))
  eta <- alpha - lambda * distance_matrix
  probability <- plogis(eta)
  use <- upper.tri(covoting_adjacency_full)

  -sum(
    covoting_adjacency_full[use] * log(probability[use] + 1e-12) +
      (1 - covoting_adjacency_full[use]) * log(1 - probability[use] + 1e-12)
  )
}

full_euclidean_lsm_gradient <- function(parameters) {
  raw_coordinates <- matrix(parameters[3:length(parameters)], ncol = 2)
  raw_norm <- sqrt(sum(raw_coordinates ^ 2))
  coordinates <- full_euclidean_lsm_coordinates_from_parameters(parameters)

  if (is.null(coordinates)) {
    return(rep(0, length(parameters)))
  }

  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  distance_matrix <- as.matrix(dist(coordinates))
  safe_distance <- pmax(distance_matrix, 1e-10)
  eta <- alpha - lambda * distance_matrix
  probability <- plogis(eta)
  residual <- probability - covoting_adjacency_full
  use <- upper.tri(covoting_adjacency_full)

  coordinate_gradient <- matrix(0, number_of_members, 2)
  pair_weight <- residual * (-lambda / safe_distance)
  diag(pair_weight) <- 0

  for (dimension_index in 1:2) {
    coordinate_difference <- outer(
      coordinates[, dimension_index], coordinates[, dimension_index], "-"
    )
    coordinate_gradient[, dimension_index] <- rowSums(pair_weight * coordinate_difference)
  }

  raw_gradient <- sqrt(number_of_members) / raw_norm * (
    coordinate_gradient[-1, , drop = FALSE] -
      raw_coordinates * sum(raw_coordinates * coordinate_gradient[-1, , drop = FALSE]) /
      raw_norm ^ 2
  )

  c(
    sum(residual[use]),
    sum(residual[use] * (-lambda * distance_matrix[use])),
    as.vector(raw_gradient)
  )
}

euclidean_lsm_start <- sweep(euclidean_coordinates, 2, euclidean_coordinates[1, ], FUN = "-")
euclidean_lsm_start_norm <- sqrt(sum(euclidean_lsm_start[-1, , drop = FALSE] ^ 2))

message("Fitting direct Euclidean LSM with analytic gradient ...")

full_euclidean_lsm_fit <- if (euclidean_lsm_start_norm < 1e-10) {
  structure(list(message = "Euclidean MDS start has zero spread."), class = "error")
} else {
  tryCatch(
    optim(
      par = c(
        qlogis(mean(covoting_adjacency_full[upper.tri(covoting_adjacency_full)])),
        log(1),
        as.vector(euclidean_lsm_start[-1, , drop = FALSE])
      ),
      fn = full_euclidean_lsm_nll,
      gr = full_euclidean_lsm_gradient,
      method = "L-BFGS-B",
      control = list(maxit = 1200, factr = 1e7, pgtol = 1e-6, trace = 1, REPORT = 20)
    ),
    error = function(error) error
  )
}

if (!inherits(full_euclidean_lsm_fit, "error") && full_euclidean_lsm_fit$convergence != 0) {
  warning("Direct Euclidean LSM did not report convergence; inspect full_euclidean_lsm_fit.")
}

full_euclidean_lsm_coordinates <- if (inherits(full_euclidean_lsm_fit, "error") ||
                                      full_euclidean_lsm_fit$convergence != 0) {
  NULL
} else {
  full_euclidean_lsm_coordinates_from_parameters(full_euclidean_lsm_fit$par)
}

#### One-dimensional Euclidean LSM -----------------------------------------

full_euclidean_1d_lsm_coordinates_from_parameters <- function(parameters) {
  raw_coordinates <- parameters[3:length(parameters)]
  raw_norm <- sqrt(sum(raw_coordinates ^ 2))
  if (!is.finite(raw_norm) || raw_norm < 1e-10) return(NULL)
  rbind(0, sqrt(number_of_members) * raw_coordinates / raw_norm)
}
full_euclidean_1d_lsm_nll <- function(parameters) {
  coordinates <- full_euclidean_1d_lsm_coordinates_from_parameters(parameters)
  if (is.null(coordinates)) return(1e12)
  alpha <- parameters[1]; lambda <- exp(parameters[2]); distance_matrix <- as.matrix(dist(coordinates)); probability <- plogis(alpha - lambda * distance_matrix); use <- upper.tri(covoting_adjacency_full)
  -sum(covoting_adjacency_full[use] * log(probability[use] + 1e-12) + (1 - covoting_adjacency_full[use]) * log(1 - probability[use] + 1e-12))
}
full_euclidean_1d_lsm_gradient <- function(parameters) {
  raw_coordinates <- parameters[3:length(parameters)]; raw_norm <- sqrt(sum(raw_coordinates ^ 2)); coordinates <- full_euclidean_1d_lsm_coordinates_from_parameters(parameters); if (is.null(coordinates)) return(rep(0, length(parameters))); alpha <- parameters[1]; lambda <- exp(parameters[2]); distance_matrix <- as.matrix(dist(coordinates)); safe_distance <- pmax(distance_matrix, 1e-10); probability <- plogis(alpha - lambda * distance_matrix); residual <- probability - covoting_adjacency_full; use <- upper.tri(covoting_adjacency_full); pair_weight <- residual * (-lambda / safe_distance); diag(pair_weight) <- 0; coordinate_difference <- outer(coordinates[, 1], coordinates[, 1], "-"); coordinate_gradient <- rowSums(pair_weight * coordinate_difference); raw_gradient <- sqrt(number_of_members) / raw_norm * (coordinate_gradient[-1] - raw_coordinates * sum(raw_coordinates * coordinate_gradient[-1]) / raw_norm ^ 2); c(sum(residual[use]), sum(residual[use] * (-lambda * distance_matrix[use])), as.vector(raw_gradient))
}
euclidean_1d_lsm_start <- euclidean_1d_coordinates - euclidean_1d_coordinates[1, 1]
full_euclidean_1d_lsm_fit <- tryCatch(optim(c(qlogis(mean(covoting_adjacency_full[upper.tri(covoting_adjacency_full)])), log(1), as.vector(euclidean_1d_lsm_start[-1, 1])), full_euclidean_1d_lsm_nll, gr = full_euclidean_1d_lsm_gradient, method = "L-BFGS-B", lower = c(-8, -5, rep(-10, number_of_members - 1)), upper = c(8, 5, rep(10, number_of_members - 1)), control = list(maxit = 3000, factr = 1e8, pgtol = 1e-5, trace = 1, REPORT = 100)), error = function(error) error)
full_euclidean_1d_lsm_coordinates <- if (inherits(full_euclidean_1d_lsm_fit, "error") || full_euclidean_1d_lsm_fit$convergence != 0) NULL else full_euclidean_1d_lsm_coordinates_from_parameters(full_euclidean_1d_lsm_fit$par)
if (!is.null(full_euclidean_1d_lsm_coordinates)) rownames(full_euclidean_1d_lsm_coordinates) <- member_ids

full_circular_lsm_nll <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  theta <- c(0, parameters[3:(number_of_members + 1)])
  eta <- alpha + lambda * cos(outer(theta, theta, "-"))
  probability <- plogis(eta)
  use <- upper.tri(covoting_adjacency_full)

  -sum(
    covoting_adjacency_full[use] * log(probability[use] + 1e-12) +
      (1 - covoting_adjacency_full[use]) *
      log(1 - probability[use] + 1e-12)
  )
}

# Analytic gradient avoids BFGS numerically perturbing roughly one parameter
# per legislator at every iteration.
full_circular_lsm_gradient <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  theta <- c(0, parameters[3:(number_of_members + 1)])
  angle_difference <- outer(theta, theta, "-")
  eta <- alpha + lambda * cos(angle_difference)
  probability <- plogis(eta)
  residual <- probability - covoting_adjacency_full
  use <- upper.tri(covoting_adjacency_full)

  gradient_alpha <- sum(residual[use])
  gradient_log_lambda <- sum(residual[use] * lambda * cos(angle_difference[use]))

  angle_gradient <- residual * lambda * sin(angle_difference)
  angle_gradient[!use] <- 0
  gradient_theta <- -rowSums(angle_gradient) + colSums(angle_gradient)

  c(gradient_alpha, gradient_log_lambda, gradient_theta[-1])
}

message("Fitting direct circular LSM with analytic gradient ...")

full_circular_lsm_fit <- tryCatch(
  optim(
    par = c(
      qlogis(mean(covoting_adjacency_full[upper.tri(covoting_adjacency_full)])),
      log(1),
      initial_angles[-1]
    ),
    fn = full_circular_lsm_nll,
    gr = full_circular_lsm_gradient,
    method = "L-BFGS-B",
    control = list(maxit = 1200, factr = 1e7, pgtol = 1e-6, trace = 1, REPORT = 20)
  ),
  error = function(error) error
)

if (!inherits(full_circular_lsm_fit, "error") && full_circular_lsm_fit$convergence != 0) {
  warning("Circular LSM did not report convergence; inspect full_circular_lsm_fit.")
}

full_circular_lsm_coordinates <- if (inherits(full_circular_lsm_fit, "error") ||
                                     full_circular_lsm_fit$convergence != 0) {
  NULL
} else {
  full_circular_angles <- c(
    0,
    full_circular_lsm_fit$par[3:(number_of_members + 1)] %% (2 * pi)
  )
  cbind(
    Dimension_1 = cos(full_circular_angles),
    Dimension_2 = sin(full_circular_angles)
  )
}

#### Circular adjacency and NOMINATE post-hoc diagnostic --------------------

if (!is.null(full_circular_lsm_coordinates) &&
    nrow(full_circular_lsm_coordinates) == nrow(members_full) &&
    all(is.finite(full_circular_lsm_coordinates))) {
  full_circular_angles <- atan2(full_circular_lsm_coordinates[, 2], full_circular_lsm_coordinates[, 1]) %% (2 * pi)
  circular_order <- order(full_circular_angles)
  next_order <- c(circular_order[-1], circular_order[1])
  circular_adjacency <- data.frame(
    position = seq_along(circular_order),
    member_id = member_ids[circular_order],
    party = members_full$party[circular_order],
    next_member_id = member_ids[next_order],
    next_party = members_full$party[next_order],
    wrap_around = ifelse(seq_along(circular_order) == length(circular_order), "yes", "no"),
    agreement_rate = agreement_full[cbind(circular_order, next_order)]
  )
  circular_adjacency$cross_party <- ifelse(circular_adjacency$party != circular_adjacency$next_party, "yes", "no")
  nominate_values <- members_full$nominate_dim1
  nominate_correlation <- if (sum(is.finite(nominate_values)) >= 3) cor(full_circular_angles, nominate_values, use = "complete.obs") else NA_real_
} else {
  message("Circular adjacency diagnostic skipped because Circular LSM coordinates are unavailable, invalid, or misaligned.")
}

#### Soft circular LSM (exploratory) ---------------------------------------

# Angles are initialized from the circular LSM; radial coordinates are then
# estimated with a penalty that keeps the configuration close to a circle.
soft_circular_lsm_lambda <- 5
soft_circular_lsm_start_radii <- pmax(sqrt(rowSums(euclidean_coordinates ^ 2)) / mean(sqrt(rowSums(euclidean_coordinates ^ 2))), 0.05)
soft_circular_lsm_dyads <- which(upper.tri(covoting_adjacency_full), arr.ind = TRUE)
soft_circular_lsm_y <- covoting_adjacency_full[cbind(soft_circular_lsm_dyads[, 1], soft_circular_lsm_dyads[, 2])]
# Use the converged direct Circular LSM angles when available. This makes the
# extension conditional on the fitted circular ordering rather than an MDS start.
soft_circular_lsm_angles <- if (is.null(full_circular_lsm_coordinates)) initial_angles else atan2(full_circular_lsm_coordinates[, 2], full_circular_lsm_coordinates[, 1])
soft_circular_lsm_angles <- as.numeric(soft_circular_lsm_angles)
if (length(soft_circular_lsm_angles) != number_of_members) stop("Soft Circular LSM angles are not aligned with the retained members.")
soft_circular_lsm_cosine <- cos(soft_circular_lsm_angles[soft_circular_lsm_dyads[, 1]] - soft_circular_lsm_angles[soft_circular_lsm_dyads[, 2]])
soft_circular_lsm_coordinates_from_parameters <- function(parameters) {
  free_log_radii <- parameters[seq.int(3, number_of_members + 1)]
  log_radii <- c(free_log_radii, -sum(free_log_radii))
  radii <- exp(log_radii)
  if (length(radii) != number_of_members || any(!is.finite(radii))) return(NULL)
  coordinates <- cbind(Dimension_1 = radii * cos(soft_circular_lsm_angles), Dimension_2 = radii * sin(soft_circular_lsm_angles))
  if (!identical(nrow(coordinates), number_of_members)) return(NULL)
  coordinates
}
soft_circular_lsm_penalized_nll <- function(parameters) {
  alpha <- parameters[1]
  distance_lambda <- exp(parameters[2])
  free_log_radii <- parameters[3:(number_of_members + 1)]
  radii <- exp(c(free_log_radii, -sum(free_log_radii)))
  if (any(!is.finite(radii))) return(.Machine$double.xmax / 100)
  radius_i <- radii[soft_circular_lsm_dyads[, 1]]
  radius_k <- radii[soft_circular_lsm_dyads[, 2]]
  fitted_distance <- sqrt(pmax(radius_i ^ 2 + radius_k ^ 2 - 2 * radius_i * radius_k * soft_circular_lsm_cosine, 1e-12))
  tie_probability <- plogis(alpha - distance_lambda * fitted_distance)
  negative_log_likelihood <- -sum(soft_circular_lsm_y * log(tie_probability + 1e-12) + (1 - soft_circular_lsm_y) * log(1 - tie_probability + 1e-12))
  radii_penalty <- soft_circular_lsm_lambda * mean((radii - mean(radii)) ^ 2)
  negative_log_likelihood + radii_penalty
}
soft_circular_lsm_gradient <- function(parameters) {
  alpha <- parameters[1]
  distance_lambda <- exp(parameters[2])
  free_log_radii <- parameters[3:(number_of_members + 1)]
  radii <- exp(c(free_log_radii, -sum(free_log_radii)))
  if (any(!is.finite(radii))) return(rep(0, length(parameters)))
  radius_i <- radii[soft_circular_lsm_dyads[, 1]]
  radius_k <- radii[soft_circular_lsm_dyads[, 2]]
  fitted_distance <- sqrt(pmax(radius_i ^ 2 + radius_k ^ 2 - 2 * radius_i * radius_k * soft_circular_lsm_cosine, 1e-12))
  tie_probability <- plogis(alpha - distance_lambda * fitted_distance)
  residual <- tie_probability - soft_circular_lsm_y
  radius_derivative_i <- distance_lambda * (radius_i - radius_k * soft_circular_lsm_cosine) / fitted_distance * radius_i
  radius_derivative_k <- distance_lambda * (radius_k - radius_i * soft_circular_lsm_cosine) / fitted_distance * radius_k
  radius_gradient <- numeric(number_of_members)
  endpoint_i <- split(radius_derivative_i * residual, soft_circular_lsm_dyads[, 1])
  endpoint_k <- split(radius_derivative_k * residual, soft_circular_lsm_dyads[, 2])
  radius_gradient[as.integer(names(endpoint_i))] <- -vapply(endpoint_i, sum, numeric(1))
  radius_gradient[as.integer(names(endpoint_k))] <- radius_gradient[as.integer(names(endpoint_k))] - vapply(endpoint_k, sum, numeric(1))
  radii_penalty_gradient <- 2 * soft_circular_lsm_lambda * (radii - mean(radii)) / number_of_members * radii
  full_log_radius_gradient <- radius_gradient + radii_penalty_gradient
  c(sum(residual), sum(-distance_lambda * fitted_distance * residual), full_log_radius_gradient[-number_of_members] - full_log_radius_gradient[number_of_members])
}
soft_circular_lsm_start_log_radii <- log(pmin(3, pmax(0.2, soft_circular_lsm_start_radii)))
soft_circular_lsm_start_log_radii <- soft_circular_lsm_start_log_radii - mean(soft_circular_lsm_start_log_radii)
soft_circular_lsm_start <- c(qlogis(pmin(0.99, pmax(0.01, mean(soft_circular_lsm_y)))), log(1), soft_circular_lsm_start_log_radii[-number_of_members])
soft_circular_lsm_fit <- tryCatch(optim(par = soft_circular_lsm_start, fn = soft_circular_lsm_penalized_nll, gr = soft_circular_lsm_gradient, method = "L-BFGS-B", lower = c(-8, -5, rep(-1, number_of_members - 1)), upper = c(8, 5, rep(1, number_of_members - 1)), control = list(maxit = 1200, factr = 1e8, pgtol = 1e-4, trace = 1, REPORT = 50)), error = function(error) error)
soft_circular_lsm_coordinates <- if (inherits(soft_circular_lsm_fit, "error") || soft_circular_lsm_fit$convergence != 0) NULL else soft_circular_lsm_coordinates_from_parameters(soft_circular_lsm_fit$par)
if (!is.null(soft_circular_lsm_coordinates) && nrow(soft_circular_lsm_coordinates) == length(member_ids)) rownames(soft_circular_lsm_coordinates) <- member_ids

full_spherical_lsm_nll <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  longitude <- parameters[3:(number_of_members + 2)]
  colatitude <- parameters[(number_of_members + 3):(2 * number_of_members + 2)]
  coordinates <- cbind(
    sin(colatitude) * cos(longitude),
    sin(colatitude) * sin(longitude),
    cos(colatitude)
  )
  eta <- alpha + lambda * tcrossprod(coordinates)
  probability <- plogis(eta)
  use <- upper.tri(covoting_adjacency_full)

  -sum(
    covoting_adjacency_full[use] * log(probability[use] + 1e-12) +
      (1 - covoting_adjacency_full[use]) *
      log(1 - probability[use] + 1e-12)
  )
}

full_spherical_lsm_gradient <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  longitude <- parameters[3:(number_of_members + 2)]
  colatitude <- parameters[(number_of_members + 3):(2 * number_of_members + 2)]
  coordinates <- cbind(
    x = sin(colatitude) * cos(longitude),
    y = sin(colatitude) * sin(longitude),
    z = cos(colatitude)
  )
  derivative_longitude <- cbind(
    x = -sin(colatitude) * sin(longitude),
    y = sin(colatitude) * cos(longitude),
    z = rep(0, number_of_members)
  )
  derivative_colatitude <- cbind(
    x = cos(colatitude) * cos(longitude),
    y = cos(colatitude) * sin(longitude),
    z = -sin(colatitude)
  )

  similarity <- tcrossprod(coordinates)
  eta <- alpha + lambda * similarity
  probability <- plogis(eta)
  residual <- probability - covoting_adjacency_full
  use <- upper.tri(covoting_adjacency_full)

  gradient_alpha <- sum(residual[use])
  gradient_log_lambda <- sum(residual[use] * lambda * similarity[use])

  pair_weight <- residual * lambda
  pair_weight[!use] <- 0
  longitude_gradient <-
    rowSums(pair_weight * (derivative_longitude %*% t(coordinates))) +
    colSums(pair_weight * (coordinates %*% t(derivative_longitude)))
  colatitude_gradient <-
    rowSums(pair_weight * (derivative_colatitude %*% t(coordinates))) +
    colSums(pair_weight * (coordinates %*% t(derivative_colatitude)))

  c(gradient_alpha, gradient_log_lambda, longitude_gradient, colatitude_gradient)
}

message("Fitting direct spherical LSM with analytic gradient ...")

full_spherical_lsm_fit <- tryCatch(
  optim(
    par = c(
      qlogis(mean(covoting_adjacency_full[upper.tri(covoting_adjacency_full)])),
      log(1),
      initial_angles,
      rep(pi / 2, number_of_members)
  ),
  fn = full_spherical_lsm_nll,
  gr = full_spherical_lsm_gradient,
  method = "L-BFGS-B",
  control = list(maxit = 1200, factr = 1e7, pgtol = 1e-6, trace = 1, REPORT = 20)
  ),
  error = function(error) error
)

if (!inherits(full_spherical_lsm_fit, "error") && full_spherical_lsm_fit$convergence != 0) {
  warning("Spherical LSM did not report convergence; inspect full_spherical_lsm_fit.")
}

full_spherical_lsm_coordinates <- if (inherits(full_spherical_lsm_fit, "error") ||
                                      full_spherical_lsm_fit$convergence != 0) {
  NULL
} else {
  full_longitude <- full_spherical_lsm_fit$par[3:(number_of_members + 2)]
  full_colatitude <- full_spherical_lsm_fit$par[(number_of_members + 3):(2 * number_of_members + 2)]
  cbind(
    Dimension_1 = sin(full_colatitude) * cos(full_longitude),
    Dimension_2 = sin(full_colatitude) * sin(full_longitude),
    Dimension_3 = cos(full_colatitude)
  )
}

#### Method-level comparison table ------------------------------------------
distance_from_coordinates <- function(coordinates, geometry = "Euclidean") {
  if (is.null(coordinates)) {
    return(NULL)
  }

  coordinates <- as.matrix(coordinates)

  if (geometry == "Spherical") {
    if (ncol(coordinates) != 3 || any(!is.finite(coordinates)) ||
        max(abs(rowSums(coordinates ^ 2) - 1)) > 1e-6) {
      return(NULL)
    }
    cosine_similarity <- tcrossprod(coordinates)
    return(acos(pmax(pmin(cosine_similarity, 1), -1)))
  }

  if (ncol(coordinates) == 1) {
    return(as.matrix(dist(coordinates[, 1, drop = FALSE])))
  }

  as.matrix(dist(coordinates[, 1:2, drop = FALSE]))
}

method_coordinates <- list(
  PCA = pca_full_coordinates,
  `Kernel PCA` = kernel_pca_full_coordinates,
  NLPCA = nlpca_full_coordinates,
  `1D Euclidean MDS` = euclidean_1d_coordinates,
  `Classical MDS` = euclidean_coordinates,
  `Metric MDS` = metric_mds_coordinates,
  `Ordinal MDS` = ordinal_mds_coordinates,
  `Circular MDS` = circular_coordinates,
  `Soft Circular MDS` = soft_circular_coordinates,
  `Spherical MDS` = spherical_coordinates,
  CA = ca_full_coordinates,
  MCA = mca_full_coordinates,
  HOMALS = homals_full_coordinates,
  NLCCA = nlcca_full_coordinates,
  PGA = pga_full_coordinates,
  `Euclidean LSM` = full_euclidean_lsm_coordinates,
  `1D Euclidean LSM` = full_euclidean_1d_lsm_coordinates,
  `Latentnet Euclidean LSM` = euclidean_lsm_full_coordinates,
  `Bilinear LSM` = bilinear_lsm_full_coordinates,
  `Circular LSM` = full_circular_lsm_coordinates,
  `Soft Circular LSM` = soft_circular_lsm_coordinates,
  `Spherical LSM` = full_spherical_lsm_coordinates
)

method_metadata <- data.frame(
  method = names(method_coordinates),
  family = c(
    rep("Dimension reduction", 15), rep("LSM", 7)
  ),
  within_family = c(
    "PCA family", "PCA family", "PCA family",
    "MDS family", "MDS family", "MDS family", "MDS family", "MDS family", "MDS family", "MDS family",
    "Correspondence / optimal scaling", "Correspondence / optimal scaling",
    "Correspondence / optimal scaling", "Correspondence / optimal scaling",
    "PGA", "LSM family", "LSM family", "Supplementary implementation", "LSM family", "LSM family", "LSM family", "LSM family"
  ),
  geometry = c(
    "Euclidean", "Euclidean", "Euclidean",
    "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Circular", "Soft circular", "Spherical",
    "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Tangent-space",
    "Euclidean", "Euclidean", "Euclidean", "Euclidean", "Circular", "Soft circular", "Spherical"
  ),
  include_geometry_comparison = c(
    FALSE, FALSE, FALSE, FALSE,
    FALSE, TRUE, FALSE, TRUE, TRUE, TRUE,
    FALSE, FALSE, FALSE, FALSE, FALSE,
    TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE
  ),
  comparison_role = c(
    rep("core", 17), "supplementary", rep("core", 4)
  ),
  stringsAsFactors = FALSE
)

parsimony_summary <- method_metadata |>
  mutate(
    intrinsic_dimension = case_when(
      method %in% c("1D Euclidean MDS", "Circular MDS", "1D Euclidean LSM", "Circular LSM") ~ 1,
      method %in% c("Spherical MDS", "Spherical LSM") ~ 2,
      TRUE ~ 2
    ),
    ambient_dimension = case_when(method %in% c("Spherical MDS", "Spherical LSM") ~ 3, TRUE ~ intrinsic_dimension),
    representation = case_when(grepl("Circular", method, ignore.case = TRUE) ~ "angular/radial", grepl("Spherical", method, ignore.case = TRUE) ~ "unit sphere", TRUE ~ "Euclidean or tangent-space")
  ) |>
  select(method, family, comparison_role, geometry, intrinsic_dimension, ambient_dimension, representation)

pca_family_methods <- c("PCA", "Kernel PCA", "NLPCA")
correspondence_family_methods <- c("CA", "MCA", "HOMALS", "NLCCA")
mds_family_methods <- c("1D Euclidean MDS", "Classical MDS", "Metric MDS", "Ordinal MDS", "Circular MDS", "Soft Circular MDS", "Spherical MDS")
lsm_family_methods <- c("1D Euclidean LSM", "Euclidean LSM", "Bilinear LSM", "Circular LSM", "Soft Circular LSM", "Spherical LSM")
core_dimension_reduction_methods <- method_metadata |>
  filter(family == "Dimension reduction", comparison_role == "core") |>
  pull(method)
core_lsm_methods <- method_metadata |>
  filter(family == "LSM", comparison_role == "core") |>
  pull(method)
core_framework_methods <- c(core_dimension_reduction_methods, core_lsm_methods)

method_metrics <- method_metadata
method_metrics$status <- "failed"
method_metrics$converged <- FALSE
method_metrics$normalized_stress <- NA_real_
method_metrics$native_information_criterion <- NA_real_
method_metrics$stress_status <- "not evaluated"

fit_converged <- c(
  PCA = TRUE,
  `Kernel PCA` = !is.null(kernel_pca_full_coordinates),
  NLPCA = !inherits(nlpca_full_fit, "error"),
  `1D Euclidean MDS` = TRUE,
  `Classical MDS` = TRUE,
  `Metric MDS` = !inherits(metric_mds_fit, "error"),
  `Ordinal MDS` = !inherits(ordinal_mds_fit, "error"),
  `Circular MDS` = circular_mds_fit$convergence == 0,
  `Soft Circular MDS` = !is.null(soft_circular_coordinates),
  `Spherical MDS` = spherical_mds_converged,
  CA = TRUE,
  MCA = !inherits(mca_full_fit, "error"),
  HOMALS = !inherits(homals_full_fit, "error"),
  NLCCA = !inherits(nlcca_full_fit, "error"),
  PGA = TRUE,
  `Euclidean LSM` = !inherits(full_euclidean_lsm_fit, "error") && full_euclidean_lsm_fit$convergence == 0,
  `1D Euclidean LSM` = !is.null(full_euclidean_1d_lsm_coordinates),
  `Latentnet Euclidean LSM` = !is.null(euclidean_lsm_full_coordinates),
  `Bilinear LSM` = !is.null(bilinear_lsm_full_coordinates),
  `Circular LSM` = !inherits(full_circular_lsm_fit, "error") && full_circular_lsm_fit$convergence == 0,
  `Soft Circular LSM` = !is.null(soft_circular_lsm_coordinates),
  `Spherical LSM` = !inherits(full_spherical_lsm_fit, "error") && full_spherical_lsm_fit$convergence == 0
)

for (method_index in seq_len(nrow(method_metrics))) {
  method_name <- method_metrics$method[method_index]
  coordinates <- method_coordinates[[method_name]]

  if (!is.null(coordinates) && isTRUE(fit_converged[method_name])) {
    fitted_distance <- distance_from_coordinates(
      coordinates,
      method_metrics$geometry[method_index]
    )

    fitted_upper_triangle <- if (is.null(fitted_distance)) numeric(0) else {
      fitted_distance[upper.tri(fitted_distance)]
    }

    if (isTRUE(length(fitted_upper_triangle) > 1) &&
        isTRUE(all(is.finite(fitted_upper_triangle))) &&
        isTRUE(sd(fitted_upper_triangle) > 1e-10)) {
      calculated_stress <- normalized_stress(disagreement_matrix, fitted_distance)
      method_metrics$normalized_stress[method_index] <- calculated_stress
      method_metrics$stress_status[method_index] <- if (calculated_stress < 1e-8) "requires manual validation" else "ok"
    }
    method_metrics$status[method_index] <- "ok"
    method_metrics$converged[method_index] <- TRUE
  }
}

if (!inherits(full_euclidean_lsm_fit, "error") && full_euclidean_lsm_fit$convergence == 0) {
  method_metrics$native_information_criterion[
    method_metrics$method == "Euclidean LSM"
  ] <- 2 * full_euclidean_lsm_fit$value +
    (2 * number_of_members - 1) * log(sum(upper.tri(covoting_adjacency_full)))
}

if (!inherits(full_euclidean_1d_lsm_fit, "error") && full_euclidean_1d_lsm_fit$convergence == 0) {
  method_metrics$native_information_criterion[method_metrics$method == "1D Euclidean LSM"] <-
    2 * full_euclidean_1d_lsm_fit$value + (number_of_members + 1) * log(sum(upper.tri(covoting_adjacency_full)))
}

if (!inherits(full_circular_lsm_fit, "error") && full_circular_lsm_fit$convergence == 0) {
  method_metrics$native_information_criterion[
    method_metrics$method == "Circular LSM"
  ] <- 2 * full_circular_lsm_fit$value +
    (number_of_members + 1) * log(sum(upper.tri(covoting_adjacency_full)))
}

if (!inherits(full_spherical_lsm_fit, "error") && full_spherical_lsm_fit$convergence == 0) {
  method_metrics$native_information_criterion[
    method_metrics$method == "Spherical LSM"
  ] <- 2 * full_spherical_lsm_fit$value +
    (2 * number_of_members - 1) * log(sum(upper.tri(covoting_adjacency_full)))
}

#### Held-out decoder comparison --------------------------------------------
set.seed(123)
observed_entries <- which(!is.na(vote_mat), arr.ind = TRUE)
entry_fold <- sample(rep(1:5, length.out = nrow(observed_entries)))

evaluate_coordinate_decoder <- function(coordinates) {
  if (is.null(coordinates)) {
    return(rep(NA_real_, 5))
  }

  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) != nrow(vote_mat) || ncol(coordinates) < 1) {
    return(rep(NA_real_, 5))
  }
  coordinates <- coordinates[, seq_len(min(2, ncol(coordinates))), drop = FALSE]

  if (any(!is.finite(coordinates))) {
    return(rep(NA_real_, 5))
  }

  fold_losses <- numeric(5)

  for (fold_index in 1:5) {
    message("  Fold ", fold_index, " of 5 ...")
    held_out_matrix <- matrix(FALSE, nrow(vote_mat), ncol(vote_mat))
    held_out_rows <- observed_entries[entry_fold == fold_index, , drop = FALSE]
    held_out_matrix[held_out_rows] <- TRUE
    loss_values <- numeric(0)

    for (roll_call_index in seq_len(ncol(vote_mat))) {
      test_members <- which(held_out_matrix[, roll_call_index])

      if (length(test_members) == 0) {
        next
      }

      training_members <- which(
        !is.na(vote_mat[, roll_call_index]) &
          !held_out_matrix[, roll_call_index]
      )
      training_votes <- vote_mat[training_members, roll_call_index]
      fallback_probability <- mean(training_votes, na.rm = TRUE)

      if (!is.finite(fallback_probability)) {
        fallback_probability <- mean(vote_mat[, roll_call_index], na.rm = TRUE)
      }

      if (length(training_members) < 10 || length(unique(training_votes)) < 2) {
        probabilities <- rep(fallback_probability, length(test_members))
      } else {
        training_design <- cbind(intercept = 1, coordinates[training_members, , drop = FALSE])
        decoder <- tryCatch(
          suppressWarnings(glm.fit(training_design, training_votes, family = binomial())),
          error = function(error) NULL
        )

        probabilities <- if (is.null(decoder) || any(!is.finite(decoder$coefficients))) {
          rep(fallback_probability, length(test_members))
        } else {
          plogis(cbind(intercept = 1, coordinates[test_members, , drop = FALSE]) %*% decoder$coefficients)
        }
      }

      probabilities[!is.finite(probabilities)] <- fallback_probability
      probabilities <- pmin(pmax(as.numeric(probabilities), 1e-6), 1 - 1e-6)
      observed_votes <- vote_mat[test_members, roll_call_index]
      loss_values <- c(
        loss_values,
        -(observed_votes * log(probabilities) +
          (1 - observed_votes) * log(1 - probabilities))
      )
    }

    fold_losses[fold_index] <- if (length(loss_values) == 0 || any(!is.finite(loss_values))) {
      NA_real_
    } else {
      mean(loss_values)
    }
  }

  fold_losses
}

cross_validation_results <- do.call(
  rbind,
  lapply(
    names(method_coordinates),
    function(method_name) {
      message("Evaluating held-out decoder for ", method_name, " ...")
      coordinates_for_decoder <- if (isTRUE(fit_converged[method_name])) {
        method_coordinates[[method_name]]
      } else {
        NULL
      }
      fold_losses <- evaluate_coordinate_decoder(coordinates_for_decoder)
      cv_status <- if (all(is.finite(fold_losses))) "ok" else "failed"
      data.frame(
        method = method_name,
        fold = 1:5,
        held_out_log_loss = fold_losses,
        cv_status = cv_status,
        evaluation_scope = "held-out decoder; coordinates fitted once",
        stringsAsFactors = FALSE
      )
    }
  )
)

#### Geometry and framework winner tables ----------------------------------
predictive_method_names <- unique(cross_validation_results$method)

predictive_summary <- data.frame(
  method = predictive_method_names,
  number_of_folds = NA_integer_,
  number_of_finite_folds = NA_integer_,
  cv_status = NA_character_,
  held_out_log_loss = NA_real_,
  held_out_log_loss_se = NA_real_,
  stringsAsFactors = FALSE
)

for (method_index in seq_along(predictive_method_names)) {
  method_name <- predictive_method_names[method_index]
  method_results <- cross_validation_results[
    cross_validation_results$method == method_name,
    ,
    drop = FALSE
  ]
  fold_losses <- as.numeric(method_results$held_out_log_loss)
  finite_fold_losses <- fold_losses[is.finite(fold_losses)]
  number_of_finite_folds <- length(finite_fold_losses)
  successful_cv <- all(method_results$cv_status == "ok") && number_of_finite_folds == 5

  predictive_summary$number_of_folds[method_index] <- length(fold_losses)
  predictive_summary$number_of_finite_folds[method_index] <- number_of_finite_folds
  predictive_summary$cv_status[method_index] <- if (successful_cv) "ok" else "failed"

  if (successful_cv) {
    mean_fold_loss <- mean(finite_fold_losses)
    predictive_summary$held_out_log_loss[method_index] <- mean_fold_loss
    predictive_summary$held_out_log_loss_se[method_index] <- sqrt(
      sum((finite_fold_losses - mean_fold_loss) ^ 2) /
        (number_of_finite_folds * (number_of_finite_folds - 1))
    )
  }
}

invalid_predictive_rows <- predictive_summary |>
  filter(
    cv_status == "ok",
    number_of_finite_folds != 5 | !is.finite(held_out_log_loss_se)
  )

if (nrow(invalid_predictive_rows) > 0) {
  print(invalid_predictive_rows)
  stop("A successful decoder comparison must have five finite fold losses and a finite standard error.")
}

method_metrics <- method_metrics |>
  left_join(predictive_summary, by = "method")

#### Clean paper comparison table ------------------------------------------

# The predictive columns are created in the preceding decoder section. This
# guard also makes the table block safe to rerun on its own after that section
# has already been executed.
if (!"held_out_log_loss" %in% names(method_metrics)) {
  if (!exists("predictive_summary")) {
    stop("Run the held-out decoder comparison first: predictive_summary is not available.")
  }
  method_metrics <- method_metrics |>
    left_join(predictive_summary, by = "method")
}

paper_table_specification <- tibble(
  method = c("PCA", "Kernel PCA", "NLPCA", "PGA", "1D Euclidean MDS", "Classical MDS", "Metric MDS", "Ordinal MDS", "Circular MDS", "Soft Circular MDS", "Spherical MDS", "CA", "MCA", "HOMALS", "NLCCA", "Euclidean LSM", "1D Euclidean LSM", "Latentnet Euclidean LSM", "Bilinear LSM", "Circular LSM", "Soft Circular LSM", "Spherical LSM"),
  comparison_family = c(rep("PCA", 4), rep("MDS", 7), rep("CA / Optimal Scaling", 4), rep("LSM", 7)),
  Family = c("PCA", rep("", 3), "MDS", rep("", 6), "CA / Optimal Scaling", rep("", 3), "LSM", rep("", 6)),
  Method = c("Linear PCA", "Kernel PCA", "Nonlinear PCA (NLPCA)", "Principal Geodesic Analysis (PGA)", "Euclidean MDS", "Classical MDS", "Metric MDS", "Ordinal MDS", "Circular MDS", "Soft Circular MDS", "Spherical MDS", "CA", "Multiple CA", "HOMALS", "Nonlinear CCA", "Euclidean LSM", "Euclidean LSM", "Latentnet Euclidean LSM", "Bilinear LSM", "Circular LSM", "Soft Circular LSM", "Spherical LSM"),
  Dimension = c(rep("2D", 4), "1D", rep("2D", 3), "1D", rep("2D", 2), rep("2D", 4), "2D", "1D", rep("2D", 2), "1D", rep("2D", 2))
)

if (!all(paper_table_specification$method %in% method_metrics$method)) stop("The paper-table specification contains a method missing from method_metrics.")

method_comparison_table <- paper_table_specification |>
  left_join(method_metrics |> select(method, held_out_log_loss, normalized_stress, native_information_criterion), by = "method") |>
  group_by(comparison_family) |>
  mutate(
    held_out_best = if (any(is.finite(held_out_log_loss))) is.finite(held_out_log_loss) & held_out_log_loss == min(held_out_log_loss[is.finite(held_out_log_loss)]) else rep(FALSE, n()),
    stress_best = if (any(is.finite(normalized_stress))) is.finite(normalized_stress) & normalized_stress == min(normalized_stress[is.finite(normalized_stress)]) else rep(FALSE, n()),
    bic_best = if (any(is.finite(native_information_criterion))) is.finite(native_information_criterion) & native_information_criterion == min(native_information_criterion[is.finite(native_information_criterion)]) else rep(FALSE, n()),
    `Held-out log-loss` = if_else(is.finite(held_out_log_loss), formatC(held_out_log_loss, format = "f", digits = 4), "NA"),
    `Normalized stress` = if_else(is.finite(normalized_stress), formatC(normalized_stress, format = "f", digits = 4), "NA"),
    BIC = if_else(is.finite(native_information_criterion), formatC(native_information_criterion, format = "f", digits = 4), "NA"),
    `Held-out log-loss` = if_else(held_out_best, paste0("**", `Held-out log-loss`, "**"), `Held-out log-loss`),
    `Normalized stress` = if_else(stress_best, paste0("**", `Normalized stress`, "**"), `Normalized stress`),
    BIC = if_else(bic_best, paste0("**", BIC, "**"), BIC)
  ) |>
  ungroup() |>
  select(Family, Method, Dimension, `Held-out log-loss`, `Normalized stress`, BIC)

print(method_comparison_table, n = Inf)
write.csv(method_comparison_table, file.path(results_directory, "method_comparison_table.csv"), row.names = FALSE, na = "NA")

within_family_comparison <- method_metrics |>
  mutate(
    criterion = case_when(
      within_family %in% c("PCA family", "Correspondence / optimal scaling", "LSM family") ~
        "Conditional held-out decoder log-loss",
      within_family == "MDS family" ~ "Normalized stress",
      TRUE ~ "No within-family comparator"
    ),
    criterion_value = case_when(
      within_family %in% c("PCA family", "Correspondence / optimal scaling", "LSM family") ~ held_out_log_loss,
      within_family == "MDS family" ~ normalized_stress,
      TRUE ~ NA_real_
    ),
    eligible_for_within_family_rank = case_when(
      within_family %in% c("PCA family", "Correspondence / optimal scaling", "LSM family") ~
        comparison_role == "core" & status == "ok" & converged & cv_status == "ok",
      within_family == "MDS family" ~
        comparison_role == "core" & status == "ok" & converged & stress_status == "ok",
      TRUE ~ FALSE
    )
  )

within_family_comparison$rank_within_family <- NA_integer_

for (family_name in unique(within_family_comparison$within_family)) {
  family_rows <- which(within_family_comparison$within_family == family_name &
                         within_family_comparison$eligible_for_within_family_rank)

  if (length(family_rows) > 0) {
    within_family_comparison$rank_within_family[family_rows] <- rank(
      within_family_comparison$criterion_value[family_rows],
      ties.method = "min"
    )
  }
}

within_family_comparison <- within_family_comparison |>
  mutate(
    comparison_complete = case_when(
      within_family == "PCA family" ~ all(pca_family_methods %in% method[eligible_for_within_family_rank]),
      within_family == "Correspondence / optimal scaling" ~
        all(correspondence_family_methods %in% method[eligible_for_within_family_rank]),
      within_family == "MDS family" ~
        all(mds_family_methods %in% method[eligible_for_within_family_rank]),
      within_family == "LSM family" ~
        all(lsm_family_methods %in% method[eligible_for_within_family_rank]),
      TRUE ~ FALSE
    )
  ) |>
  arrange(within_family, rank_within_family, method)

mds_geometry_comparison <- method_metrics |>
  filter(method %in% c("1D Euclidean MDS", "Circular MDS", "Metric MDS", "Soft Circular MDS", "Spherical MDS"),
         comparison_role == "core", status == "ok", converged, stress_status == "ok") |>
  mutate(
    comparison_group = if_else(method %in% c("1D Euclidean MDS", "Circular MDS"), "MDS 1D geometry", "MDS 2D geometry"),
    criterion = "Normalized stress",
    criterion_value = normalized_stress
  ) |>
  group_by(comparison_group) |>
  mutate(comparison_complete = if_else(comparison_group == "MDS 1D geometry", n() == 2, n() == 3)) |>
  ungroup() |>
  arrange(criterion_value)

lsm_geometry_comparison <- method_metrics |>
  filter(method %in% c("1D Euclidean LSM", "Circular LSM", "Euclidean LSM", "Soft Circular LSM", "Spherical LSM"),
         comparison_role == "core", status == "ok", converged, !is.na(native_information_criterion)) |>
  mutate(
    comparison_group = if_else(method %in% c("1D Euclidean LSM", "Circular LSM"), "LSM 1D geometry", "LSM 2D geometry"),
    criterion = "BIC",
    criterion_value = native_information_criterion
  ) |>
  group_by(comparison_group) |>
  mutate(comparison_complete = if_else(comparison_group == "LSM 1D geometry", n() == 2, n() == 3)) |>
  ungroup() |>
  arrange(criterion_value)

geometry_comparison_full <- bind_rows(
  mds_geometry_comparison,
  lsm_geometry_comparison
)

framework_comparison <- method_metrics |>
  filter(comparison_role == "core", status == "ok", converged, cv_status == "ok") |>
  group_by(family) |>
  mutate(rank_within_family = rank(held_out_log_loss, ties.method = "min")) |>
  ungroup() |>
  arrange(held_out_log_loss)

framework_family_winners <- framework_comparison |>
  group_by(family) |>
  slice_min(held_out_log_loss, n = 1, with_ties = FALSE) |>
  ungroup()

dimension_reduction_vs_lsm <- framework_family_winners |>
  transmute(Comparison = "Dimension reduction versus LSM", Family = family, Method = method, `Held-out log-loss` = held_out_log_loss)

one_vs_two_dimensional <- bind_rows(
  mds_geometry_comparison |>
    transmute(Comparison = comparison_group, Method = method, Geometry = geometry, `Intrinsic dimension` = if_else(comparison_group == "MDS 1D geometry", 1L, 2L), Criterion = criterion, Value = criterion_value),
  lsm_geometry_comparison |>
    transmute(Comparison = comparison_group, Method = method, Geometry = geometry, `Intrinsic dimension` = if_else(comparison_group == "LSM 1D geometry", 1L, 2L), Criterion = criterion, Value = criterion_value)
)

geometry_plot_full <- ggplot(
  mds_geometry_comparison,
  aes(
    x = normalized_stress,
    y = reorder(method, normalized_stress),
    colour = geometry
  )
) +
  geom_point(size = 3) +
  labs(
    title = "Full-data MDS geometry comparison",
    subtitle = "Lower normalized stress indicates better distance fit",
    x = "Normalized stress",
    y = NULL,
    colour = "Geometry"
  ) +
  theme_minimal()

ggsave(
  file.path(results_directory, "geometry_comparison.png"),
  geometry_plot_full,
  width = 9,
  height = 5,
  dpi = 220
)

best_mds_1d <- mds_geometry_comparison |> filter(comparison_group == "MDS 1D geometry", comparison_complete) |> slice_min(criterion_value, n = 1, with_ties = FALSE)
best_mds_2d <- mds_geometry_comparison |> filter(comparison_group == "MDS 2D geometry", comparison_complete) |> slice_min(criterion_value, n = 1, with_ties = FALSE)

best_lsm_1d <- lsm_geometry_comparison |> filter(comparison_group == "LSM 1D geometry", comparison_complete) |> slice_min(criterion_value, n = 1, with_ties = FALSE)
best_lsm_2d <- lsm_geometry_comparison |> filter(comparison_group == "LSM 2D geometry", comparison_complete) |> slice_min(criterion_value, n = 1, with_ties = FALSE)

framework_comparison_complete <- all(core_framework_methods %in% framework_comparison$method)

best_framework <- if (framework_comparison_complete) {
  framework_comparison[1, , drop = FALSE]
} else {
  data.frame(method = "Dimension-reduction versus LSM comparison incomplete")
}

winner_lines <- c(
  "Unified Members' Votes method comparison",
  "",
  paste("Best 1D MDS geometry by normalized stress:", if (nrow(best_mds_1d) > 0) best_mds_1d$method[1] else "comparison incomplete"),
  paste("Best 2D MDS geometry by normalized stress:", if (nrow(best_mds_2d) > 0) best_mds_2d$method[1] else "comparison incomplete"),
  paste("Best 1D LSM geometry by BIC:", if (nrow(best_lsm_1d) > 0) best_lsm_1d$method[1] else "comparison incomplete"),
  paste("Best 2D LSM geometry by BIC:", if (nrow(best_lsm_2d) > 0) best_lsm_2d$method[1] else "comparison incomplete"),
  paste("Best framework by held-out decoder log-loss:", best_framework$method[1]),
  "",
  "All methods originate from the same cleaned Members' Votes sample.",
  "Dimension reduction uses member-by-roll-call vote profiles.",
  "MDS uses pairwise voting disagreement, and LSM uses the derived co-voting network.",
  "H119_rollcalls.json is used only for metadata auditing and is not a model feature source.",
  "The framework table reports the best available method within each family.",
  "An overall framework winner is reported only after every core method has five finite fold losses.",
  "Latentnet Euclidean LSM is a supplementary implementation check and is not ranked with direct LSM fits.",
  "Soft Circular LSM is exploratory; its penalized objective is not treated as an ordinary BIC.",
  "Important: held-out decoder scores use coordinates fitted once on the full matrix.",
  "They are a predictive diagnostic and should not be described as fully nested five-fold re-estimation."
)

writeLines(
  winner_lines,
  file.path(results_directory, "winner_summary.txt")
)

print(method_metrics)
print(geometry_comparison_full)
print(framework_comparison)

