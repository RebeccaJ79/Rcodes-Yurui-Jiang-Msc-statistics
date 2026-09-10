#### Voteview graphical interpretation ---------------------------------------
# This script uses a small reproducible subset of the 119th U.S. House data.

# Packages and working directory
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(FactoMineR)
library(smacof)
library(vegan)
library(Gifi)
library(latentnet)
library(network)

setwd("E:/MSc_dissertation/data")

set.seed(123)

party_colours <- c(
  "Democrat" = "#1f77b4",
  "Republican" = "#d62728",
  "Other" = "grey55"
)


# Read Voteview data and recode votes
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


# Select the reproducible demonstration subset
selected_members <- members_119 |>
  filter(
    party %in% c("Democrat", "Republican"),
    !is.na(nominate_dim1),
    !is.na(nominate_dim2)
  ) |>
  group_by(
    party
  ) |>
  slice_sample(
    n = 30
  ) |>
  ungroup()

votes_selected_members <- votes_119 |>
  filter(
    icpsr %in% selected_members$icpsr
  )

vote_summary <- votes_selected_members |>
  group_by(
    rollnumber
  ) |>
  summarise(
    number_observed = sum(!is.na(vote_binary)),
    yea_share = mean(vote_binary, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(
    number_observed >= 0.90 * nrow(selected_members),
    yea_share >= 0.20,
    yea_share <= 0.80
  )

selected_votes <- vote_summary |>
  slice_sample(
    n = min(30, nrow(vote_summary))
  )

votes_demo <- votes_119 |>
  filter(
    icpsr %in% selected_members$icpsr,
    rollnumber %in% selected_votes$rollnumber
  )


# Construct the legislator-by-roll-call vote matrix 
vote_matrix_demo <- votes_demo |>
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

# Convert the tibble into a numeric matrix
vote_mat <- as.matrix(tibble::column_to_rownames(vote_matrix_demo, "icpsr"))

# Keep legislators with at least 80% observed votes
member_completeness <- rowMeans(!is.na(vote_mat))
vote_mat <- vote_mat[member_completeness >= 0.80, , drop = FALSE]

# Keep roll calls observed for at least 80% of retained legislators
vote_completeness <- colMeans(!is.na(vote_mat))
vote_mat <- vote_mat[, vote_completeness >= 0.80, drop = FALSE]

# Store retained legislator IDs and matching metadata
member_ids <- rownames(vote_mat)
members_demo <- selected_members |>
  filter(
    icpsr %in% as.numeric(member_ids)
  ) |>
  arrange(
    match(icpsr, as.numeric(member_ids))
  )

# Check the final demonstration dataset
dim(vote_mat)
mean(is.na(vote_mat))

#### Plot the vote matrix and DW-NOMINATE benchmark -------------------------
vote_heatmap_data <- as.data.frame(as.table(vote_mat))
names(vote_heatmap_data) <- c("icpsr", "roll_call", "vote")

# Vote heatmap
vote_heatmap <- ggplot(
  vote_heatmap_data,
  aes(
    x = roll_call,
    y = icpsr,
    fill = vote
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "black",
    na.value = "grey85",
    name = "Vote"
  ) +
  labs(
    title = "",
    subtitle = "White = Nay, black = Yea, grey = missing",
    x = "Roll call",
    y = "Legislator"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
print(vote_heatmap)

# Supplied DW-NOMINATE benchmark
nominate_plot <- ggplot(
  members_demo,
  aes(
    x = nominate_dim1,
    y = nominate_dim2,
    colour = party
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  scale_colour_manual(
    values = party_colours
  ) +
  labs(
    title = "",
    x = "NOMINATE Dimension 1",
    y = "NOMINATE Dimension 2",
    colour = "Party"
  ) +
  theme_minimal()
print(nominate_plot)

ggsave(
  "vote_matrix.png", vote_heatmap, width = 6.8, height = 5, units = "in", dpi = 300
)
ggsave(
  "DW_NOMINATE.png", nominate_plot, width = 6.8, height = 5, units = "in", dpi = 300
)

#### Static W-NOMINATE benchmark -------------------------------------------
if (requireNamespace("wnominate", quietly = TRUE)) {
  message(
    "Package 'wnominate' is installed. Check its local version-specific " ,
    "roll-call object interface before fitting a static W-NOMINATE model."
  )
} else {
  message(
    "Package 'wnominate' is not installed; supplied DW-NOMINATE remains " ,
    "an external benchmark only."
  )
}

# Impute missing votes
vote_mat_imputed <- vote_mat
for (vote_index in seq_len(ncol(vote_mat_imputed))) {
  observed_yea_rate <- mean(
    vote_mat_imputed[, vote_index],
    na.rm = TRUE
  )
  vote_mat_imputed[
    is.na(vote_mat_imputed[, vote_index]),
    vote_index
  ] <- observed_yea_rate
}


# Construct pairwise voting-disagreement matrix
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

# Replace unavailable pairwise distances by the median observed disagreement
median_disagreement <- median(disagreement_matrix[upper.tri(disagreement_matrix)], na.rm = TRUE)

disagreement_matrix[is.na(disagreement_matrix)] <- median_disagreement
diag(disagreement_matrix) <- 0


#### PCA ---------------------------------------------------------------------
pca_fit <- prcomp(vote_mat_imputed, center = TRUE, scale. = TRUE)
pca_explained_variance <- pca_fit$sdev ^ 2 / sum(pca_fit$sdev ^ 2)
pca_scree_data <- data.frame(
  component = seq_along(pca_explained_variance),
  explained_variance = pca_explained_variance,
  cumulative_variance = cumsum(pca_explained_variance)
)

pca_scree_plot <- ggplot(
  pca_scree_data,
  aes(
    x = component,
    y = explained_variance
  )
) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = " ",
    x = "Principal component",
    y = "Explained variance"
  ) +
  theme_minimal()
print(pca_scree_plot)

pca_scores <- data.frame(icpsr = as.numeric(rownames(pca_fit$x)), 
                         Dimension_1 = pca_fit$x[, 1], Dimension_2 = pca_fit$x[, 2]) |>
  left_join(
    members_demo |>
      select(icpsr, party),
    by = "icpsr"
  )

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = Dimension_1,
    y = Dimension_2,
    colour = party
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  scale_colour_manual(values = party_colours) +
  labs(
    title = " ",
    x = "PC1",
    y = "PC2",
    colour = "Party"
  ) +
  theme_minimal()
print(pca_plot)
ggsave(
  "PCA.png", pca_scree_plot + pca_plot + patchwork::plot_layout(ncol = 2),
  width = 9, height = 4.2, units = "in", dpi = 300
)

#### Classical, metric, and ordinal MDS -------------------------------------
classical_mds <- cmdscale(as.dist(disagreement_matrix), k = 2, eig = TRUE)
metric_mds <- smacofSym(as.dist(disagreement_matrix), ndim = 2, type = "ratio")
ordinal_mds <- smacofSym(as.dist(disagreement_matrix), ndim = 2,
  type = "ordinal", ties = "primary")

make_score_plot <- function(
    coordinates,
    title,
    x_label,
    y_label,
    show_party_legend = TRUE
) {
  plot_data <- data.frame(
    icpsr = as.numeric(rownames(coordinates)),
    Dimension_1 = coordinates[, 1],
    Dimension_2 = coordinates[, 2]
  ) |>
    left_join(
      members_demo |>
        select(icpsr, party),
      by = "icpsr"
    )

  score_plot <- ggplot(
    plot_data,
    aes(
      x = Dimension_1,
      y = Dimension_2,
      colour = party
    )
  ) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = party_colours) +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      colour = "Party"
    ) +
    theme_minimal()

  if (!show_party_legend) {
    score_plot <- score_plot +
      theme(
        legend.position = "none"
      )
  }
  score_plot
}

classical_mds_plot <- make_score_plot(
  classical_mds$points,
  "Classical MDS",
  "MDS Dimension 1",
  "MDS Dimension 2",
  show_party_legend = FALSE
)
classical_mds_plot

metric_mds_plot <- make_score_plot(
  metric_mds$conf,
  "Metric MDS",
  "Metric Dimension 1",
  "Metric Dimension 2",
  show_party_legend = FALSE
)
metric_mds_plot

ordinal_mds_plot <- make_score_plot(
  ordinal_mds$conf,
  "Ordinal MDS",
  "Ordinal Dimension 1",
  "Ordinal Dimension 2",
  show_party_legend = TRUE
)
ordinal_mds_plot


#### Circular MDS ------------------------------------------------------------
# Rescale observed disagreement to the angular interval [0, pi]
angular_target <- pi * disagreement_matrix / max(disagreement_matrix)

# Initialise angles using the classical MDS configuration
initial_angles <- atan2(classical_mds$points[, 2], classical_mds$points[, 1]) %% (2 * pi)

# The first legislator is fixed at zero to remove rotational ambiguity.
circular_stress <- function(free_angles) {
  theta <- c(0, free_angles %% (2 * pi))
  angular_difference <- abs(outer(theta, theta, "-"))
  fitted_distance <- pmin(
    angular_difference,
    2 * pi - angular_difference
  )

  sum(
    (
      fitted_distance[upper.tri(fitted_distance)] -
        angular_target[upper.tri(angular_target)]
    ) ^ 2
  )
}

circular_fit <- optim(
  par = initial_angles[-1],
  fn = circular_stress,
  method = "BFGS",
  control = list(maxit = 3000)
)

circular_angles <- c(0, circular_fit$par %% (2 * pi))
circular_coordinates <- cbind(
  Dimension_1 = cos(circular_angles),
  Dimension_2 = sin(circular_angles)
)
rownames(circular_coordinates) <- member_ids

# Coordinates used to draw the boundary of the unit circle.
circle_outline <- data.frame(
  theta = seq(
    0,
    2 * pi,
    length.out = 500
  )
) |>
  mutate(
    x = cos(theta),
    y = sin(theta)
  )

circular_plot <- make_score_plot(
  circular_coordinates,
  "Circular MDS",
  "cos(theta)",
  "sin(theta)"
) +
  geom_path(
    data = circle_outline,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.4
  ) +
  coord_equal(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  )
circular_plot

#### Circular projection and circularity index ------------------------------
classical_mds_centre <- colMeans(classical_mds$points)
classical_mds_centered <- sweep(classical_mds$points, 2, classical_mds_centre, FUN = "-")

# A low circularity index means the Euclidean representation is close to a circle.
classical_mds_radius <- sqrt(rowSums(classical_mds_centered ^ 2))

circularity_index <- sd(classical_mds_radius) / mean(classical_mds_radius)
print(circularity_index)

# Project the Euclidean MDS configuration onto the unit circle.
classical_mds_projection <- classical_mds_centered / classical_mds_radius
rownames(classical_mds_projection) <- member_ids

projection_plot <- make_score_plot(
  classical_mds_projection,
  "Circular projection",
  "Projected Dimension 1",
  "Projected Dimension 2",
  show_party_legend = FALSE
) +
  geom_path(
    data = circle_outline,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.4
  ) +
  coord_equal() +
  labs(
    subtitle = paste0(
      "Circularity index = ",
      round(circularity_index, 2)
    )
  )
projection_plot
ggsave("MDS.png",
       ((classical_mds_plot + metric_mds_plot) / ordinal_mds_plot /
          (projection_plot + circular_plot) +
          patchwork::plot_layout(heights = c(1, 1, 1))
       ) & theme(legend.position = "bottom"),
       width = 6, height = 10.5, units = "in", dpi = 300
)

#### Soft circular extension -----------------------------------------------
soft_circular_lambda <- 5
soft_circular_start <- pmin(3, pmax(0.2, classical_mds_radius / mean(classical_mds_radius)))
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

if (inherits(soft_circular_fit, "error")) stop("Soft circular MDS failed: ", conditionMessage(soft_circular_fit))
soft_circular_radii <- exp(soft_circular_fit$par)
soft_circular_coordinates <- cbind(Dimension_1 = soft_circular_radii * cos(circular_angles), Dimension_2 = soft_circular_radii * sin(circular_angles))
rownames(soft_circular_coordinates) <- member_ids

soft_circular_outline <- data.frame(theta = seq(0, 2 * pi, length.out = 500)) |>
  mutate(x = mean(soft_circular_radii) * cos(theta), y = mean(soft_circular_radii) * sin(theta))

soft_circular_plot <- 
  make_score_plot(soft_circular_coordinates, "", "Soft circular Dimension 1", "Soft circular Dimension 2") +
  geom_path(data = soft_circular_outline, aes(x = x, y = y), inherit.aes = FALSE, colour = "black", linetype = "dashed", linewidth = 0.4) +
  coord_equal() + labs(subtitle = paste0("Radial variation = ", round(sd(soft_circular_radii) / mean(soft_circular_radii), 2)))
print(soft_circular_plot)

ggsave("Soft_Circular.png", soft_circular_plot, width = 6.8, height = 4.8, units = "in", dpi = 300)

#### Correspondence analysis -----------------------------------------------
# Convert every roll call into a categorical Yea/Nay variable for CA/MCA.
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

ca_fit <- FactoMineR::CA(vote_mat_imputed,graph = FALSE)
ca_row_coordinates <- as.data.frame(ca_fit$row$coord[, 1:2, drop = FALSE])
names(ca_row_coordinates) <- c("Dimension_1", "Dimension_2")
ca_row_coordinates$icpsr <- rownames(ca_fit$row$coord)
ca_row_coordinates <- ca_row_coordinates |>
  left_join(
    members_demo |>
      transmute(
        icpsr = as.character(icpsr),
        party
      ),
    by = "icpsr"
  )

ca_plot <- ggplot(
  ca_row_coordinates,
  aes(Dimension_1, Dimension_2, colour = party)
) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = party_colours, drop = FALSE) +
  coord_fixed(ratio = 4) +
  theme_minimal() +
  labs(
    title = "",
    x = "CA Dimension 1",
    y = "CA Dimension 2",
    colour = "Party"
  )
print(ca_plot)

ca_scree_data <- data.frame(dimension = seq_along(ca_fit$eig[, 1]), inertia = ca_fit$eig[, 1])
ca_scree_plot <- ggplot(
  ca_scree_data,
  aes(dimension, inertia)
) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    title = "",
    x = "Dimension",
    y = "Inertia"
  )
ca_scree_plot

ggsave(
  "CA.png", ca_scree_plot + ca_plot + patchwork::plot_layout(ncol = 2),
  width = 8, height = 4.2, units = "in", dpi = 300
)

#### Multiple correspondence analysis --------------------------------------

mca_fit <- FactoMineR::MCA(vote_factor_data, ncp = 5, graph = FALSE)

mca_row_coordinates <- as.data.frame(mca_fit$ind$coord[, 1:2, drop = FALSE])

names(mca_row_coordinates) <- c("Dimension_1", "Dimension_2")
mca_row_coordinates$icpsr <- rownames(mca_fit$ind$coord)

mca_row_coordinates <- mca_row_coordinates |>
  left_join(
    members_demo |>
      transmute(
        icpsr = as.character(icpsr),
        party
      ),
    by = "icpsr"
  )

mca_plot <- ggplot(
  mca_row_coordinates,
  aes(Dimension_1, Dimension_2, colour = party)
) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = party_colours, drop = FALSE) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Multiple correspondence analysis",
    subtitle = "Yea/Nay categories treated as nominal",
    x = "MCA Dimension 1",
    y = "MCA Dimension 2",
    colour = "Party"
  )

mca_scree_data <- data.frame(
  dimension = seq_len(nrow(mca_fit$eig)),
  eigenvalue = mca_fit$eig[, 1]
)

mca_scree_plot <- ggplot(
  mca_scree_data,
  aes(dimension, eigenvalue)
) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    title = "MCA eigenvalues",
    x = "Dimension",
    y = "Eigenvalue"
  )
print(mca_scree_plot)

#### Constrained canonical correspondence analysis -------------------------

# Party and supplied NOMINATE coordinates are excluded from the constraints.
# State and district are non-party member covariates used only for this
# graphical constrained-ordination adaptation.
cca_members <- members_demo |>
  transmute(
    state = factor(ifelse(is.na(state_icpsr), "Unknown", state_icpsr)),
    district = factor(ifelse(is.na(district_code), "Unknown", district_code))
  )

cca_fit <- vegan::cca(vote_mat_imputed ~ state + district, data = cca_members)
cca_site_scores <- as.data.frame(vegan::scores(cca_fit, display = "sites", choices = 1:2))

names(cca_site_scores) <- c("Dimension_1", "Dimension_2")
cca_site_scores$party <- members_demo$party

cca_plot <- ggplot(
  cca_site_scores,
  aes(Dimension_1, Dimension_2, colour = party)
) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = party_colours, drop = FALSE) +
  scale_x_continuous(
    expand = expansion(mult = c(0.20, 0.20))
  ) +
  coord_fixed(ratio = 0.70) +
  theme_minimal() +
  labs(
    title = "Constrained CA",
    x = "CCA Dimension 1",
    y = "CCA Dimension 2",
    colour = "Party"
  )
print(cca_plot)

#### Principal geodesic analysis --------------------------------------------

# Normalize centered voting profiles to the unit sphere and perform PCA in
# the tangent space at their normalized mean direction.
pga_center <- colMeans(vote_mat_imputed)
pga_centered <- sweep(vote_mat_imputed, 2, pga_center, FUN = "-")
pga_norms <- sqrt(rowSums(pga_centered ^ 2))
pga_unit_profiles <- pga_centered / pmax(pga_norms, 1e-12)
pga_frechet_mean <- colMeans(pga_unit_profiles)
pga_frechet_mean <- pga_frechet_mean /
  sqrt(sum(pga_frechet_mean ^ 2))

pga_projection_coefficients <- as.vector(pga_unit_profiles %*% pga_frechet_mean)

pga_tangent_scores <- pga_unit_profiles - pga_projection_coefficients * matrix(pga_frechet_mean, nrow = nrow(pga_unit_profiles), ncol = ncol(pga_unit_profiles), byrow = TRUE)

pga_fit <- prcomp(pga_tangent_scores, center = TRUE, scale. = FALSE)

pga_coordinates <- pga_fit$x[, 1:2, drop = FALSE]
rownames(pga_coordinates) <- member_ids

pga_plot <- make_score_plot(
  pga_coordinates,
  "Principal geodesic analysis",
  "Tangent Dimension 1",
  "Tangent Dimension 2",
  show_party_legend = TRUE
) 

#### Kernel PCA --------------------------------------------------------------

if (requireNamespace("kernlab", quietly = TRUE)) {
  kernel_pca_fit <- kernlab::kpca(
    as.matrix(vote_mat_imputed),
    kernel = "rbfdot",
    kpar = list(sigma = 0.05),
    features = 2
  )

  kernel_pca_coordinates <- as.matrix(kernlab::rotated(kernel_pca_fit))
  rownames(kernel_pca_coordinates) <- member_ids

  kernel_pca_plot <- make_score_plot(
    kernel_pca_coordinates,
    "Kernel PCA",
    "Kernel Component 1",
    "Kernel Component 2",
    show_party_legend = FALSE
  )
} else {
  message("Package 'kernlab' is not installed; Kernel PCA plot skipped.")
  kernel_pca_plot <- ggplot() +
    theme_void() +
    labs(title = "Kernel PCA unavailable")
}

ggsave(
  "MCA.png", mca_scree_plot + mca_plot + patchwork::plot_layout(ncol = 2),
  width = 9, height = 4.2, units = "in", dpi = 300
)
ggsave("CCA.png", cca_plot, width = 4, height = 4.2, units = "in", dpi = 300)
ggsave("PGA.png", pga_plot, width = 6.8, height = 4.2, units = "in", dpi = 300)
ggsave("Kernel_PCA.png", kernel_pca_plot, width = 6.8, height = 4.2, units = "in", dpi = 300)

#### Nonlinear PCA -----------------------------------------------------------
# Convert every roll call into an ordinal Yea/Nay variable.
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
        levels = c("Nay", "Yea"),
        ordered = TRUE
      )
    }
  )
)

rownames(vote_factor_data) <- member_ids

nlpca_fit <- princals(
  vote_factor_data,
  ndim = 2,
  ordinal = TRUE
)

nlpca_scores <- nlpca_fit$objectscores[, 1:2, drop = FALSE]
rownames(nlpca_scores) <- member_ids

nlpca_plot <- make_score_plot(
  nlpca_scores,
  "Nonlinear PCA",
  "Nonlinear Component 1",
  "Nonlinear Component 2",
  show_party_legend = FALSE
)
print(nlpca_plot)
ggsave(
  "NLPCA.png", nlpca_plot, width = 6.8, height = 4.2, units = "in", dpi = 300
)

#### HOMALS / nonlinear MCA --------------------------------------------------
homals_fit <- Gifi::homals(vote_factor_data, ndim = 2)
homals_scores <- homals_fit$objectscores[, 1:2, drop = FALSE]
rownames(homals_scores) <- member_ids

homals_plot <- make_score_plot(
  homals_scores,
  "HOMALS",
  "HOMALS Dimension 1",
  "HOMALS Dimension 2",
  show_party_legend = TRUE
)
print(homals_plot)

homals_eigenvalue_data <- data.frame(
  dimension = seq_along(homals_fit$evals),
  eigenvalue = homals_fit$evals
)

homals_scree_plot <- ggplot(
  homals_eigenvalue_data,
  aes(
    x = dimension,
    y = eigenvalue
  )
) +
  geom_line() +
  geom_point() +
  labs(
    title = "HOMALS eigenvalues",
    x = "Dimension",
    y = "Eigenvalue"
  ) +
  theme_minimal()
print(homals_scree_plot)

ggsave(
  "HOMALS.png", homals_scree_plot + homals_plot + patchwork::plot_layout(ncol = 2),
  width = 9, height = 4.2, units = "in", dpi = 300
)

#### Nonlinear canonical correlation analysis --------------------------------
# Split roll calls into odd and even numbered sets.
odd_roll_calls <- vote_factor_data[seq(1, ncol(vote_factor_data), by = 2)]
even_roll_calls <- vote_factor_data[seq(2, ncol(vote_factor_data), by = 2)]

# corAspect requires unique variable names across the two variable sets.
names(odd_roll_calls) <- paste0("odd_roll_", seq_len(ncol(odd_roll_calls)))
names(even_roll_calls) <- paste0("even_roll_", seq_len(ncol(even_roll_calls))) 

# Optimal-scale each roll-call set separately.
odd_nlpca <- princals(odd_roll_calls, ndim = 2, ordinal = TRUE)
even_nlpca <- princals(even_roll_calls, ndim = 2, ordinal = TRUE)

# Use the optimally scaled object scores as the two variable sets for CCA.
nlcca_x <- odd_nlpca$objectscores[, 1:2]
nlcca_y <- even_nlpca$objectscores[, 1:2]

nlcca_cancor <- cancor(nlcca_x, nlcca_y)

nlcca_u <- nlcca_x %*% nlcca_cancor$xcoef[, 1:2, drop = FALSE]
nlcca_v <- nlcca_y %*% nlcca_cancor$ycoef[, 1:2, drop = FALSE]

# Average aligned canonical scores to obtain one two-dimensional member map.
nlcca_scores <- (scale(nlcca_u) + scale(nlcca_v)) / 2
rownames(nlcca_scores) <- member_ids

nlcca_plot <- make_score_plot(
  nlcca_scores,
  "Nonlinear CCA",
  "NLCCA Dimension 1",
  "NLCCA Dimension 2",
  show_party_legend = TRUE
)
print(nlcca_plot)
ggsave("NLCCA.png", nlcca_plot, width = 6.8, height = 4.2, units = "in", dpi = 300)

#### Horseshoe diagnostic ----------------------------------------------------
horseshoe_fit <- lm(
  ordinal_mds$conf[, 2] ~
    ordinal_mds$conf[, 1] +
    I(ordinal_mds$conf[, 1] ^ 2)
)

horseshoe_data <- data.frame(
  Dimension_1 = ordinal_mds$conf[, 1],
  Dimension_2 = ordinal_mds$conf[, 2],
  party = members_demo$party
)

horseshoe_plot <- ggplot(
  horseshoe_data,
  aes(
    x = Dimension_1,
    y = Dimension_2,
    colour = party
  )
) +
  geom_point(size = 2.5) +
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x ^ 2),
    se = FALSE,
    colour = "black"
  ) +
  scale_colour_manual(values = party_colours) +
  labs(
    title = "Horseshoe diagnostic for ordinal MDS",
    subtitle = paste0(
      "Quadratic R-squared = ",
      round(summary(horseshoe_fit)$r.squared, 3)
    ),
    x = "Ordinal Dimension 1",
    y = "Ordinal Dimension 2",
    colour = "Party"
  ) +
  theme_minimal()

ggsave(
  "Horseshoe.png",
  horseshoe_plot,
  width = 5.5,
  height = 4.2,
  units = "in",
  dpi = 300
)

#### Latent space models -----------------------------------------------------
# Construct an undirected strong co-voting network.
agreement_matrix <- 1 - disagreement_matrix
agreement_threshold <- quantile(
  agreement_matrix[upper.tri(agreement_matrix)],
  probs = 0.75,
  na.rm = TRUE
)
covoting_adjacency <- ifelse(agreement_matrix >= agreement_threshold, 1, 0)
diag(covoting_adjacency) <- 0
covoting_network <- network(
  covoting_adjacency,
  directed = FALSE,
  matrix.type = "adjacency"
)
network.vertex.names(covoting_network) <- members_demo$bioname
set.vertex.attribute(
  covoting_network,
  "party",
  members_demo$party
)

#### Euclidean latent-distance LSM ------------------------------------------
set.seed(123)
euclidean_lsm <- ergmm(covoting_network ~ euclidean(d = 2), family = "Bernoulli")

# Posterior mean Euclidean latent coordinates.
euclidean_lsm_coordinates <- apply(euclidean_lsm$sample$Z, c(2, 3), mean)
rownames(euclidean_lsm_coordinates) <- member_ids

euclidean_lsm_plot <- make_score_plot(
  euclidean_lsm_coordinates,
  "Euclidean LSM",
  "Latent Dimension 1",
  "Latent Dimension 2",
  show_party_legend = FALSE
)
euclidean_lsm_plot

#### Bilinear LSM ------------------------------------------------------------
set.seed(123)
bilinear_lsm <- ergmm(covoting_network ~ bilinear(d = 2), family = "Bernoulli")

# Posterior mean bilinear latent coordinates.
bilinear_coordinates <- apply(bilinear_lsm$sample$Z, c(2, 3), mean)
rownames(bilinear_coordinates) <- member_ids

bilinear_lsm_plot <- make_score_plot(
  bilinear_coordinates,
  "Bilinear LSM",
  "Latent Dimension 1",
  "Latent Dimension 2",
  show_party_legend = FALSE
) +
  geom_segment(
    data = data.frame(
      x = 0,
      y = 0,
      xend = bilinear_coordinates[, 1],
      yend = bilinear_coordinates[, 2]
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = "grey60",
    arrow = grid::arrow(length = grid::unit(0.08, "inches"))
  )

#### Circular projection of bilinear LSM ------------------------------------
bilinear_radius <- sqrt(rowSums(bilinear_coordinates ^ 2))
bilinear_circular_projection <- bilinear_coordinates / bilinear_radius
rownames(bilinear_circular_projection) <- member_ids

bilinear_projection_plot <- make_score_plot(
  bilinear_circular_projection,
  "Circular projection of bilinear LSM",
  "Projected Dimension 1",
  "Projected Dimension 2",
  show_party_legend = FALSE
) +
  geom_path(
    data = circle_outline,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.4
  ) +
  coord_equal(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) 
print(bilinear_projection_plot)

#### Direct circular LSM -----------------------------------------------------
circular_lsm_negative_log_likelihood <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  theta <- c(0, parameters[3:(number_of_members + 1)])

  linear_predictor <- matrix(
    alpha,
    nrow = number_of_members,
    ncol = number_of_members
  ) +
    lambda * cos(outer(theta, theta, "-"))

  tie_probability <- plogis(linear_predictor)
  # Each undirected legislator pair is used once, not twice.
  use_dyad <- upper.tri(covoting_adjacency)

  -sum(
    covoting_adjacency[use_dyad] *
      log(tie_probability[use_dyad] + 1e-12) +
      (1 - covoting_adjacency[use_dyad]) *
      log(1 - tie_probability[use_dyad] + 1e-12)
  )
}

set.seed(123)
circular_lsm_start <- c(
  qlogis(mean(covoting_adjacency[upper.tri(covoting_adjacency)])),
  log(1),
  runif(number_of_members - 1, -pi, pi)
)

circular_lsm_fit <- optim(
  par = circular_lsm_start,
  fn = circular_lsm_negative_log_likelihood,
  method = "BFGS",
  control = list(maxit = 5000)
)

if (circular_lsm_fit$convergence != 0) {
  warning("Circular LSM did not report convergence; inspect circular_lsm_fit.")
}

circular_lsm_angles <- c(0, circular_lsm_fit$par[3:(number_of_members + 1)] %% (2 * pi))
circular_lsm_coordinates <- cbind(
  Dimension_1 = cos(circular_lsm_angles),
  Dimension_2 = sin(circular_lsm_angles)
)
rownames(circular_lsm_coordinates) <- member_ids

circular_lsm_plot <- make_score_plot(
  circular_lsm_coordinates,
  "Circular LSM",
  "cos(theta)",
  "sin(theta)",
  show_party_legend = TRUE
) +
  geom_path(
    data = circle_outline,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.4
  ) +
  coord_equal(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  )
circular_lsm_plot

#### Soft circular LSM ------------------------------------------------------
soft_circular_lsm_penalty <- 5
soft_circular_lsm_negative_log_likelihood <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  theta <- c(0, parameters[3:(number_of_members + 1)])
  log_radii <- c(0, parameters[(number_of_members + 2):(2 * number_of_members)])
  log_radii <- log_radii - mean(log_radii)
  radii <- exp(log_radii)
  soft_circular_lsm_coordinates <- cbind(radii * cos(theta), radii * sin(theta))
  fitted_distance <- as.matrix(dist(soft_circular_lsm_coordinates))
  tie_probability <- plogis(alpha - lambda * fitted_distance)
  use_dyad <- upper.tri(covoting_adjacency)
  radius_penalty <- soft_circular_lsm_penalty * sum((radii - mean(radii)) ^ 2)
  -sum(covoting_adjacency[use_dyad] * log(tie_probability[use_dyad] + 1e-12) + (1 - covoting_adjacency[use_dyad]) * log(1 - tie_probability[use_dyad] + 1e-12)) + radius_penalty
}

soft_circular_lsm_start <- c(
  circular_lsm_fit$par[1],
  circular_lsm_fit$par[2],
  circular_lsm_angles[-1],
  rep(0, number_of_members - 1)
)

soft_circular_lsm_fit <- optim(
  par = soft_circular_lsm_start,
  fn = soft_circular_lsm_negative_log_likelihood,
  method = "BFGS",
  control = list(maxit = 2000)
)

if (soft_circular_lsm_fit$convergence != 0) {
  warning("Soft circular LSM did not report convergence; inspect soft_circular_lsm_fit.")
}

soft_circular_lsm_angles <- c(0, soft_circular_lsm_fit$par[3:(number_of_members + 1)] %% (2 * pi))
soft_circular_lsm_log_radii <- c(0, soft_circular_lsm_fit$par[(number_of_members + 2):(2 * number_of_members)])
soft_circular_lsm_log_radii <- soft_circular_lsm_log_radii - mean(soft_circular_lsm_log_radii)
soft_circular_lsm_radii <- exp(soft_circular_lsm_log_radii)
soft_circular_lsm_coordinates <- cbind(Dimension_1 = soft_circular_lsm_radii * cos(soft_circular_lsm_angles), Dimension_2 = soft_circular_lsm_radii * sin(soft_circular_lsm_angles))
rownames(soft_circular_lsm_coordinates) <- member_ids

soft_circular_lsm_outline <- data.frame(theta = seq(0, 2 * pi, length.out = 500)) |>
  mutate(x = mean(soft_circular_lsm_radii) * cos(theta), y = mean(soft_circular_lsm_radii) * sin(theta))

soft_circular_lsm_plot <- make_score_plot(
  soft_circular_lsm_coordinates,
  "",
  "Soft circular Dimension 1",
  "Soft circular Dimension 2",
  show_party_legend = FALSE
) +
  geom_path(data = soft_circular_lsm_outline, aes(x = x, y = y), inherit.aes = FALSE, colour = "black", linetype = "dashed", linewidth = 0.4) +
  coord_equal() +
  labs(subtitle = paste0("Radial variation = ", round(sd(soft_circular_lsm_radii) / mean(soft_circular_lsm_radii), 2)))
print(soft_circular_lsm_plot)
ggsave("Soft_Circular_LSM.png", soft_circular_lsm_plot, width = 6.8, height = 4.8, units = "in", dpi = 300)

#### Direct spherical LSM ----------------------------------------------------
spherical_lsm_dyads <- which(upper.tri(covoting_adjacency), arr.ind = TRUE)
spherical_lsm_outcomes <- covoting_adjacency[cbind(spherical_lsm_dyads[, 1], spherical_lsm_dyads[, 2])]
spherical_lsm_negative_log_likelihood <- function(parameters) {
  alpha <- parameters[1]
  lambda <- exp(parameters[2])
  longitude <- c(0, parameters[3:(number_of_members + 1)])
  colatitude <- parameters[(number_of_members + 2):(2 * number_of_members + 1)]
  spherical_coordinates <- cbind(sin(colatitude) * cos(longitude), sin(colatitude) * sin(longitude), cos(colatitude))
  similarity <- rowSums(spherical_coordinates[spherical_lsm_dyads[, 1], , drop = FALSE] * spherical_coordinates[spherical_lsm_dyads[, 2], , drop = FALSE])
  tie_probability <- plogis(alpha + lambda * similarity)
  -sum(spherical_lsm_outcomes * log(tie_probability + 1e-12) + (1 - spherical_lsm_outcomes) * log(1 - tie_probability + 1e-12))
}

set.seed(123)
spherical_initial <- cmdscale(as.dist(disagreement_matrix), k = 3)
spherical_initial <- sweep(
  spherical_initial,
  1,
  pmax(sqrt(rowSums(spherical_initial ^ 2)), 1e-8),
  FUN = "/"
)
spherical_initial_longitude <- atan2(spherical_initial[, 2], spherical_initial[, 1])
spherical_initial_colatitude <- acos(pmin(1, pmax(-1, spherical_initial[, 3])))
spherical_lsm_start <- c(qlogis(pmin(0.99, pmax(0.01, mean(spherical_lsm_outcomes)))), log(1), spherical_initial_longitude[-1], spherical_initial_colatitude)
spherical_lsm_fit <- optim(par = spherical_lsm_start, fn = spherical_lsm_negative_log_likelihood, method = "L-BFGS-B", lower = c(-8, -5, rep(-pi, number_of_members - 1), rep(1e-3, number_of_members)), upper = c(8, 5, rep(pi, number_of_members - 1), rep(pi - 1e-3, number_of_members)), control = list(maxit = 800, factr = 1e7, pgtol = 1e-5))

if (spherical_lsm_fit$convergence != 0) {
  warning("Spherical LSM did not report convergence; inspect spherical_lsm_fit.")
}

spherical_longitude <- c(0, spherical_lsm_fit$par[3:(number_of_members + 1)])
spherical_colatitude <- spherical_lsm_fit$par[(number_of_members + 2):(2 * number_of_members + 1)]
spherical_lsm_coordinates <- cbind(
  Dimension_1 = sin(spherical_colatitude) * cos(spherical_longitude),
  Dimension_2 = sin(spherical_colatitude) * sin(spherical_longitude),
  Dimension_3 = cos(spherical_colatitude)
)
rownames(spherical_lsm_coordinates) <- member_ids

spherical_lsm_plot <- make_score_plot(
  spherical_lsm_coordinates[, 1:2],
  "Spherical LSM",
  "Sphere coordinate 1",
  "Sphere coordinate 2",
  show_party_legend = FALSE
) +
  coord_equal()
spherical_lsm_plot

bilinear_lsm_cell <- bilinear_lsm_plot + bilinear_projection_plot + patchwork::plot_layout(ncol = 2)
ggsave(
  "LSM.png",
  ((euclidean_lsm_plot + soft_circular_lsm_plot) /
      bilinear_lsm_cell /
      (circular_lsm_plot + spherical_lsm_plot) +
      patchwork::plot_layout(guides = "collect", heights = c(1, 1, 1))
  ) & theme(legend.position = "bottom"),
  width = 8, height = 13, units = "in", dpi = 300
)