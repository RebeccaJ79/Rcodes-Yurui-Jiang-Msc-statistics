#### 119th House Voteview dataset description -------------------------------
# This tutorial-style script describes the cleaned Members' Votes data.
# Run the sections from top to bottom in RStudio.

#### Packages and working directory ------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

setwd("E:/MSc_dissertation/data")

#### Read and clean the 119th House voting records ---------------------------
votes_119 <- read_csv("H119_votes.csv", show_col_types = FALSE) |>
  filter(chamber == "House") |>
  mutate(vote_binary = case_when(cast_code %in% c(1, 2, 3) ~ 1, cast_code %in% c(4, 5, 6) ~ 0, TRUE ~ NA_real_))

members_119 <- read_csv("H119_members.csv", show_col_types = FALSE) |>
  filter(chamber == "House") |>
  mutate(party = case_when(party_code == 100 ~ "Democrat", party_code == 200 ~ "Republican", TRUE ~ "Other"))

roll_call_summary <- votes_119 |>
  group_by(rollnumber) |>
  summarise(number_observed = sum(!is.na(vote_binary)), yea_share = mean(vote_binary, na.rm = TRUE), .groups = "drop") |>
  filter(number_observed >= 0.80 * max(number_observed), yea_share >= 0.05, yea_share <= 0.95)

vote_matrix <- votes_119 |>
  filter(icpsr %in% members_119$icpsr, rollnumber %in% roll_call_summary$rollnumber) |>
  select(icpsr, rollnumber, vote_binary) |>
  distinct() |>
  pivot_wider(names_from = rollnumber, values_from = vote_binary, names_prefix = "vote_") |>
  column_to_rownames("icpsr") |>
  as.matrix()

#### Apply the common 80% completeness rule -------------------------------
repeat {
  old_dimensions <- dim(vote_matrix)
  member_completeness <- rowMeans(!is.na(vote_matrix))
  vote_matrix <- vote_matrix[member_completeness >= 0.80, , drop = FALSE]
  roll_call_completeness <- colMeans(!is.na(vote_matrix))
  vote_matrix <- vote_matrix[, roll_call_completeness >= 0.80, drop = FALSE]
  if (identical(dim(vote_matrix), old_dimensions)) break
}

member_ids <- rownames(vote_matrix)
members_retained <- members_119 |>
  filter(icpsr %in% as.numeric(member_ids)) |>
  arrange(match(icpsr, as.numeric(member_ids)))

#### Basic dataset checks ----------------------------------------------------
message("Retained ", nrow(vote_matrix), " legislators and ", ncol(vote_matrix), " roll calls.")
message("Overall missing proportion: ", round(mean(is.na(vote_matrix)), 3))
print(table(members_retained$party))

#### Party composition and total Yea/Nay counts -----------------------------
party_counts <- members_retained |>
  count(party)

vote_counts <- data.frame(
  vote = c("Yea", "Nay"),
  count = c(sum(vote_matrix == 1, na.rm = TRUE), sum(vote_matrix == 0, na.rm = TRUE))
)

party_plot <- ggplot(party_counts, aes(x = party, y = n, fill = party)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  labs(title = "Party composition", x = NULL, y = "Number of legislators") +
  theme_minimal()

vote_count_plot <- ggplot(vote_counts, aes(x = vote, y = count, fill = vote)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = c(Yea = "#222222", Nay = "#BDBDBD")) +
  labs(title = "Observed Yea and Nay votes", x = NULL, y = "Number of votes") +
  theme_minimal()

basic_distribution_plot <- party_plot + vote_count_plot + patchwork::plot_layout(ncol = 2)
print(basic_distribution_plot)
ggsave("voteview_basic_distributions.png", basic_distribution_plot, width = 8, height = 4.2, units = "in", dpi = 300)

#### Pairwise voting disagreement ------------------------------------------
number_of_members <- nrow(vote_matrix)
disagreement_matrix <- matrix(0, nrow = number_of_members, ncol = number_of_members, dimnames = list(member_ids, member_ids))

for (member_i in seq_len(number_of_members - 1)) {
  for (member_k in (member_i + 1):number_of_members) {
    jointly_observed <- !is.na(vote_matrix[member_i, ]) & !is.na(vote_matrix[member_k, ])
    if (sum(jointly_observed) >= 5) {
      disagreement_matrix[member_i, member_k] <- 1 - mean(vote_matrix[member_i, jointly_observed] == vote_matrix[member_k, jointly_observed])
      disagreement_matrix[member_k, member_i] <- disagreement_matrix[member_i, member_k]
    } else {
      disagreement_matrix[member_i, member_k] <- NA_real_
      disagreement_matrix[member_k, member_i] <- NA_real_
    }
  }
}

disagreement_values <- disagreement_matrix[upper.tri(disagreement_matrix)]
disagreement_values <- disagreement_values[is.finite(disagreement_values)]
disagreement_summary <- c(mean = mean(disagreement_values), sd = sd(disagreement_values), quantile(disagreement_values, c(0.25, 0.5, 0.75)))
print(round(disagreement_summary, 3))

pairwise_plot <- ggplot(data.frame(disagreement = disagreement_values), aes(x = disagreement)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "#4C78A8", colour = "white") +
  geom_density(colour = "#D62728", linewidth = 0.8) +
  geom_vline(xintercept = mean(disagreement_values), linetype = "dashed", colour = "black") +
  labs(title = "Pairwise voting disagreement", subtitle = paste0("Mean = ", round(mean(disagreement_values), 3), "; median = ", round(median(disagreement_values), 3)), x = "Disagreement", y = "Density") +
  theme_minimal()
print(pairwise_plot)
ggsave("voteview_pairwise_disagreement.png", pairwise_plot, width = 6.8, height = 4.8, units = "in", dpi = 300)
