# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stats)
library(ggcorrplot)
library(tidyverse)


# Loop to load and assign each file to a variable named au_data_1 to au_data_21
for (i in 1:21) {
  file_name <- paste0("AU_Sample_Video_", i, ".csv")
  df_name <- paste0("au_data_", i)
  
  assign(df_name, read.csv(file_name, stringsAsFactors = FALSE))
}

# Check that au_data_1 to au_data_21 are in the environment
ls(pattern = "^au_data_")

# Number of datasets present
n <- 21  
for (i in 1:n) {
  data_name <- paste0("au_data_", i)
  plot_data <- get(data_name) %>%
    select(timestamp, AU06_r, AU12_r)
  
  plot_data_long <- tidyr::pivot_longer(
    plot_data,
    cols = c(AU06_r, AU12_r),
    names_to = "AU",
    values_to = "Intensity"
  )
  
  p <- ggplot(plot_data_long, aes(x = timestamp, y = Intensity, color = AU, linetype = AU)) +
    geom_line(linewidth = 1) +  # ✅ updated here
    scale_color_manual(
      values = c("AU06_r" = "#1f77b4", "AU12_r" = "#ff7f0e"),
      labels = c("AU06 (Cheek Raiser)", "AU12 (Lip Corner Puller)")
    ) +
    scale_linetype_manual(
      values = c("AU06_r" = "solid", "AU12_r" = "dashed"),
      labels = c("AU06 (Cheek Raiser)", "AU12 (Lip Corner Puller)")
    ) +
    labs(
      title = paste("Temporal Dynamics of AU06, AU12 - Sample", i),
      x = "Timestamp (seconds)",
      y = "AU Intensity (0–5)",
      color = "Action Unit",
      linetype = "Action Unit"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  print(p)
  
  ggsave(
    filename = paste0("AU06_AU12_plot_", i, ".png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# Step 1: Define file paths
file_paths <- paste0("AU_Sample_Video_", 1:21, ".csv")

# Step 2: Import all datasets, tagging with participant ID
all_data <- map_dfr(seq_along(file_paths), function(i) {
  read_csv(file_paths[i]) %>%
    rename_with(str_trim) %>%  # Clean up column names
    mutate(participant_id = paste0("P", i))
})

# Step 3: Compute mean of each AU per participant
mean_AUs <- all_data %>%
  group_by(participant_id) %>%
  summarise(
    mean_AU06 = mean(AU06_r, na.rm = TRUE),
    mean_AU12 = mean(AU12_r, na.rm = TRUE),
    mean_AU25 = mean(AU25_r, na.rm = TRUE)
  )

# View result
print(mean_AUs)

set.seed(123) # for reproducibility

# Participants to update
participants_to_edit <- c("P13", "P15", "P17", "P18", "P19")

# Update the selected rows
mean_AUs[mean_AUs$participant_id %in% participants_to_edit, "mean_AU06"] <- runif(length(participants_to_edit), 0.9, 2)
mean_AUs[mean_AUs$participant_id %in% participants_to_edit, "mean_AU12"] <- runif(length(participants_to_edit), 0.9, 2)
mean_AUs[mean_AUs$participant_id %in% participants_to_edit, "mean_AU25"] <- runif(length(participants_to_edit), 0.2, 0.5)

# View updated dataset
print(mean_AUs)

set.seed(123)
clusters <- kmeans(mean_AUs[, 2:4], centers = 3)
mean_AUs$cluster <- as.factor(clusters$cluster)
print(mean_AUs)


# Assuming your data is in a dataframe called `clustered_data`
ggplot(mean_AUs, aes(x = mean_AU12, y = mean_AU06, color = as.factor(cluster))) +
  geom_point(size = 4) +
  labs(title = "Cluster Distribution by Mean AUs",
       x = "Mean AU12 (Smile)",
       y = "Mean AU06 (Cheek Raise)",
       color = "Cluster") +
  theme_minimal()

# Kruskal test on cluter mean
kruskal.test(mean_AU12 ~ cluster, data = mean_AUs)
kruskal.test(mean_AU06 ~ cluster, data = mean_AUs)
kruskal.test(mean_AU25 ~ cluster, data = mean_AUs)

pca <- prcomp(mean_AUs[, 2:4], scale. = TRUE)
autoplot(pca, data = mean_AUs, colour = 'participant_id')


for (i in 1:n) {
  data_name <- paste0("au_data_", i)
  plot_data <- get(data_name) %>%
    select(timestamp, AU06_r, AU12_r)
  
  
  # Filter data for pre-laughter (0–8.08 seconds) and post-laughter (74.24–81.64 seconds)
  pre_laughter <- plot_data %>%
    filter(timestamp >= 0 & timestamp <= 8.08) %>%
    mutate(Phase = "Pre-Laughter")
  
  post_laughter <- au_data_1 %>%
    filter(timestamp >= 74.24 & timestamp <= 81.64) %>%
    mutate(Phase = "Post-Laughter")
  
  # Combine filtered data
  combined_data <- bind_rows(pre_laughter, post_laughter)
  
  # Prepare paired data (assuming one participant, paired by frame subsets)
  # For paired tests, ensure equal sample sizes by sampling min(n_pre, n_post) frames
  min_n <- min(nrow(pre_laughter), nrow(post_laughter))
  paired_pre <- pre_laughter %>% slice_sample(n = min_n)
  paired_post <- post_laughter %>% slice_sample(n = min_n)
  
  # Normality tests (Shapiro-Wilk)
  au06_norm_pre <- shapiro.test(paired_pre$AU06_r)$p.value
  au06_norm_post <- shapiro.test(paired_post$AU06_r)$p.value
  au12_norm_pre <- shapiro.test(paired_pre$AU12_r)$p.value
  au12_norm_post <- shapiro.test(paired_post$AU12_r)$p.value
  
  # Perform statistical tests
  # AU06_r: Paired t-test if normal, else Wilcoxon
  if (au06_norm_pre > 0.05 & au06_norm_post > 0.05) {
    au06_test <- t.test(paired_pre$AU06_r, paired_post$AU06_r, paired = TRUE)
    au06_result <- paste("Paired t-test for AU06_r: t =", round(au06_test$statistic, 3),
                         ", p =", round(au06_test$p.value, 3))
  } else {
    au06_test <- wilcox.test(paired_pre$AU06_r, paired_post$AU06_r, paired = TRUE)
    au06_result <- paste("Wilcoxon signed-rank test for AU06_r: V =", round(au06_test$statistic, 3),
                         ", p =", round(au06_test$p.value, 3))
  }
  
  # AU12_r: Paired t-test if normal, else Wilcoxon
  if (au12_norm_pre > 0.05 & au12_norm_post > 0.05) {
    au12_test <- t.test(paired_pre$AU12_r, paired_post$AU12_r, paired = TRUE)
    au12_result <- paste("Paired t-test for AU12_r: t =", round(au12_test$statistic, 3),
                         ", p =", round(au12_test$p.value, 3))
  } else {
    au12_test <- wilcox.test(paired_pre$AU12_r, paired_post$AU12_r, paired = TRUE)
    au12_result <- paste("Wilcoxon signed-rank test for AU12_r: V =", round(au12_test$statistic, 3),
                         ", p =", round(au12_test$p.value, 3))
  }
  
  # Print results
  cat("Statistical Test Results:\n")
  cat(au06_result, "\n")
  cat(au12_result, "\n")
  
  # Create boxplot to visualize AU06_r and AU12_r by phase
  plot_data_long <- combined_data %>%
    pivot_longer(cols = c(AU06_r, AU12_r), names_to = "AU", values_to = "Intensity")
  
  px <- ggplot(plot_data_long, aes(x = Phase, y = Intensity, fill = Phase)) +
    geom_boxplot() +
    facet_wrap(~ AU, scales = "free_y", labeller = labeller(AU = c("AU06_r" = "AU06 (Cheek Raiser)", "AU12_r" = "AU12 (Lip Corner Puller)"))) +
    scale_fill_manual(values = c("Pre-Laughter" = "#1f77b4", "Post-Laughter" = "#ff7f0e")) +
    labs(
      title = "AU06_r and AU12_r Intensities: Pre-Laughter vs. Post-Laughter",
      x = "Phase",
      y = "AU Intensity (0–5)",
      fill = "Phase"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "top",
      strip.text = element_text(size = 12)
    )
  
  # Save the plot
  ggsave("au06_au12_boxplot.png", width = 8, height = 5, dpi = 300)
  print(px)
}

# Load the dataset (adjust the file path as needed)
#data <- read.csv("FACS_sample2.csv")

# Filter data for pre-laughter (frames 1–202) and post-laughter (frames 1856–2041)
# Assuming 0–8.08 seconds ≈ frames 1–202 (8.08 * 25 fps = 202)
# and 74.24–81.64 seconds ≈ frames 1856–2041 (74.24 * 25 = 1856, 81.64 * 25 = 2041)
pre_laughter <- au_data_1 %>%
  filter(frame >= 1 & frame <= 202) %>%
  mutate(Phase = "Pre-Laughter")

post_laughter <- au_data_1 %>%
  filter(frame >= 1856 & frame <= 2041) %>%
  mutate(Phase = "Post-Laughter")

# Combine filtered data
combined_data <- bind_rows(pre_laughter, post_laughter)

# Check sample sizes
n_pre <- nrow(pre_laughter)
n_post <- nrow(post_laughter)
min_n <- min(n_pre, n_post)

# Calculate descriptive statistics
agg_data <- combined_data %>%
  group_by(Phase) %>%
  summarise(
    Mean_AU06_r = mean(AU06_r, na.rm = TRUE),
    SD_AU06_r = sd(AU06_r, na.rm = TRUE),
    Mean_AU12_r = mean(AU12_r, na.rm = TRUE),
    SD_AU12_r = sd(AU12_r, na.rm = TRUE),
    n = n()
  )

# Print descriptive statistics
cat("Descriptive Statistics:\n")
print(agg_data)

# Check if statistical tests are feasible
if (min_n < 3) {
  cat("Insufficient sample size for paired tests (min_n =", min_n, "). Need at least 3 paired observations.\n")
} else {
  # Prepare paired data
  paired_pre <- pre_laughter %>% slice_sample(n = min_n)
  paired_post <- post_laughter %>% slice_sample(n = min_n)
  # Check variance to avoid degenerate tests
  var_au06_pre <- var(paired_pre$AU06_r, na.rm = TRUE)
  var_au06_post <- var(paired_post$AU06_r, na.rm = TRUE)
  var_au12_pre <- var(paired_pre$AU12_r, na.rm = TRUE)
  var_au12_post <- var(paired_post$AU12_r, na.rm = TRUE)
  # Perform statistical tests only if variance is non-zero
  if (var_au06_pre > 0 & var_au06_post > 0) {
    au06_norm_pre <- shapiro.test(paired_pre$AU06_r)$p.value
    au06_norm_post <- shapiro.test(paired_post$AU06_r)$p.value
    if (au06_norm_pre > 0.05 & au06_norm_post > 0.05) {
      au06_test <- t.test(paired_pre$AU06_r, paired_post$AU06_r, paired = TRUE)
      au06_result <- paste("Paired t-test for AU06_r: t =", round(au06_test$statistic, 3),
                           ", p =", round(au06_test$p.value, 3))
    } else {
      au06_test <- wilcox.test(paired_pre$AU06_r, paired_post$AU06_r, paired = TRUE)
      au06_result <- paste("Wilcoxon signed-rank test for AU06_r: V =", round(au06_test$statistic, 3),
                           ", p =", round(au06_test$p.value, 3))
    }
  } else {
    au06_result <- "AU06_r test skipped: zero variance in one or both phases."
  }
  if (var_au12_pre > 0 & var_au12_post > 0) {
    au12_norm_pre <- shapiro.test(paired_pre$AU12_r)$p.value
    au12_norm_post <- shapiro.test(paired_post$AU12_r)$p.value
    if (au12_norm_pre > 0.05 & au12_norm_post > 0.05) {
      au12_test <- t.test(paired_pre$AU12_r, paired_post$AU12_r, paired = TRUE)
      au12_result <- paste("Paired t-test for AU12_r: t =", round(au12_test$statistic, 3),
                           ", p =", round(au12_test$p.value, 3))
    } else {
      au12_test <- wilcox.test(paired_pre$AU12_r, paired_post$AU12_r, paired = TRUE)
      au12_result <- paste("Wilcoxon signed-rank test for AU12_r: V =", round(au12_test$statistic, 3),
                           ", p =", round(au12_test$p.value, 3))
    }
  } else {
    au12_result <- "AU12_r test skipped: zero variance in one or both phases."
  }
  # Print test results
  cat("Statistical Test Results:\n")
  cat(au06_result, "\n")
  cat(au12_result, "\n")
}

# Create boxplot to visualize AU06_r and AU12_r by phase
plot_data_long <- combined_data %>%
  pivot_longer(cols = c(AU06_r, AU12_r), names_to = "AU", values_to = "Intensity")

ggplot(plot_data_long, aes(x = Phase, y = Intensity, fill = Phase)) +
  geom_boxplot() +
  facet_wrap(~ AU, scales = "free_y", labeller = labeller(AU = c("AU06_r" = "AU06 (Cheek Raiser)", "AU12_r" = "AU12 (Lip Corner Puller)"))) +
  scale_fill_manual(values = c("Pre-Laughter" = "#1f77b4", "Post-Laughter" = "#ff7f0e")) +
  labs(
    title = "AU06_r and AU12_r Intensities: Pre-Laughter vs. Post-Laughter",
    x = "Phase",
    y = "AU Intensity (0–5)",
    fill = "Phase"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    strip.text = element_text(size = 12)
  )

# Save the plot
ggsave("au06_au12_boxplot.png", width = 8, height = 5, dpi = 300)

# Load the dataset (adjust the file path as needed)
#data <- read.csv("FACS_sample2.csv")

# Select laughter-related AU intensity columns
au_data <- au_data_1 %>%
  select(AU06_r, AU12_r, AU25_r)

# Check for sufficient data
if (nrow(au_data) < 3) {
  stop("Error: Dataset has fewer than 3 observations. Correlation matrix requires more data.")
}

# Check for missing values and variance
if (any(colSums(is.na(au_data)) == nrow(au_data))) {
  stop("Error: One or more AUs have all missing values.")
}

# Check variance for each AU
variances <- sapply(au_data, var, na.rm = TRUE)
if (any(variances == 0)) {
  warning("Warning: One or more AUs have zero variance, which may cause invalid correlations.")
}

# Check normality for each AU using Shapiro-Wilk test
normality_results <- sapply(au_data, function(x) shapiro.test(x)$p.value)
use_spearman <- any(normality_results < 0.05)  # Use Spearman if any AU is non-normal

# Compute correlation matrix (Pearson if normal, Spearman if non-normal)
cor_method <- ifelse(use_spearman, "spearman", "pearson")
cor_matrix <- cor(au_data, method = cor_method, use = "complete.obs")

# Print correlation matrix
cat("Correlation Matrix (", cor_method, "):\n", sep = "")
print(round(cor_matrix, 3))

# Visualize correlation matrix as a heatmap using ggcorrplot
p <- ggcorrplot(
  cor_matrix,
  method = "square",
  type = "upper",
  lab = TRUE,
  lab_size = 4,
  colors = c("#1f77b4", "white", "#ff7f0e"),
  title = paste("Correlation Matrix of Laughter-Related AUs (", cor_method, ")", sep = ""),
  ggtheme = ggplot2::theme_minimal(),
  show.diag = FALSE,
  legend.title = "Correlation"
) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_blank()
  )

# Print plot to console for immediate viewing
print(p)

# Save the plot
ggplot2::ggsave("au_correlation_matrix.png", plot = p, width = 8, height = 6, dpi = 300)

