# scripts/analysis_smk.R
#
# Combined analysis script: connectivity + top-CC extraction + DICE plots.
# DICE datasets are driven entirely by snakemake@input — no hardcoding.

library(tidyverse)
library(ggridges)
library(scales)

# ── Inputs / outputs / params ─────────────────────────────────────────────────

connectivity_csv  <- snakemake@input[["connectivity"]]
out_dir           <- snakemake@output[["out_dir"]]

K_REFINED         <- as.integer(snakemake@params[["k_refined"]])
K_COMPARE_SCATTER <- as.integer(snakemake@params[["k_compare"]])
K_COMPARE_TOP     <- as.integer(snakemake@params[["k_top"]])
TOP_N             <- as.integer(snakemake@params[["top_n"]])

# All inputs except "connectivity" are dice CSVs, keyed by dataset slug.
# Slugs use "_" as word separator and "__" as dataset/other_dir separator.
# e.g. "Fundus-AVSeg" -> "Fundus_AVSeg"
#      "leuven-haifa/train" -> "leuven_haifa__train"
dice_inputs <- snakemake@input[names(snakemake@input) != "connectivity"]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — CONNECTIVITY
# ══════════════════════════════════════════════════════════════════════════════

PALETTE <- c("0" = "#888888", "3" = "#d9f0a3", "4" = "#addd8e",
             "5" = "#78c679", "6" = "#41ab5d", "7" = "#238443", "8" = "#005a32")

df_conn <- read_csv(connectivity_csv, show_col_types = FALSE) |>
  distinct() |>
  mutate(
    k        = factor(k, levels = c(0, K_REFINED)),
    k_label  = if_else(k == "0", "Unrefined\n(k=0)", paste0("k=", k)),
    k_label  = factor(k_label, levels = c("Unrefined\n(k=0)", paste0("k=", K_REFINED))),
    ds_label = if_else(is.na(other_dir), dataset, paste0(dataset, "/", other_dir)),
    refined  = k != "0"
  )

cat("Connectivity rows after dedup:", nrow(df_conn), "\n")

# 1. Density — all datasets pooled
p1 <- df_conn |>
  ggplot(aes(x = num_components, fill = k, colour = k)) +
  geom_density(alpha = 0.45, linewidth = 0.5) +
  scale_fill_manual(values = PALETTE, name = "k") +
  scale_colour_manual(values = PALETTE, name = "k") +
  labs(title    = "Distribution of connected components — all datasets pooled",
       subtitle = "k=0 is unrefined",
       x = "Number of connected components", y = "Density") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "01_density_cc_all.pdf"), p1, width = 9, height = 5)

# 2. Density — faceted by dataset
p2 <- df_conn |>
  ggplot(aes(x = num_components, fill = k, colour = k)) +
  geom_density(alpha = 0.45, linewidth = 0.4) +
  facet_wrap(~ds_label, scales = "free_y") +
  scale_fill_manual(values = PALETTE, name = "k") +
  scale_colour_manual(values = PALETTE, name = "k") +
  labs(title = "Connected components distribution — per dataset",
       x = "Number of connected components", y = "Density") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "02_density_cc_per_dataset.pdf"), p2, width = 14, height = 9)

# 3. Ridge plot
p3 <- df_conn |>
  ggplot(aes(x = num_components, y = fct_rev(k_label), fill = k)) +
  geom_density_ridges(alpha = 0.75, scale = 1.2,
                      quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  labs(title    = "Connected components — ridge plot across k values",
       subtitle = "Vertical line = median",
       x = "Number of connected components", y = NULL)
ggsave(file.path(out_dir, "03_ridge_cc.pdf"), p3, width = 8, height = 6)

# 4. Violin — all datasets pooled
p4 <- df_conn |>
  ggplot(aes(x = k_label, y = num_components, fill = k)) +
  geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", linewidth = 0.5) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  labs(title = "Connected components before and after refinement — all datasets",
       x = NULL, y = "Number of connected components") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "04_violin_cc_all.pdf"), p4, width = 10, height = 6)

# 5. Violin — faceted by dataset
p5 <- df_conn |>
  ggplot(aes(x = k_label, y = num_components, fill = k)) +
  geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  facet_wrap(~ds_label, scales = "free_y") +
  scale_fill_manual(values = PALETTE, guide = "none") +
  labs(title = "Connected components — per dataset",
       x = NULL, y = "Number of connected components") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "05_violin_cc_per_dataset.pdf"), p5, width = 15, height = 10)

# 6. Mean CC curve per dataset (±1 SE ribbon)
summary_k <- df_conn |>
  group_by(ds_label, k, k_label) |>
  summarise(mean_cc = mean(num_components, na.rm = TRUE),
            se_cc   = sd(num_components,   na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

p6 <- summary_k |>
  ggplot(aes(x = k_label, y = mean_cc, colour = ds_label, group = ds_label)) +
  geom_ribbon(aes(ymin = mean_cc - se_cc, ymax = mean_cc + se_cc, fill = ds_label),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(title    = "Mean connected components vs refinement level",
       subtitle = "Ribbon = ±1 SE",
       x = NULL, y = "Mean number of connected components",
       colour = "Dataset", fill = "Dataset") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "06_mean_cc_curve.pdf"), p6, width = 10, height = 6)

# 7. Violin — OD connectivity
p7 <- df_conn |>
  ggplot(aes(x = k_label, y = proportion_length_connected, fill = k)) +
  geom_violin(trim = TRUE, alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Proportion of vessel length connected to optic disc",
       x = NULL, y = "Proportion of total length") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "07_violin_od_connectivity.pdf"), p7, width = 10, height = 6)

# 8. Scatter — delta CC vs delta OD connectivity (k=0 → K_COMPARE_SCATTER)
paired_scatter <- inner_join(
  df_conn |> filter(k == "0") |>
    select(ds_label, image_id, cc_before = num_components,
           od_before = proportion_length_connected),
  df_conn |> filter(k == as.character(K_COMPARE_SCATTER)) |>
    select(ds_label, image_id, cc_after = num_components,
           od_after = proportion_length_connected),
  by = c("ds_label", "image_id")
) |>
  mutate(delta_cc = cc_after - cc_before,
         delta_od = od_after  - od_before)

p8 <- paired_scatter |>
  ggplot(aes(x = delta_cc, y = delta_od, colour = ds_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(labels = percent_format()) +
  labs(title    = paste0("Per-image changes: k=0 → k=", K_COMPARE_SCATTER),
       subtitle = "Bottom-left = fewer CCs AND better OD connectivity",
       x = "Change in number of connected components",
       y = "Change in proportion of length connected to OD",
       colour = "Dataset") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "08_scatter_delta_cc_vs_od.pdf"), p8, width = 9, height = 7)

# 9. Summary table
summary_table <- df_conn |>
  group_by(ds_label, k) |>
  summarise(
    n                        = n(),
    num_components_mean      = mean(num_components, na.rm = TRUE),
    num_components_median    = median(num_components, na.rm = TRUE),
    num_components_sd        = sd(num_components, na.rm = TRUE),
    prop_length_od_mean      = mean(proportion_length_connected, na.rm = TRUE),
    prop_length_od_sd        = sd(proportion_length_connected, na.rm = TRUE),
    largest_cc_prop_len_mean = mean(largest_component_proportion_length, na.rm = TRUE),
    largest_cc_prop_len_sd   = sd(largest_component_proportion_length, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(ds_label, k)
write_csv(summary_table, file.path(out_dir, "summary_connectivity.csv"))

# 10. Paired Wilcoxon — each k vs k=0, per dataset
wilcoxon_results <- df_conn |>
  group_by(ds_label) |>
  group_modify(function(data, key) {
    map_dfr(K_REFINED, function(k_val) {
      matched <- inner_join(
        data |> filter(k == "0")                 |> select(image_id, before = num_components),
        data |> filter(k == as.character(k_val)) |> select(image_id, after  = num_components),
        by = "image_id"
      )
      if (nrow(matched) < 5) return(tibble())
      w <- wilcox.test(matched$before, matched$after, paired = TRUE)
      tibble(
        k             = k_val,
        n_pairs       = nrow(matched),
        median_before = median(matched$before),
        median_after  = median(matched$after),
        delta_median  = median(matched$after) - median(matched$before),
        p_value       = w$p.value,
        signif        = case_when(
          w$p.value < 0.001 ~ "***",
          w$p.value < 0.01  ~ "**",
          w$p.value < 0.05  ~ "*",
          TRUE              ~ "ns"
        )
      )
    })
  }) |>
  ungroup()

write_csv(wilcoxon_results, file.path(out_dir, "wilcoxon_num_components.csv"))

# 4b. Violin with Wilcoxon annotations
wilcox_labels <- wilcoxon_results |>
  filter(ds_label == wilcoxon_results$ds_label[1]) |>
  mutate(k_label = paste0("k=", k),
         y_pos   = max(df_conn$num_components, na.rm = TRUE) * 1.05)

p4b <- df_conn |>
  ggplot(aes(x = k_label, y = num_components, fill = k)) +
  geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_text(data = wilcox_labels,
            aes(x = k_label, y = y_pos, label = signif),
            inherit.aes = FALSE, size = 5, vjust = -0.5) +
  scale_fill_manual(values = PALETTE, guide = "none") +
  labs(title    = "Connected components before and after refinement",
       subtitle = "Stars = significance vs k=0 (paired Wilcoxon)",
       x = NULL, y = "Number of connected components") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "04b_violin_cc_all_wilcox.pdf"), p4b, width = 10, height = 6)

cat("Section 1 (connectivity) done.\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — TOP CC REDUCTIONS
# ══════════════════════════════════════════════════════════════════════════════

paired_top <- df_conn |>
  filter(k %in% c("0", as.character(K_COMPARE_TOP))) |>
  select(dataset, other_dir, image_id, k, num_components) |>
  pivot_wider(names_from = k, values_from = num_components, names_prefix = "cc_k") |>
  drop_na() |>
  mutate(
    cc_reduction      = cc_k0 - .data[[paste0("cc_k", K_COMPARE_TOP)]],
    percent_reduction = (cc_reduction / cc_k0) * 100,
    ds_label          = if_else(is.na(other_dir), dataset, paste0(dataset, "/", other_dir))
  )

top_images <- paired_top |>
  group_by(ds_label) |>
  slice_max(order_by = cc_reduction, n = TOP_N, with_ties = FALSE) |>
  ungroup() |>
  arrange(ds_label, desc(cc_reduction))

write_csv(top_images, file.path(out_dir, "top_cc_reduction.csv"))
cat("Section 2 (top CC reductions) done.\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — DICE PLOTS (fully dynamic)
# ══════════════════════════════════════════════════════════════════════════════

color_map  <- c("macro" = "#9467BD", "A" = "#E84040", "V" = "#4080E8")
color_labs <- c("macro" = "Macro avg", "A" = "Artery", "V" = "Vein")

load_dice <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(image, matches("^dice_(unref|k\\d+)_vs_gt_(macro|A|V)$")) %>%
    pivot_longer(cols = -image,
                 names_to     = c("iter", "class"),
                 names_pattern = "dice_(.+)_vs_gt_(.+)") %>%
    mutate(
      k       = case_when(iter == "unref" ~ 0L,
                          TRUE ~ as.integer(str_extract(iter, "\\d+"))),
      k_label = ifelse(k == 0, "unref", paste0("k=", k)),
      k_label = factor(k_label, levels = unique(k_label[order(k)]))
    )
}

k_scale <- function(ks) {
  scale_x_continuous(breaks = ks,
                     labels = ifelse(ks == 0, "unref", paste0("k=", ks)))
}

plot_line <- function(df_long, title) {
  df_sum <- df_long %>%
    group_by(k, class) %>%
    summarise(mean_dice = mean(value, na.rm = TRUE),
              sd_dice   = sd(value,   na.rm = TRUE), .groups = "drop")
  ggplot(df_sum, aes(x = k, y = mean_dice, color = class, group = class)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_dice - sd_dice, ymax = mean_dice + sd_dice),
                  width = 0.2, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
    annotate("text", x = 0.15, y = min(df_sum$mean_dice),
             label = "unrefined", hjust = 0, size = 3, color = "grey40") +
    k_scale(unique(df_sum$k)) +
    scale_color_manual(values = color_map, labels = color_labs) +
    labs(title = title, x = "Refinement iteration (k)",
         y = "Mean DICE ± SD", color = "Class") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
}

plot_box <- function(df_long, title) {
  ggplot(df_long, aes(x = k_label, y = value, fill = class)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, linewidth = 0.4) +
    scale_fill_manual(values = color_map, labels = color_labs) +
    labs(title = title, x = "Refinement iteration (k)",
         y = "DICE score", fill = "Class") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
}

plot_violin <- function(df_long, title) {
  ggplot(df_long, aes(x = k_label, y = value, fill = class)) +
    geom_violin(trim = TRUE, linewidth = 0.3, alpha = 0.8) +
    geom_boxplot(width = 0.08, outlier.shape = NA, linewidth = 0.4,
                 position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = color_map, labels = color_labs) +
    labs(title = title, x = "Refinement iteration (k)",
         y = "DICE score", fill = "Class") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
}

for (slug in names(dice_inputs)) {
  path     <- dice_inputs[[slug]]
  # "__" separates dataset from other_dir in the slug; "-" used within names
  ds_title <- gsub("__", "/", gsub("_", "-", slug))

  cat(sprintf("Processing DICE for: %s\n", ds_title))
  df_long <- load_dice(path)

  if (nrow(df_long) == 0 || !any(!is.na(df_long$value))) {
    cat(sprintf("  Skipping %s: no valid DICE values.\n", ds_title))
    next
  }

  plots <- list(
    line   = plot_line  (df_long, sprintf("%s — DICE vs refinement iteration (mean ± SD)", ds_title)),
    box    = plot_box   (df_long, sprintf("%s — DICE distribution per iteration (boxplot)",  ds_title)),
    violin = plot_violin(df_long, sprintf("%s — DICE distribution per iteration (violin)",   ds_title))
  )

  for (plot_type in names(plots)) {
    ggsave(file.path(out_dir, sprintf("%s_dice_%s.png", slug, plot_type)),
           plot = plots[[plot_type]], width = 10, height = 5, dpi = 150)
  }
}

cat("All done. Outputs in:", out_dir, "\n")
