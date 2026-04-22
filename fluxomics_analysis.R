# ============================================================
# Metabolomics Flux Analysis — 13C6-Glucose Tracing
# ============================================================
# Author : Agnieszka Karaś
# GitHub : github.com/agkaras
# ============================================================
#
# Experimental design:
#   Murine thoracic aorta 
#   2 h 13C6-glucose labelling, 24 h total incubation
#   Factors: Age (young / old) × Treatment (CTR / IL-1β)
#
# Input : Excel file from FIMM metabolomics platform
# Output: Per-metabolite PNG plots + multi-panel PNG figure
# ============================================================

# ============================================================
# CONFIGURATION — edit these before running
# ============================================================
DATA_FILE  <- "your_file.xlsx"          # <-- replace with your filename
OUTPUT_DIR <- "plots_individual"
PANEL_PNG  <- "metabolomics_flux_panel.png"
PLOT_W     <- 10    # panel figure width, inches
PLOT_H     <- 8     # panel figure height, inches

# ---- 1. INSTALL / LOAD PACKAGES ----
# install.packages(c("readxl", "tidyverse", "ggplot2", "ggbeeswarm",
#                    "rstatix", "ggpubr", "patchwork"))

library(readxl)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
library(patchwork)

# ---- 2. LOAD DATA ----
raw <- read_excel(DATA_FILE,
                  sheet     = "Results normalized TIC",
                  skip      = 14,
                  col_names = TRUE)

# ---- 3. DEFINE SAMPLE GROUPS ----
# Column names must match exactly what is in your Excel file.
# Run colnames(raw) after loading to verify.
# The names below correspond to the dataset used in Karaś et al. (2026)
# Free Radical Biology and Medicine, doi:10.1016/j.freeradbiomed.2025.12.046
samples_0h         <- c("32_0h_young", "33_0h_young", "34_0h_old")

samples_ctr_young  <- c("1_young_CTR",  "3_young_CTR",  "5_young_CTR",
                         "14_young_CTR", "16_young_CTR", "18_young_CTR")

samples_ctr_old    <- c("2_old_CTR",  "4_old_CTR",  "6_old_CTR",
                         "15_old_CTR", "17_old_CTR", "19_old_CTR")

samples_il1b_young <- c("7_young_IL1b",  "9_young_IL1b",  "11_young_IL1b",
                          "20_young_IL1b", "22_young_IL1b", "24_young_IL1b")

samples_il1b_old   <- c("8_old_IL1b",  "10_old_IL1b",  "12_old_IL1b",
                          "21_old_IL1b", "23_old_IL1b", "25_old_IL1b")

# ---- 4. CLEAN & TIDY DATA ----
data_clean <- raw %>%
  rename(Compound = Compounds,
         Flux     = FLUX) %>%
  filter(!is.na(Compound), Compound != "Compounds") %>%
  mutate(across(all_of(c(samples_0h, samples_ctr_young, samples_ctr_old,
                          samples_il1b_young, samples_il1b_old)),
                ~ suppressWarnings(as.numeric(.x))))

make_long <- function(df, sample_cols, group_label) {
  df %>%
    select(Pathway, Compound, Flux, all_of(sample_cols)) %>%
    pivot_longer(cols      = all_of(sample_cols),
                 names_to  = "SampleID",
                 values_to = "Intensity") %>%
    mutate(Group = group_label)
}

GROUP_LEVELS <- c("CTR young", "IL-1\u03b2 young", "CTR old", "IL-1\u03b2 old")

df_long <- bind_rows(
  make_long(data_clean, samples_0h,         "0h baseline"),
  make_long(data_clean, samples_ctr_young,  "CTR young"),
  make_long(data_clean, samples_ctr_old,    "CTR old"),
  make_long(data_clean, samples_il1b_young, "IL-1\u03b2 young"),
  make_long(data_clean, samples_il1b_old,   "IL-1\u03b2 old")
) %>%
  mutate(Group = factor(Group, levels = c("0h baseline", GROUP_LEVELS)))

# ---- 5. COLOUR PALETTE ----
group_colours <- c(
  "0h baseline"       = "#6aada0",
  "CTR young"         = "#458B74",
  "IL-1\u03b2 young"  = "#00688B",
  "CTR old"           = "#5D478B",
  "IL-1\u03b2 old"    = "#8B3A3A"
)

# ---- 6. THEME ----
theme_metabolomics <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "grey40", fill = NA),
    axis.text         = element_text(colour = "black", size = 9),
    axis.title        = element_text(colour = "black", size = 10),
    axis.title.y      = element_text(margin = margin(r = 6)),
    plot.title        = element_text(colour = "black", size = 11,
                                     face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(colour = "grey50", size = 7,
                                     hjust = 0.5),
    plot.tag          = element_text(face = "bold", size = 13),
    legend.position   = "none"
  )
}

# ---- 7. PLOTTING FUNCTION ----
# Single panel — 4 boxplots (CTR young / IL-1β young / CTR old / IL-1β old).
#
# X-axis labels: "young\nCTR", "young\nIL-1β", "old\nCTR", "old\nIL-1β"
#   set via scale_x_discrete(labels = x_labs).
#
# Boxplot: median, IQR box, whiskers with end-caps (stat_boxplot geom="errorbar").
# Fill uses alpha = 0.35 so individual points stay visible against the box.
# Border, whiskers and points all use the full group colour.
#
# Stats: one-way ANOVA Šidák post-hoc contrasts (CTR vs IL-1β within each age
#        group). Pre-computed FDRs are passed in via stat_data (results_full).
#        Bracket y-position based on max raw Intensity so it clears all points.
#        Group NAMES used for xmin/xmax — unambiguous on a discrete axis.

plot_metabolite <- function(df_long,
                            metabolite_name,
                            y_label     = "Peak Area normalized to TIC",
                            groups_plot = GROUP_LEVELS,
                            add_stats   = TRUE,
                            stat_data   = NULL,   # pass results_full here
                            min_n_stats = 2) {

  if (!metabolite_name %in% unique(df_long$Compound))
    stop(paste0("Metabolite not found: '", metabolite_name, "'"))

  # ---- Subset ----
  df_sub <- df_long %>%
    filter(Compound == metabolite_name, Group %in% groups_plot) %>%
    mutate(
      Group = factor(Group, levels = groups_plot),
      treatment = factor(
        ifelse(grepl("CTR", Group), "CTR", "IL-1\u03b2"),
        levels = c("CTR", "IL-1\u03b2")
      ),
      age = factor(
        ifelse(grepl("young", Group), "young", "old"),
        levels = c("young", "old")
      )
    )

  if (nrow(df_sub) == 0 || all(is.na(df_sub$Intensity))) {
    message("Skipping '", metabolite_name, "': no data in target groups.")
    return(invisible(NULL))
  }

  # ---- Summary (for n check) ----
  df_sum <- df_sub %>%
    group_by(Group, treatment, age) %>%
    summarise(
      n_obs = sum(!is.na(Intensity)),
      .groups = "drop"
    )

  # ---- x-axis labels ----
  x_labs <- c(
    "CTR young"        = "young\nCTR",
    "IL-1\u03b2 young" = "young\nIL-1\u03b2",
    "CTR old"          = "old\nCTR",
    "IL-1\u03b2 old"   = "old\nIL-1\u03b2"
  )

  # ---- Build plot ----
  p <- ggplot(df_sub, aes(x = Group, y = Intensity,
                           fill = Group, colour = Group)) +

    stat_boxplot(
      geom = "errorbar", width = 0.28, linewidth = 0.5, na.rm = TRUE
    ) +

    geom_boxplot(
      width = 0.55, linewidth = 0.5,
      outlier.shape = NA, alpha = 0.35, na.rm = TRUE
    ) +

    geom_quasirandom(
      size = 2.3, alpha = 0.9, width = 0.17, na.rm = TRUE
    ) +

    scale_fill_manual(values   = group_colours, drop = FALSE) +
    scale_colour_manual(values = group_colours, drop = FALSE) +
    scale_x_discrete(labels = x_labs) +

    labs(title = metabolite_name, x = NULL, y = y_label) +
    theme_metabolomics()

  # ---- Statistics: ANOVA Šidák post-hoc from pre-computed results_full ----
  if (add_stats && !is.null(stat_data)) {

    met_row <- stat_data %>%
      filter(Compound == metabolite_name)

    if (nrow(met_row) == 1) {

      # 3 contrasts — raw Šidák p-values (corrected within metabolite).
      # Use p_* columns (not FDR_*) to match manual pairwise comparisons.
      all_brackets <- tibble(
        p_col  = c("p_young: CTR vs IL-1b",
                   "p_old: CTR vs IL-1b",
                   "p_CTR: young vs old"),
        group1 = c("CTR young",        "CTR old",   "CTR young"),
        group2 = c("IL-1\u03b2 young", "IL-1\u03b2 old", "CTR old")
      ) %>%
        rowwise() %>%
        mutate(pval = met_row[[p_col]]) %>%
        ungroup() %>%
        filter(!is.na(pval) & pval < 0.05) %>%
        mutate(
          p.signif = case_when(
            pval < 0.001 ~ "***",
            pval < 0.01  ~ "**",
            TRUE         ~ "*"
          )
        )

      if (nrow(all_brackets) > 0) {

        # Base y: highest data point in the whole plot
        y_max <- max(df_sub$Intensity, na.rm = TRUE)
        step  <- y_max * 0.12   # vertical gap between stacked brackets

        # Sort brackets by "span width" (wider brackets go higher)
        # so short brackets sit close to the data and long ones are above.
        group_pos <- setNames(seq_along(GROUP_LEVELS), GROUP_LEVELS)
        all_brackets <- all_brackets %>%
          mutate(
            x1   = group_pos[group1],
            x2   = group_pos[group2],
            span = abs(x2 - x1)
          ) %>%
          arrange(span) %>%
          mutate(y.position = y_max * 1.10 + (row_number() - 1) * step) %>%
          select(-x1, -x2, -span, -p_col, -pval)

        p <- p + stat_pvalue_manual(
          all_brackets,
          label      = "p.signif",
          tip.length = 0.01,
          colour     = "grey30",
          size       = 3.5
        )

        # Expand top margin to fit all brackets — added here, once, after stats
        p <- p + scale_y_continuous(
          expand = expansion(mult = c(0.02, 0.15 + nrow(all_brackets) * 0.07))
        )

      } else {
        # No significant brackets — add default scale
        p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
      }
    } else {
      p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
    }
  } else {
    # add_stats = FALSE or no stat_data
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
  }

  return(p)
}


# ============================================================
# ---- 8. ONE-WAY ANOVA + SIDAK POST-HOC ----
# ============================================================
#
# Computed here so results are available for both the boxplot
# significance brackets (sections 9-10) and the heatmaps (12-13).
#
# Model: Intensity ~ Group  (3 groups as a single factor)
# Post-hoc: 3 a priori contrasts with Šidák correction:
#   1. young CTR vs young IL-1β
#   2. old CTR   vs old IL-1β
#   3. young CTR vs old CTR
# FDR (BH) applied across metabolites per test.
# ============================================================

# install.packages("emmeans")
library(emmeans)

# ---- Helper: one-way ANOVA + Sidak post-hoc for one metabolite ----
run_anova <- function(compound, df) {

  dat <- df %>%
    filter(Compound == compound, Group %in% GROUP_LEVELS,
           !is.na(Intensity)) %>%
    mutate(Group = factor(Group, levels = GROUP_LEVELS))

  # Need >= 2 observations per group
  grp_n <- dat %>% group_by(Group) %>% summarise(n = n(), .groups = "drop")
  if (any(grp_n$n < 2)) return(NULL)

  fit <- tryCatch(
    aov(Intensity ~ Group, data = dat),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # ---- Omnibus F and p ----
  av <- summary(fit)[[1]]
  anova_row <- tibble(
    Compound = compound,
    F_value  = av["Group", "F value"],
    p_value  = av["Group", "Pr(>F)"]
  )

  # ---- Post-hoc: 3 contrasts with Šidák correction ----
  # emmeans follows GROUP_LEVELS order:
  #   [1] CTR young, [2] IL-1b young, [3] CTR old, [4] IL-1b old
  emm <- tryCatch(emmeans(fit, ~ Group), error = function(e) NULL)
  if (is.null(emm)) return(list(anova = anova_row, contrasts = NULL))

  contrasts_res <- tryCatch({
    contrast(emm,
             method = list(
               "young: CTR vs IL-1b" = c( 1, -1,  0,  0),
               "old: CTR vs IL-1b"   = c( 0,  0,  1, -1),
               "CTR: young vs old"   = c( 1,  0, -1,  0)
             ),
             adjust = "sidak") %>%
      as_tibble() %>%
      select(contrast, p.value) %>%    # p.value = Šidák-corrected within metabolite
      mutate(Compound = compound)
  }, error = function(e) NULL)

  list(anova = anova_row, contrasts = contrasts_res)
}

# ---- Run for all metabolites ----
cat("\nRunning one-way ANOVA for all metabolites...\n")
compounds_all <- unique(df_long$Compound[df_long$Group %in% GROUP_LEVELS])

anova_list <- lapply(compounds_all, function(cmp) {
  tryCatch(run_anova(cmp, df_long), error = function(e) NULL)
})
names(anova_list) <- compounds_all

# ---- Combine ANOVA results + FDR across metabolites ----
anova_df <- bind_rows(lapply(anova_list, `[[`, "anova")) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

# ---- Combine contrast results: raw Šidák p + FDR across metabolites ----
contrast_raw <- bind_rows(lapply(anova_list, `[[`, "contrasts")) %>%
  # Raw Šidák p (corrected within metabolite across 3 contrasts)
  pivot_wider(names_from  = contrast,
              values_from = p.value,
              names_prefix = "p_")

contrast_fdr <- bind_rows(lapply(anova_list, `[[`, "contrasts")) %>%
  group_by(contrast) %>%
  mutate(FDR_contrast = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  select(Compound, contrast, FDR_contrast) %>%
  pivot_wider(names_from  = contrast,
              values_from = FDR_contrast,
              names_prefix = "FDR_")

# ---- Merge and save full results CSV ----
# Columns: F, p (omnibus), FDR (omnibus),
#          p_* (Šidák within metabolite — use for plot brackets),
#          FDR_* (BH across metabolites — use for heatmap)
results_full <- left_join(anova_df, contrast_raw,  by = "Compound") %>%
  left_join(contrast_fdr, by = "Compound") %>%
  arrange(FDR)
write.csv(results_full, "anova_results.csv", row.names = FALSE)
cat("Saved: anova_results.csv\n")

# ---- Prepare significance heatmap ----
FDR_THRESHOLD <- 0.05
LOG_CAP        <- 10

hm_data <- results_full %>%
  transmute(
    Compound,
    "ANOVA"                     = FDR,
    "young\nCTR vs IL-1\u03b2" = `FDR_young: CTR vs IL-1b`,
    "old\nCTR vs IL-1\u03b2"   = `FDR_old: CTR vs IL-1b`,
    "CTR\nyoung vs old"         = `FDR_CTR: young vs old`
  ) %>%
  filter(if_any(-Compound, ~ !is.na(.) & . < FDR_THRESHOLD))

col_order <- c("ANOVA",
               "young\nCTR vs IL-1\u03b2", "old\nCTR vs IL-1\u03b2",
               "CTR\nyoung vs old")

row_order <- hm_data %>% arrange(ANOVA) %>% pull(Compound)

hm_long <- hm_data %>%
  pivot_longer(-Compound, names_to = "Effect", values_to = "FDR") %>%
  mutate(
    Effect    = factor(Effect, levels = col_order),
    Compound  = factor(Compound, levels = rev(row_order)),
    log_fdr   = pmin(-log10(replace_na(FDR, 1)), LOG_CAP),
    sig_label = case_when(
      !is.na(FDR) & FDR < 0.001 ~ "***",
      !is.na(FDR) & FDR < 0.01  ~ "**",
      !is.na(FDR) & FDR < 0.05  ~ "*",
      TRUE                       ~ ""
    )
  )

heatmap_sig <- ggplot(hm_long,
                      aes(x = Effect, y = Compound, fill = log_fdr)) +

  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label),
            size = 2.8, colour = "white", fontface = "bold", vjust = 0.75) +

  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +

  scale_fill_gradient2(
    low = "white", mid = "#4f93d4", high = "#8B1A1A",
    midpoint = LOG_CAP / 2, limits = c(0, LOG_CAP),
    name = expression(-log[10](FDR)), na.value = "grey92"
  ) +
  scale_x_discrete(position = "top") +

  labs(
    title    = "One-way ANOVA: 4 groups",
    subtitle = paste0("FDR < ", FDR_THRESHOLD,
                      "  |  post-hoc: Sid\u00e1k, FDR across metabolites"),
    x = NULL, y = NULL
  ) +

  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border    = element_blank(),
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 8, colour = "black",
                                   angle = 45, hjust = 0, vjust = 0),
    axis.text.y     = element_text(size = 7,  colour = "black"),
    axis.ticks      = element_blank(),
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle   = element_text(size = 8, colour = "grey50", hjust = 0.5),
    plot.margin     = margin(30, 10, 10, 10),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7)
  )

ggsave("heatmap_anova.png",
       plot   = heatmap_sig,
       width  = 7,
       height = max(5, nrow(hm_data) * 0.22 + 2.5),
       dpi    = 300, bg = "white")
cat("Saved: heatmap_anova.png (", nrow(hm_data), "metabolites)\n")

# ---- 9. EXAMPLE: 4-PANEL FIGURE ----
# Replace with compound names from your dataset.
# Run unique(data_clean$Compound) to see all available names.
metabolites_to_plot <- c(
  "Compound_1",   # e.g. "009_Pyruvate +3"
  "Compound_2",
  "Compound_3",
  "Compound_4"
)

available <- intersect(metabolites_to_plot, unique(data_clean$Compound))
cat("Metabolites found:", paste(available, collapse = ", "), "\n")

missing <- setdiff(metabolites_to_plot, unique(data_clean$Compound))
if (length(missing) > 0)
  warning("NOT found (check spelling): ", paste(missing, collapse = ", "))

panel_list <- lapply(available, function(m) {
  tryCatch(
    plot_metabolite(df_long, m, stat_data = results_full),
    error = function(e) { message("Error plotting ", m, ": ", e$message); NULL }
  )
})
panel_list <- Filter(Negate(is.null), panel_list)

if (length(panel_list) >= 2) {
  combined_plot <- wrap_plots(panel_list, ncol = 2) +
    plot_annotation(
      title      = "13C6-Glucose Flux \u2013 Murine Thoracic Aorta",
      subtitle   = "Young vs Old \u00d7 CTR vs IL-1\u03b2  |  2 h labelling, 24 h total",
      tag_levels = "A",
      theme = theme(
        plot.title    = element_text(colour = "black", size = 14,
                                     face = "bold", hjust = 0.5),
        plot.subtitle = element_text(colour = "grey50", size = 10,
                                     hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    ) &
    theme(plot.tag = element_text(face = "bold", size = 13))

  ggsave(PANEL_PNG, plot = combined_plot,
         width = PLOT_W, height = PLOT_H, dpi = 300, bg = "white")
  cat("Saved:", PANEL_PNG, "\n")
}

# ---- 10. LOOP: SAVE EVERY METABOLITE AS INDIVIDUAL PNG ----
dir.create(OUTPUT_DIR, showWarnings = FALSE)

for (met in unique(data_clean$Compound)) {
  p <- tryCatch(
    plot_metabolite(df_long, met, stat_data = results_full),
    error = function(e) { message("Error for '", met, "': ", e$message); NULL }
  )
  if (!is.null(p)) {
    tryCatch(
      ggsave(file.path(OUTPUT_DIR,
                       paste0(gsub("[^A-Za-z0-9_-]", "_", met), ".png")),
             plot = p, width = 5, height = 4, dpi = 300, bg = "white"),
      error = function(e) message("Save failed '", met, "': ", e$message)
    )
  }
}
cat("Individual plots saved to:", OUTPUT_DIR, "\n")

# ---- 11. QUICK DATA EXPLORATION ----
cat("\nPathways in dataset:\n")
print(table(data_clean$Pathway))

cat("\nFirst 20 compounds:\n")
print(head(unique(data_clean$Compound), 20))

if (length(available) > 0) {
  cat("\nSummary for:", available[1], "\n")
  df_long %>%
    filter(Compound == available[1], Group %in% GROUP_LEVELS) %>%
    group_by(Group) %>%
    summarise(
      n    = sum(!is.na(Intensity)),
      Mean = round(mean(Intensity, na.rm = TRUE), 0),
      SD   = round(sd(Intensity,   na.rm = TRUE), 0),
      SE   = round(sd(Intensity,   na.rm = TRUE) / sqrt(sum(!is.na(Intensity))), 0),
      .groups = "drop"
    ) %>%
    print()
}


# ---- 12. CLUSTERED HEATMAP — TOP 50 METABOLITES ----
# ============================================================
#
# Rows    : top 50 metabolites by minimum FDR across all ANOVA effects
# Columns : individual samples (all 4 groups)
# Values  : z-score per row (mean = 0, SD = 1 across samples)
# Clustering: hierarchical on both axes (euclidean + complete linkage)
# Annotation: top bar coloured by experimental group
# ============================================================

# install.packages("pheatmap")
library(pheatmap)

# ---- Select top 50 metabolites ----
# Top 50 by minimum FDR across omnibus test and all 3 contrasts
top50 <- results_full %>%
  rowwise() %>%
  mutate(min_FDR = min(FDR,
                       `FDR_young: CTR vs IL-1b`,
                       `FDR_old: CTR vs IL-1b`,
                       `FDR_CTR: young vs old`,
                       na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(min_FDR)) %>%
  arrange(min_FDR) %>%
  slice_head(n = 50) %>%
  pull(Compound)

# ---- Build sample-level matrix ----
# Use raw (non-summarised) intensities from df_long, all 4 target groups
all_samples <- c(samples_ctr_young, samples_il1b_young,
                 samples_ctr_old,   samples_il1b_old)

mat <- data_clean %>%
  filter(Compound %in% top50) %>%
  select(Compound, all_of(all_samples)) %>%
  column_to_rownames("Compound") %>%
  as.matrix()

# Keep only rows with at least 1 non-NA value; replace remaining NAs with row mean
mat <- mat[rowSums(!is.na(mat)) > 0, ]
mat <- apply(mat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}) %>% t()

# Row-wise z-score
mat_z <- t(scale(t(mat)))

# Preserve top50 ordering as initial sort (pheatmap will recluster)
mat_z <- mat_z[intersect(top50, rownames(mat_z)), ]

# ---- Shorten row labels for readability ----
rownames(mat_z) <- gsub("_", " ", rownames(mat_z))

# ---- Column annotation (group membership) ----
ann_col <- data.frame(
  Group = c(
    rep("young CTR",  length(samples_ctr_young)),
    rep("young IL-1\u03b2", length(samples_il1b_young)),
    rep("old CTR",    length(samples_ctr_old)),
    rep("old IL-1\u03b2",   length(samples_il1b_old))
  ),
  row.names = all_samples
)

ann_colours <- list(
  Group = c(
    "young CTR"        = group_colours[["CTR young"]],
    "young IL-1\u03b2" = group_colours[["IL-1\u03b2 young"]],
    "old CTR"          = group_colours[["CTR old"]],
    "old IL-1\u03b2"   = group_colours[["IL-1\u03b2 old"]]
  )
)

# ---- Colour scale: blue — white — red, capped at ±3 SD ----
breaks     <- seq(-3, 3, length.out = 101)
hm_colours <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

# ---- Draw and save ----
pheatmap(
  mat_z,
  color            = hm_colours,
  breaks           = breaks,
  annotation_col   = ann_col,
  annotation_colors= ann_colours,
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_colnames    = FALSE,   # sample IDs are not informative
  fontsize_row     = 7,
  border_color     = NA,
  main             = "Top 50 metabolites by ANOVA significance\n(z-score per metabolite)",
  filename         = "heatmap_top50.png",
  width            = 8,
  height           = 12,
  dpi              = 300
)
cat("Saved: heatmap_top50.png\n")

cat("\nSession info:\n")
sessionInfo()
