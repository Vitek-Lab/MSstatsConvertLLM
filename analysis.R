library(tidyverse)
library(data.table)

summary_results <- fread("results/summary_20260527_102710.csv")


summary_results %>% filter(tool=="spectronaut") %>%
    ggplot() + 
    geom_bar(aes(x = model, y = exact_accuracy, fill=prompt), stat = "identity", position = "dodge") + 
    theme_bw(base_size = 28)

summary_results %>% filter(tool=="spectronaut") %>%
    ggplot() + 
    geom_bar(aes(x = model, y = acceptable_accuracy, fill=prompt), stat = "identity", position = "dodge") + 
    theme_bw(base_size = 28)


summary_results %>%
    ggplot() + 
    geom_point(aes(x=acceptable_accuracy, y=mean_conf, color=model)) + theme_bw(base_size = 28)

field_scores <- fread("results/field_scores_20260527_102710.csv")    

field_scores %>% filter(expected != "" & model == "llama3.2-3b") %>%
    ggplot() + 
    geom_boxplot(aes(x=exact_match, y=confidence, fill=prompt)) + theme_bw(base_size = 28)



library(ggplot2)
library(cowplot)
library(dplyr)

# Shared palette — colorblind-friendly
prompt_colors <- c(
  "lean"                = "#E69F00",
  "constrained"         = "#56B4E9",
  "filter_aware"        = "#009E73",
  "constrained_filter"  = "#CC79A7"
)

prompt_labels <- c(
  "lean"                = "Lean",
  "constrained"         = "Constrained",
  "filter_aware"        = "Filter-Aware",
  "constrained_filter"  = "Constrained +\nFilter"
)

tool_labels <- c(
  "spectronaut"         = "Spectronaut (DIA)",
  "proteome_discoverer" = "Proteome Discoverer (DDA)",
  "metamorpheus"        = "MetaMorpheus (DDA)"
)

summary_results = summary_results %>% filter(prompt !=  "filter_aware")

# Clean up model labels for display
summary_plot <- summary_results %>%
  mutate(
    model_label = case_when(
      grepl("3b", model)       ~ "Llama 3.2\n3B",
      grepl("8b", model)       ~ "Llama 3.1\n8B",
      grepl("deepseek", model) ~ "DeepSeek-R1\n14B",
      TRUE ~ model
    ),
    model_label = factor(model_label, levels = c("Llama 3.2\n3B", 
                                                  "Llama 3.1\n8B", 
                                                  "DeepSeek-R1\n14B")),
    prompt_label = prompt_labels[prompt],
    prompt_label = factor(prompt_label, levels = prompt_labels),
    tool_label = tool_labels[tool],
    tool_label = factor(tool_label, levels = tool_labels)
  )

# --- Top row: Exact accuracy ---
p_exact <- ggplot(summary_plot, 
                  aes(x = model_label, y = exact_accuracy, fill = prompt)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", exact_accuracy * 100)),
            position = position_dodge(0.8), vjust = -0.4, size = 4.75) +
  facet_wrap(~ tool_label, nrow = 1) +
  scale_fill_manual(values = prompt_colors, labels = prompt_labels,
                    name = "Prompt") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                     labels = scales::percent) +
  labs(y = "Exact Match Accuracy") +
  theme_bw(base_size = 22) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 28),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 0, 10)
  )

# --- Bottom row: Acceptable accuracy ---
p_acceptable <- ggplot(summary_plot, 
                       aes(x = model_label, y = acceptable_accuracy, fill = prompt)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", acceptable_accuracy * 100)),
            position = position_dodge(0.8), vjust = -0.4, size = 4.75) +
  facet_wrap(~ tool_label, nrow = 1) +
  scale_fill_manual(values = prompt_colors, labels = prompt_labels,
                    name = "Prompt") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                     labels = scales::percent) +
  labs(y = "Acceptable Match Accuracy") +
  theme_bw(base_size = 22) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 22),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 10, 10, 10)
  )

# --- Combine ---
p_panel <- plot_grid(
  p_exact, p_acceptable,
  ncol = 1, rel_heights = c(0.47, 0.53)
  # labels = c("A", "B"), label_size = 22, label_fontface = "bold"
)

ggsave("results/accuracy_panel.pdf", p_panel, 
       width = 18, height = 10, dpi = 300)
ggsave("results/accuracy_panel.png", p_panel, 
       width = 18, height = 10, dpi = 300)

cat("Saved to results/accuracy_panel.pdf\n")

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

filters <- fread("results/filters_20260527_102710.csv")

# ============================================================
# 1. Classify filters by quality
# ============================================================

# Known good filters per tool (from MSstatsConvert source)
KNOWN_GOOD <- list(
  spectronaut = c("F.ExcludedFromQuantification", "F.FrgLossType", "F.PossibleInterference"),
  proteome_discoverer = c("XProteins"),
  metamorpheus = c("Decoy.Peptide", "Random.RT")
)

# Known harmful filters — would corrupt the analysis
KNOWN_HARMFUL <- list(
  spectronaut = c("R.Condition", "EG.Library", "EG.iRTPredicted"),  # filtering on condition removes experimental groups
  proteome_discoverer = c(),
  metamorpheus = c()
)

classify_filter <- function(column, tool) {
  good    <- KNOWN_GOOD[[tool]]
  harmful <- KNOWN_HARMFUL[[tool]]
  
  if (column %in% good)    return("Correct")
  if (column %in% harmful) return("Harmful")
  return("Novel")  # not in our ground truth — needs expert review
}

filters[, quality := mapply(classify_filter, column, tool)]
filters[, quality := factor(quality, levels = c("Correct", "Novel", "Harmful"))]

# ============================================================
# 2. Summary: how many filters per condition, and what quality?
# ============================================================

# Deduplicate across reps (take rep 1 as representative since they're mostly identical)
filters_dedup <- filters[rep == 1]

cat("=== Filter Discovery Summary ===\n\n")

# Count by model × tool × prompt × quality
filter_counts <- filters_dedup[, .N, by = .(model, tool, prompt, quality)]
filter_counts_wide <- dcast(filter_counts, model + tool + prompt ~ quality, 
                             value.var = "N", fill = 0)
print(filter_counts_wide)

cat("\n=== Most common filters per tool ===\n\n")
common <- filters_dedup[, .N, by = .(tool, column, quality)]
common <- common[order(tool, -N)]
print(common)

# ============================================================
# 3. Reproducibility across reps
# ============================================================
cat("\n=== Filter Reproducibility (across 3 reps) ===\n")

repro <- filters[, .(
  n_reps_present = uniqueN(rep)
), by = .(model, tool, prompt, column)]

repro_summary <- repro[, .(
  always_present = sum(n_reps_present == 3),
  sometimes      = sum(n_reps_present < 3 & n_reps_present > 0),
  total_unique   = .N
), by = .(model, tool, prompt)]

print(repro_summary)

# ============================================================
# 4. Plots
# ============================================================

# --- Plot A: Filter count by quality, per model × tool ---
# Aggregate: one bar per model × tool, stacked by quality
plot_data_a <- filters_dedup[, .N, by = .(model, tool, prompt, quality)]

# Clean labels
plot_data_a[, model_label := fcase(
  grepl("3b", model),       "Llama 3.2\n3B",
  grepl("8b", model),       "Llama 3.1\n8B",
  grepl("deepseek", model), "DeepSeek-R1\n14B"
)]
plot_data_a[, model_label := factor(model_label, 
  levels = c("Llama 3.2\n3B", "Llama 3.1\n8B", "DeepSeek-R1\n14B"))]

tool_labels <- c(
  "spectronaut"         = "Spectronaut (DIA)",
  "proteome_discoverer" = "Proteome Discoverer (DDA)",
  "metamorpheus"        = "MetaMorpheus (DDA)"
)
plot_data_a[, tool_label := tool_labels[tool]]
plot_data_a[, tool_label := factor(tool_label, levels = tool_labels)]

prompt_labels <- c("filter_aware" = "Filter-Aware", 
                   "constrained_filter" = "Constrained + Filter")
plot_data_a[, prompt_label := prompt_labels[prompt]]

quality_colors <- c("Correct" = "#2CA02C", "Novel" = "#1F77B4", "Harmful" = "#D62728")

p_counts <- ggplot(plot_data_a, 
                   aes(x = model_label, y = N, fill = quality)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), 
            size = 4, color = "white", fontface = "bold") +
  facet_grid(prompt_label ~ tool_label) +
  scale_fill_manual(values = quality_colors, name = "Filter Quality") +
  labs(y = "Number of Filters Suggested", x = NULL) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11)
  )

# --- Plot B: Heatmap of which filters each model discovers ---
# Show unique column × tool combinations, which models find them
heatmap_data <- filters_dedup[, .(
  discovered = TRUE
), by = .(model, tool, column, quality)]

# Make a complete grid
all_combos <- CJ(
  model  = unique(filters_dedup$model),
  tool   = unique(filters_dedup$tool),
  column = unique(filters_dedup$column)
)
# Only keep column × tool combos that actually appeared
valid_pairs <- unique(filters_dedup[, .(tool, column, quality)])
all_combos <- merge(all_combos, valid_pairs, by = c("tool", "column"), 
                     allow.cartesian = TRUE)
heatmap_data <- merge(all_combos, heatmap_data, 
                       by = c("model", "tool", "column", "quality"), 
                       all.x = TRUE)
heatmap_data[is.na(discovered), discovered := FALSE]

heatmap_data[, model_label := fcase(
  grepl("3b", model),       "Llama 3.2 3B",
  grepl("8b", model),       "Llama 3.1 8B",
  grepl("deepseek", model), "DeepSeek-R1 14B"
)]
heatmap_data[, tool_label := tool_labels[tool]]
heatmap_data[, tool_label := factor(tool_label, levels = tool_labels)]

# Order columns: correct first, then novel, then harmful
heatmap_data[, column := factor(column, 
  levels = unique(heatmap_data[order(quality, column)]$column))]

p_heatmap <- ggplot(heatmap_data, 
                    aes(x = model_label, y = column, fill = interaction(discovered, quality))) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ tool_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(
      "TRUE.Correct"  = "#2CA02C",
      "TRUE.Novel"    = "#1F77B4", 
      "TRUE.Harmful"  = "#D62728",
      "FALSE.Correct" = "grey90",
      "FALSE.Novel"   = "grey90",
      "FALSE.Harmful" = "grey90"
    ),
    breaks = c("TRUE.Correct", "TRUE.Novel", "TRUE.Harmful", "FALSE.Correct"),
    labels = c("Discovered (Correct)", "Discovered (Novel)", 
               "Discovered (Harmful)", "Not discovered"),
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 22) +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 15),
    axis.text.y = element_text(size = 12, face = "italic"),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

ggsave(file.path(outdir, "filter_heatmap.pdf"), p_heatmap, 
       width = 14, height = 6, dpi = 300)
ggsave(file.path(outdir, "filter_heatmap.png"), p_heatmap, 
       width = 14, height = 6, dpi = 300)

# --- Plot C: Does the 3B model over-suggest? Filter count distribution ---
n_filters <- filters_dedup[, .N, by = .(model, tool, prompt)]
n_filters[, model_label := fcase(
  grepl("3b", model),       "Llama 3.2 3B",
  grepl("8b", model),       "Llama 3.1 8B",
  grepl("deepseek", model), "DeepSeek-R1 14B"
)]
n_filters[, model_label := factor(model_label, 
  levels = c("Llama 3.2 3B", "Llama 3.1 8B", "DeepSeek-R1 14B"))]
n_filters[, tool_label := tool_labels[tool]]
n_filters[, prompt_label := prompt_labels[prompt]]

p_nfilters <- ggplot(n_filters, aes(x = model_label, y = N, fill = prompt_label)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = N), position = position_dodge(0.7), vjust = -0.3, size = 4) +
  facet_wrap(~ tool_label) +
  scale_fill_manual(values = c("Filter-Aware" = "#009E73", 
                                "Constrained + Filter" = "#CC79A7"),
                    name = "Prompt") +
  labs(y = "Number of Filters Suggested", x = NULL) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# --- Save ---
outdir <- "results"

ggsave(file.path(outdir, "filter_quality_stacked.pdf"), p_counts, 
       width = 14, height = 7, dpi = 300)
ggsave(file.path(outdir, "filter_quality_stacked.png"), p_counts, 
       width = 14, height = 7, dpi = 300)

ggsave(file.path(outdir, "filter_heatmap.pdf"), p_heatmap, 
       width = 14, height = 6, dpi = 300)
ggsave(file.path(outdir, "filter_heatmap.png"), p_heatmap, 
       width = 14, height = 6, dpi = 300)

ggsave(file.path(outdir, "filter_count_by_model.pdf"), p_nfilters, 
       width = 14, height = 5, dpi = 300)
ggsave(file.path(outdir, "filter_count_by_model.png"), p_nfilters, 
       width = 14, height = 5, dpi = 300)

cat("\nPlots saved to results/\n")