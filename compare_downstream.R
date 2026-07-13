# ============================================================
# compare_downstream.R — Full statistical comparison
#
# Runs both converters -> dataProcess -> groupComparison,
# then compares log2FC and adjusted p-values on scatter plots.
# ============================================================

library(MSstatsConvertLLM)     # run_trial(), score_trial(), LLMtoMSstatsFormat()

library(MSstats)
library(ggplot2)
library(data.table)

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
MODEL      <- "llama3.2-3b"
TOOL       <- "proteome_discoverer"
PROMPT     <- "filter_aware"
ANNOT_PATH <- system.file("tinytest/annotations/annot_pd.csv",
                          package = "MSstats")

# Contrast matrix — set to NULL to auto-generate all pairwise
# Or define manually, e.g.:
# CONTRAST_MATRIX <- matrix(c(1, -1, 0, 0), nrow = 1)
# rownames(CONTRAST_MATRIX) <- c("C2-C1")
CONTRAST_MATRIX <- NULL
# ------------------------------------------------------------------

cat("=== Downstream Statistical Comparison ===\n")
cat("Model:", MODEL, "| Tool:", TOOL, "| Prompt:", PROMPT, "\n\n")

# ==================================================================
# Step 1: Get LLM mapping
# ==================================================================
cat("--- Step 1: Running LLM inference ---\n")
trial <- run_trial(MODEL, TOOL, PROMPT)
scores <- score_trial(trial)
print_scorecard(scores, trial)

llm_mapping <- trial$mapping

# ==================================================================
# Step 2: Run both converters
# ==================================================================
cat("\n--- Step 2a: Running LLM converter ---\n")
raw_data   <- load_test_dataset(TOOL)
annotation <- data.table::fread(ANNOT_PATH)

llm_converted <- tryCatch(
  LLMtoMSstatsFormat(
    input        = raw_data,
    annotation   = annotation,
    llm_mapping  = llm_mapping,
    verbose      = TRUE,
    use_log_file = FALSE
  ),
  error = function(e) {
    message("LLM converter FAILED: ", e$message)
    NULL
  }
)

# Extract and display diagnostics
llm_diagnostics <- NULL
if (!is.null(llm_converted)) {
  llm_diagnostics <- attr(llm_converted, "diagnostics")
  if (!is.null(llm_diagnostics)) {
    cat("\n--- LLM Converter Diagnostics ---\n")
    if (nrow(llm_diagnostics$hallucinated_columns) > 0) {
      cat("  HALLUCINATED COLUMNS:\n")
      for (i in seq_len(nrow(llm_diagnostics$hallucinated_columns))) {
        r <- llm_diagnostics$hallucinated_columns[i]
        cat(sprintf("    %s -> '%s' (does not exist in data)\n",
                    r$field, r$hallucinated_name))
      }
    } else {
      cat("  No hallucinated columns\n")
    }
    if (nrow(llm_diagnostics$skipped_filters) > 0) {
      cat("  SKIPPED FILTERS:\n")
      for (i in seq_len(nrow(llm_diagnostics$skipped_filters))) {
        r <- llm_diagnostics$skipped_filters[i]
        cat(sprintf("    %s %s %s — reason: %s\n",
                    r$column, r$operation, r$value, r$reason))
      }
    }
    if (nrow(llm_diagnostics$applied_filters) > 0) {
      cat("  APPLIED FILTERS:\n")
      for (i in seq_len(nrow(llm_diagnostics$applied_filters))) {
        r <- llm_diagnostics$applied_filters[i]
        cat(sprintf("    %s %s %s -> removed %d rows\n",
                    r$column, r$operation, r$value, r$rows_removed))
      }
    }
    cat("\n")
  }
}

cat("\n--- Step 2b: Running hand-written converter ---\n")
raw_paths <- list(
  spectronaut         = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                                    package = "MSstatsConvert"),
  proteome_discoverer = system.file("tinytest/raw_data/PD/pd_input.csv",
                                    package = "MSstatsConvert"),
  metamorpheus        = system.file("tinytest/raw_data/Metamorpheus/QuantifiedPeaks.tsv",
                                    package = "MSstatsConvert")
)

raw_original <- data.table::fread(raw_paths[[TOOL]])

hand_converted <- tryCatch({
  if (TOOL == "spectronaut") {
    SpectronauttoMSstatsFormat(raw_original, annotation, use_log_file = FALSE)
  } else if (TOOL == "proteome_discoverer") {
    PDtoMSstatsFormat(raw_original, annotation, use_log_file = FALSE)
  } else if (TOOL == "metamorpheus") {
    MetamorpheusToMSstatsFormat(raw_original, annotation, use_log_file = FALSE)
  } else {
    stop("No hand-written converter for: ", TOOL)
  }
}, error = function(e) {
  message("Hand-written converter FAILED: ", e$message)
  NULL
})

if (is.null(llm_converted) || is.null(hand_converted)) {
  stop("One or both converters failed. Cannot proceed with comparison.")
}

# Coerce to data.table
if (!is.data.table(hand_converted)) hand_converted <- as.data.table(hand_converted)
if (!is.data.table(llm_converted))  llm_converted  <- as.data.table(llm_converted)

# ==================================================================
# Step 3: dataProcess on both
# ==================================================================
cat("\n--- Step 3: Running dataProcess ---\n")

cat("  dataProcess (hand-written)...\n")
hand_processed <- tryCatch(
  dataProcess(hand_converted, use_log_file = FALSE),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

cat("  dataProcess (LLM)...\n")
llm_processed <- tryCatch(
  dataProcess(llm_converted, use_log_file = FALSE),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

if (is.null(hand_processed) || is.null(llm_processed)) {
  stop("dataProcess failed for one or both. Check converter output format.")
}

# ==================================================================
# Step 4: Build contrast matrix and run groupComparison
# ==================================================================
cat("\n--- Step 4: Running groupComparison ---\n")

# Auto-generate all pairwise contrasts if not specified
if (is.null(CONTRAST_MATRIX)) {
  conditions <- levels(hand_processed$ProteinLevelData$GROUP)
  if (is.null(conditions)) {
    conditions <- sort(unique(as.character(hand_processed$ProteinLevelData$GROUP)))
  }
  n_cond <- length(conditions)
  
  pairs <- combn(n_cond, 2)
  CONTRAST_MATRIX <- matrix(0, nrow = ncol(pairs), ncol = n_cond)
  colnames(CONTRAST_MATRIX) <- conditions
  contrast_names <- character(ncol(pairs))
  
  for (i in seq_len(ncol(pairs))) {
    CONTRAST_MATRIX[i, pairs[1, i]] <-  1
    CONTRAST_MATRIX[i, pairs[2, i]] <- -1
    contrast_names[i] <- paste0(conditions[pairs[1, i]], "-", 
                                 conditions[pairs[2, i]])
  }
  rownames(CONTRAST_MATRIX) <- contrast_names
  
  cat("  Auto-generated contrasts:", paste(contrast_names, collapse = ", "), "\n")
}

cat("  groupComparison (hand-written)...\n")
hand_comparison <- tryCatch(
  groupComparison(CONTRAST_MATRIX, hand_processed, use_log_file = FALSE),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

cat("  groupComparison (LLM)...\n")
llm_comparison <- tryCatch(
  groupComparison(CONTRAST_MATRIX, llm_processed, use_log_file = FALSE),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

if (is.null(hand_comparison) || is.null(llm_comparison)) {
  stop("groupComparison failed. Check that both dataProcess outputs have matching conditions.")
}

# ==================================================================
# Step 5: Merge results and build comparison table
# ==================================================================
cat("\n--- Step 5: Merging results ---\n")

hand_res <- as.data.table(hand_comparison$ComparisonResult)
llm_res  <- as.data.table(llm_comparison$ComparisonResult)

# Merge on Protein + Label (comparison name)
merged <- merge(
  hand_res[, .(Protein, Label, 
               log2FC_hand = log2FC, 
               pvalue_hand = pvalue,
               adj.pvalue_hand = adj.pvalue,
               SE_hand = SE,
               DF_hand = DF)],
  llm_res[, .(Protein, Label, 
              log2FC_llm = log2FC,
              pvalue_llm = pvalue,
              adj.pvalue_llm = adj.pvalue,
              SE_llm = SE,
              DF_llm = DF)],
  by = c("Protein", "Label"),
  all = FALSE  # inner join — only shared proteins
)

cat(sprintf("  Proteins in hand-written: %d\n", uniqueN(hand_res$Protein)))
cat(sprintf("  Proteins in LLM:         %d\n", uniqueN(llm_res$Protein)))
cat(sprintf("  Shared (for comparison):  %d\n", uniqueN(merged$Protein)))
cat(sprintf("  Total comparison rows:    %d\n", nrow(merged)))

# Summary stats
valid_fc <- complete.cases(merged[, .(log2FC_hand, log2FC_llm)])
if (sum(valid_fc) > 2) {
  r_fc <- cor(merged$log2FC_hand[valid_fc], merged$log2FC_llm[valid_fc])
  cat(sprintf("\n  log2FC correlation:    r = %.4f (n = %d)\n", r_fc, sum(valid_fc)))
}

valid_pv <- complete.cases(merged[, .(adj.pvalue_hand, adj.pvalue_llm)])
if (sum(valid_pv) > 2) {
  r_pv <- cor(-log10(merged$adj.pvalue_hand[valid_pv]), 
              -log10(merged$adj.pvalue_llm[valid_pv]))
  cat(sprintf("  -log10(adj.pvalue) corr: r = %.4f (n = %d)\n", r_pv, sum(valid_pv)))
}

# ==================================================================
# Step 5b: Merge protein-level summarized values from dataProcess
# ==================================================================
cat("\n--- Step 5b: Merging summarized protein-level data ---\n")

hand_prot <- as.data.table(hand_processed$ProteinLevelData)
llm_prot  <- as.data.table(llm_processed$ProteinLevelData)

merged_summary <- merge(
  hand_prot[, .(Protein, GROUP, SUBJECT, RUN, LogIntensities_hand = LogIntensities)],
  llm_prot[, .(Protein, GROUP, SUBJECT, RUN, LogIntensities_llm = LogIntensities)],
  by = c("Protein", "GROUP", "SUBJECT", "RUN"),
  all = FALSE
)

# If RUN doesn't match (different Run column choices), try Protein + GROUP + SUBJECT only
if (nrow(merged_summary) == 0) {
  cat("  RUN-level merge empty, trying Protein + GROUP + SUBJECT...\n")
  merged_summary <- merge(
    hand_prot[, .(Protein, GROUP, SUBJECT, LogIntensities_hand = LogIntensities)],
    llm_prot[, .(Protein, GROUP, SUBJECT, LogIntensities_llm = LogIntensities)],
    by = c("Protein", "GROUP", "SUBJECT"),
    all = FALSE
  )
}

valid_summary <- complete.cases(merged_summary[, .(LogIntensities_hand, LogIntensities_llm)])
if (sum(valid_summary) > 2) {
  r_sum <- cor(merged_summary$LogIntensities_hand[valid_summary],
               merged_summary$LogIntensities_llm[valid_summary])
  cat(sprintf("  Summarized intensity correlation: r = %.4f (n = %d)\n",
              r_sum, sum(valid_summary)))
}

# ==================================================================
# Step 6: Scatter plots
# ==================================================================
cat("\n--- Step 6: Generating plots ---\n")

outdir <- "results"
dir.create(outdir, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Plot 1: Summarized protein-level intensities ---
p_summary <- ggplot(merged_summary[valid_summary], 
                    aes(x = LogIntensities_hand, y = LogIntensities_llm)) +
  geom_point(alpha = 0.4, size = 1.2, color = "#2C8A5F") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "grey40", linewidth = 0.7) +
  labs(
    x = "log2 Intensity (Hand-written)",
    y = "log2 Intensity (LLM)",
    title = sprintf("Summarized Intensities: %s", TOOL),
    subtitle = sprintf("r = %.4f, n = %d observations",
                       cor(merged_summary$LogIntensities_hand[valid_summary],
                           merged_summary$LogIntensities_llm[valid_summary]),
                       sum(valid_summary))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

sum_path <- file.path(outdir, paste0("summarized_scatter_", TOOL, "_", MODEL, "_",
                                      timestamp, ".pdf"))
ggsave(sum_path, p_summary, width = 6, height = 6, dpi = 300)
cat("  Saved:", sum_path, "\n")

# --- Plot 2: log2FC ---
p_fc <- ggplot(merged[valid_fc], aes(x = log2FC_hand, y = log2FC_llm)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#2C5F8A") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "grey40", linewidth = 0.7) +
  labs(
    x = "log2 Fold Change (Hand-written converter)",
    y = "log2 Fold Change (LLM converter)",
    title = sprintf("log2FC Comparison: %s (%s, %s prompt)",
                    TOOL, MODEL, PROMPT),
    subtitle = sprintf("r = %.4f, n = %d proteins",
                       cor(merged$log2FC_hand[valid_fc], 
                           merged$log2FC_llm[valid_fc]),
                       uniqueN(merged$Protein[valid_fc]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

fc_path <- file.path(outdir, paste0("log2FC_scatter_", TOOL, "_", MODEL, "_",
                                     timestamp, ".pdf"))
ggsave(fc_path, p_fc, width = 6, height = 6, dpi = 300)
cat("  Saved:", fc_path, "\n")

# --- Plot 3: Adjusted p-value (-log10) ---
pv_data <- merged[valid_pv & adj.pvalue_hand > 0 & adj.pvalue_llm > 0]

p_pv <- ggplot(pv_data, aes(x = -log10(adj.pvalue_hand), 
                              y = -log10(adj.pvalue_llm))) +
  geom_point(alpha = 0.5, size = 1.5, color = "#8A2C2C") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "grey40", linewidth = 0.7) +
  labs(
    x = "-log10(adj. p-value) Hand-written",
    y = "-log10(adj. p-value) LLM",
    title = sprintf("Adjusted P-value Comparison: %s (%s, %s prompt)",
                    TOOL, MODEL, PROMPT),
    subtitle = sprintf("r = %.4f, n = %d proteins",
                       cor(-log10(pv_data$adj.pvalue_hand), 
                           -log10(pv_data$adj.pvalue_llm)),
                       uniqueN(pv_data$Protein))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

pv_path <- file.path(outdir, paste0("adjpval_scatter_", TOOL, "_", MODEL, "_",
                                     timestamp, ".pdf"))
ggsave(pv_path, p_pv, width = 6, height = 6, dpi = 300)
cat("  Saved:", pv_path, "\n")

# --- Plot 4: Combined panel (for poster) ---
p_combined <- cowplot::plot_grid(p_summary, p_fc, p_pv, ncol = 3, 
                                 labels = c("A", "B", "C"))

combined_path <- file.path(outdir, paste0("comparison_panel_", TOOL, "_", MODEL, "_",
                                           timestamp, ".pdf"))
ggsave(combined_path, p_combined, width = 18, height = 6, dpi = 300)
cat("  Saved:", combined_path, "\n")

# ==================================================================
# Step 7: Save all results
# ==================================================================
saveRDS(
  list(
    trial           = trial,
    scores          = scores,
    merged_results  = merged,
    merged_summary  = merged_summary,
    diagnostics     = llm_diagnostics,
    hand_comparison = hand_comparison$ComparisonResult,
    llm_comparison  = llm_comparison$ComparisonResult,
    config = list(model = MODEL, tool = TOOL, prompt = PROMPT)
  ),
  file.path(outdir, paste0("full_comparison_", TOOL, "_", MODEL, "_",
                            timestamp, ".rds"))
)

cat("\n=== Done ===\n")
cat(sprintf("Shared proteins compared: %d\n", uniqueN(merged$Protein)))
cat(sprintf("log2FC correlation: %.4f\n", 
            cor(merged$log2FC_hand, merged$log2FC_llm, use = "complete.obs")))