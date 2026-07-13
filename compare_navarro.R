# ============================================================
# compare_downstream.R — Full statistical comparison
#
# Runs both converters -> dataProcess -> groupComparison,
# then compares summarized values, log2FC, adjusted p-values,
# and empirical FDR using species-labeled benchmark data.
# ============================================================

library(MSstatsConvertLLM)     # run_trial(), score_trial(), LLMtoMSstatsFormat()

library(MSstats)
library(ggplot2)
library(data.table)
library(stringr)
library(jsonlite)  # toJSON()/fromJSON() used directly below

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
MODEL      <- "deepseek-r1-14b"
TOOL       <- "spectronaut"
PROMPT     <- "constrained_filter"
ACQUISITION <- "DIA"  # "DIA", "DDA", or NULL to let LLM guess

# --- Data paths ---
# For Navarro 2016 benchmark: set these to your local paths
RAW_PATH   <- "data/Navarro2016_DIA_Spectronaut_input.xls"
ANNOT_PATH <- "data/Navarro2016_DIA_Spectronaut_annotation.csv"

# --- Contrast matrix ---
# For Navarro: 2 conditions (A vs B)
CONTRAST_MATRIX <- matrix(c(1, -1), nrow = 1)
rownames(CONTRAST_MATRIX) <- c("A-B")
colnames(CONTRAST_MATRIX) <- c("A", "B")

# Set to NULL for auto-generated all pairwise:
# CONTRAST_MATRIX <- NULL

# --- Species labels for empirical FDR ---
# Proteins matching TRUE_NEGATIVE_PATTERN should NOT be DE (FP if significant)
# Everything else is a true positive (ecoli, yeast in Navarro)
# Set to NULL to skip species-based FDR analysis
TRUE_NEGATIVE_PATTERN <- "HUMAN"

# --- dataProcess options (match your benchmark settings) ---
DP_FEATURE_SUBSET <- "topN"
DP_N_TOP_FEATURE  <- 10
DP_NORMALIZATION  <- FALSE
DP_MBIMPUTE       <- TRUE

# ------------------------------------------------------------------

cat("=== Downstream Statistical Comparison ===\n")
cat("Model:", MODEL, "| Tool:", TOOL, "| Prompt:", PROMPT, "\n\n")

# ==================================================================
# Step 1: Get LLM mapping
# ==================================================================
cat("--- Step 1: Running LLM inference ---\n")

# Load raw data and feed to LLM (normalized headers)
raw_data <- normalize_headers(data.table::fread(RAW_PATH))
preview <- make_json_preview(raw_data, n_rows = 3)

system_prompt <- PROMPT_VERSIONS[[PROMPT]]
user_prompt   <- build_user_prompt(preview, acquisition = ACQUISITION)

chat <- create_chat(MODEL)
chat$set_system_prompt(system_prompt)

# Pick structured type based on prompt (internal helper; builds an ellmer type)
result_type <- MSstatsConvertLLM:::.select_mapping_type(
  has_filters      = PROMPT %in% c("filter_aware", "constrained_filter"),
  allow_transforms = FALSE
)

t0 <- Sys.time()
llm_mapping <- tryCatch(
  chat$chat_structured(user_prompt, type = result_type, convert = TRUE),
  error = function(e) {
    message("  [WARN] chat_structured failed: ", e$message)
    chat2 <- create_chat(MODEL)
    chat2$set_system_prompt(system_prompt)
    raw_text <- chat2$chat(user_prompt)
    fromJSON(raw_text, simplifyVector = FALSE)
  }
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("  LLM inference: %.1fs\n", elapsed))
cat("  Mapping:\n")
cat(toJSON(llm_mapping, auto_unbox = TRUE, pretty = TRUE), "\n\n")

# ==================================================================
# Step 2: Run both converters
# ==================================================================
cat("--- Step 2a: Running LLM converter ---\n")
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
raw_original <- data.table::fread(RAW_PATH)

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

# MSstats needs this column (may be absent from LLM output)
# if (!"AnomalyScores" %in% colnames(llm_converted))  llm_converted[, AnomalyScores := NA]
# if (!"AnomalyScores" %in% colnames(hand_converted)) hand_converted[, AnomalyScores := NA]

cat("  dataProcess (hand-written)...\n")
hand_processed <- tryCatch(
  dataProcess(hand_converted, 
              featureSubset = DP_FEATURE_SUBSET,
              n_top_feature = DP_N_TOP_FEATURE,
              normalization = DP_NORMALIZATION,
              MBimpute = DP_MBIMPUTE,
              use_log_file = FALSE,
              numberOfCores=12),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

cat("  dataProcess (LLM)...\n")
llm_processed <- tryCatch(
  dataProcess(llm_converted,
              featureSubset = DP_FEATURE_SUBSET,
              n_top_feature = DP_N_TOP_FEATURE,
              normalization = DP_NORMALIZATION,
              MBimpute = DP_MBIMPUTE,
              use_log_file = FALSE,
              numberOfCores=12),
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
  groupComparison(CONTRAST_MATRIX, hand_processed, use_log_file = FALSE, numberOfCores=12),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

cat("  groupComparison (LLM)...\n")
llm_comparison <- tryCatch(
  groupComparison(CONTRAST_MATRIX, llm_processed, use_log_file = FALSE, numberOfCores=12),
  error = function(e) { message("  FAILED: ", e$message); NULL }
)

if (is.null(hand_comparison) || is.null(llm_comparison)) {
  stop("groupComparison failed. Check that both dataProcess outputs have matching conditions.")
}

# ==================================================================
# Step 5a: Merge groupComparison results
# ==================================================================
cat("\n--- Step 5a: Merging groupComparison results ---\n")

hand_res <- as.data.table(hand_comparison$ComparisonResult)
llm_res  <- as.data.table(llm_comparison$ComparisonResult)

merged <- merge(
  hand_res[, .(Protein, Label, 
               log2FC_hand = log2FC, pvalue_hand = pvalue,
               adj.pvalue_hand = adj.pvalue, SE_hand = SE, DF_hand = DF)],
  llm_res[, .(Protein, Label, 
              log2FC_llm = log2FC, pvalue_llm = pvalue,
              adj.pvalue_llm = adj.pvalue, SE_llm = SE, DF_llm = DF)],
  by = c("Protein", "Label"),
  all = FALSE
)

cat(sprintf("  Proteins in hand-written: %d\n", uniqueN(hand_res$Protein)))
cat(sprintf("  Proteins in LLM:         %d\n", uniqueN(llm_res$Protein)))
cat(sprintf("  Shared (for comparison):  %d\n", uniqueN(merged$Protein)))

valid_fc <- complete.cases(merged[, .(log2FC_hand, log2FC_llm)])
if (sum(valid_fc) > 2) {
  r_fc <- cor(merged$log2FC_hand[valid_fc], merged$log2FC_llm[valid_fc])
  cat(sprintf("  log2FC correlation:      r = %.4f (n = %d)\n", r_fc, sum(valid_fc)))
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
# Step 5c: Species-based empirical FDR analysis
# ==================================================================
compute_fdr_stats <- function(comparison_result, tn_pattern, alpha = 0.05) {
  dt <- as.data.table(comparison_result)
  dt[, species := ifelse(str_detect(Protein, tn_pattern), "TN", "TP")]
  
  sig <- dt[adj.pvalue < alpha]
  
  n_total_sig <- nrow(sig)
  n_fp        <- sig[species == "TN", .N]
  n_tp        <- sig[species == "TP", .N]
  empirical_fdr <- if (n_total_sig > 0) n_fp / n_total_sig else NA_real_
  
  # Also compute across thresholds for ROC-like curve
  thresholds <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2)
  fdr_curve <- rbindlist(lapply(thresholds, function(a) {
    s <- dt[adj.pvalue < a]
    n_s <- nrow(s)
    data.table(
      threshold = a,
      n_sig     = n_s,
      n_fp      = s[species == "TN", .N],
      n_tp      = s[species == "TP", .N],
      fdr       = if (n_s > 0) s[species == "TN", .N] / n_s else NA_real_
    )
  }))
  
  list(
    n_proteins    = uniqueN(dt$Protein),
    n_tn          = uniqueN(dt[species == "TN"]$Protein),
    n_tp          = uniqueN(dt[species == "TP"]$Protein),
    n_sig_005     = n_total_sig,
    n_fp_005      = n_fp,
    n_tp_005      = n_tp,
    empirical_fdr = empirical_fdr,
    fdr_curve     = fdr_curve,
    per_protein   = dt
  )
}

fdr_results <- NULL
if (!is.null(TRUE_NEGATIVE_PATTERN)) {
  cat("\n--- Step 5c: Species-based empirical FDR ---\n")
  cat(sprintf("  True negative pattern: '%s'\n", TRUE_NEGATIVE_PATTERN))
  
  hand_fdr <- compute_fdr_stats(hand_comparison$ComparisonResult, TRUE_NEGATIVE_PATTERN)
  llm_fdr  <- compute_fdr_stats(llm_comparison$ComparisonResult, TRUE_NEGATIVE_PATTERN)
  
  cat(sprintf("\n  Hand-written converter:\n"))
  cat(sprintf("    Total proteins: %d (TN: %d, TP: %d)\n",
              hand_fdr$n_proteins, hand_fdr$n_tn, hand_fdr$n_tp))
  cat(sprintf("    Significant (adj.p < 0.05): %d (FP: %d, TP: %d)\n",
              hand_fdr$n_sig_005, hand_fdr$n_fp_005, hand_fdr$n_tp_005))
  cat(sprintf("    Empirical FDR: %.4f\n", hand_fdr$empirical_fdr))
  
  cat(sprintf("\n  LLM converter:\n"))
  cat(sprintf("    Total proteins: %d (TN: %d, TP: %d)\n",
              llm_fdr$n_proteins, llm_fdr$n_tn, llm_fdr$n_tp))
  cat(sprintf("    Significant (adj.p < 0.05): %d (FP: %d, TP: %d)\n",
              llm_fdr$n_sig_005, llm_fdr$n_fp_005, llm_fdr$n_tp_005))
  cat(sprintf("    Empirical FDR: %.4f\n", llm_fdr$empirical_fdr))
  
  fdr_results <- list(hand = hand_fdr, llm = llm_fdr)
}

# ==================================================================
# Step 6: Plots
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
    title = "Summarized Intensities",
    subtitle = sprintf("r = %.4f, n = %d",
                       cor(merged_summary$LogIntensities_hand[valid_summary],
                           merged_summary$LogIntensities_llm[valid_summary]),
                       sum(valid_summary))
  ) +
  theme_minimal(base_size = 24) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank()) +
  coord_fixed()

ggsave(file.path(outdir, paste0("summarized_scatter_", timestamp, ".pdf")),
       p_summary, width = 6, height = 6, dpi = 300)

# --- Plot 2: log2FC ---
p_fc <- ggplot(merged[valid_fc], aes(x = log2FC_hand, y = log2FC_llm)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#2C5F8A") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "grey40", linewidth = 0.7) +
  labs(
    x = "log2FC (Hand-written)",
    y = "log2FC (LLM)",
    title = "log2 Fold Change",
    subtitle = sprintf("r = %.4f, n = %d proteins",
                       cor(merged$log2FC_hand[valid_fc], merged$log2FC_llm[valid_fc]),
                       uniqueN(merged$Protein[valid_fc]))
  ) +
  theme_minimal(base_size = 24) +
  theme(plot.title = element_text(face = "bold", size = 24),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank()) +
  coord_fixed()

ggsave(file.path(outdir, paste0("log2FC_scatter_", timestamp, ".pdf")),
       p_fc, width = 6, height = 6, dpi = 300)

# --- Plot 3: Adjusted p-value ---
pv_data <- merged[valid_pv & adj.pvalue_hand > 0 & adj.pvalue_llm > 0]

p_pv <- ggplot(pv_data, aes(x = -log10(adj.pvalue_hand), 
                              y = -log10(adj.pvalue_llm))) +
  geom_point(alpha = 0.5, size = 1.5, color = "#8A2C2C") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "grey40", linewidth = 0.7) +
  labs(
    x = "-log10(adj. p-value) Hand-written",
    y = "-log10(adj. p-value) LLM",
    title = "Adjusted P-values",
    subtitle = sprintf("r = %.4f, n = %d proteins",
                       cor(-log10(pv_data$adj.pvalue_hand), 
                           -log10(pv_data$adj.pvalue_llm)),
                       uniqueN(pv_data$Protein))
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank()) +
  coord_fixed()

ggsave(file.path(outdir, paste0("adjpval_scatter_", timestamp, ".pdf")),
       p_pv, width = 6, height = 6, dpi = 300)

# --- Plot 4: Standard error ---
valid_se <- complete.cases(merged[, .(SE_hand, SE_llm)])
se_data <- merged[valid_se & SE_hand > 0 & SE_llm > 0]

if (sum(valid_se) > 2) {
  r_se <- cor(merged$SE_hand[valid_se], merged$SE_llm[valid_se])
  cat(sprintf("  SE correlation:          r = %.4f (n = %d)\n", r_se, sum(valid_se)))
}

p_se <- ggplot(se_data, aes(x = SE_hand, y = SE_llm)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#5F2C8A") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.7) +
  labs(
    x = "Standard Error (Hand-written)",
    y = "Standard Error (LLM)",
    title = "Standard Error",
    subtitle = sprintf("r = %.4f, n = %d proteins",
                       cor(se_data$SE_hand, se_data$SE_llm),
                       uniqueN(se_data$Protein))
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank()) +
  coord_fixed()

ggsave(file.path(outdir, paste0("SE_scatter_", timestamp, ".pdf")),
       p_se, width = 6, height = 6, dpi = 300)

# --- Plot 5: Empirical FDR comparison (if species labels available) ---
p_fdr <- NULL
if (!is.null(fdr_results)) {
  # Bar chart: hand vs LLM empirical FDR at 0.05
  fdr_bar_data <- data.table(
    Converter = c("Hand-written", "LLM"),
    FDR = c(fdr_results$hand$empirical_fdr, fdr_results$llm$empirical_fdr),
    TP  = c(fdr_results$hand$n_tp_005, fdr_results$llm$n_tp_005),
    FP  = c(fdr_results$hand$n_fp_005, fdr_results$llm$n_fp_005)
  )
  
  p_fdr_bar <- ggplot(fdr_bar_data, aes(x = Converter, y = FDR, fill = Converter)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = sprintf("%.3f\n(TP:%d, FP:%d)", FDR, TP, FP)), 
              vjust = -0.3, size = 3.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = c("Hand-written" = "#2C5F8A", "LLM" = "#8A5F2C")) +
    scale_y_continuous(limits = c(0, max(fdr_bar_data$FDR, 0.05) * 1.4)) +
    labs(x = NULL, y = "Empirical FDR",
         title = "Empirical FDR at adj.p < 0.05",
         subtitle = sprintf("HUMAN = true negative, ECOLI/YEAST = true positive")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey40"),
          legend.position = "none",
          panel.grid.minor = element_blank())
  
  ggsave(file.path(outdir, paste0("fdr_comparison_", timestamp, ".pdf")),
         p_fdr_bar, width = 5, height = 5, dpi = 300)
  
  # FDR across thresholds
  fdr_curve <- rbind(
    fdr_results$hand$fdr_curve[, Converter := "Hand-written"],
    fdr_results$llm$fdr_curve[, Converter := "LLM"]
  )
  
  p_fdr_curve <- ggplot(fdr_curve, aes(x = threshold, y = fdr, 
                                        color = Converter, shape = Converter)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "grey60", linewidth = 0.5) +
    scale_color_manual(values = c("Hand-written" = "#2C5F8A", "LLM" = "#8A5F2C")) +
    scale_x_continuous(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2)) +
    labs(x = "Adj. p-value threshold", y = "Empirical FDR",
         title = "FDR Calibration",
         subtitle = "Dashed line = nominal FDR") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  ggsave(file.path(outdir, paste0("fdr_calibration_", timestamp, ".pdf")),
         p_fdr_curve, width = 6, height = 5, dpi = 300)
  
  p_fdr <- p_fdr_bar
}

# --- Plot 6: Venn diagram of significant proteins (adj.p < 0.05) ---
ALPHA <- 0.05
hand_sig <- unique(hand_res[!is.na(adj.pvalue) & adj.pvalue < ALPHA]$Protein)
llm_sig  <- unique(llm_res[!is.na(adj.pvalue) & adj.pvalue < ALPHA]$Protein)

n_both     <- length(intersect(hand_sig, llm_sig))
n_hand_only <- length(setdiff(hand_sig, llm_sig))
n_llm_only  <- length(setdiff(llm_sig, hand_sig))

cat(sprintf("\n  Differential proteins (adj.p < %.2f):\n", ALPHA))
cat(sprintf("    Hand-written: %d  |  LLM: %d  |  Shared: %d\n",
            length(hand_sig), length(llm_sig), n_both))

# Area-proportional Venn (Euler diagram): circle area ∝ set size,
# overlap area ∝ intersection size. Solves circle-circle intersection
# geometry to position the two circles correctly.
make_manual_venn <- function(n_left, n_both, n_right,
                             label_left = "Hand-written",
                             label_right = "LLM",
                             color_left = "#2C5F8A",
                             color_right = "#8A5F2C",
                             title = NULL, subtitle = NULL) {
  size_left  <- n_left  + n_both
  size_right <- n_right + n_both

  # Scale so areas are proportional to set sizes
  unit_area <- max(size_left, size_right)
  r1 <- sqrt(size_left  / unit_area)
  r2 <- sqrt(size_right / unit_area)
  target_overlap_area <- pi * (n_both / unit_area)

  # Solve for distance d between centers s.t. lens area = target_overlap_area
  lens_area <- function(d) {
    if (d >= r1 + r2) return(0)
    if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
    a <- (d^2 + r1^2 - r2^2) / (2 * d * r1)
    b <- (d^2 + r2^2 - r1^2) / (2 * d * r2)
    a <- max(-1, min(1, a)); b <- max(-1, min(1, b))
    r1^2 * acos(a) + r2^2 * acos(b) -
      0.5 * sqrt(max(0, (-d + r1 + r2) * (d + r1 - r2) *
                         (d - r1 + r2) * (d + r1 + r2)))
  }

  d <- if (n_both == 0) {
    r1 + r2 + 0.1
  } else if (n_both >= min(size_left, size_right)) {
    abs(r1 - r2)
  } else {
    tryCatch(
      uniroot(function(d) lens_area(d) - target_overlap_area,
              interval = c(abs(r1 - r2) + 1e-6, r1 + r2 - 1e-6))$root,
      error = function(e) (r1 + r2) * 0.6
    )
  }

  cx_left  <- -d / 2
  cx_right <-  d / 2
  theta <- seq(0, 2 * pi, length.out = 200)
  circles <- rbind(
    data.frame(set = label_left,
               x = cx_left  + r1 * cos(theta), y = r1 * sin(theta)),
    data.frame(set = label_right,
               x = cx_right + r2 * cos(theta), y = r2 * sin(theta))
  )

  # Label x-positions: centroids of each region
  x_left_only  <- cx_left  - r1 * 0.55
  x_right_only <- cx_right + r2 * 0.55
  x_overlap    <- (cx_left + cx_right) / 2
  y_top        <- max(r1, r2)

  ggplot(circles, aes(x = x, y = y, fill = set, group = set)) +
    geom_polygon(alpha = 0.45, color = "grey20", linewidth = 0.5) +
    annotate("text", x = x_left_only,  y = 0, label = n_left,
             size = 5, fontface = "bold") +
    annotate("text", x = x_overlap,    y = 0, label = n_both,
             size = 5, fontface = "bold") +
    annotate("text", x = x_right_only, y = 0, label = n_right,
             size = 5, fontface = "bold") +
    annotate("text", x = cx_left,  y = y_top + 0.15,
             label = label_left,  size = 4.2, fontface = "bold",
             color = color_left) +
    annotate("text", x = cx_right, y = y_top + 0.15,
             label = label_right, size = 4.2, fontface = "bold",
             color = color_right) +
    scale_fill_manual(values = setNames(c(color_left, color_right),
                                        c(label_left, label_right))) +
    coord_fixed(clip = "off") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(color = "grey40", hjust = 0.5),
          legend.position = "none",
          plot.margin = margin(20, 10, 10, 10))
}

p_venn <- tryCatch({
  if (!requireNamespace("eulerr", quietly = TRUE)) stop("eulerr not installed")
  if (!requireNamespace("ggplotify", quietly = TRUE)) stop("ggplotify not installed")
  fit <- eulerr::euler(c(
    `Hand-written`        = n_hand_only,
    `LLM`                 = n_llm_only,
    `Hand-written&LLM`    = n_both
  ))
  ggplotify::as.ggplot(plot(fit,
    fills = list(fill = c("#2C5F8A", "#8A5F2C"), alpha = 0.5),
    edges = list(col = "grey20", lwd = 1),
    labels = list(font = 2, cex = 1.1),
    quantities = list(cex = 1.1, font = 2),
    main = sprintf("Differential proteins (adj.p < %.2f)", ALPHA)
  ))
}, error = function(e) {
  message("  eulerr unavailable (", e$message,
          ") — drawing area-proportional Venn manually")
  make_manual_venn(
    n_left = n_hand_only, n_both = n_both, n_right = n_llm_only,
    title    = sprintf("Differential proteins (adj.p < %.2f)", ALPHA),
    subtitle = sprintf("Hand: %d  |  LLM: %d  |  Shared: %d",
                       length(hand_sig), length(llm_sig), n_both)
  )
})

ggsave(file.path(outdir, paste0("venn_diffproteins_", timestamp, ".pdf")),
       p_venn, width = 5, height = 5, dpi = 300)

# --- Plot 7: Combined poster panel ---
if (!is.null(p_fdr)) {
  p_combined <- cowplot::plot_grid(
    p_summary, p_fc, p_pv, p_se, p_fdr, p_venn,
    ncol = 6, labels = c("A", "B", "C", "D", "E", "F")
  )
  combined_path <- file.path(outdir, paste0("comparison_panel_", timestamp, ".pdf"))
  ggsave(combined_path, p_combined, width = 32, height = 6, dpi = 300)
} else {
  p_combined <- cowplot::plot_grid(
    p_summary, p_fc, p_pv, p_se, p_venn,
    ncol = 5, labels = c("A", "B", "C", "D", "E")
  )
  combined_path <- file.path(outdir, paste0("comparison_panel_", timestamp, ".pdf"))
  ggsave(combined_path, p_combined, width = 27, height = 6, dpi = 300)
}
cat("  Saved:", combined_path, "\n")

# ==================================================================
# Step 7: Save all results
# ==================================================================
saveRDS(
  list(
    trial           = list(model = MODEL, tool = TOOL, prompt = PROMPT, 
                           elapsed_sec = elapsed),
    llm_mapping     = llm_mapping,
    diagnostics     = llm_diagnostics,
    merged_results  = merged,
    merged_summary  = merged_summary,
    fdr_results     = fdr_results,
    diff_proteins   = list(alpha = ALPHA, hand = hand_sig, llm = llm_sig),
    hand_comparison = hand_comparison$ComparisonResult,
    llm_comparison  = llm_comparison$ComparisonResult
  ),
  file.path(outdir, paste0("full_comparison_", timestamp, ".rds"))
)

cat("\n=== Done ===\n")
cat(sprintf("Shared proteins compared: %d\n", uniqueN(merged$Protein)))
cat(sprintf("log2FC correlation: %.4f\n", 
            cor(merged$log2FC_hand, merged$log2FC_llm, use = "complete.obs")))
if (!is.null(fdr_results)) {
  cat(sprintf("Empirical FDR — Hand: %.4 f, LLM: %.4f\n",
              fdr_results$hand$empirical_fdr, fdr_results$llm$empirical_fdr))
}