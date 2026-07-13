# ============================================================
# run_all.R — Full benchmark matrix
#
# Usage:
#   Rscript run_all.R                  # run everything
#   Rscript run_all.R --models llama3.1-8b --tools spectronaut
#
# Install the package first (devtools::install(".")) or, during
# development, replace library() with devtools::load_all(".").
# ============================================================

library(MSstatsConvertLLM)
library(data.table)  # for rbindlist/fwrite/data.table used in this script

# ------------------------------------------------------------------
# CLI overrides (optional)
# ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_flag <- function(flag, default) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  val <- args[idx + 1]
  strsplit(val, ",")[[1]]
}

selected_models  <- parse_flag("--models",  names(MODEL_REGISTRY))
selected_tools   <- parse_flag("--tools",   names(GROUND_TRUTH))
selected_prompts <- parse_flag("--prompts", names(PROMPT_VERSIONS))
n_reps           <- as.integer(parse_flag("--reps", "5"))

cat("========================================\n")
cat("LLM-MSstats Schema Inference Benchmark\n")
cat("========================================\n")
cat("Models:  ", paste(selected_models, collapse = ", "), "\n")
cat("Tools:   ", paste(selected_tools, collapse = ", "), "\n")
cat("Prompts: ", paste(selected_prompts, collapse = ", "), "\n")
cat("Reps:    ", n_reps, "\n")
cat("Total trials:", length(selected_models) * length(selected_tools) *
      length(selected_prompts) * n_reps, "\n")
cat("========================================\n\n")

# ------------------------------------------------------------------
# Run the matrix
# ------------------------------------------------------------------
all_results <- list()
all_scores  <- list()
idx <- 0

for (model_key in selected_models) {
  for (tool_name in selected_tools) {
    for (prompt_key in selected_prompts) {
      for (rep in seq_len(n_reps)) {
        idx <- idx + 1
        cat(sprintf("[%d] %s | %s | %s | rep %d ... ",
                    idx, model_key, tool_name, prompt_key, rep))
        
        trial <- tryCatch(
          run_trial(model_key, tool_name, prompt_key,
                    acquisition = TOOL_ACQUISITION[[tool_name]],
                    allow_transforms = isTRUE(TOOL_TRANSFORMS[[tool_name]])),
          error = function(e) {
            message("FAILED: ", e$message)
            list(
              model = model_key, tool = tool_name, prompt = prompt_key,
              mapping = list(mappings = list()), elapsed_sec = NA,
              timestamp = Sys.time(), error = e$message
            )
          }
        )
        
        cat(sprintf("%.1fs\n", trial$elapsed_sec))
        
        scores <- tryCatch(
          score_trial(trial),
          error = function(e) {
            message("  Scoring failed: ", e$message)
            data.table(
              field = MSSTATS_FIELDS, predicted = NA, expected = NA,
              n_acceptable = 0L, exact_match = FALSE, 
              acceptable_match = FALSE, candidate_hit = FALSE,
              confidence = NA, has_candidates = FALSE, n_candidates = 0L
            )
          }
        )
        
        scores[, `:=`(model = model_key, tool = tool_name,
                      prompt = prompt_key, rep = rep,
                      elapsed_sec = trial$elapsed_sec)]
        
        print_scorecard(scores, trial)
        
        all_results[[idx]] <- c(trial, list(rep = rep))
        all_scores[[idx]]  <- scores
      }
    }
  }
}

# ------------------------------------------------------------------
# Aggregate and save
# ------------------------------------------------------------------
results_dt <- rbindlist(all_scores, fill = TRUE)

# Summary table
summary_dt <- results_dt[, .(
  exact_accuracy      = mean(exact_match, na.rm = TRUE),
  acceptable_accuracy = mean(acceptable_match, na.rm = TRUE),
  candidate_rate      = mean(candidate_hit, na.rm = TRUE),
  mean_conf           = mean(confidence, na.rm = TRUE),
  mean_time_sec       = mean(elapsed_sec, na.rm = TRUE),
  n_trials            = .N / length(MSSTATS_FIELDS)
), by = .(model, tool, prompt)]

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")
print(summary_dt)

# Save outputs
outdir <- "results"
dir.create(outdir, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
fwrite(results_dt, file.path(outdir, paste0("field_scores_", timestamp, ".csv")))
fwrite(summary_dt, file.path(outdir, paste0("summary_", timestamp, ".csv")))

# Extract filter results from trials that used filter-aware prompts
filter_rows <- list()
for (trial in all_results) {
  if (!is.null(trial$mapping$filters)) {
    fl <- trial$mapping$filters
    trial_rep <- if (!is.null(trial$rep)) trial$rep else NA_integer_
    if (is.data.frame(fl) && nrow(fl) > 0) {
      fl$model  <- trial$model
      fl$tool   <- trial$tool
      fl$prompt <- trial$prompt
      fl$rep    <- trial_rep
      filter_rows[[length(filter_rows) + 1]] <- fl
    } else if (is.list(fl) && length(fl) > 0) {
      for (f in fl) {
        filter_rows[[length(filter_rows) + 1]] <- data.table(
          model = trial$model, tool = trial$tool, prompt = trial$prompt,
          rep = trial_rep,
          column = f[["column"]], dtype = f[["dtype"]],
          operation = f[["operation"]], value = as.character(f[["value"]]),
          description = f[["description"]],
          confidence = as.numeric(f[["confidence"]])
        )
      }
    }
  }
}
if (length(filter_rows) > 0) {
  filters_dt <- rbindlist(filter_rows, fill = TRUE)
  fwrite(filters_dt, file.path(outdir, paste0("filters_", timestamp, ".csv")))
  cat("  filters_",     timestamp, ".csv  (LLM-suggested filters)\n")
}

# Save raw LLM outputs for auditing
saveRDS(all_results, file.path(outdir, paste0("raw_trials_", timestamp, ".rds")))

cat("\nResults saved to: ", outdir, "/\n")
cat("  field_scores_", timestamp, ".csv  (per-field detail)\n")
cat("  summary_",      timestamp, ".csv  (aggregate)\n")
cat("  raw_trials_",   timestamp, ".rds  (full LLM outputs)\n")