# ============================================================
# run_all.R — Full benchmark matrix
#
# Usage:
#   Rscript run_all.R                      # fixtures, local models
#   Rscript run_all.R --source benchmark   # deposited datasets
#   Rscript run_all.R --models llama3.1-8b --tools spectronaut
#   Rscript run_all.R --models claude-sonnet
#   Rscript run_all.R --tools spectronaut,proteome_discoverer,metamorpheus
#
# A hosted model, such as claude-sonnet, runs only when --models names it.
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

# Default scope is the three formats reported in the paper. diann is still
# available via --tools but is not run by default.
PAPER_TOOLS <- c("spectronaut", "proteome_discoverer", "metamorpheus")

# Only the local models run by default. A hosted model sends the preview rows
# of each tool output to its provider, so it runs only when --models names it.
is_hosted     <- function(key) MODEL_REGISTRY[[key]]$provider != "ollama"
provider_of   <- function(key) tools::toTitleCase(MODEL_REGISTRY[[key]]$provider)
HOSTED_MODELS <- Filter(is_hosted, names(MODEL_REGISTRY))

selected_models  <- parse_flag("--models",  setdiff(names(MODEL_REGISTRY), HOSTED_MODELS))
selected_tools   <- parse_flag("--tools",   PAPER_TOOLS)
selected_prompts <- parse_flag("--prompts", names(PROMPT_VERSIONS))
n_reps           <- as.integer(parse_flag("--reps", "3"))

# "fixture" uses the example files bundled with MSstatsConvert; "benchmark"
# uses the deposited datasets in benchmark_data/.
selected_source  <- parse_flag("--source", "fixture")

unknown <- setdiff(selected_models, names(MODEL_REGISTRY))
if (length(unknown) > 0) {
  stop("Unknown model: ", paste(unknown, collapse = ", "),
       "\nAvailable: ", paste(names(MODEL_REGISTRY), collapse = ", "),
       call. = FALSE)
}
describe_model <- function(key) {
  sprintf("%s (%s)", key,
          if (is_hosted(key)) paste("hosted by", provider_of(key)) else "local")
}

cat("========================================\n")
cat("LLM-MSstats Schema Inference Benchmark\n")
cat("========================================\n")
cat("Models:  ", paste(vapply(selected_models, describe_model, ""),
                       collapse = ", "), "\n")
if (!"--models" %in% args && length(HOSTED_MODELS) > 0) {
  cat("Not run: ", paste(HOSTED_MODELS, collapse = ", "),
      "(hosted; name it with --models to run it)\n")
}
cat("Tools:   ", paste(selected_tools, collapse = ", "), "\n")
cat("Prompts: ", paste(selected_prompts, collapse = ", "), "\n")
cat("Source:  ", selected_source, "\n")
cat("Reps:    ", n_reps, "\n")
cat("Total trials:", length(selected_models) * length(selected_tools) *
      length(selected_prompts) * n_reps, "\n")
cat("========================================\n\n")

n_preview <- formals(run_benchmark_trial)$n_preview
for (key in Filter(is_hosted, selected_models)) {
  warning(sprintf(paste("%s is a hosted model: its calls send the column names",
                        "and the first %d rows of each tool output to %s."),
                  key, n_preview, provider_of(key)),
          call. = FALSE, immediate. = TRUE)
}

# ------------------------------------------------------------------
# Output paths, opened before the run
#
# Scores are appended after every trial rather than held in memory until the
# end. A full matrix is hours long, and a dropped tunnel used to cost the
# entire run.
#
# The source is part of the file name so a benchmark run never overwrites a
# fixture run, and the two can be told apart later.
# ------------------------------------------------------------------
outdir <- "results"
dir.create(outdir, showWarnings = FALSE)

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag         <- paste0(selected_source, "_", timestamp)
scores_path <- file.path(outdir, paste0("field_scores_", tag, ".csv"))

cat("Appending scores to:", scores_path, "\n\n")

# ------------------------------------------------------------------
# Run the matrix
# ------------------------------------------------------------------
all_results <- list()
all_scores  <- list()
idx <- 0

for (model_key in selected_models) {

  # Warm-up: the model loads into VRAM on its first call, so without this the
  # first cell of each model block runs cold and the rest run warm. Results are
  # discarded.
  cat(sprintf("[warmup] %s ... ", model_key))
  warm <- tryCatch({
    run_benchmark_trial(model_key, selected_tools[1], selected_prompts[1],
                        acquisition = TOOL_ACQUISITION[[selected_tools[1]]],
                        allow_transforms = isTRUE(TOOL_TRANSFORMS[[selected_tools[1]]]),
                        source = selected_source)
  }, error = function(e) NULL)
  cat(if (is.null(warm)) "failed\n\n" else sprintf("%.1fs\n\n", warm$elapsed_sec))

  for (tool_name in selected_tools) {
    for (prompt_key in selected_prompts) {

      # The first call after the prompt prefix changes does a full prefill;
      # later identical calls reuse the cached prefix. The two paths are not
      # bit-identical, which is enough to flip a near-tie between two candidate
      # columns. Without this, rep 1 of every block is not comparable to the
      # rest. Discarded.
      tryCatch(
        run_benchmark_trial(model_key, tool_name, prompt_key,
                            acquisition = TOOL_ACQUISITION[[tool_name]],
                            allow_transforms = isTRUE(TOOL_TRANSFORMS[[tool_name]]),
                            source = selected_source),
        error = function(e) NULL)

      for (rep in seq_len(n_reps)) {
        idx <- idx + 1
        cat(sprintf("[%d] %s | %s | %s | rep %d ... ",
                    idx, model_key, tool_name, prompt_key, rep))
        
        trial <- tryCatch(
          run_benchmark_trial(model_key, tool_name, prompt_key,
                              acquisition = TOOL_ACQUISITION[[tool_name]],
                              allow_transforms = isTRUE(TOOL_TRANSFORMS[[tool_name]]),
                              source = selected_source),
          error = function(e) {
            # On one line and without colour codes, so that it fits in one
            # field of the CSV
            msg <- gsub("\\s*\n\\s*", " ", cli::ansi_strip(conditionMessage(e)))
            message("FAILED: ", msg)
            list(model = model_key, tool = tool_name, prompt = prompt_key,
                 elapsed_sec = NA, timestamp = Sys.time(), error = msg)
          }
        )
        failed <- !is.null(trial$error)

        cat(if (failed) "failed\n" else sprintf("%.1fs\n", trial$elapsed_sec))

        # A trial whose model call failed has no answer to score. One row records
        # the error in place of the field rows, and the trial is left out of the
        # summary and of the saved trials.
        scores <- if (failed) {
          data.table(
            field = NA_character_, predicted = NA_character_,
            expected = NA_character_, n_acceptable = NA_integer_,
            exact_match = NA, acceptable_match = NA, candidate_hit = NA,
            confidence = NA_real_, has_candidates = NA, n_candidates = NA_integer_
          )
        } else {
          tryCatch(
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
        }
        
        scores[, `:=`(model = model_key, tool = tool_name,
                      prompt = prompt_key, rep = rep,
                      source = selected_source,
                      elapsed_sec = trial$elapsed_sec,
                      status = if (failed) "error" else "scored",
                      error = if (failed) trial$error else NA_character_)]
        
        if (!failed) print_scorecard(scores, trial)

        # Crash-safe: on disk before the next trial starts.
        fwrite(scores, scores_path, append = file.exists(scores_path))
        
        if (!failed) {
          all_results[[length(all_results) + 1]] <-
            c(trial, list(rep = rep, source = selected_source))
        }
        all_scores[[idx]]  <- scores
      }
    }
  }
}

# ------------------------------------------------------------------
# Aggregate and save
# ------------------------------------------------------------------
results_dt <- rbindlist(all_scores, fill = TRUE)

# Summary table, over the trials that were scored
summary_dt <- results_dt[status == "scored", .(
  exact_accuracy      = mean(exact_match, na.rm = TRUE),
  acceptable_accuracy = mean(acceptable_match, na.rm = TRUE),
  candidate_rate      = mean(candidate_hit, na.rm = TRUE),
  mean_conf           = mean(confidence, na.rm = TRUE),
  mean_time_sec       = mean(elapsed_sec, na.rm = TRUE),
  n_trials            = .N / length(MSSTATS_FIELDS)
), by = .(model, tool, prompt, source)]

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")
print(summary_dt)

n_failed <- results_dt[status == "error", .N]
if (n_failed > 0) {
  cat(sprintf(paste("\n%d of %d trials failed and are left out of the summary.",
                    "Their errors are in %s.\n"), n_failed, idx, scores_path))
}

fwrite(summary_dt, file.path(outdir, paste0("summary_", tag, ".csv")))

# Save raw LLM outputs first. Filter parsing below depends on the model
# returning well-formed objects, which it does not always do; if that throws,
# the raw trials are already on disk.
saveRDS(all_results, file.path(outdir, paste0("raw_trials_", tag, ".rds")))

# Extract filter results from trials that used filter-aware prompts
filter_rows <- list()

# Models omit keys. Missing becomes NA rather than an error that ends the run.
fget <- function(f, key) {
  v <- tryCatch(f[[key]], error = function(e) NULL)
  if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1]
}

for (trial in all_results) {
  if (!is.null(trial$mapping$filters)) {
    fl <- trial$mapping$filters
    trial_rep <- if (!is.null(trial$rep)) trial$rep else NA_integer_
    if (is.data.frame(fl) && nrow(fl) > 0) {
      fl$model  <- trial$model
      fl$tool   <- trial$tool
      fl$prompt <- trial$prompt
      fl$rep    <- trial_rep
      fl$source <- selected_source
      filter_rows[[length(filter_rows) + 1]] <- fl
    } else if (is.list(fl) && length(fl) > 0) {
      for (f in fl) {
        filter_rows[[length(filter_rows) + 1]] <- data.table(
          model = trial$model, tool = trial$tool, prompt = trial$prompt,
          rep = trial_rep, source = selected_source,
          column = fget(f, "column"), dtype = fget(f, "dtype"),
          operation = fget(f, "operation"), value = fget(f, "value"),
          description = fget(f, "description"),
          confidence = suppressWarnings(as.numeric(fget(f, "confidence")))
        )
      }
    }
  }
}
if (length(filter_rows) > 0) {
  filters_dt <- rbindlist(filter_rows, fill = TRUE)
  fwrite(filters_dt, file.path(outdir, paste0("filters_", tag, ".csv")))
  cat("  filters_",     tag, ".csv  (LLM-suggested filters)\n")
}

cat("\nResults saved to: ", outdir, "/\n")
cat("  field_scores_", tag, ".csv  (per-field detail)\n")
cat("  summary_",      tag, ".csv  (aggregate)\n")
cat("  raw_trials_",   tag, ".rds  (full LLM outputs)\n")