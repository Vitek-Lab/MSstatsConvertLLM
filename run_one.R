# ============================================================
# run_one.R — Quick interactive test of a single combination
#
# Use this while iterating on prompts or testing a new model.
# Install the package first (devtools::install(".")) or, during
# development, run devtools::load_all(".") instead of library().
# ============================================================

library(MSstatsConvertLLM)
library(jsonlite)  # for toJSON() in the summary below

# ----- Change these to test different combos -----
MODEL   <- "llama3.2-3b" # llama3.2-3b, llama3.1-8b, deepseek-r1-14b, claude-sonnet
TOOL    <- "diann"    # spectronaut, proteome_discoverer, metamorpheus, diann
PROMPT  <- "constrained_filter"    # "lean" or "constrained" or "constrained_filter"
ACQUISITION <- TOOL_ACQUISITION[[TOOL]]  # "DIA", "DDA", or NULL
TRANSFORMS  <- TOOL_TRANSFORMS[[TOOL]]   # TRUE for DIA-NN, FALSE for others
# --------------------------------------------------
 
cat("Running:", MODEL, "|", TOOL, "|", PROMPT, "|", ACQUISITION, 
    "| transforms:", TRANSFORMS, "\n\n")
 
trial  <- run_trial(MODEL, TOOL, PROMPT, 
                    acquisition = ACQUISITION,
                    allow_transforms = TRANSFORMS)
scores <- score_trial(trial)
print_scorecard(scores, trial)
 
# Inspect the raw mapping
cat("\nRaw LLM mapping output:\n")
cat(toJSON(trial$mapping, auto_unbox = TRUE, pretty = TRUE), "\n")