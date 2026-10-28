# make_new_format_mapping.R - infer the mapping the new-format vignette reads
#
# Needs Ollama running and deepseek-r1-14b pulled. Run from the repository
# root, after inst/scripts/make_new_format_data.R has written the subset:
#
#   Rscript inst/scripts/make_new_format_mapping.R
#
# Point MSSTATS_OLLAMA_URL at the server if it is not on localhost:11434.

suppressPackageStartupMessages({
  library(data.table)
  library(MSstatsConvertLLM)
})

MODEL  <- "deepseek-r1-14b"
PROMPT <- "filter_aware"
INPUT  <- "inst/extdata/AlphaPept_subset.csv"
OUTPUT <- "inst/extdata/new_format_mapping.rds"

stopifnot(file.exists(INPUT))
input <- fread(INPUT)

cat(sprintf("Ollama at %s, model %s, prompt %s\n", OLLAMA_URL, MODEL, PROMPT))

# The first call after a change of prompt is the odd one out: the model has the
# prompt to digest as well as the data, and tends to answer differently from
# every call after it. Make it and throw it away, so what gets saved is a
# settled answer.
cat("Warm-up call, discarded ...\n")
invisible(infer_mapping(input, MODEL, PROMPT, acquisition = "DDA"))

cat("Inferring the mapping ...\n")
mapping <- infer_mapping(input, MODEL, PROMPT, acquisition = "DDA")

cat(sprintf("Model took %.1f s\n", attr(mapping, "elapsed_sec")))
print(as.data.frame(mapping$mappings)[, c("field", "from")])
if (NROW(mapping$filters) > 0) {
  cat("\nProposed filters:\n")
  print(as.data.frame(mapping$filters)[, c("column", "operation", "value")])
} else {
  cat("\nNo filters proposed.\n")
}

validation <- validate_mapping(mapping)
cat("\nvalidate_mapping passed:", validation$passed, "\n")
if (length(validation$issues)) cat(paste0("  ", validation$issues, collapse = "\n"), "\n")

saveRDS(mapping, OUTPUT)
cat("\nWritten:", OUTPUT, "\n")
