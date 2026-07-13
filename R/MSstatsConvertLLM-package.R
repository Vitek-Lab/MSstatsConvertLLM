#' MSstatsConvertLLM: LLM-driven schema inference and conversion for MSstats
#'
#' Uses locally-hosted large language models to map arbitrary proteomics tool
#' output onto the MSstats analysis schema, then feeds that mapping into the
#' MSstatsConvert pipeline via [LLMtoMSstatsFormat()]. Also provides a
#' benchmarking harness ([run_trial()], [score_trial()], [print_scorecard()])
#' for evaluating model / prompt combinations against curated ground truth.
#'
#' @keywords internal
#' @import data.table
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom utils head
"_PACKAGE"

# data.table uses non-standard evaluation, so the column names referenced bare
# inside `[...]` would otherwise trip R CMD check's "no visible binding" note.
utils::globalVariables(c(
  "Intensity", "PrecursorCharge", "ProductCharge", "Qvalue",
  "PeptideSequence", "ProteinName", "N"
))
