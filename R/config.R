# ============================================================
# config.R — Models, schemas, ground truth, and helper functions
# ============================================================

# ------------------------------------------------------------------
# 1. Model registry — add new models here
# ------------------------------------------------------------------

#' Registry of benchmark models
#'
#' A named list of model specifications used by [create_chat()] and
#' [run_trial()]. Each entry has `provider` (`"ollama"` or `"anthropic"`),
#' `model` (the provider-specific model id), and a human-readable `label`.
#' Edit this object in `R/config.R` (and reinstall / `devtools::load_all()`) to
#' add models to the benchmark.
#'
#' @export
MODEL_REGISTRY <- list(
  "llama3.2-3b" = list(
    provider = "ollama",
    model    = "llama3.2:3b-instruct-q4_K_M",
    label    = "Llama 3.2 3B"
  ),
  "llama3.1-8b" = list(
    provider = "ollama",
    model    = "llama3.1:8b",
    label    = "Llama 3.1 8B"
  ),
  "deepseek-r1-14b" = list(
    provider = "ollama",
    model    = "deepseek-r1:14b",
    label    = "DeepSeek-R1 14B"
  )
)

# ------------------------------------------------------------------
# 2. MSstats target schema
# ------------------------------------------------------------------

#' The eight MSstats schema fields, in canonical order
#' @export
MSSTATS_FIELDS <- c(
  "ProteinName", "PeptideSequence", "PrecursorCharge",
  "FragmentIon", "ProductCharge", "Run", "Intensity", "Qvalue"
)

# ------------------------------------------------------------------
# 3. Ground truth mappings (from MSstatsConvert converter logic)
#    Each entry: MSstats field -> expected source column name
#    NULL = field not present / must be inferred or NA
# ------------------------------------------------------------------

#' Ground-truth column mappings per tool
#'
#' A named list keyed by tool. Each value is a list mapping every MSstats field
#' to the expected source column name(s): a single string, a character vector
#' (first element is the primary answer, the rest are accepted alternates), or
#' `NULL` when the field is not present in that tool's output. Used by
#' [score_trial()].
#'
#' @export
GROUND_TRUTH <- list(
  spectronaut = list(
    ProteinName     = c("PG.ProteinGroups", "PG.ProteinAccessions"),
    PeptideSequence = c("EG.ModifiedSequence", "PEP.StrippedSequence"),
    PrecursorCharge = "FG.Charge",
    FragmentIon     = "F.FrgIon",
    ProductCharge   = "F.Charge",
    Run             = "R.FileName",
    Intensity       = c("F.PeakArea", "F.NormalizedPeakArea", 
                        "F.PeakHeight", "F.NormalizedPeakHeight"),
    Qvalue          = c("EG.Qvalue", "PG.Qvalue")
  ),
  proteome_discoverer = list(
    ProteinName     = "Protein.Group.Accessions",
    PeptideSequence = "Sequence",
    PrecursorCharge = "Charge",
    FragmentIon     = NULL,
    ProductCharge   = NULL,
    Run             = "Spectrum.File",
    Intensity       = "Intensity",
    Qvalue          = NULL
  ),
  metamorpheus = list(
    ProteinName     = "Protein Group",
    PeptideSequence = c("Full Sequence", "Base Sequence"),
    PrecursorCharge = c("Precursor Charge", "Peak Charge"),
    FragmentIon     = NULL,
    ProductCharge   = NULL,
    Run             = "File Name",
    Intensity       = "Peak intensity",
    Qvalue          = "PIP Q-Value"
  ),
  # ---------------------------------------------------------------
  # Add new tools here, e.g.:
  # diann_v18 = list(...)
  #
  # Example with alternates:
  #   ProteinName = c("Protein.Group", "Master.Protein", "Protein.Accessions")
  # ---------------------------------------------------------------
  diann = list(
    ProteinName     = c("Protein.Names", "Protein.Group", "Protein.Ids", "Genes"),
    PeptideSequence = c("Modified.Sequence", "Stripped.Sequence"),
    PrecursorCharge = c("Precursor.Charge"),
    FragmentIon     = NULL,
    ProductCharge   = NULL,
    Run             = "File.Name",
    Intensity       = c("Fragment.Quant.Raw", "Fragment.Quant.Corrected"),
    Qvalue          = c("Q.Value", "Protein.Q.Value", "Global.PG.Q.Value", "Global.Q.Value")
  )
)

# ------------------------------------------------------------------
# 3b. Default acquisition type per tool
#     Used to inform the LLM which fields to expect
# ------------------------------------------------------------------

#' Default acquisition type (`"DIA"` / `"DDA"`) per tool
#' @export
TOOL_ACQUISITION <- list(
  spectronaut         = "DIA",
  proteome_discoverer = "DDA",
  metamorpheus        = "DDA",
  diann               = "DIA"
)

# ------------------------------------------------------------------
# 3c. Whether to enable packed-column transform detection per tool
# ------------------------------------------------------------------

#' Whether to enable packed-column transform detection per tool
#' @export
TOOL_TRANSFORMS <- list(
  spectronaut         = FALSE,
  proteome_discoverer = FALSE,
  metamorpheus        = FALSE,
  diann               = TRUE
)

# ------------------------------------------------------------------
# 4. Test dataset loader
# ------------------------------------------------------------------

#' Load a registered benchmark test dataset
#'
#' Reads a proteomics tool's example output (bundled with `MSstatsConvert`) and
#' normalizes its column headers.
#'
#' @param tool_name Character key, one of the tools registered in the function
#'   body (e.g. `"spectronaut"`, `"proteome_discoverer"`, `"metamorpheus"`,
#'   `"diann"`).
#' @return A `data.table` of the raw tool output with normalized headers.
#' @export
load_test_dataset <- function(tool_name) {
  paths <- list(
    spectronaut = system.file(
      "tinytest/raw_data/Spectronaut/spectronaut_input.csv",
      package = "MSstatsConvert"
    ),
    proteome_discoverer = system.file(
      "tinytest/raw_data/PD/pd_input.csv",
      package = "MSstatsConvert"
    ),
    metamorpheus = system.file(
      "tinytest/raw_data/Metamorpheus/QuantifiedPeaks.tsv",
      package = "MSstatsConvert"
    ),
    diann = system.file(
      "tinytest/raw_data/DIANN/diann_input.tsv",
      package = "MSstatsConvert"
    )
  )
  
  path <- paths[[tool_name]]
  if (is.null(path) || !nzchar(path)) {
    stop("No test data registered for tool: ", tool_name,
         "\nAvailable: ", paste(names(paths), collapse = ", "))
  }
  
  normalize_headers(data.table::fread(path))
}

# ------------------------------------------------------------------
# 5. Utility functions
# ------------------------------------------------------------------

#' Canonical column-name normalizer
#'
#' Used for both raw data-frame headers and [GROUND_TRUTH] values so the LLM,
#' the data, and the expected answers all speak the same dialect.
#'
#' @param x Character vector of column names.
#' @return Normalized character vector.
#' @export
normalize_colname <- function(x) {
  x <- gsub("[^A-Za-z0-9._]+", ".", x)
  x <- gsub("\\.+", ".", x)
  x <- sub("^\\.|\\.$", "", x)
  x
}

#' Normalize the headers of a data frame in place
#' @param df A data frame / data.table.
#' @return `df` with normalized column names.
#' @export
normalize_headers <- function(df) {
  names(df) <- normalize_colname(names(df))
  df
}

# Run GROUND_TRUTH values through the same normalizer so expected names
# match the headers the LLM will see.
GROUND_TRUTH <- lapply(GROUND_TRUTH, function(tool_truth) {
  lapply(tool_truth, function(v) if (is.null(v)) NULL else normalize_colname(v))
})

#' Build a compact JSON preview of a data frame's columns
#'
#' Produces one entry per column with a few example values, used as the
#' `DATA HEADER` block injected into the LLM prompt.
#'
#' @param df A data frame / data.table.
#' @param n_rows Number of example rows/values per column.
#' @param max_chars Truncate each example value to this many characters.
#' @return A JSON string.
#' @export
make_json_preview <- function(df, n_rows = 3, max_chars = 50) {
  df <- head(df, n_rows)
  profiles <- lapply(names(df), function(col) {
    vals <- as.character(df[[col]])
    vals <- substr(vals, 1L, max_chars)
    vals <- vals[!is.na(vals) & vals != ""]
    vals <- unique(vals)
    vals <- head(vals, n_rows)
    list(column = col, examples = vals)
  })
  jsonlite::toJSON(profiles, auto_unbox = TRUE, pretty = TRUE)
}

# ------------------------------------------------------------------
# 6. Structured output types (for ellmer chat_structured)
#
# These are built lazily by the functions below so the package can be
# installed and loaded without the (optional) 'ellmer' dependency. They are
# only evaluated from within run_trial(), which checks for ellmer first.
# ------------------------------------------------------------------

#' @keywords internal
.mapping_entry_type <- function() {
  ellmer::type_object(
    field      = ellmer::type_string(),
    from       = ellmer::type_string(),
    confidence = ellmer::type_number(),
    candidates = ellmer::type_array(
      ellmer::type_object(from = ellmer::type_string(),
                          score = ellmer::type_number())
    )
  )
}

#' @keywords internal
.filter_entry_type <- function() {
  ellmer::type_object(
    column      = ellmer::type_string(),   # QC column name
    dtype       = ellmer::type_string(),   # "boolean", "string", "numeric"
    operation   = ellmer::type_string(),   # "equals", "not_equals", "less_than", etc.
    value       = ellmer::type_string(),   # string representation of filter value
    description = ellmer::type_string(),
    confidence  = ellmer::type_number()
  )
}

# Transform type for packed columns (e.g., DIA-NN Fragment.Quant.Corrected)
#' @keywords internal
.transform_entry_type <- function() {
  ellmer::type_object(
    field         = ellmer::type_string(),  # MSstats field produced (e.g., "Intensity")
    source_column = ellmer::type_string(),  # packed column name in data
    action        = ellmer::type_string(),  # "split_long"
    delimiter     = ellmer::type_string(),  # separator character
    id_field      = ellmer::type_string(),  # MSstats field for generated IDs
    id_prefix     = ellmer::type_string()   # prefix for generated IDs (e.g., "Frag")
  )
}

#' @keywords internal
.mapping_result_type <- function() {
  ellmer::type_object(
    "LLM mapping proposal for MSstats",
    mappings   = ellmer::type_array(.mapping_entry_type()),
    confidence = ellmer::type_number(),
    notes      = ellmer::type_array(ellmer::type_string()),
    warnings   = ellmer::type_array(ellmer::type_string())
  )
}

#' @keywords internal
.mapping_with_filters_type <- function() {
  ellmer::type_object(
    "LLM mapping proposal for MSstats with filter discovery",
    mappings   = ellmer::type_array(.mapping_entry_type()),
    filters    = ellmer::type_array(.filter_entry_type()),
    confidence = ellmer::type_number(),
    notes      = ellmer::type_array(ellmer::type_string()),
    warnings   = ellmer::type_array(ellmer::type_string())
  )
}

# Full type: mappings + filters + transforms
#' @keywords internal
.mapping_full_type <- function() {
  ellmer::type_object(
    "LLM mapping proposal for MSstats with filters and transforms",
    mappings   = ellmer::type_array(.mapping_entry_type()),
    filters    = ellmer::type_array(.filter_entry_type()),
    transforms = ellmer::type_array(.transform_entry_type()),
    confidence = ellmer::type_number(),
    notes      = ellmer::type_array(ellmer::type_string()),
    warnings   = ellmer::type_array(ellmer::type_string())
  )
}

#' Select the appropriate structured-output type for a trial
#' @param has_filters logical; TRUE for filter-aware prompts.
#' @param allow_transforms logical; TRUE to include packed-column transforms.
#' @return An `ellmer` type object.
#' @keywords internal
.select_mapping_type <- function(has_filters, allow_transforms) {
  if (allow_transforms) {
    # transforms (with or without filters) use the full type; filters may be empty
    .mapping_full_type()
  } else if (has_filters) {
    .mapping_with_filters_type()
  } else {
    .mapping_result_type()
  }
}