# ============================================================
# LLMtoMSstatsFormat.R — Generic converter driven by LLM schema mapping
#
# The LLM supplies the rename map and any discovered filters. Everything
# after that is the same MSstatsConvert pipeline the hand-written
# converters use.
# ============================================================

#' Convert arbitrary proteomics output to MSstats format using LLM mapping
#'
#' @param input data.table of raw tool output
#' @param annotation data.table with Run, Condition, BioReplicate columns, or
#'   NULL to extract them from `input`
#' @param llm_mapping list — parsed LLM output with $mappings and optionally $filters
#' @param useUniquePeptide logical, remove shared peptides (default TRUE)
#' @param removeFewMeasurements logical (default TRUE)
#' @param removeProtein_with1Peptide logical (default FALSE)
#' @param summaryforMultipleRows function for summarizing multiple PSMs (default max)
#' @param filter_with_Qvalue logical (default FALSE, as in the converters). TRUE
#'   fills intensities above `qvalue_cutoff` with NA, treated as censored.
#' @param qvalue_cutoff numeric (default 0.01, as in the converters)
#' @param use_log_file logical (default TRUE)
#' @param append logical (default FALSE)
#' @param verbose logical (default TRUE)
#' @param log_file_path character or NULL
#'
#' @return data.table in MSstats format, ready for dataProcess()
#' @export
LLMtoMSstatsFormat <- function(
    input,
    annotation = NULL,
    llm_mapping,
    useUniquePeptide = TRUE,
    removeFewMeasurements = TRUE,
    removeProtein_with1Peptide = FALSE,
    summaryforMultipleRows = max,
    filter_with_Qvalue = FALSE,
    qvalue_cutoff = 0.01,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path)
  
  input <- data.table::copy(input)
  
  diagnostics <- list(
    hallucinated_columns = data.table(
      field = character(0), 
      hallucinated_name = character(0)
    ),
    skipped_filters = data.table(
      column = character(0), 
      operation = character(0), 
      value = character(0), 
      reason = character(0)
    ),
    applied_filters = data.table(
      column = character(0), 
      operation = character(0), 
      value = character(0), 
      kind = character(0)
    )
  )
  
  # ---------------------------------------------------------------
  # 1. Parse the LLM mapping into a column rename map
  # ---------------------------------------------------------------
  mapping <- .parseLLMMapping(llm_mapping)
  
  if (verbose) {
    message("== LLM Column Mapping ==")
    for (field in names(mapping$rename_map)) {
      message(sprintf("  %s <- %s", field, mapping$rename_map[[field]]))
    }
    if (length(mapping$fill_columns) > 0) {
      message("  Columns to fill with NA: ", 
              paste(mapping$fill_columns, collapse = ", "))
    }
  }
  
  # ---------------------------------------------------------------
  # 2. Translate LLM filters into MSstatsPreprocess configs
  #    Handed over with behavior = "fill", not applied here: a failing
  #    row keeps its place with a blank intensity, so MSstats treats it
  #    as censored rather than absent.
  # ---------------------------------------------------------------
  filter_cfg <- .buildFilterConfigs(input, llm_mapping, verbose)
  diagnostics$skipped_filters <- filter_cfg$skipped
  diagnostics$applied_filters <- filter_cfg$applied
  
  # ---------------------------------------------------------------
  # 3. Select and rename columns to MSstats names
  # ---------------------------------------------------------------
  source_cols <- unlist(mapping$rename_map)
  missing_in_data <- setdiff(source_cols, colnames(input))
  if (length(missing_in_data) > 0) {
    if (verbose) {
      message(sprintf("  WARNING: LLM hallucinated columns not in data: %s",
                      paste(missing_in_data, collapse = ", ")))
      message("  -> These fields will be filled with NA instead.")
    }
    for (field in names(mapping$rename_map)) {
      if (mapping$rename_map[[field]] %in% missing_in_data) {
        diagnostics$hallucinated_columns <- rbind(
          diagnostics$hallucinated_columns,
          data.table(field = field, 
                     hallucinated_name = mapping$rename_map[[field]])
        )
        mapping$fill_columns <- c(mapping$fill_columns, field)
      }
    }
    mapping$rename_map <- mapping$rename_map[
      !mapping$rename_map %in% missing_in_data
    ]
    source_cols <- unlist(mapping$rename_map)
  }
  
  # Filter columns must survive the subset; MSstatsPreprocess drops them.
  annot_candidates <- grep("Condition|Replicate", colnames(input),
                           value = TRUE, ignore.case = TRUE)
  keep_cols <- unique(c(source_cols, annot_candidates, filter_cfg$keep_cols))
  keep_cols <- intersect(keep_cols, colnames(input))
  input <- input[, keep_cols, with = FALSE]
  
  for (msstats_name in names(mapping$rename_map)) {
    src <- mapping$rename_map[[msstats_name]]
    if (src %in% colnames(input)) {
      data.table::setnames(input, src, msstats_name, skip_absent = TRUE)
    }
  }
  
  # ---------------------------------------------------------------
  # 4. Fill missing MSstats columns
  #    FragmentIon and ProductCharge go via columns_to_fill, as in the
  #    converters; filling them here drops them from the output.
  # ---------------------------------------------------------------
  preprocess_fill <- intersect(mapping$fill_columns,
                               c("FragmentIon", "ProductCharge"))
  for (col in setdiff(mapping$fill_columns, preprocess_fill)) {
    input[[col]] <- NA
  }
  
  if (!"Run" %in% colnames(input)) {
    stop("LLM mapping did not produce a 'Run' column.")
  }
  
  # ---------------------------------------------------------------
  # 5. Type coercion. Intensity == 0 -> NA matches the converters; the
  #    previous < 1 also nulled values the converters keep.
  # ---------------------------------------------------------------
  if ("Intensity" %in% colnames(input)) {
    suppressWarnings(
      input[, Intensity := as.numeric(as.character(Intensity))]
    )
    input[, Intensity := ifelse(Intensity == 0, NA_real_, Intensity)]
  }
  if ("PrecursorCharge" %in% colnames(input)) {
    suppressWarnings(
      input[, PrecursorCharge := as.integer(as.character(PrecursorCharge))]
    )
  }
  if ("ProductCharge" %in% colnames(input) && !all(is.na(input$ProductCharge))) {
    suppressWarnings(
      input[, ProductCharge := as.integer(as.character(ProductCharge))]
    )
  }
  if ("Qvalue" %in% colnames(input)) {
    suppressWarnings(
      input[, Qvalue := as.numeric(as.character(Qvalue))]
    )
  }
  
  # ---------------------------------------------------------------
  # 6. Annotation
  # ---------------------------------------------------------------
  if (is.null(annotation)) {
    cond_col <- grep("Condition", colnames(input), value = TRUE,
                     ignore.case = TRUE)[1]
    brep_col <- grep("Replicate", colnames(input), value = TRUE,
                     ignore.case = TRUE)[1]
    if (is.na(cond_col) || is.na(brep_col)) {
      stop("annotation is NULL and no Condition / BioReplicate columns ",
           "were found in input.")
    }
    annotation <- unique(input[, .(Run = as.character(Run),
                                   Condition = get(cond_col),
                                   BioReplicate = get(brep_col))])
    if (verbose) {
      message(sprintf("  Annotation extracted from input: %s, %s",
                      cond_col, brep_col))
    }
    input[, (setdiff(c(cond_col, brep_col), "Run")) := NULL]
  } else {
    if (is.character(annotation) || is.factor(annotation)) {
      annotation <- data.table::fread(annotation)
    }
    annotation <- data.table::as.data.table(annotation)
  }
  
  # Standardizes run labels, stripping '.' and '%'. Every converter calls it.
  annotation <- MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
  
  # ---------------------------------------------------------------
  # 7. Hand off to the shared MSstatsConvert pipeline. Shared peptide
  #    removal, single-feature removal, PSM summarization and the
  #    few-measurements filter were reimplemented by hand here, which is
  #    why removeFewMeasurements was accepted but never used.
  # ---------------------------------------------------------------
  feature_columns <- c("PeptideSequence", "PrecursorCharge",
                       "FragmentIon", "ProductCharge")
  
  # A field the tool does not report is not a feature.
  feature_columns <- feature_columns[
    vapply(feature_columns, function(col) {
      col %in% colnames(input) && !all(is.na(input[[col]]))
    }, logical(1))
  ]
  
  preprocess_feature_columns <- if ("IsotopeLabelType" %in% colnames(input)) {
    c(feature_columns, "IsotopeLabelType")
  } else {
    feature_columns
  }
  columns_to_fill <- stats::setNames(
    rep(list(NA), length(preprocess_fill)), preprocess_fill)
  if (!"IsotopeLabelType" %in% colnames(input)) {
    columns_to_fill[["IsotopeLabelType"]] <- "L"
  }
  
  score_filtering <- filter_cfg$score
  if ("Qvalue" %in% colnames(input)) {
    score_filtering$llm_qvalue <- list(
      score_column    = "Qvalue",
      score_threshold = qvalue_cutoff,
      direction       = "smaller",
      behavior        = "fill",
      handle_na       = "keep",
      fill_value      = NA_real_,
      filter          = filter_with_Qvalue,
      drop_column     = TRUE
    )
  }
  
  input <- MSstatsConvert::MSstatsPreprocess(
    input,
    annotation,
    preprocess_feature_columns,
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = removeProtein_with1Peptide,
    feature_cleaning = list(
      remove_features_with_few_measurements = removeFewMeasurements,
      summarize_multiple_psms = summaryforMultipleRows
    ),
    score_filtering = score_filtering,
    exact_filtering = filter_cfg$exact,
    columns_to_fill = columns_to_fill
  )
  
  input <- MSstatsConvert::MSstatsBalancedDesign(
    input, feature_columns,
    remove_few = removeFewMeasurements
  )
  
  if (verbose) {
    message(sprintf("== LLM converter finished: %d rows, %d proteins, %d runs ==",
                    nrow(input),
                    data.table::uniqueN(input$ProteinName),
                    data.table::uniqueN(input$Run)))
    if (nrow(diagnostics$hallucinated_columns) > 0) {
      message(sprintf("  Hallucinated columns: %d", 
                      nrow(diagnostics$hallucinated_columns)))
    }
    if (nrow(diagnostics$skipped_filters) > 0) {
      message(sprintf("  Skipped filters: %d", 
                      nrow(diagnostics$skipped_filters)))
    }
  }
  
  attr(input, "diagnostics") <- diagnostics
  input
}


# ==================================================================
# Internal helpers
# ==================================================================

#' Parse LLM mapping JSON into a rename map and list of fill columns
#' @param llm_mapping list with $mappings (array of field/from objects)
#' @return list(rename_map = named list, fill_columns = character vector)
#' @keywords internal
.parseLLMMapping <- function(llm_mapping) {
  mappings <- llm_mapping$mappings
  
  if (is.data.frame(mappings)) {
    mappings <- lapply(seq_len(nrow(mappings)), function(i) as.list(mappings[i, ]))
  }
  
  rename_map    <- list()
  fill_columns  <- character(0)
  
  for (m in mappings) {
    field <- m[["field"]]
    from  <- m[["from"]]
    
    if (is.null(field)) next
    
    if (is.null(from) || is.na(from) || from == "" || from == "null") {
      fill_columns <- c(fill_columns, field)
    } else {
      rename_map[[field]] <- from
    }
  }
  
  list(rename_map = rename_map, fill_columns = fill_columns)
}


#' Translate LLM-discovered filters into MSstatsPreprocess filter configs
#'
#' Numeric comparisons become `score_filtering`, value matches become
#' `exact_filtering`, both with `behavior = "fill"`.
#'
#' `equals` on a logical column is inverted to `not_equals`, since
#' `filter_symbols` names the values to remove. On any other column type it
#' cannot be expressed without enumerating the column, so it is skipped.
#'
#' @param input data.table of raw tool output, before renaming
#' @param llm_mapping list with optional $filters array
#' @param verbose logical
#' @return list(score, exact, skipped, applied, keep_cols)
#' @keywords internal
.buildFilterConfigs <- function(input, llm_mapping, verbose = TRUE) {
  skipped <- data.table::data.table(
    column = character(0), operation = character(0),
    value = character(0), reason = character(0))
  applied <- data.table::data.table(
    column = character(0), operation = character(0),
    value = character(0), kind = character(0))
  
  score <- list()
  exact <- list()
  keep_cols <- character(0)
  
  filters <- llm_mapping$filters
  if (is.null(filters)) {
    return(list(score = score, exact = exact, skipped = skipped,
                applied = applied, keep_cols = keep_cols))
  }
  
  if (is.data.frame(filters)) {
    filters <- lapply(seq_len(nrow(filters)), function(i) as.list(filters[i, ]))
  }
  
  .skip <- function(col, op, val, reason) {
    skipped <<- rbind(skipped, data.table::data.table(
      column = as.character(col), operation = as.character(op),
      value = as.character(val), reason = reason))
    if (verbose) message(sprintf("  Filter skipped: %s", reason))
  }
  
  i <- 0L
  for (f in filters) {
    col <- f[["column"]]
    op  <- f[["operation"]]
    val <- as.character(f[["value"]])
    
    if (is.null(col) || !col %in% colnames(input)) {
      .skip(col, op, val,
            if (is.null(col)) "null column name"
            else paste0("column '", col, "' not in data"))
      next
    }
    
    if (is.null(val) || is.na(val) || val == "" || 
        tolower(val) %in% c("na", "null", "none")) {
      .skip(col, op, val, "no valid threshold value")
      next
    }
    
    i <- i + 1L
    key <- paste0("llm_", i)
    col_vals <- input[[col]]
  
    if (op %in% c("less_than", "less_than_or_equals",
                  "greater_than", "greater_than_or_equals")) {
      threshold <- suppressWarnings(as.numeric(val))
      if (is.na(threshold)) {
        .skip(col, op, val, "non-numeric threshold for a numeric comparison")
        next
      }
      score[[key]] <- list(
        score_column    = col,
        score_threshold = threshold,
        direction       = if (grepl("^less", op)) "smaller" else "greater",
        behavior        = "fill",
        handle_na       = "keep",
        fill_value      = NA_real_,
        filter          = TRUE,
        drop_column     = TRUE
      )
      keep_cols <- c(keep_cols, col)
      applied <- rbind(applied, data.table::data.table(
        column = col, operation = as.character(op),
        value = val, kind = "score"))
  
    } else if (op == "not_equals") {
      exact[[key]] <- list(
        col_name       = col,
        filter_symbols = .coerceFilterValue(val, col_vals),
        behavior       = "fill",
        fill_value     = NA_real_,
        filter         = TRUE,
        drop_column    = TRUE
      )
      keep_cols <- c(keep_cols, col)
      applied <- rbind(applied, data.table::data.table(
        column = col, operation = as.character(op),
        value = val, kind = "exact"))
  
    } else if (op == "equals") {
      if (!is.logical(col_vals)) {
        .skip(col, op, val,
              "equals on a non-logical column cannot be expressed as a removal set")
        next
      }
      target <- suppressWarnings(as.logical(val))
      if (is.na(target)) {
        .skip(col, op, val, "value is not a logical")
        next
      }
      exact[[key]] <- list(
        col_name       = col,
        filter_symbols = !target,
        behavior       = "fill",
        fill_value     = NA_real_,
        filter         = TRUE,
        drop_column    = TRUE
      )
      keep_cols <- c(keep_cols, col)
      applied <- rbind(applied, data.table::data.table(
        column = col, operation = "equals (inverted)",
        value = val, kind = "exact"))
  
    } else {
      .skip(col, op, val,
            paste0("operation '", op, "' has no MSstatsPreprocess equivalent"))
      next
    }
  
    if (verbose) {
      message(sprintf("  Filter registered: %s %s %s", col, op, val))
    }
  }
  
  list(score = score, exact = exact, skipped = skipped,
       applied = applied, keep_cols = unique(keep_cols))
}
  
  
#' Coerce a filter value to the column's type
#' @keywords internal
.coerceFilterValue <- function(val, col_vals) {
  if (is.logical(col_vals)) {
    out <- suppressWarnings(as.logical(val))
    if (!is.na(out)) return(out)
  }
  if (is.numeric(col_vals)) {
    out <- suppressWarnings(as.numeric(val))
    if (!is.na(out)) return(out)
  }
  val
}