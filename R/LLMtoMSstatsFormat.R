# ============================================================
# LLMtoMSstatsFormat.R — Generic converter driven by LLM schema mapping
#
# Uses the LLM's mapping + filter output to transform arbitrary
# proteomics tool output into MSstats format, leveraging the
# existing MSstatsConvert pipeline for all downstream processing.
# ============================================================

#' Convert arbitrary proteomics output to MSstats format using LLM mapping
#'
#' @param input data.table of raw tool output
#' @param annotation data.table with Run, Condition, BioReplicate columns
#' @param llm_mapping list — parsed LLM output with $mappings and optionally $filters
#' @param useUniquePeptide logical, remove shared peptides (default TRUE)
#' @param removeFewMeasurements logical (default TRUE)
#' @param removeProtein_with1Peptide logical (default FALSE)
#' @param summaryforMultipleRows function for summarizing multiple PSMs (default max)
#' @param qvalue_cutoff numeric, keep only rows with Qvalue < cutoff (default 0.05).
#'   Set to NULL to disable. Rows with NA Qvalue are retained.
#' @param use_log_file logical (default TRUE)
#' @param append logical (default FALSE)
#' @param verbose logical (default TRUE)
#' @param log_file_path character or NULL
#'
#' @return data.table in MSstats format, ready for dataProcess()
#' @export
LLMtoMSstatsFormat <- function(
    input,
    annotation,
    llm_mapping,
    useUniquePeptide = TRUE,
    removeFewMeasurements = TRUE,
    removeProtein_with1Peptide = FALSE,
    summaryforMultipleRows = max,
    qvalue_cutoff = 0.05,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path)
  
  input <- data.table::copy(input)
  
  # Track issues for auditing / poster
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
      rows_removed = integer(0)
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
  # 2. Apply LLM-discovered filters BEFORE column renaming
  #    (filter column names are in the original tool namespace)
  # ---------------------------------------------------------------
  n_before <- nrow(input)
  filter_result <- .applyLLMFilters(input, llm_mapping, verbose)
  input <- filter_result$data
  diagnostics$skipped_filters <- filter_result$skipped
  diagnostics$applied_filters <- filter_result$applied
  n_after <- nrow(input)
  
  if (verbose) {
    message(sprintf("== LLM Filters: %d -> %d rows (removed %d) ==",
                    n_before, n_after, n_before - n_after))
  }
  
  # ---------------------------------------------------------------
  # 3. Select and rename columns to MSstats names
  # ---------------------------------------------------------------
  # Keep only the columns we need — demote hallucinated columns to NA fill
  source_cols <- unlist(mapping$rename_map)
  missing_in_data <- setdiff(source_cols, colnames(input))
  if (length(missing_in_data) > 0) {
    if (verbose) {
      message(sprintf("  WARNING: LLM hallucinated columns not in data: %s",
                      paste(missing_in_data, collapse = ", ")))
      message("  -> These fields will be filled with NA instead.")
    }
    # Record and demote hallucinated mappings to fill list
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
  
  # Subset to mapped columns (plus any extras needed for annotation merge)
  keep_cols <- unique(c(source_cols, 
                        intersect(c("Condition", "BioReplicate"), colnames(input))))
  input <- input[, keep_cols, with = FALSE]
  
  # Rename: source name -> MSstats name
  for (msstats_name in names(mapping$rename_map)) {
    src <- mapping$rename_map[[msstats_name]]
    if (src %in% colnames(input)) {
      data.table::setnames(input, src, msstats_name, skip_absent = TRUE)
    }
  }
  
  # ---------------------------------------------------------------
  # 4. Fill missing MSstats columns with NA
  # ---------------------------------------------------------------
  for (col in mapping$fill_columns) {
    input[[col]] <- NA
  }
  
  # IsotopeLabelType is always required
  if (!"IsotopeLabelType" %in% colnames(input)) {
    input[["IsotopeLabelType"]] <- "L"
  }
  
  # ---------------------------------------------------------------
  # 5. Type coercion for MSstats expectations
  # ---------------------------------------------------------------
  if ("Intensity" %in% colnames(input)) {
    suppressWarnings(
      input[, Intensity := as.numeric(as.character(Intensity))]
    )
    n_below <- sum(!is.na(input$Intensity) & input$Intensity < 1)
    if (n_below > 0) {
      input[Intensity < 1, Intensity := NA_real_]
      if (verbose) {
        message(sprintf("  Set %d Intensity values < 1 to NA", n_below))
      }
    }
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
    if (!is.null(qvalue_cutoff)) {
      n_before <- nrow(input)
      input <- input[is.na(Qvalue) | Qvalue < qvalue_cutoff]
      if (verbose) {
        message(sprintf("  Qvalue filter (< %s): removed %d rows",
                        qvalue_cutoff, n_before - nrow(input)))
      }
    }
  } else if (!is.null(qvalue_cutoff) && verbose) {
    message("  Qvalue filter skipped: no Qvalue column mapped")
  }
  
  # ---------------------------------------------------------------
  # 6. Merge annotation
  # ---------------------------------------------------------------
  if (!"Run" %in% colnames(input)) {
    stop("LLM mapping did not produce a 'Run' column.")
  }
  
  if (is.character(annotation) || is.factor(annotation)) {
    annotation <- data.table::fread(annotation)
  }
  annotation <- data.table::as.data.table(annotation)
  
  # Find the annotation column whose values match the Run values in input.
  # This avoids depending on annotation column naming conventions.
  input_runs <- unique(as.character(input$Run))
  annot_run_col <- NULL
  
  for (acol in colnames(annotation)) {
    annot_vals <- unique(as.character(annotation[[acol]]))
    overlap <- length(intersect(input_runs, annot_vals))
    if (overlap > 0 && overlap >= length(input_runs) * 0.5) {
      annot_run_col <- acol
      break
    }
  }
  
  if (is.null(annot_run_col)) {
    warning("Could not find annotation column matching Run values. ",
            "Condition and BioReplicate will be NA.\n",
            "  Input runs: ", paste(head(input_runs, 3), collapse = ", "), "\n",
            "  Annotation cols: ", paste(colnames(annotation), collapse = ", "))
  } else {
    if (verbose) {
      message(sprintf("  Annotation merge: input$Run <-> annotation$%s", annot_run_col))
    }
    # Rename annotation column to Run for merge
    if (annot_run_col != "Run") {
      data.table::setnames(annotation, annot_run_col, "Run")
    }
    input <- merge(input, annotation, by = "Run", all.x = TRUE,
                   suffixes = c("", ".annot"))
  }
  
  # ---------------------------------------------------------------
  # 7. Determine feature columns based on what's available
  # ---------------------------------------------------------------
  possible_features <- c("PeptideSequence", "PrecursorCharge", 
                         "FragmentIon", "ProductCharge")
  feature_columns <- intersect(possible_features, colnames(input))
  # Only use feature columns that aren't entirely NA
  feature_columns <- feature_columns[
    sapply(feature_columns, function(col) !all(is.na(input[[col]])))
  ]
  
  # At minimum, PeptideSequence must be a feature column
  if (!"PeptideSequence" %in% feature_columns) {
    stop("PeptideSequence column is missing or entirely NA after mapping.")
  }
  
  # ---------------------------------------------------------------
  # 8. Standard MSstats cleanup
  # ---------------------------------------------------------------
  # Remove shared peptides if requested
  if (useUniquePeptide && "ProteinName" %in% colnames(input)) {
    # Simple shared peptide removal: drop peptides mapping to multiple proteins
    pep_prot <- unique(input[, .(PeptideSequence, ProteinName)])
    shared <- pep_prot[, .N, by = PeptideSequence][N > 1]$PeptideSequence
    if (length(shared) > 0) {
      input <- input[!PeptideSequence %in% shared]
      if (verbose) {
        message(sprintf("  Removed %d shared peptides", length(shared)))
      }
    }
  }
  
  # Remove proteins with single feature
  if (removeProtein_with1Peptide && "ProteinName" %in% colnames(input)) {
    feat_count <- unique(input[, c("ProteinName", feature_columns), with = FALSE])
    feat_count <- feat_count[, .N, by = ProteinName]
    single_feat <- feat_count[N <= 1]$ProteinName
    if (length(single_feat) > 0) {
      input <- input[!ProteinName %in% single_feat]
      if (verbose) {
        message(sprintf("  Removed %d single-feature proteins", length(single_feat)))
      }
    }
  }
  
  # Summarize multiple PSMs per feature (e.g., take max intensity)
  key_cols <- c("ProteinName", feature_columns, "Run",
                "Condition", "BioReplicate", "IsotopeLabelType")
  key_cols <- intersect(key_cols, colnames(input))
  if (anyDuplicated(input, by = key_cols)) {
    input <- input[, .(Intensity = summaryforMultipleRows(Intensity, na.rm = TRUE)),
                   by = key_cols]
    if (verbose) message("  Summarized multiple PSMs per feature")
  }
  
  # ---------------------------------------------------------------
  # 9. Ensure required MSstats columns are present
  # ---------------------------------------------------------------
  required_final <- c("ProteinName", "PeptideSequence", "PrecursorCharge",
                      "FragmentIon", "ProductCharge", "IsotopeLabelType",
                      "Condition", "BioReplicate", "Run", "Intensity")
  
  for (col in required_final) {
    if (!col %in% colnames(input)) {
      input[[col]] <- NA
    }
  }
  
  # Select final columns in MSstats order
  input <- input[, required_final, with = FALSE]
  
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
  
  # Normalize data.frame to list-of-lists
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


#' Apply LLM-discovered filters to raw data
#' @param input data.table
#' @param llm_mapping list with optional $filters array
#' @param verbose logical
#' @return filtered data.table
#' @keywords internal
.applyLLMFilters <- function(input, llm_mapping, verbose = TRUE) {
  skipped <- data.table(column = character(0), operation = character(0),
                        value = character(0), reason = character(0))
  applied <- data.table(column = character(0), operation = character(0),
                        value = character(0), rows_removed = integer(0))
  
  # browser()
  filters <- llm_mapping$filters
  if (is.null(filters)) return(list(data = input, skipped = skipped, applied = applied))
  
  # Normalize
  if (is.data.frame(filters)) {
    filters <- lapply(seq_len(nrow(filters)), function(i) as.list(filters[i, ]))
  }
  
  for (f in filters) {
    col <- f[["column"]]
    op  <- f[["operation"]]
    val <- as.character(f[["value"]])
    
    if (is.null(col) || !col %in% colnames(input)) {
      reason <- if (is.null(col)) "null column name" else paste0("column '", col, "' not in data")
      skipped <- rbind(skipped, data.table(
        column = as.character(col), operation = as.character(op),
        value = val, reason = reason))
      if (verbose) message(sprintf("  Filter skipped: %s", reason))
      next
    }
    
    # Skip filters where the LLM couldn't determine a threshold
    if (is.null(val) || is.na(val) || val == "" || 
        tolower(val) %in% c("na", "null", "none")) {
      skipped <- rbind(skipped, data.table(
        column = col, operation = as.character(op),
        value = val, reason = "no valid threshold value"))
      if (verbose) {
        message(sprintf("  Filter skipped: '%s %s' has no valid threshold value",
                        col, op))
      }
      next
    }
    
    n_before <- nrow(input)
    
    input <- tryCatch({
      col_vals <- input[[col]]
      
      # Coerce value to match column type
      if (is.numeric(col_vals)) {
        val <- as.numeric(val)
      }
      
      keep <- switch(op,
        "equals"               = col_vals == val,
        "not_equals"           = col_vals != val,
        "less_than"            = as.numeric(col_vals) < as.numeric(val),
        "greater_than"         = as.numeric(col_vals) > as.numeric(val),
        "less_than_or_equals"  = as.numeric(col_vals) <= as.numeric(val),
        "greater_than_or_equals" = as.numeric(col_vals) >= as.numeric(val),
        "contains"             = grepl(val, col_vals, fixed = TRUE),
        "not_contains"         = !grepl(val, col_vals, fixed = TRUE),
        {
          skipped <<- rbind(skipped, data.table(
            column = col, operation = as.character(op),
            value = as.character(val), reason = paste0("unknown operation '", op, "'")))
          if (verbose) message(sprintf("  Filter skipped: unknown operation '%s'", op))
          rep(TRUE, nrow(input))
        }
      )
      
      # Handle NAs in the keep vector — don't drop rows with NA filter values
      keep[is.na(keep)] <- TRUE
      input[keep]
    }, error = function(e) {
      skipped <<- rbind(skipped, data.table(
        column = col, operation = as.character(op),
        value = as.character(val), reason = paste0("error: ", e$message)))
      if (verbose) message(sprintf("  Filter failed on '%s': %s", col, e$message))
      input
    })
    
    rows_removed <- n_before - nrow(input)
    applied <- rbind(applied, data.table(
      column = col, operation = as.character(op),
      value = as.character(val), rows_removed = as.integer(rows_removed)))
    
    if (verbose) {
      message(sprintf("  Filter: %s %s %s -> removed %d rows",
                      col, op, val, rows_removed))
    }
  }
  
  list(data = input, skipped = skipped, applied = applied)
}