# ============================================================
# benchmark.R — Run a single benchmark trial
# ============================================================

# ------------------------------------------------------------------
# Create a chat object for a given model key
# ------------------------------------------------------------------

#' Create an `ellmer` chat object for a registered model
#'
#' @param model_key Key into [MODEL_REGISTRY].
#' @return An `ellmer` chat object configured for the model's provider.
#' @details Requires the optional `ellmer` package.
#' @export
create_chat <- function(model_key) {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop("The 'ellmer' package is required for live inference. ",
         "Install it with install.packages('ellmer').", call. = FALSE)
  }

  spec <- MODEL_REGISTRY[[model_key]]
  if (is.null(spec)) {
    stop("Unknown model: ", model_key,
         "\nAvailable: ", paste(names(MODEL_REGISTRY), collapse = ", "))
  }

  if (spec$provider == "ollama") {
    ellmer::chat_ollama(
      model = spec$model,
      base_url = OLLAMA_URL,
      params = ellmer::params(
        format      = "json",
        temperature = 0,
        top_p       = 1,
        seed        = 1,
        num_ctx     = 4096,
        num_predict = -1
      )
    )
  } else if (spec$provider == "anthropic") {
    ellmer::chat_anthropic(
      model = spec$model,
      params = ellmer::params(
        temperature = 0,
        max_tokens  = 4096
      )
    )
  } else {
    stop("Unknown provider: ", spec$provider)
  }
}

# ------------------------------------------------------------------
# Run one benchmark trial
# ------------------------------------------------------------------

#' Run one schema-inference trial
#'
#' Loads a tool's test dataset, builds the prompt, queries a local model, and
#' returns the parsed mapping together with validation and timing metadata.
#' Requires the optional `ellmer` package and a running model backend
#' (e.g. Ollama).
#'
#' @param model_key Key into [MODEL_REGISTRY] (e.g. `"llama3.1-8b"`).
#' @param tool_name Key into [GROUND_TRUTH] / [load_test_dataset()].
#' @param prompt_key One of `names(PROMPT_VERSIONS)` (e.g. `"lean"`,
#'   `"constrained"`, `"constrained_filter"`).
#' @param acquisition `"DIA"`, `"DDA"`, or `NULL` to let the model guess.
#' @param allow_transforms Enable packed-column transform detection.
#' @param n_preview Number of rows to include in the data preview.
#' @param use_structured Use `chat_structured()` (`TRUE`) or raw `chat()`
#'   (`FALSE`).
#' @return A list with `model`, `tool`, `prompt`, `mapping` (parsed),
#'   `validation`, `elapsed_sec`, and `timestamp`.
#' @export
run_trial <- function(model_key,
                      tool_name,
                      prompt_key = "constrained",
                      acquisition = NULL,
                      allow_transforms = FALSE,
                      n_preview  = 3,
                      use_structured = TRUE) {
  
  # Load data and build prompt
  df <- load_test_dataset(tool_name)
  preview <- make_json_preview(df, n_rows = n_preview)
  system_prompt <- PROMPT_VERSIONS[[prompt_key]]
  user_prompt   <- build_user_prompt(preview, acquisition = acquisition,
                                     allow_transforms = allow_transforms)
  
  # Create fresh chat with system prompt
  chat <- create_chat(model_key)
  chat$set_system_prompt(system_prompt)
  
  # Run inference with timing
  t0 <- Sys.time()
  
  # Pick the right structured type based on prompt and transforms
  has_filters <- prompt_key %in% c("filter_aware", "constrained_filter")
  result_type <- .select_mapping_type(has_filters, allow_transforms)
  
  if (use_structured) {
    result <- tryCatch(
      chat$chat_structured(user_prompt, type = result_type, convert = TRUE),
      error = function(e) {
        message("  [WARN] chat_structured failed: ", e$message)
        message("  [INFO] Falling back to raw chat()")
        NULL
      }
    )
    
    # Fallback: raw chat + manual parse
    if (is.null(result)) {
      chat2 <- create_chat(model_key)
      chat2$set_system_prompt(system_prompt)
      raw_text <- chat2$chat(user_prompt)
      result <- tryCatch(
        jsonlite::fromJSON(raw_text, simplifyVector = FALSE),
        error = function(e) {
          message("  [WARN] JSON parse failed: ", e$message)
          list(mappings = list(), raw_response = raw_text)
        }
      )
    }
  } else {
    raw_text <- chat$chat(user_prompt)
    result <- tryCatch(
      jsonlite::fromJSON(raw_text, simplifyVector = FALSE),
      error = function(e) {
        message("  [WARN] JSON parse failed: ", e$message)
        list(mappings = list(), raw_response = raw_text)
      }
    )
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # Post-hoc validation
  validation <- validate_mapping(result)
  
  list(
    model       = model_key,
    tool        = tool_name,
    prompt      = prompt_key,
    mapping     = result,
    validation  = validation,
    elapsed_sec = round(elapsed, 2),
    timestamp   = Sys.time()
  )
}

# ------------------------------------------------------------------
# Post-hoc validation of LLM mapping output
# ------------------------------------------------------------------

#' Post-hoc validation of an LLM mapping result
#'
#' Catches issues the LLM gets wrong regardless of prompt: duplicate column
#' assignments, missing required fields, and out-of-range confidences.
#'
#' @param mapping_result A parsed mapping (with a `$mappings` element).
#' @return A list with `passed` (logical) and `issues` (character vector).
#' @export
validate_mapping <- function(mapping_result) {
  issues <- character(0)
  mappings <- mapping_result$mappings
  
  # Normalize to list-of-lists
  if (is.data.frame(mappings)) {
    mappings <- lapply(seq_len(nrow(mappings)), function(i) as.list(mappings[i, ]))
  }
  
  # --- Check 1: No duplicate column assignments ---
  assigned <- vapply(mappings, function(m) {
    val <- m[["from"]]
    if (is.null(val) || is.na(val)) NA_character_ else val
  }, character(1))
  
  names(assigned) <- vapply(mappings, function(m) m[["field"]], character(1))
  non_na <- assigned[!is.na(assigned)]
  
  dupes <- non_na[duplicated(non_na)]
  if (length(dupes) > 0) {
    for (d in unique(dupes)) {
      fields_using <- names(non_na[non_na == d])
      issues <- c(issues, sprintf(
        "DUPLICATE: Column '%s' assigned to multiple fields: %s",
        d, paste(fields_using, collapse = ", ")
      ))
    }
  }
  
  # --- Check 2: All 8 required fields present ---
  returned_fields <- vapply(mappings, function(m) m[["field"]], character(1))
  missing <- setdiff(MSSTATS_FIELDS, returned_fields)
  if (length(missing) > 0) {
    issues <- c(issues, sprintf(
      "MISSING FIELDS: %s", paste(missing, collapse = ", ")
    ))
  }
  
  # --- Check 3: No invented columns ---
  # (requires column names from the input data — store in mapping_result 
  #  if you want this check. Skipped for now.)
  
  # --- Check 4: Confidence in valid range ---
  confs <- vapply(mappings, function(m) {
    val <- m[["confidence"]]
    if (is.null(val)) NA_real_ else as.numeric(val)
  }, numeric(1))
  
  bad_conf <- which(confs < 0 | confs > 1)
  if (length(bad_conf) > 0) {
    issues <- c(issues, sprintf(
      "INVALID CONFIDENCE: fields %s have confidence outside [0,1]",
      paste(returned_fields[bad_conf], collapse = ", ")
    ))
  }
  
  list(passed = length(issues) == 0, issues = issues)
}

# ------------------------------------------------------------------
# Score a single trial against ground truth
# ------------------------------------------------------------------

#' Score a trial against ground truth
#'
#' @param trial_result A list returned by [run_trial()].
#' @return A `data.table` with one row per MSstats field: `field`, `predicted`,
#'   `expected`, `n_acceptable`, `exact_match`, `acceptable_match`,
#'   `candidate_hit`, `confidence`, `has_candidates`, `n_candidates`.
#' @export
score_trial <- function(trial_result) {
  tool <- trial_result$tool
  truth <- GROUND_TRUTH[[tool]]
  
  if (is.null(truth)) {
    stop("No ground truth for tool: ", tool)
  }
  
  mappings_raw <- trial_result$mapping$mappings
  
  # Normalize to list-of-lists regardless of what chat_structured returned
  if (is.data.frame(mappings_raw)) {
    mappings <- lapply(seq_len(nrow(mappings_raw)), function(i) as.list(mappings_raw[i, ]))
  } else if (is.list(mappings_raw) && !is.null(names(mappings_raw))) {
    # Single mapping returned as a named list instead of list-of-lists
    mappings <- list(mappings_raw)
  } else {
    mappings <- mappings_raw
  }
  
  rows <- lapply(MSSTATS_FIELDS, function(field) {
    # Find this field in the LLM output
    entry <- NULL
    for (m in mappings) {
      if (identical(m$field, field) || identical(m[["field"]], field)) {
        entry <- m
        break
      }
    }
    
    raw_from   <- entry[["from"]]
    if (length(raw_from) > 1) raw_from <- raw_from[1]
    predicted  <- if (is.null(raw_from) || is.na(raw_from) || 
                      raw_from == "" || tolower(raw_from) == "null") {
      NA_character_
    } else {
      raw_from
    }
    # Ground truth may be NULL, a single string, or a character vector
    # (first element = primary, rest = acceptable alternates)
    accepted   <- truth[[field]]
    accepted   <- if (is.null(accepted)) character(0) else accepted
    expected   <- if (length(accepted) == 0) NA_character_ else accepted[1]
    confidence <- if (!is.null(entry[["confidence"]])) entry[["confidence"]] else NA_real_

    norm <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))

    # Exact = matches primary expected name
    exact <- if (is.na(predicted) && is.na(expected)) {
      TRUE
    } else if (is.na(predicted) || is.na(expected)) {
      FALSE
    } else {
      norm(predicted) == norm(expected)
    }

    # Acceptable = matches any element of the ground-truth vector
    acceptable <- if (is.na(predicted) && length(accepted) == 0) {
      TRUE
    } else if (is.na(predicted) || length(accepted) == 0) {
      FALSE
    } else {
      norm(predicted) %in% norm(accepted)
    }
    
    # Candidates may be a data.frame or list-of-lists
    cands <- entry[["candidates"]]
    n_candidates <- if (is.null(cands)) 0L
                    else if (is.data.frame(cands)) nrow(cands)
                    else length(cands)
    has_candidates <- n_candidates > 0L
    
    candidate_hit <- FALSE
    if (!acceptable && has_candidates && length(accepted) > 0) {
      if (is.data.frame(cands)) {
        cand_names <- cands$from
      } else {
        cand_names <- sapply(cands, function(c) c[["from"]])
      }
      candidate_hit <- any(norm(cand_names) %in% norm(accepted))
    }
    
    data.table(
      field            = field,
      predicted        = predicted,
      expected         = expected,
      n_acceptable     = length(accepted),
      exact_match      = exact,
      acceptable_match = acceptable,
      candidate_hit    = candidate_hit,
      confidence       = confidence,
      has_candidates   = has_candidates,
      n_candidates     = n_candidates
    )
  })
  
  rbindlist(rows)
}

# ------------------------------------------------------------------
# Pretty-print a scored trial
# ------------------------------------------------------------------

#' Pretty-print a scored trial to the console
#'
#' @param scores A `data.table` returned by [score_trial()].
#' @param trial_result The list returned by [run_trial()].
#' @return `scores`, invisibly.
#' @export
print_scorecard <- function(scores, trial_result) {
  cat("\n", strrep("=", 60), "\n")
  cat(sprintf("Model: %-25s  Prompt: %s\n",
              trial_result$model, trial_result$prompt))
  cat(sprintf("Tool:  %-25s  Time:   %.1fs\n",
              trial_result$tool, trial_result$elapsed_sec))
  cat(strrep("-", 60), "\n")
  
  for (i in seq_len(nrow(scores))) {
    r <- scores[i]
    status <- if (r$exact_match) "\u2705"
              else if (r$acceptable_match) "\U0001f7e2"
              else if (r$candidate_hit) "\U0001f7e1"
              else "\u274c"
    expected_str <- if (is.na(r$expected)) {
      "(null)"
    } else if (r$n_acceptable > 1) {
      sprintf("%s (+%d alt)", r$expected, r$n_acceptable - 1L)
    } else {
      r$expected
    }
    cat(sprintf("  %s %-17s  predicted=%-25s expected=%s\n",
                status, r$field,
                ifelse(is.na(r$predicted), "(null)", r$predicted),
                expected_str))
  }
  
  cat(strrep("-", 60), "\n")
  cat(sprintf("  Strict:     %d/%d (%.0f%%)   primary-only matches\n",
              sum(scores$exact_match), nrow(scores),
              mean(scores$exact_match) * 100))
  cat(sprintf("  Acceptable: %d/%d (%.0f%%)   primary or alternate\n",
              sum(scores$acceptable_match), nrow(scores),
              mean(scores$acceptable_match) * 100))
  cat(sprintf("  Candidate recovery: %d\n", sum(scores$candidate_hit)))
  
  # Show validation issues
  v <- trial_result$validation
  if (!is.null(v) && !v$passed) {
    cat(sprintf("  \u26a0\ufe0f  Validation issues (%d):\n", length(v$issues)))
    for (issue in v$issues) {
      cat("     ", issue, "\n")
    }
  } else if (!is.null(v) && v$passed) {
    cat("  \u2705 All post-hoc checks passed\n")
  }
  
  # Show filter results if this was a filter-aware prompt
  if (trial_result$prompt %in% c("filter_aware", "constrained_filter")) {
    predicted_filters <- trial_result$mapping$filters
    if (is.data.frame(predicted_filters)) {
      predicted_filters <- lapply(
        seq_len(nrow(predicted_filters)),
        function(i) as.list(predicted_filters[i, ])
      )
    }
    if (is.null(predicted_filters)) predicted_filters <- list()
    
    cat(strrep("-", 60), "\n")
    cat(sprintf("  Discovered filters (%d):\n", length(predicted_filters)))
    if (length(predicted_filters) == 0) {
      cat("    (none)\n")
    } else {
      for (f in predicted_filters) {
        cat(sprintf("    %s %s %s  (conf=%.2f)\n",
                    f[["column"]], f[["operation"]], f[["value"]],
                    as.numeric(f[["confidence"]])))
        cat(sprintf("      %s\n", f[["description"]]))
      }
    }
  }
  cat(strrep("=", 60), "\n\n")
  
  invisible(scores)
}