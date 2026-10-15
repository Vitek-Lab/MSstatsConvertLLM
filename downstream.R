# ============================================================
# downstream.R: shared engine for the downstream comparisons
#
# Sourced by compare_solivais.R and compare_pxd005642.R. Each dataset script
# builds a configuration and calls run_downstream(), which sends one raw file
# through the hand-written converter and the LLM-driven converter, runs both
# through dataProcess and groupComparison with identical settings, and
# compares the results on a controlled mixture with known truth.
#
# Environment variables:
#   LLM_MODEL, LLM_PROMPT  the saved mapping to test. LLM_MODEL=reference uses
#                          the converter's own default columns instead
#   FORCE=1                recompute cached steps, needed after converter edits
#   N_CORES                cores for dataProcess and groupComparison
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MSstatsConvert)
  library(MSstats)
})

# Load the working tree, so an older installed copy of the package cannot
# supply the pre-parity converter.
if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(MSstatsConvertLLM)
}
if (!"filter_with_Qvalue" %in% names(formals(LLMtoMSstatsFormat))) {
  stop("LLMtoMSstatsFormat predates the parity fix. Check out a branch that includes it.")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
say <- function(template, ...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(template, ...), "\n", sep = "")
fmt <- function(n) format(n, big.mark = ",", scientific = FALSE)

TRIALS_RDS <- "results/raw_trials_benchmark_20261007_163702.rds"
FORCE      <- nzchar(Sys.getenv("FORCE"))
n_detected <- parallel::detectCores()
if (is.na(n_detected)) n_detected <- 2L
N_CORES <- as.integer(Sys.getenv("N_CORES", as.character(max(1L, min(8L, n_detected - 1L)))))

COL_HAND  <- "#2C5F8A"
COL_LLM   <- "#8A5F2C"
COL_POS   <- "#2C8A5F"
COL_NEG   <- "grey65"
COL_OTHER <- "#C9A227"

cached <- function(name, dir, expr) {
  path <- file.path(dir, paste0(name, ".rds"))
  if (!FORCE && file.exists(path)) {
    say("  %s: from cache", name)
    return(readRDS(path))
  }
  t0 <- Sys.time()
  value <- expr
  saveRDS(value, path)
  say("  %s: %.1f min", name, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  value
}

# ------------------------------------------------------------------
# Mappings
# ------------------------------------------------------------------
mapping_vector <- function(m) {
  mp <- m$mappings
  if (is.data.frame(mp)) return(setNames(as.character(mp$from), as.character(mp$field)))
  setNames(vapply(mp, function(x) as.character(x[["from"]] %||% NA_character_)[1], ""),
           vapply(mp, function(x) as.character(x[["field"]])[1], ""))
}

saved_mapping <- function(model, tool, prompt) {
  if (!file.exists(TRIALS_RDS)) stop("Missing ", TRIALS_RDS)
  trials <- readRDS(TRIALS_RDS)
  picked <- Filter(function(t) identical(t$model, model) && identical(t$tool, tool) &&
                     identical(t$prompt, prompt), trials)
  if (!length(picked)) stop("No saved trials for ", model, " / ", tool, " / ", prompt)
  vecs  <- lapply(picked, function(t) { v <- mapping_vector(t$mapping); v[order(names(v))] })
  agree <- all(vapply(vecs, identical, logical(1), vecs[[1]]))
  say("  %d saved repeats, mapping identical across repeats: %s", length(picked), agree)
  if (!agree) warning("Repeats disagree on the mapping; using the first.")
  picked[[1]]$mapping
}

# The converter's own default columns, in the structure of a saved mapping.
# Testing it separates what the generic pipeline cannot express from anything
# a model chose.
reference_mapping <- function(template, columns, filters = NULL) {
  m <- template
  if (is.data.frame(m$mappings)) {
    m$mappings$from <- unname(columns[as.character(m$mappings$field)])
  } else {
    m$mappings <- lapply(m$mappings, function(x) { x$from <- unname(columns[x$field]); x })
  }
  # The converter's own built-in filters, where it has any; otherwise none
  if (!is.null(filters)) {
    m$filters <- filters
  } else if (is.data.frame(m$filters)) {
    m$filters <- m$filters[0, , drop = FALSE]
  } else if (!is.null(m$filters)) {
    m$filters <- list()
  }
  m
}

print_mapping <- function(m) {
  mv <- mapping_vector(m)
  for (f in names(mv)) say("    %-16s <- %s", f, mv[[f]])
  fl <- m$filters
  if (is.data.frame(fl) && nrow(fl)) {
    for (i in seq_len(nrow(fl))) {
      say("    proposed filter: %s %s %s", as.character(fl$column[i]),
          as.character(fl$operation[i]), as.character(unlist(fl$value[i])))
    }
  } else if (is.list(fl) && length(fl)) {
    for (f in fl) {
      say("    proposed filter: %s %s %s", as.character(f$column %||% "?"),
          as.character(f$operation %||% "?"), as.character(f$value %||% "?"))
    }
  } else {
    say("    no filters proposed")
  }
  invisible(mv)
}

# Reads LLM_MODEL and LLM_PROMPT. The reference option borrows the structure of
# a filter-free saved mapping and replaces every column with the converter's.
# REF_SWAP=Field=Column changes exactly one reference column, so a single
# choice can be tested with everything else held at the converter's defaults.
choose_mapping <- function(tool, reference_columns, reference_filters = NULL,
                           template = c(model = "llama3.1-8b", prompt = "constrained_v3")) {
  model  <- Sys.getenv("LLM_MODEL",  "deepseek-r1-14b")
  prompt <- Sys.getenv("LLM_PROMPT", "filter_aware")
  if (identical(model, "reference")) {
    cols  <- reference_columns
    label <- "reference"
    description <- "reference mapping (the converter's default columns)"
    swap <- Sys.getenv("REF_SWAP", "")
    if (nzchar(swap)) {
      kv <- strsplit(swap, "=", fixed = TRUE)[[1]]
      if (length(kv) != 2 || !kv[1] %in% names(cols)) {
        stop("REF_SWAP must look like Field=Column, with Field one of: ", paste(names(cols), collapse = ", "))
      }
      cols[[kv[1]]] <- kv[2]
      label <- sprintf("reference_%s_%s", kv[1], kv[2])
      description <- sprintf("reference mapping with %s from %s", kv[1], kv[2])
    }
    say("1. Mapping: %s", description)
    m <- reference_mapping(saved_mapping(template[["model"]], tool, template[["prompt"]]),
                           cols, reference_filters)
  } else {
    say("1. Mapping: %s, %s prompt, from %s", model, prompt, TRIALS_RDS)
    m <- saved_mapping(model, tool, prompt)
    label <- paste(model, prompt, sep = "_")
    description <- sprintf("%s, %s prompt", model, prompt)
  }
  print_mapping(m)
  list(mapping = m, label = label, description = description)
}

print_diagnostics <- function(dg) {
  if (is.null(dg)) return(invisible())
  if (NROW(dg$hallucinated_columns)) {
    say("  hallucinated columns: %s", paste(dg$hallucinated_columns$hallucinated_name, collapse = ", "))
  }
  for (i in seq_len(NROW(dg$applied_filters))) {
    a <- dg$applied_filters[i]
    say("  filter applied: %s %s %s (%s)", as.character(a$column), as.character(a$operation),
        as.character(unlist(a$value)), as.character(a$kind %||% ""))
  }
  for (i in seq_len(NROW(dg$skipped_filters))) {
    s <- dg$skipped_filters[i]
    say("  filter skipped: %s %s %s (%s)", as.character(s$column), as.character(s$operation),
        as.character(unlist(s$value)), as.character(s$reason %||% ""))
  }
  invisible()
}

# ------------------------------------------------------------------
# Comparisons
# ------------------------------------------------------------------
# sequence_key rewrites PeptideSequence before rows are matched, for converters
# that encode it differently from a single source column.
compare_converted <- function(a, b, sequence_key = identity, protein_key = identity) {
  prep <- function(x) {
    x <- copy(x)
    for (j in names(x)) if (is.factor(x[[j]])) set(x, j = j, value = as.character(x[[j]]))
    if ("ProteinName" %in% names(x)) set(x, j = "ProteinName", value = protein_key(x$ProteinName))
    # Some converters write the peptide column as PeptideModifiedSequence
    if (!"PeptideSequence" %in% names(x) && "PeptideModifiedSequence" %in% names(x)) {
      setnames(x, "PeptideModifiedSequence", "PeptideSequence")
    }
    if ("PeptideSequence" %in% names(x)) set(x, j = "PeptideSequence", value = sequence_key(x$PeptideSequence))
    x
  }
  a <- prep(a); b <- prep(b)
  cols <- intersect(names(a), names(b))
  a <- a[, ..cols]; b <- b[, ..cols]
  key <- intersect(c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                     "ProductCharge", "IsotopeLabelType", "Run"), cols)
  key <- key[!vapply(key, function(k) all(is.na(a[[k]])) && all(is.na(b[[k]])), logical(1))]
  feat <- setdiff(key, "Run")
  n_feat <- c(hand = uniqueN(a, by = feat), llm = uniqueN(b, by = feat))
  setorderv(a, cols, na.last = TRUE); setorderv(b, cols, na.last = TRUE)
  if (isTRUE(all.equal(a, b, check.attributes = FALSE))) {
    return(list(identical = TRUE, rows_hand = nrow(a), rows_llm = nrow(b),
                features_hand = n_feat[["hand"]], features_llm = n_feat[["llm"]]))
  }
  # Repeated keys would make the joins below explode. Report them with an
  # example, then match on the first row of each key.
  dup_hand <- sum(duplicated(a, by = key)); dup_llm <- sum(duplicated(b, by = key))
  for (side in c("hand", "llm")) {
    x <- if (side == "hand") a else b
    if (sum(duplicated(x, by = key))) {
      say("    %s path: %s rows repeat a key (%s); example:", side,
          fmt(sum(duplicated(x, by = key))), paste(key, collapse = ", "))
      ex <- x[duplicated(x, by = key) | duplicated(x, by = key, fromLast = TRUE)]
      setorderv(ex, key)
      print(head(ex, 4))
    }
  }
  if (dup_hand) a <- unique(a, by = key)
  if (dup_llm)  b <- unique(b, by = key)

  # Columns other than the key and Intensity that differ, if any
  shared_rows <- merge(a, b, by = key, suffixes = c("_hand", "_llm"))
  other <- setdiff(cols, c(key, "Intensity"))
  other_diff <- other[vapply(other, function(j) {
    x <- shared_rows[[paste0(j, "_hand")]]; y <- shared_rows[[paste0(j, "_llm")]]
    !isTRUE(all.equal(x, y, check.attributes = FALSE))
  }, logical(1))]
  ka <- a[, c(key, "Intensity"), with = FALSE][, in_hand := TRUE]
  kb <- b[, c(key, "Intensity"), with = FALSE][, in_llm := TRUE]
  m  <- merge(ka, kb, by = key, all = TRUE, suffixes = c("_hand", "_llm"))
  both <- m[in_hand %in% TRUE & in_llm %in% TRUE]
  ih <- both$Intensity_hand; il <- both$Intensity_llm
  same <- (is.na(ih) & is.na(il)) |
          (!is.na(ih) & !is.na(il) & abs(ih - il) <= 1e-8 * pmax(1, abs(ih)))
  list(identical = FALSE, rows_hand = nrow(a), rows_llm = nrow(b),
       features_hand = n_feat[["hand"]], features_llm = n_feat[["llm"]],
       only_hand = m[is.na(in_llm), .N], only_llm = m[is.na(in_hand), .N],
       intensity_differs = sum(!same),
       blank_in_llm_only = sum(!is.na(ih) & is.na(il)),
       blank_in_hand_only = sum(is.na(ih) & !is.na(il)),
       other_columns_differ = other_diff,
       duplicate_keys_hand = dup_hand, duplicate_keys_llm = dup_llm)
}

levels_of <- function(p) {
  g <- p$ProteinLevelData$GROUP
  if (is.factor(g)) levels(g) else sort(unique(as.character(g)))
}

thin <- function(d, n = 50000) if (nrow(d) > n) { set.seed(1); d[sample(.N, n)] } else d

venn_plot <- function(n_hand, n_both, n_llm, title) {
  a <- n_hand + n_both; b <- n_llm + n_both
  if (a + b == 0) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "No significant proteins") +
             theme_void() + labs(title = title))
  }
  u  <- max(a, b)
  r1 <- sqrt(a / u); r2 <- sqrt(b / u)
  lens <- function(d) {
    if (d >= r1 + r2) return(0)
    if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
    x <- max(-1, min(1, (d^2 + r1^2 - r2^2) / (2 * d * r1)))
    y <- max(-1, min(1, (d^2 + r2^2 - r1^2) / (2 * d * r2)))
    r1^2 * acos(x) + r2^2 * acos(y) -
      0.5 * sqrt(max(0, (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)))
  }
  d <- if (n_both == 0) {
    r1 + r2 + 0.1
  } else if (n_both >= min(a, b)) {
    abs(r1 - r2)
  } else {
    tryCatch(uniroot(function(d) lens(d) - pi * n_both / u,
                     c(abs(r1 - r2) + 1e-6, r1 + r2 - 1e-6))$root,
             error = function(e) (r1 + r2) * 0.6)
  }
  th <- seq(0, 2 * pi, length.out = 200)
  circ <- rbind(data.table(set = "Hand-written", x = -d / 2 + r1 * cos(th), y = r1 * sin(th)),
                data.table(set = "LLM",          x =  d / 2 + r2 * cos(th), y = r2 * sin(th)))
  ggplot(circ, aes(x, y, fill = set, group = set)) +
    geom_polygon(alpha = 0.45, colour = "grey20", linewidth = 0.4) +
    annotate("text", x = -d / 2 - r1 * 0.5, y = 0, label = fmt(n_hand), fontface = "bold") +
    annotate("text", x = 0, y = 0, label = fmt(n_both), fontface = "bold") +
    annotate("text", x = d / 2 + r2 * 0.5, y = 0, label = fmt(n_llm), fontface = "bold") +
    scale_fill_manual(values = c("Hand-written" = COL_HAND, "LLM" = COL_LLM)) +
    coord_fixed(clip = "off") + labs(title = title) + theme_void(base_size = 11) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

pkg_source <- function(p) {
  d <- packageDescription(p)
  src <- d$RemoteType %||% (if (!is.null(d$git_branch)) paste("Bioconductor", d$git_branch)) %||%
    d$Repository %||% "source unknown"
  sha <- d$RemoteSha %||% d$git_last_commit
  sprintf("%s %s, %s%s", p, d$Version, src, if (!is.null(sha)) paste0(" ", substr(sha, 1, 7)) else "")
}

# ------------------------------------------------------------------
# The pipeline
#
# cfg fields: dataset, slug, label, description, mapping, raw (original
# headers), annotation, run_col, hand_convert(raw, annotation),
# species_of(protein) returning positive / negative / "OTHER" (NA stops the
# run), positive, negative, dp_args, build_contrast(levels),
# designed_ratio(label) or NULL, alpha.
# ------------------------------------------------------------------
run_downstream <- function(cfg) {
  pos <- cfg$positive; neg <- cfg$negative; alpha <- cfg$alpha %||% 0.05
  # Maps every protein identifier a mapping can produce to one canonical key,
  # so two identifier vocabularies for the same proteins compare as equal
  protein_key <- cfg$protein_key %||% identity
  run_name    <- sprintf("downstream_%s_%s", cfg$slug, cfg$label)
  out_dir     <- file.path("results", run_name)                                  # not tracked
  llm_cache   <- file.path(out_dir, "cache")
  hand_cache  <- file.path("results", sprintf("downstream_%s_hand", cfg$slug), "cache")
  summary_csv <- file.path("results", paste0(run_name, "_summary.csv"))           # tracked
  for (d in c(llm_cache, hand_cache)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  raw <- cfg$raw; annotation <- cfg$annotation

  runs  <- unique(raw[[cfg$run_col]])
  miss  <- setdiff(annotation$Run, runs)
  extra <- setdiff(runs, annotation$Run)
  if (length(miss) || length(extra)) {
    stop("Annotation and ", cfg$run_col, " do not match.",
         "\n  only in annotation: ", paste(head(miss, 3), collapse = ", "),
         "\n  only in data: ", paste(head(extra, 3), collapse = ", "))
  }
  say("  %d runs, all matched to the annotation", length(runs))

  # 3. Convert with both paths -------------------------------------------
  say("3. Converting")
  hand_converted <- cached("hand_converted", hand_cache,
                           as.data.table(cfg$hand_convert(copy(raw), copy(annotation))))
  llm_out <- cached("llm_converted", llm_cache, {
    out <- LLMtoMSstatsFormat(input = normalize_headers(copy(raw)), annotation = copy(annotation),
                              llm_mapping = cfg$mapping, use_log_file = FALSE)
    list(data = as.data.table(out), diagnostics = attr(out, "diagnostics"))
  })
  llm_converted <- llm_out$data
  print_diagnostics(llm_out$diagnostics)

  # Fail before the slow step if the truth labels did not reach ProteinName
  # The raw table is no longer needed; free it before the memory-heavy steps
  rm(raw); cfg$raw <- NULL; invisible(gc())

  for (which_path in c("hand", "llm")) {
    x  <- if (which_path == "hand") hand_converted else llm_converted
    sp <- cfg$species_of(protein_key(unique(as.character(x$ProteinName))))
    if (anyNA(sp)) {
      stop(which_path, " path: ", sum(is.na(sp)), " proteins have no truth label; check the ProteinName mapping.")
    }
    if (!any(sp == pos) || !any(sp == neg)) {
      stop(which_path, " path: no ", if (!any(sp == pos)) pos else neg,
           " proteins found; check the ProteinName mapping.")
    }
  }

  cc <- compare_converted(hand_converted, llm_converted, cfg$sequence_key %||% identity, protein_key)
  if (cc$identical) {
    say("  converted data identical: %s rows in each path", fmt(cc$rows_hand))
  } else {
    say("  converted data differ: hand %s rows, LLM %s rows", fmt(cc$rows_hand), fmt(cc$rows_llm))
    say("    rows only in hand %s, only in LLM %s", fmt(cc$only_hand), fmt(cc$only_llm))
    say("    shared rows with a different intensity %s (blank in LLM only %s, blank in hand only %s)",
        fmt(cc$intensity_differs), fmt(cc$blank_in_llm_only), fmt(cc$blank_in_hand_only))
    say("    other columns that differ on shared rows: %s",
        if (length(cc$other_columns_differ)) paste(cc$other_columns_differ, collapse = ", ") else "none")
    if (cc$duplicate_keys_hand || cc$duplicate_keys_llm) {
      say("    repeated keys collapsed before matching: hand %s, LLM %s",
          fmt(cc$duplicate_keys_hand), fmt(cc$duplicate_keys_llm))
    }
  }
  say("  features: hand %s, LLM %s", fmt(cc$features_hand), fmt(cc$features_llm))

  # 4. dataProcess, same settings on both paths --------------------------
  say("4. dataProcess; this is the slow step")
  # The data is passed by name, so a failure never deparses the whole table
  run_dp <- function(x) {
    do.call("dataProcess", c(list(quote(x)), cfg$dp_args,
                             list(use_log_file = FALSE, numberOfCores = N_CORES)),
            envir = environment())
  }
  hand_processed <- cached("hand_processed", hand_cache, run_dp(hand_converted))
  llm_processed  <- cached("llm_processed",  llm_cache,  run_dp(llm_converted))

  # 5. groupComparison ---------------------------------------------------
  levs <- levels_of(hand_processed)
  if (!identical(levs, levels_of(llm_processed))) stop("The two paths produced different conditions.")
  contrast <- cfg$build_contrast(levs)
  say("5. groupComparison: %s", paste(rownames(contrast), collapse = ", "))
  run_gc <- function(p) {
    groupComparison(contrast.matrix = contrast, data = p, use_log_file = FALSE,
                    numberOfCores = N_CORES)
  }
  hand_cmp <- cached("hand_comparison", hand_cache, run_gc(hand_processed))
  llm_cmp  <- cached("llm_comparison",  llm_cache,  run_gc(llm_processed))

  # 6. Compare the results -----------------------------------------------
  say("6. Comparing")
  prep_res <- function(cmp) {
    r <- as.data.table(cmp$ComparisonResult)
    r[, Protein := protein_key(as.character(Protein))][, Label := as.character(Label)]
    if (!"issue" %in% names(r)) r[, issue := NA_character_]
    r[, species := cfg$species_of(Protein)]
    if (anyNA(r$species)) stop("Some proteins in the comparison result have no truth label.")
    # Proteins MSstats flags with an issue (a condition entirely missing) have
    # no usable test; they are counted separately rather than called significant.
    r[, testable := is.na(issue)]
    r[, sig := testable & !is.na(adj.pvalue) & adj.pvalue < alpha]
    r
  }
  rh <- prep_res(hand_cmp)
  rl <- prep_res(llm_cmp)

  per_path <- function(r, prefix) {
    s <- r[, .(tested = sum(testable), with_issue = sum(!testable), significant = sum(sig),
               tp = sum(sig & species == pos), fp = sum(sig & species == neg),
               pos_median = median(log2FC[testable & species == pos & is.finite(log2FC)]),
               neg_median = median(log2FC[testable & species == neg & is.finite(log2FC)])),
           by = Label]
    setnames(s, c("pos_median", "neg_median"), paste0(tolower(c(pos, neg)), "_median_log2fc"))
    s[, fdr := fifelse(tp + fp > 0, fp / (tp + fp), NA_real_)]
    setnames(s, setdiff(names(s), "Label"), paste0(prefix, setdiff(names(s), "Label")))
    s
  }

  overlap <- rbindlist(lapply(rownames(contrast), function(lab) {
    sig_h <- rh[Label == lab & sig, unique(Protein)]
    sig_l <- rl[Label == lab & sig, unique(Protein)]
    ht <- rh[Label == lab & testable, unique(Protein)]
    lt <- rl[Label == lab & testable, unique(Protein)]
    fc <- merge(rh[Label == lab & testable, .(Protein, h = log2FC)],
                rl[Label == lab & testable, .(Protein, l = log2FC)], by = "Protein")
    fc <- fc[is.finite(h) & is.finite(l)]
    data.table(Label = lab,
               sig_shared = length(intersect(sig_h, sig_l)),
               sig_only_hand = length(setdiff(sig_h, sig_l)), sig_only_llm = length(setdiff(sig_l, sig_h)),
               tested_only_hand = length(setdiff(ht, lt)), tested_only_llm = length(setdiff(lt, ht)),
               log2fc_r = if (nrow(fc) > 2) cor(fc$h, fc$l) else NA_real_, log2fc_n = nrow(fc))
  }))

  summ <- Reduce(function(x, y) merge(x, y, by = "Label"),
                 list(per_path(rh, "hand_"), per_path(rl, "llm_"), overlap))
  has_design <- !is.null(cfg$designed_ratio)
  if (has_design) {
    summ[, designed_ratio := cfg$designed_ratio(Label)]
    summ[, designed_log2fc := log2(designed_ratio)]
  }
  summ <- summ[match(rownames(contrast), Label)]
  if (has_design) setcolorder(summ, c("Label", "designed_ratio", "designed_log2fc"))

  ph <- as.data.table(hand_processed$ProteinLevelData)
  pl <- as.data.table(llm_processed$ProteinLevelData)
  run_c <- if ("originalRUN" %in% names(ph) && "originalRUN" %in% names(pl)) "originalRUN" else "RUN"
  ints <- merge(ph[, .(Protein = protein_key(as.character(Protein)), run = as.character(get(run_c)), hand = LogIntensities)],
                pl[, .(Protein = protein_key(as.character(Protein)), run = as.character(get(run_c)), llm = LogIntensities)],
                by = c("Protein", "run"))
  ints  <- ints[is.finite(hand) & is.finite(llm)]
  int_r <- if (nrow(ints) > 2) cor(ints$hand, ints$llm) else NA_real_

  paired <- merge(rh[testable == TRUE, .(Protein, Label, species, h = log2FC, hp = adj.pvalue)],
                  rl[testable == TRUE, .(Protein, Label, l = log2FC, lp = adj.pvalue)],
                  by = c("Protein", "Label"))
  pf <- paired[is.finite(h) & is.finite(l)]
  pq <- paired[!is.na(hp) & !is.na(lp)]
  max_fc <- if (nrow(pf)) max(abs(pf$h - pf$l)) else NA_real_
  max_p  <- if (nrow(pq)) max(abs(pq$hp - pq$lp)) else NA_real_

  fwrite(summ, summary_csv)

  # 7. Figures -----------------------------------------------------------
  say("7. Figures to %s", out_dir)
  save_plot <- function(p, file, w, h) ggsave(file.path(out_dir, file), p, width = w, height = h)
  lab_levels <- rownames(contrast)
  g_cols <- if (length(lab_levels) <= 4) length(lab_levels) else 3
  g_rows <- ceiling(length(lab_levels) / g_cols)

  venns <- lapply(seq_len(nrow(summ)), function(i) {
    title <- if (has_design) sprintf("%s (%gx)", summ$Label[i], summ$designed_ratio[i]) else summ$Label[i]
    venn_plot(summ$sig_only_hand[i], summ$sig_shared[i], summ$sig_only_llm[i], title)
  })
  if (requireNamespace("cowplot", quietly = TRUE)) {
    save_plot(cowplot::plot_grid(plotlist = venns, ncol = g_cols), "venn_by_contrast.pdf",
              4 * g_cols, 4.5 * g_rows)
  } else {
    for (i in seq_along(venns)) save_plot(venns[[i]], sprintf("venn_%s.pdf", summ$Label[i]), 4.5, 4.5)
  }

  paired[, Label := factor(Label, levels = lab_levels)]
  species_colours <- scale_colour_manual(values = setNames(c(COL_POS, COL_NEG, COL_OTHER),
                                                           c(pos, neg, "OTHER")))
  p_fc <- ggplot(thin(paired[is.finite(h) & is.finite(l)]), aes(h, l, colour = species)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.6, alpha = 0.4) + facet_wrap(~ Label, ncol = g_cols) + species_colours +
    labs(x = "log2FC, hand-written", y = "log2FC, LLM", colour = NULL,
         title = "log2 fold change by contrast") +
    theme_minimal(base_size = 11) + coord_fixed()
  save_plot(p_fc, "log2fc_scatter.pdf", 4 * g_cols, 4.8 * g_rows)

  p_pv <- ggplot(thin(paired[!is.na(hp) & !is.na(lp) & hp > 0 & lp > 0]),
                 aes(-log10(hp), -log10(lp), colour = species)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.6, alpha = 0.4) + facet_wrap(~ Label, ncol = g_cols) + species_colours +
    labs(x = "-log10 adj. p, hand-written", y = "-log10 adj. p, LLM", colour = NULL,
         title = "Adjusted p-values by contrast") +
    theme_minimal(base_size = 11) + coord_fixed()
  save_plot(p_pv, "adjp_scatter.pdf", 4 * g_cols, 4.8 * g_rows)

  p_int <- ggplot(thin(ints), aes(hand, llm)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.4, alpha = 0.3, colour = COL_HAND) +
    labs(x = "log2 intensity, hand-written", y = "log2 intensity, LLM",
         title = "Summarized protein intensities",
         subtitle = sprintf("r = %.4f over %s protein-run values%s", int_r, fmt(nrow(ints)),
                            if (nrow(ints) > 50000) ", 50,000 shown" else "")) +
    theme_minimal(base_size = 11) + coord_fixed()
  save_plot(p_int, "intensity_scatter.pdf", 6, 6)

  pos_cols <- paste0(c("hand_", "llm_"), tolower(pos), "_median_log2fc")
  if (has_design) {
    rr <- melt(summ[, c("Label", "designed_log2fc", pos_cols), with = FALSE],
               id.vars = c("Label", "designed_log2fc"), variable.name = "path", value.name = "observed")
    rr[, path := fifelse(startsWith(as.character(path), "hand_"), "Hand-written", "LLM")]
    p_rr <- ggplot(rr, aes(designed_log2fc, observed, colour = path)) +
      geom_abline(linetype = "dashed", colour = "grey50") +
      geom_point(size = 3, position = position_dodge(width = 0.04)) +
      scale_colour_manual(values = c("Hand-written" = COL_HAND, LLM = COL_LLM)) +
      labs(x = "Designed log2 ratio", y = sprintf("Median %s log2FC", pos), colour = NULL,
           title = "Recovered spike-in ratios") +
      theme_minimal(base_size = 11)
    save_plot(p_rr, "recovered_ratios.pdf", 6, 5)
  }

  fd <- melt(summ[, .(Label, `Hand-written` = hand_fdr, LLM = llm_fdr)],
             id.vars = "Label", variable.name = "path", value.name = "fdr")
  p_fdr <- ggplot(fd, aes(factor(Label, levels = lab_levels), fdr, fill = path)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = c("Hand-written" = COL_HAND, LLM = COL_LLM)) +
    labs(x = NULL, y = "Empirical FDR", fill = NULL,
         title = sprintf("Empirical FDR at adj. p < %.2f", alpha)) +
    theme_minimal(base_size = 11)
  save_plot(p_fdr, "fdr_by_contrast.pdf", max(7, 1.3 * length(lab_levels)), 5)

  versions <- c(R.version.string, pkg_source("MSstats"), pkg_source("MSstatsConvert"))
  saveRDS(list(config = cfg[c("dataset", "label", "description", "positive", "negative", "dp_args")],
               contrast = contrast, cores = N_CORES, mapping = cfg$mapping, converted = cc,
               summary = summ, intensity_r = int_r, max_log2fc_diff = max_fc, max_adjp_diff = max_p,
               hand = rh, llm = rl, versions = versions),
          file.path(out_dir, "full_results.rds"))

  # Results --------------------------------------------------------------
  cat("\n================ RESULTS ================\n")
  cat(sprintf("%s, %s\n", cfg$dataset, cfg$description))
  cat(sprintf("Converted data identical: %s\n", cc$identical))
  if (!cc$identical) {
    cat(sprintf("  hand %s rows, LLM %s rows; only in hand %s, only in LLM %s\n",
                fmt(cc$rows_hand), fmt(cc$rows_llm), fmt(cc$only_hand), fmt(cc$only_llm)))
    cat(sprintf("  shared rows with a different intensity %s (blank in LLM only %s, blank in hand only %s)\n",
                fmt(cc$intensity_differs), fmt(cc$blank_in_llm_only), fmt(cc$blank_in_hand_only)))
    cat(sprintf("  other columns that differ: %s\n",
                if (length(cc$other_columns_differ)) paste(cc$other_columns_differ, collapse = ", ") else "none"))
    if (cc$duplicate_keys_hand || cc$duplicate_keys_llm) {
      cat(sprintf("  repeated keys collapsed before matching: hand %s, LLM %s\n",
                  fmt(cc$duplicate_keys_hand), fmt(cc$duplicate_keys_llm)))
    }
  }
  cat(sprintf("Features: hand %s, LLM %s\n", fmt(cc$features_hand), fmt(cc$features_llm)))
  cat(sprintf("Summarized intensities: r = %.4f over %s values\n", int_r, fmt(nrow(ints))))
  cat(sprintf("Largest difference over shared tested proteins: log2FC %.1e, adj. p %.1e\n", max_fc, max_p))
  cat(sprintf("Significance: adj. p < %.2f, excluding proteins MSstats flags with an issue\n", alpha))
  cat(sprintf("FDR counts %s proteins as true positives and %s proteins as true negatives\n\n", pos, neg))
  for (i in seq_len(nrow(summ))) {
    s <- summ[i]
    n_union <- s$sig_shared + s$sig_only_hand + s$sig_only_llm
    if (has_design) {
      cat(sprintf("%s  designed %gx (log2 %.2f)\n", s$Label, s$designed_ratio, s$designed_log2fc))
    } else {
      cat(s$Label, "\n", sep = "")
    }
    cat(sprintf("  significant      hand %s   LLM %s   shared %s   only hand %s   only LLM %s   agreement %s\n",
                fmt(s$hand_significant), fmt(s$llm_significant), fmt(s$sig_shared),
                fmt(s$sig_only_hand), fmt(s$sig_only_llm),
                if (n_union > 0) sprintf("%.1f%%", 100 * s$sig_shared / n_union) else "NA"))
    cat(sprintf("  empirical FDR    hand %.3f   LLM %.3f\n", s$hand_fdr, s$llm_fdr))
    cat(sprintf("  %s median   hand %.3f   LLM %.3f%s\n", pos, s[[pos_cols[1]]], s[[pos_cols[2]]],
                if (has_design) sprintf("   designed %.3f", s$designed_log2fc) else ""))
    cat(sprintf("  tested           hand %s   LLM %s   only hand %s   only LLM %s\n",
                fmt(s$hand_tested), fmt(s$llm_tested), fmt(s$tested_only_hand), fmt(s$tested_only_llm)))
    cat(sprintf("  with an issue    hand %s   LLM %s\n", fmt(s$hand_with_issue), fmt(s$llm_with_issue)))
    cat(sprintf("  log2FC r = %.4f over %s proteins\n\n", s$log2fc_r, fmt(s$log2fc_n)))
  }
  cat(sprintf("Versions: %s\n", paste(versions, collapse = "; ")))
  cat(sprintf("Summary table: %s\nFigures: %s\n", summary_csv, out_dir))

  invisible(list(summary = summ, hand = rh, llm = rl, converted = cc, contrast = contrast))
}
