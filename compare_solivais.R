# ============================================================
# compare_solivais.R: downstream comparison on Solivais2024 (MetaMorpheus)
#
# One raw file goes through the hand-written converter and the LLM-driven
# converter, then dataProcess and groupComparison, and the two results are
# compared. Solivais2024 is a controlled mixture: E. coli is spiked into a
# constant human background at designed ratios, so a significant human
# protein is a false positive.
#
# The LLM mapping comes from the saved benchmark run rather than a live
# model call, so the result uses exactly the mapping behind the accuracy
# table and needs no GPU.
#
# Run from the repository root:
#   caffeinate -i Rscript compare_solivais.R
# Slow steps are cached under results/. Set FORCE=1 to recompute, which is
# required after any change to the converter code.
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

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
MODEL  <- Sys.getenv("LLM_MODEL",  "deepseek-r1-14b")
PROMPT <- Sys.getenv("LLM_PROMPT", "filter_aware")
TOOL   <- "metamorpheus"

DATA_DIR   <- Sys.getenv("SOLIVAIS_DATA_DIR", "data/Solivais2024")
PEAKS_TSV  <- file.path(DATA_DIR, "QuantifiedPeaks.tsv")
PROT_TSV   <- file.path(DATA_DIR, "QuantifiedProteins.tsv")
ANNOT_CSV  <- file.path(DATA_DIR, "annotation.csv")
TRIALS_RDS <- "results/raw_trials_benchmark_20261007_163702.rds"

ECOLI_ORG <- "Escherichia coli (strain K12)"
HUMAN_ORG <- "Homo sapiens"

# Designed E. coli ratio of each condition against A, from code-tony.R
REFERENCE   <- "A"
SPIKE_RATIO <- c(E = 3, D = 2.5, C = 2, B = 1.5)
ALPHA       <- 0.05

n_detected <- parallel::detectCores()
if (is.na(n_detected)) n_detected <- 2L
N_CORES <- as.integer(Sys.getenv("N_CORES", as.character(max(1L, min(8L, n_detected - 1L)))))
FORCE   <- nzchar(Sys.getenv("FORCE"))

RUN_NAME    <- sprintf("downstream_solivais_%s_%s", MODEL, PROMPT)
OUT_DIR     <- file.path("results", RUN_NAME)                       # not tracked
LLM_CACHE   <- file.path(OUT_DIR, "cache")
HAND_CACHE  <- file.path("results", "downstream_solivais_hand", "cache")  # shared across configs
SUMMARY_CSV <- file.path("results", paste0(RUN_NAME, "_summary.csv"))    # tracked
for (d in c(LLM_CACHE, HAND_CACHE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

COL_HAND  <- "#2C5F8A"
COL_LLM   <- "#8A5F2C"
COL_ECOLI <- "#2C8A5F"
COL_HUMAN <- "grey65"

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

missing_inputs <- Filter(Negate(file.exists), c(PEAKS_TSV, PROT_TSV, ANNOT_CSV, TRIALS_RDS))
if (length(missing_inputs)) {
  stop("Missing input files:\n  ", paste(missing_inputs, collapse = "\n  "))
}

say("Downstream comparison: Solivais2024, %s, %s prompt, %d cores", MODEL, PROMPT, N_CORES)

# ==================================================================
# 1. LLM mapping from the saved benchmark run
# ==================================================================
say("1. Mapping from %s", TRIALS_RDS)

mapping_vector <- function(m) {
  mp <- m$mappings
  if (is.data.frame(mp)) return(setNames(as.character(mp$from), as.character(mp$field)))
  setNames(vapply(mp, function(x) as.character(x[["from"]] %||% NA_character_)[1], ""),
           vapply(mp, function(x) as.character(x[["field"]])[1], ""))
}

trials <- readRDS(TRIALS_RDS)
picked <- Filter(function(t) identical(t$model, MODEL) && identical(t$tool, TOOL) &&
                   identical(t$prompt, PROMPT), trials)
rm(trials)
if (!length(picked)) stop("No saved trials for ", MODEL, " / ", TOOL, " / ", PROMPT)

vecs  <- lapply(picked, function(t) { v <- mapping_vector(t$mapping); v[order(names(v))] })
agree <- all(vapply(vecs, identical, logical(1), vecs[[1]]))
say("  %d saved repeats, mapping identical across repeats: %s", length(picked), agree)
if (!agree) warning("Repeats disagree on the mapping; using the first.")

llm_mapping <- picked[[1]]$mapping
mv <- mapping_vector(llm_mapping)
for (f in names(mv)) say("    %-16s <- %s", f, mv[[f]])

fl <- llm_mapping$filters
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

# The species label is appended to Protein Group, so it only reaches the
# LLM path's ProteinName if the model mapped ProteinName to that column.
if (!identical(unname(mv["ProteinName"]), "Protein.Group")) {
  warning("ProteinName maps to ", mv["ProteinName"], " rather than Protein.Group.")
}

# ==================================================================
# 2. Read and preprocess, identically for both paths (as code-tony.R)
# ==================================================================
say("2. Reading %s", DATA_DIR)
peaks      <- fread(PEAKS_TSV, showProgress = FALSE)
proteins   <- fread(PROT_TSV,  showProgress = FALSE)
annotation <- fread(ANNOT_CSV, showProgress = FALSE)
say("  QuantifiedPeaks: %s rows, %d columns", fmt(nrow(peaks)), ncol(peaks))

# Whether match-between-runs peaks are present decides whether the
# MetaMorpheus converter's two-column MBR filter has anything to act on.
if ("Peak Detection Type" %in% names(peaks)) {
  tb <- peaks[, .N, by = "Peak Detection Type"]
  say("  Peak Detection Type: %s", paste0(tb[[1]], " ", fmt(tb$N), collapse = ", "))
}
if ("PIP Q-Value" %in% names(peaks)) {
  pq <- peaks$`PIP Q-Value`
  say("  PIP Q-Value present on %s of %s rows",
      fmt(sum(!is.na(pq) & nzchar(as.character(pq)))), fmt(nrow(peaks)))
}

preprocess <- function(peaks, proteins) {
  rows <- c(raw = nrow(peaks))
  peaks <- peaks[!grepl(";", `Protein Group`, fixed = TRUE)]
  rows["single protein group"] <- nrow(peaks)
  peaks <- peaks[!grepl("DECOY", `Protein Group`, fixed = TRUE)]
  rows["decoys removed"] <- nrow(peaks)

  ecoli <- unique(proteins[Organism == ECOLI_ORG, `Protein Groups`])
  human <- unique(proteins[Organism == HUMAN_ORG, `Protein Groups`])
  if (length(intersect(ecoli, human))) stop("Some protein groups are listed under both organisms.")

  # Drops the contaminant species and any group absent from QuantifiedProteins
  peaks <- peaks[`Protein Group` %chin% c(ecoli, human)]
  rows["E. coli or human only"] <- nrow(peaks)

  peaks[, `Protein Group` := paste0(`Protein Group`,
                                    fifelse(`Protein Group` %chin% ecoli, "|ECOLI", "|HUMAN"))]
  list(data = peaks, rows = rows)
}

pp  <- preprocess(peaks, proteins)
pre <- pp$data
rm(peaks); invisible(gc())
for (s in names(pp$rows)) say("  %-24s %s rows", s, fmt(pp$rows[[s]]))
sp <- pre[, .(groups = uniqueN(`Protein Group`)), by = .(species = sub(".*\\|", "", `Protein Group`))]
say("  protein groups kept: %s", paste0(sp$species, " ", fmt(sp$groups), collapse = ", "))

runs  <- unique(pre$`File Name`)
miss  <- setdiff(annotation$Run, runs)
extra <- setdiff(runs, annotation$Run)
if (length(miss) || length(extra)) {
  stop("Annotation and File Name do not match.",
       "\n  only in annotation: ", paste(head(miss, 3), collapse = ", "),
       "\n  only in data: ", paste(head(extra, 3), collapse = ", "))
}
say("  %d runs, all matched to the annotation", length(runs))

# ==================================================================
# 3. Convert with both paths
# ==================================================================
say("3. Converting")
hand_converted <- cached("hand_converted", HAND_CACHE, as.data.table(
  MetamorpheusToMSstatsFormat(copy(pre), annotation = copy(annotation), use_log_file = FALSE)))

llm_out <- cached("llm_converted", LLM_CACHE, {
  out <- LLMtoMSstatsFormat(input = normalize_headers(copy(pre)), annotation = copy(annotation),
                            llm_mapping = llm_mapping, use_log_file = FALSE)
  list(data = as.data.table(out), diagnostics = attr(out, "diagnostics"))
})
llm_converted <- llm_out$data

dg <- llm_out$diagnostics
if (!is.null(dg)) {
  if (nrow(dg$hallucinated_columns)) {
    say("  hallucinated columns: %s", paste(dg$hallucinated_columns$hallucinated_name, collapse = ", "))
  }
  for (i in seq_len(nrow(dg$applied_filters))) {
    a <- dg$applied_filters[i]
    say("  filter applied: %s %s %s (%s)", as.character(a$column), as.character(a$operation),
        as.character(unlist(a$value)), as.character(a$kind %||% ""))
  }
  for (i in seq_len(nrow(dg$skipped_filters))) {
    s <- dg$skipped_filters[i]
    say("  filter skipped: %s %s %s (%s)", as.character(s$column), as.character(s$operation),
        as.character(unlist(s$value)), as.character(s$reason %||% ""))
  }
}

# Fail before the slow step if the species label did not reach ProteinName
for (which_path in c("hand", "llm")) {
  x <- if (which_path == "hand") hand_converted else llm_converted
  bad <- sum(!grepl("\\|(ECOLI|HUMAN)$", unique(as.character(x$ProteinName))))
  if (bad) stop(which_path, " path: ", bad, " protein names have no species label.")
}

compare_converted <- function(a, b) {
  prep <- function(x) {
    x <- copy(x)
    for (j in names(x)) if (is.factor(x[[j]])) set(x, j = j, value = as.character(x[[j]]))
    x
  }
  a <- prep(a); b <- prep(b)
  cols <- intersect(names(a), names(b))
  a <- a[, ..cols]; b <- b[, ..cols]
  setorderv(a, cols, na.last = TRUE); setorderv(b, cols, na.last = TRUE)
  if (isTRUE(all.equal(a, b, check.attributes = FALSE))) {
    return(list(identical = TRUE, rows_hand = nrow(a), rows_llm = nrow(b)))
  }

  key <- intersect(c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                     "ProductCharge", "IsotopeLabelType", "Run"), cols)
  key <- key[!vapply(key, function(k) all(is.na(a[[k]])) && all(is.na(b[[k]])), logical(1))]
  ka <- a[, c(key, "Intensity"), with = FALSE][, in_hand := TRUE]
  kb <- b[, c(key, "Intensity"), with = FALSE][, in_llm := TRUE]
  m  <- merge(ka, kb, by = key, all = TRUE, suffixes = c("_hand", "_llm"))
  both <- m[in_hand %in% TRUE & in_llm %in% TRUE]
  ih <- both$Intensity_hand; il <- both$Intensity_llm
  same <- (is.na(ih) & is.na(il)) |
          (!is.na(ih) & !is.na(il) & abs(ih - il) <= 1e-8 * pmax(1, abs(ih)))
  list(identical = FALSE, rows_hand = nrow(a), rows_llm = nrow(b),
       only_hand = m[is.na(in_llm), .N], only_llm = m[is.na(in_hand), .N],
       intensity_differs = sum(!same),
       blank_in_llm_only = sum(!is.na(ih) & is.na(il)),
       blank_in_hand_only = sum(is.na(ih) & !is.na(il)))
}

cc <- compare_converted(hand_converted, llm_converted)
if (cc$identical) {
  say("  converted data identical: %s rows in each path", fmt(cc$rows_hand))
} else {
  say("  converted data differ: hand %s rows, LLM %s rows", fmt(cc$rows_hand), fmt(cc$rows_llm))
  say("    rows only in hand %s, only in LLM %s", fmt(cc$only_hand), fmt(cc$only_llm))
  say("    shared rows with a different intensity %s (blank in LLM only %s, blank in hand only %s)",
      fmt(cc$intensity_differs), fmt(cc$blank_in_llm_only), fmt(cc$blank_in_hand_only))
}

# ==================================================================
# 4. dataProcess, same settings on both paths (as code-tony.R)
# ==================================================================
say("4. dataProcess; this is the slow step")
run_dp <- function(x) {
  dataProcess(x, normalization = FALSE, use_log_file = FALSE, numberOfCores = N_CORES)
}
hand_processed <- cached("hand_processed", HAND_CACHE, run_dp(hand_converted))
llm_processed  <- cached("llm_processed",  LLM_CACHE,  run_dp(llm_converted))

# ==================================================================
# 5. groupComparison, each spiked condition against A
# ==================================================================
levels_of <- function(p) {
  g <- p$ProteinLevelData$GROUP
  if (is.factor(g)) levels(g) else sort(unique(as.character(g)))
}
levs <- levels_of(hand_processed)
if (!identical(levs, levels_of(llm_processed))) stop("The two paths produced different conditions.")
if (!setequal(levs, c(REFERENCE, names(SPIKE_RATIO)))) {
  stop("Unexpected conditions: ", paste(levs, collapse = ", "))
}

contrast <- matrix(0, nrow = length(SPIKE_RATIO), ncol = length(levs),
                   dimnames = list(paste0(names(SPIKE_RATIO), "-", REFERENCE), levs))
for (g in names(SPIKE_RATIO)) contrast[paste0(g, "-", REFERENCE), c(g, REFERENCE)] <- c(1, -1)

say("5. groupComparison: %s", paste(rownames(contrast), collapse = ", "))
run_gc <- function(p) {
  groupComparison(contrast.matrix = contrast, data = p, use_log_file = FALSE,
                  numberOfCores = N_CORES)
}
hand_cmp <- cached("hand_comparison", HAND_CACHE, run_gc(hand_processed))
llm_cmp  <- cached("llm_comparison",  LLM_CACHE,  run_gc(llm_processed))

# ==================================================================
# 6. Compare the results
# ==================================================================
say("6. Comparing")

prep_res <- function(cmp) {
  r <- as.data.table(cmp$ComparisonResult)
  r[, Protein := as.character(Protein)][, Label := as.character(Label)]
  if (!"issue" %in% names(r)) r[, issue := NA_character_]
  r[, species := fcase(grepl("\\|ECOLI$", Protein), "ECOLI", grepl("\\|HUMAN$", Protein), "HUMAN")]
  if (anyNA(r$species)) stop("Some proteins in the comparison result have no species label.")
  # Proteins MSstats flags with an issue (a condition entirely missing) have
  # no usable test; they are counted separately rather than called significant.
  r[, testable := is.na(issue)]
  r[, sig := testable & !is.na(adj.pvalue) & adj.pvalue < ALPHA]
  r
}
rh <- prep_res(hand_cmp)
rl <- prep_res(llm_cmp)

per_path <- function(r, prefix) {
  s <- r[, .(tested = sum(testable), with_issue = sum(!testable), significant = sum(sig),
             tp = sum(sig & species == "ECOLI"), fp = sum(sig & species == "HUMAN"),
             ecoli_median_log2fc = median(log2FC[testable & species == "ECOLI" & is.finite(log2FC)]),
             human_median_log2fc = median(log2FC[testable & species == "HUMAN" & is.finite(log2FC)])),
         by = Label]
  s[, fdr := fifelse(significant > 0, fp / significant, NA_real_)]
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
summ[, designed_ratio := unname(SPIKE_RATIO[sub(paste0("-", REFERENCE, "$"), "", Label)])]
summ[, designed_log2fc := log2(designed_ratio)]
summ <- summ[match(rownames(contrast), Label)]
setcolorder(summ, c("Label", "designed_ratio", "designed_log2fc"))

ph <- as.data.table(hand_processed$ProteinLevelData)
pl <- as.data.table(llm_processed$ProteinLevelData)
run_col <- if ("originalRUN" %in% names(ph) && "originalRUN" %in% names(pl)) "originalRUN" else "RUN"
ints <- merge(ph[, .(Protein = as.character(Protein), run = as.character(get(run_col)), hand = LogIntensities)],
              pl[, .(Protein = as.character(Protein), run = as.character(get(run_col)), llm = LogIntensities)],
              by = c("Protein", "run"))
ints  <- ints[is.finite(hand) & is.finite(llm)]
int_r <- if (nrow(ints) > 2) cor(ints$hand, ints$llm) else NA_real_

fwrite(summ, SUMMARY_CSV)

# ==================================================================
# 7. Figures
# ==================================================================
say("7. Figures to %s", OUT_DIR)
thin <- function(d, n = 50000) if (nrow(d) > n) { set.seed(1); d[sample(.N, n)] } else d
save_plot <- function(p, file, w, h) ggsave(file.path(OUT_DIR, file), p, width = w, height = h)
lab_levels <- rownames(contrast)

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

venns <- lapply(seq_len(nrow(summ)), function(i) {
  venn_plot(summ$sig_only_hand[i], summ$sig_shared[i], summ$sig_only_llm[i],
            sprintf("%s (%gx)", summ$Label[i], summ$designed_ratio[i]))
})
if (requireNamespace("cowplot", quietly = TRUE)) {
  save_plot(cowplot::plot_grid(plotlist = venns, nrow = 1), "venn_by_contrast.pdf", 16, 4.5)
} else {
  for (i in seq_along(venns)) save_plot(venns[[i]], sprintf("venn_%s.pdf", summ$Label[i]), 4.5, 4.5)
}

paired <- merge(rh[testable == TRUE, .(Protein, Label, species, h = log2FC, hp = adj.pvalue)],
               rl[testable == TRUE, .(Protein, Label, l = log2FC, lp = adj.pvalue)],
               by = c("Protein", "Label"))
paired[, Label := factor(Label, levels = lab_levels)]
species_colours <- scale_colour_manual(values = c(ECOLI = COL_ECOLI, HUMAN = COL_HUMAN))

p_fc <- ggplot(thin(paired[is.finite(h) & is.finite(l)]), aes(h, l, colour = species)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.6, alpha = 0.4) + facet_wrap(~ Label, nrow = 1) + species_colours +
  labs(x = "log2FC, hand-written", y = "log2FC, LLM", colour = NULL,
       title = "log2 fold change by contrast") +
  theme_minimal(base_size = 11) + coord_fixed()
save_plot(p_fc, "log2fc_scatter.pdf", 16, 4.8)

p_pv <- ggplot(thin(paired[!is.na(hp) & !is.na(lp) & hp > 0 & lp > 0]),
               aes(-log10(hp), -log10(lp), colour = species)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.6, alpha = 0.4) + facet_wrap(~ Label, nrow = 1) + species_colours +
  labs(x = "-log10 adj. p, hand-written", y = "-log10 adj. p, LLM", colour = NULL,
       title = "Adjusted p-values by contrast") +
  theme_minimal(base_size = 11) + coord_fixed()
save_plot(p_pv, "adjp_scatter.pdf", 16, 4.8)

p_int <- ggplot(thin(ints), aes(hand, llm)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.3, colour = COL_HAND) +
  labs(x = "log2 intensity, hand-written", y = "log2 intensity, LLM",
       title = "Summarized protein intensities",
       subtitle = sprintf("r = %.4f over %s protein-run values%s", int_r, fmt(nrow(ints)),
                          if (nrow(ints) > 50000) ", 50,000 shown" else "")) +
  theme_minimal(base_size = 11) + coord_fixed()
save_plot(p_int, "intensity_scatter.pdf", 6, 6)

rr <- melt(summ[, .(Label, designed_log2fc,
                    `Hand-written` = hand_ecoli_median_log2fc, LLM = llm_ecoli_median_log2fc)],
           id.vars = c("Label", "designed_log2fc"), variable.name = "path", value.name = "observed")
p_rr <- ggplot(rr, aes(designed_log2fc, observed, colour = path)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.04)) +
  scale_colour_manual(values = c("Hand-written" = COL_HAND, LLM = COL_LLM)) +
  labs(x = "Designed log2 ratio", y = "Median E. coli log2FC", colour = NULL,
       title = "Recovered spike-in ratios") +
  theme_minimal(base_size = 11)
save_plot(p_rr, "recovered_ratios.pdf", 6, 5)

fd <- melt(summ[, .(Label, `Hand-written` = hand_fdr, LLM = llm_fdr)],
           id.vars = "Label", variable.name = "path", value.name = "fdr")
p_fdr <- ggplot(fd, aes(factor(Label, levels = lab_levels), fdr, fill = path)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("Hand-written" = COL_HAND, LLM = COL_LLM)) +
  labs(x = NULL, y = "Empirical FDR", fill = NULL,
       title = sprintf("Empirical FDR at adj. p < %.2f", ALPHA)) +
  theme_minimal(base_size = 11)
save_plot(p_fdr, "fdr_by_contrast.pdf", 7, 5)

versions <- list(R = R.version.string,
                 MSstats = as.character(packageVersion("MSstats")),
                 MSstatsConvert = as.character(packageVersion("MSstatsConvert")))
saveRDS(list(config = list(model = MODEL, prompt = PROMPT, alpha = ALPHA, contrast = contrast,
                           cores = N_CORES),
             mapping = llm_mapping, preprocessing = pp$rows, converted = cc,
             summary = summ, intensity_r = int_r, hand = rh, llm = rl, versions = versions),
        file.path(OUT_DIR, "full_results.rds"))

# ==================================================================
# Results
# ==================================================================
cat("\n================ RESULTS ================\n")
cat(sprintf("Solivais2024, %s, %s prompt\n", MODEL, PROMPT))
cat(sprintf("Converted data identical: %s\n", cc$identical))
cat(sprintf("Summarized intensities: r = %.4f over %s values\n", int_r, fmt(nrow(ints))))
cat(sprintf("Significance: adj. p < %.2f, excluding proteins MSstats flags with an issue\n\n", ALPHA))
for (i in seq_len(nrow(summ))) {
  s <- summ[i]
  cat(sprintf("%s  designed %gx (log2 %.2f)\n", s$Label, s$designed_ratio, s$designed_log2fc))
  cat(sprintf("  significant      hand %s   LLM %s   shared %s   only hand %s   only LLM %s\n",
              fmt(s$hand_significant), fmt(s$llm_significant), fmt(s$sig_shared),
              fmt(s$sig_only_hand), fmt(s$sig_only_llm)))
  cat(sprintf("  empirical FDR    hand %.3f   LLM %.3f\n", s$hand_fdr, s$llm_fdr))
  cat(sprintf("  E. coli median   hand %.3f   LLM %.3f   designed %.3f\n",
              s$hand_ecoli_median_log2fc, s$llm_ecoli_median_log2fc, s$designed_log2fc))
  cat(sprintf("  tested           hand %s   LLM %s   only hand %s   only LLM %s\n",
              fmt(s$hand_tested), fmt(s$llm_tested), fmt(s$tested_only_hand), fmt(s$tested_only_llm)))
  cat(sprintf("  with an issue    hand %s   LLM %s\n", fmt(s$hand_with_issue), fmt(s$llm_with_issue)))
  cat(sprintf("  log2FC r = %.4f over %s proteins\n\n", s$log2fc_r, fmt(s$log2fc_n)))
}
cat(sprintf("Versions: %s, MSstats %s, MSstatsConvert %s\n",
            versions$R, versions$MSstats, versions$MSstatsConvert))
cat(sprintf("Summary table: %s\nFigures: %s\n", SUMMARY_CSV, OUT_DIR))
