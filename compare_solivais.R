# ============================================================
# compare_solivais.R: downstream comparison on Solivais2024 (MetaMorpheus)
#
# E. coli is spiked into a constant human background at designed ratios, so a
# significant human protein is a false positive. Both paths use dataProcess
# with medians equalized, the setting used on all three datasets, which assumes
# most proteins do not change in abundance. The shared machinery is in
# downstream.R.
#
# Run from the repository root:
#   caffeinate -i Rscript compare_solivais.R
#   LLM_MODEL=llama3.1-8b LLM_PROMPT=lean caffeinate -i Rscript compare_solivais.R
# ============================================================

source("downstream.R")

DATA_DIR  <- Sys.getenv("SOLIVAIS_DATA_DIR", "data/Solivais2024")
ECOLI_ORG <- "Escherichia coli (strain K12)"
HUMAN_ORG <- "Homo sapiens"

# Designed E. coli ratio of each condition against A, from code-tony.R
REFERENCE   <- "A"
SPIKE_RATIO <- c(E = 3, D = 2.5, C = 2, B = 1.5)

say("Downstream comparison: Solivais2024 (MetaMorpheus), %d cores", N_CORES)

chosen <- choose_mapping("metamorpheus", c(
  ProteinName = "Protein.Group", PeptideSequence = "Full.Sequence",
  PrecursorCharge = "Precursor.Charge", FragmentIon = "null", ProductCharge = "null",
  Run = "File.Name", Intensity = "Peak.intensity", Qvalue = "PIP.Q.Value"))

# The species label is appended to Protein Group, so it only reaches the LLM
# path's ProteinName if the mapping used that column.
if (!identical(unname(mapping_vector(chosen$mapping)["ProteinName"]), "Protein.Group")) {
  warning("ProteinName maps to a column other than Protein.Group; species labels may not carry through.")
}

files <- file.path(DATA_DIR, c("QuantifiedPeaks.tsv", "QuantifiedProteins.tsv", "annotation.csv"))
missing_inputs <- files[!file.exists(files)]
if (length(missing_inputs)) stop("Missing input files:\n  ", paste(missing_inputs, collapse = "\n  "))

say("2. Reading %s", DATA_DIR)
peaks      <- fread(files[1], showProgress = FALSE)
proteins   <- fread(files[2], showProgress = FALSE)
annotation <- fread(files[3], showProgress = FALSE)
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

run_downstream(list(
  dataset = "Solivais2024", slug = "solivais",
  label = chosen$label, description = chosen$description, mapping = chosen$mapping,
  raw = pre, annotation = annotation, run_col = "File Name",
  hand_convert = function(x, a) MetamorpheusToMSstatsFormat(x, annotation = a, use_log_file = FALSE),
  species_of = function(p) fcase(grepl("\\|ECOLI$", p), "ECOLI", grepl("\\|HUMAN$", p), "HUMAN"),
  positive = "ECOLI", negative = "HUMAN",
  dp_args = list(normalization = "equalizeMedians"),
  build_contrast = function(levs) {
    if (!setequal(levs, c(REFERENCE, names(SPIKE_RATIO)))) {
      stop("Unexpected conditions: ", paste(levs, collapse = ", "))
    }
    cm <- matrix(0, nrow = length(SPIKE_RATIO), ncol = length(levs),
                 dimnames = list(paste0(names(SPIKE_RATIO), "-", REFERENCE), levs))
    for (g in names(SPIKE_RATIO)) cm[paste0(g, "-", REFERENCE), c(g, REFERENCE)] <- c(1, -1)
    cm
  },
  designed_ratio = function(lab) unname(SPIKE_RATIO[sub(paste0("-", REFERENCE, "$"), "", lab)]),
  alpha = 0.05))
