# ============================================================
# compare_pxd005642.R: downstream comparison on PXD005642 (Proteome Discoverer)
#
# Thirty commercial proteins in three subsets are spiked into a constant
# E. coli background across five mixtures, in triplicate. The design is from
# the lab's MassIVE reanalysis of this dataset, RMSV000000263. Spike-ins carry
# an SSn_ prefix on their accession, so they are labelled from the protein
# name itself.
#
# Only contrasts in which all three subsets change amount are tested. In
# those, every spike-in is a true positive and every E. coli protein a true
# negative, so the false discovery rate is exact without knowing each
# protein's subset.
#
# dataProcess follows that reanalysis: top-3 features, otherwise the MSstats
# defaults. The shared machinery is in downstream.R.
#
# Run from the repository root:
#   LLM_MODEL=reference caffeinate -i Rscript compare_pxd005642.R
#   caffeinate -i Rscript compare_pxd005642.R
#   LLM_MODEL=llama3.1-8b LLM_PROMPT=constrained_v3 caffeinate -i Rscript compare_pxd005642.R
# ============================================================

source("downstream.R")

DATA_DIR <- Sys.getenv("PXD005642_DATA_DIR", "data/PXD005642")

# Relative amount of each subset in mixtures 1 to 5 (100, 200 and 400 fmol/ug)
DESIGN <- rbind(SS1 = c(1, 2, 4, 2, 1),
                SS2 = c(4, 1, 2, 1, 2),
                SS3 = c(2, 4, 1, 2, 1))
# The six pairs of mixtures in which every subset changes amount
PAIRS <- list(c(1, 2), c(1, 3), c(2, 3), c(2, 5), c(3, 4), c(4, 5))
stopifnot(all(vapply(PAIRS, function(p) all(DESIGN[, p[1]] != DESIGN[, p[2]]), logical(1))))

say("Downstream comparison: PXD005642 (Proteome Discoverer), %d cores", N_CORES)

# The reference columns are PDtoMSstatsFormat's own defaults, which.proteinid,
# which.sequence and which.quantification, with its charge and run columns.
chosen <- choose_mapping("proteome_discoverer", c(
  ProteinName = "Protein.Group.Accessions", PeptideSequence = "Sequence",
  PrecursorCharge = "Charge", FragmentIon = "null", ProductCharge = "null",
  Run = "Spectrum.File", Intensity = "Precursor.Area", Qvalue = "null"))

files <- file.path(DATA_DIR, c("ControlMixture_DDA_ProteomeDiscoverer_input.csv", "annotation.csv"))
missing_inputs <- files[!file.exists(files)]
if (length(missing_inputs)) stop("Missing input files:\n  ", paste(missing_inputs, collapse = "\n  "))

say("2. Reading %s", DATA_DIR)
raw        <- fread(files[1], showProgress = FALSE)
annotation <- fread(files[2], showProgress = FALSE)
say("  input: %s rows, %d columns", fmt(nrow(raw)), ncol(raw))

# The hand-written converter combines Modifications with Sequence. A mapping
# names one column per field, so this count bounds what the LLM path can match.
mods <- as.character(raw$Modifications)
say("  rows carrying a modification: %s of %s", fmt(sum(!is.na(mods) & nzchar(mods))), fmt(nrow(raw)))

# SPIKED when every accession in the group is a spike-in, ECOLI when none is,
# OTHER for a group mixing the two or an empty name
species_of <- function(p) {
  p <- as.character(p)
  parts <- strsplit(p, ";", fixed = TRUE)
  vapply(seq_along(p), function(i) {
    if (is.na(p[i])) return("OTHER")
    x <- trimws(parts[[i]])
    x <- x[nzchar(x)]
    if (!length(x)) return("OTHER")
    s <- grepl("^SS[0-9]+_", x)
    if (all(s)) "SPIKED" else if (!any(s)) "ECOLI" else "OTHER"
  }, character(1))
}

res <- run_downstream(list(
  dataset = "PXD005642", slug = "pxd005642",
  label = chosen$label, description = chosen$description, mapping = chosen$mapping,
  raw = raw, annotation = annotation, run_col = "Spectrum File",
  hand_convert = function(x, a) PDtoMSstatsFormat(x, annotation = a, use_log_file = FALSE),
  species_of = species_of,
  # The converter writes Sequence_Modifications, and an unmodified peptide as
  # Sequence_ with nothing after it. Dropping that bare underscore lets
  # unmodified rows match; the modified ones are what the LLM path merges.
  sequence_key = function(s) sub("_$", "", s),
  positive = "SPIKED", negative = "ECOLI",
  dp_args = list(normalization = "equalizeMedians", featureSubset = "topN", n_top_feature = 3),
  build_contrast = function(levs) {
    need <- paste0("Condition", 1:5)
    if (!setequal(levs, need)) stop("Unexpected conditions: ", paste(levs, collapse = ", "))
    labs <- vapply(PAIRS, function(p) sprintf("Mix%d-Mix%d", p[1], p[2]), "")
    cm <- matrix(0, nrow = length(PAIRS), ncol = length(levs), dimnames = list(labs, levs))
    for (k in seq_along(PAIRS)) {
      cm[k, paste0("Condition", PAIRS[[k]][1])] <- 1
      cm[k, paste0("Condition", PAIRS[[k]][2])] <- -1
    }
    cm
  },
  designed_ratio = NULL,
  alpha = 0.05))

# Check the design table and the subset labels against the data: each
# subset's spike-ins should sit near that subset's designed log2 ratio.
# SS4 is not in the published design, so its row shows which pattern it follows.
sp <- res$hand[species == "SPIKED" & testable & is.finite(log2FC)]
sp[, subset := sub("^(SS[0-9]+)_.*$", "\\1", Protein)]
observed <- dcast(sp[, .(m = round(median(log2FC), 2)), by = .(subset, Label)],
                  subset ~ Label, value.var = "m")
designed <- dcast(rbindlist(lapply(rownames(DESIGN), function(s) {
  data.table(subset = s, Label = rownames(res$contrast),
             m = round(vapply(PAIRS, function(p) log2(DESIGN[s, p[1]] / DESIGN[s, p[2]]), 1), 2))
})), subset ~ Label, value.var = "m")
cat("\nMedian log2FC of spike-ins by subset, hand-written path:\n"); print(observed)
cat("\nDesigned log2 ratio by subset:\n"); print(designed)
cat("\nSpike-ins tested per subset:\n")
print(sp[Label == rownames(res$contrast)[1], .N, by = subset][order(subset)])
