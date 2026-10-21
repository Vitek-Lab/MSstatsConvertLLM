# ============================================================
# compare_puyvelde2022.R: downstream comparison on Puyvelde2022 (Spectronaut, DIA)
#
# Human, yeast and E. coli mixed in two conditions. Human is constant, so a
# significant human protein is a false positive; yeast and E. coli change
# between conditions and are true positives. Three biological samples per
# condition, each measured three times.
#
# Converter defaults, dataProcess with medians equalized (the setting used on
# all three datasets), one A against B contrast. The report is
# first cut on Explorer to the columns the converters read. The shared
# machinery is in downstream.R.
#
# Run from the repository root, in this order:
#   LLM_MODEL=reference N_CORES=2 caffeinate -i Rscript compare_puyvelde2022.R
#   N_CORES=2 caffeinate -i Rscript compare_puyvelde2022.R
#   LLM_MODEL=llama3.1-8b LLM_PROMPT=constrained_v3 N_CORES=2 caffeinate -i Rscript compare_puyvelde2022.R
#   LLM_MODEL=reference REF_SWAP=Intensity=F.NormalizedPeakArea N_CORES=2 caffeinate -i Rscript compare_puyvelde2022.R
#   LLM_MODEL=reference REF_SWAP=ProteinName=PG.ProteinNames N_CORES=2 caffeinate -i Rscript compare_puyvelde2022.R
#   LLM_MODEL=reference REF_FILTERS=none N_CORES=2 caffeinate -i Rscript compare_puyvelde2022.R
# ============================================================

source("downstream.R")

DATA_DIR <- Sys.getenv("PUYVELDE_DATA_DIR", "data/Puyvelde2022")
REPORT   <- file.path(DATA_DIR, "report_cut.tsv")
ANNOT    <- file.path(DATA_DIR, "annot.csv")

# Expected amount in A relative to B, from this benchmark's human, yeast and
# E. coli design (65/65, 30/15 and 5/20 percent). Checked against the data below.
DESIGNED_RATIO <- c(HUMAN = 1, YEAST = 2, ECOLI = 0.25)

say("Downstream comparison: Puyvelde2022 (Spectronaut, DIA), %d cores", N_CORES)

# The reference is the converter's defaults (intensity "PeakArea", peptide
# EG.ModifiedSequence, protein PG.ProteinGroups) plus its two built-in rules:
# blank fragments excluded from quantification, keep only "noloss" fragments.
chosen <- choose_mapping("spectronaut", c(
  ProteinName = "PG.ProteinGroups", PeptideSequence = "EG.ModifiedSequence",
  PrecursorCharge = "FG.Charge", FragmentIon = "F.FrgIon", ProductCharge = "F.Charge",
  Run = "R.FileName", Intensity = "F.PeakArea", Qvalue = "EG.Qvalue"),
  reference_filters = data.frame(
    column      = c("F.ExcludedFromQuantification", "F.FrgLossType"),
    dtype       = c("boolean", "string"),
    operation   = c("equals", "equals"),
    value       = c("FALSE", "noloss"),
    description = c("the converter's excludedFromQuantificationFilter",
                    "the converter keeps only noloss fragments"),
    confidence  = 1, stringsAsFactors = FALSE))

missing_inputs <- c(REPORT, ANNOT)[!file.exists(c(REPORT, ANNOT))]
if (length(missing_inputs)) stop("Missing input files:\n  ", paste(missing_inputs, collapse = "\n  "))

say("2. Reading %s", DATA_DIR)

# Species and a canonical protein key for every identifier a mapping can use.
# The key is PG.ProteinGroups; species comes from the UniProt mnemonics.
ids <- unique(fread(REPORT, select = c("PG.ProteinGroups", "PG.ProteinAccessions", "PG.ProteinNames"),
                    showProgress = FALSE))
name_species <- function(n) {
  vapply(strsplit(as.character(n), ";", fixed = TRUE), function(x) {
    s <- unique(sub("^.*_", "", trimws(x[nzchar(x)])))
    if (length(s) == 1 && s %in% names(DESIGNED_RATIO)) s else "OTHER"
  }, character(1))
}
ids[, species3 := name_species(PG.ProteinNames)]
key_map <- unique(rbind(ids[, .(id = PG.ProteinGroups,     key = PG.ProteinGroups)],
                        ids[, .(id = PG.ProteinAccessions, key = PG.ProteinGroups)],
                        ids[, .(id = PG.ProteinNames,      key = PG.ProteinGroups)]))
n_amb <- key_map[, .N, by = id][N > 1, .N]
if (n_amb) warning(n_amb, " identifiers map to more than one protein group; the first is used.")
key_map <- key_map[!duplicated(id)]
sp_map  <- ids[!duplicated(PG.ProteinGroups), .(key = PG.ProteinGroups, species3)]
say("  protein groups: %s", paste0(names(table(ids$species3)), " ", fmt(as.integer(table(ids$species3))), collapse = ", "))
rm(ids); invisible(gc())

protein_key <- function(p) {
  p <- as.character(p)
  k <- key_map$key[match(p, key_map$id)]
  fifelse(is.na(k), p, k)
}
species3_of <- function(p) {
  s <- sp_map$species3[match(protein_key(p), sp_map$key)]
  fifelse(is.na(s), "OTHER", s)
}
species_of <- function(p) {
  s <- species3_of(p)
  fcase(s == "HUMAN", "HUMAN", s %in% c("YEAST", "ECOLI"), "SPIKED", default = "OTHER")
}

res <- run_downstream(list(
  dataset = "Puyvelde2022", slug = "puyvelde2022",
  label = chosen$label, description = chosen$description, mapping = chosen$mapping,
  # Read here so the list holds the only copy, which the engine frees after converting
  raw = fread(REPORT, showProgress = FALSE),
  annotation = fread(ANNOT, showProgress = FALSE), run_col = "R.FileName",
  hand_convert = function(x, a) SpectronauttoMSstatsFormat(x, annotation = a, use_log_file = FALSE),
  species_of = species_of, protein_key = protein_key,
  positive = "SPIKED", negative = "HUMAN",
  dp_args = list(normalization = "equalizeMedians"),
  build_contrast = function(levs) {
    if (!setequal(levs, c("A", "B"))) stop("Unexpected conditions: ", paste(levs, collapse = ", "))
    matrix(c(1, -1)[match(levs, c("A", "B"))], nrow = 1, dimnames = list("A-B", levs))
  },
  designed_ratio = NULL,
  alpha = 0.05))

# Each species against its designed log2 ratio, on both paths. Tony's script
# records an FDR of 0.21 for this dataset, a check on the hand-written path.
by_species <- function(r, path) {
  d <- r[testable & is.finite(log2FC)]
  d[, species3 := species3_of(Protein)]
  d[species3 %in% names(DESIGNED_RATIO),
    .(path = path, proteins = .N, significant = sum(sig), median_log2fc = round(median(log2FC), 3)),
    by = species3]
}
sp <- rbind(by_species(res$hand, "hand"), by_species(res$llm, "LLM"))
sp[, designed_log2fc := round(log2(DESIGNED_RATIO[species3]), 3)]
cat("\nBy species, contrast A-B:\n"); print(sp[order(species3, path)])
