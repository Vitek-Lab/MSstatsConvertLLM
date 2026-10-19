# diagnose_pd_intensity.R
# Why did Intensity (llama3.1-8b mapping) lose most of the power on PXD005642,
# where the converter reads Precursor Area? Significance rests on t = log2FC / SE,
# so compare both parts for the spiked proteins, plus how many features each
# path measured. Reads saved results only.
#
# Run from the repo root:  Rscript diagnose_pd_intensity.R

suppressPackageStartupMessages(library(data.table))

r   <- readRDS("results/downstream_pxd005642_llama3.1-8b_constrained_v3/full_results.rds")
pos <- r$config$positive
h   <- as.data.table(r$hand)   # Precursor Area, hand-written converter
l   <- as.data.table(r$llm)    # Intensity, LLM path

m <- merge(h[species == pos & is.finite(log2FC), .(Protein, Label, fc_area = log2FC, se_area = SE)],
           l[species == pos & is.finite(log2FC), .(Protein, Label, fc_int  = log2FC, se_int  = SE)],
           by = c("Protein", "Label"))
med <- function(v) median(v, na.rm = TRUE)

cat("=== Spiked proteins, tested in both paths ===\n")
cat(sprintf("Protein-comparison pairs: %d\n", nrow(m)))
cat(sprintf("Median |log2FC|: Precursor Area %.2f, Intensity %.2f\n", med(abs(m$fc_area)), med(abs(m$fc_int))))
cat(sprintf("Median SE:       Precursor Area %.2f, Intensity %.2f\n", med(m$se_area), med(m$se_int)))
cat(sprintf("Median |t|:      Precursor Area %.2f, Intensity %.2f\n",
            med(abs(m$fc_area / m$se_area)), med(abs(m$fc_int / m$se_int))))
fit <- lm(fc_int ~ fc_area, data = m)
cat(sprintf("Intensity log2FC regressed on Precursor Area log2FC: slope %.2f, r = %.2f\n",
            coef(fit)[2], cor(m$fc_area, m$fc_int)))

cat("\n=== By comparison ===\n")
print(m[, .(pairs = .N,
            fc_area = round(med(fc_area), 2), fc_int = round(med(fc_int), 2),
            se_area = round(med(se_area), 2), se_int = round(med(se_int), 2)), by = Label])

cat("\n=== Proteins tested in each path (all species) ===\n")
cat(sprintf("Precursor Area: %d proteins with a finite log2FC; Intensity: %d\n",
            uniqueN(h[is.finite(log2FC), Protein]), uniqueN(l[is.finite(log2FC), Protein])))

# Feature counts, completeness and replicate CV from the conversion caches.
# The converter and the LLM path may name the peptide column differently, so it
# is found by name, and features are also paired across the two paths on the
# plain peptide sequence and charge, so the CVs compare the same peptides.
pick_table <- function(x) {
  if (is.data.frame(x)) return(as.data.table(x))
  if (is.list(x)) for (e in x) if (is.data.frame(e) && "Intensity" %in% names(e)) return(as.data.table(e))
  stop("No converted table with an Intensity column found.")
}
peptide_col <- function(dt) {
  hit <- intersect(c("PeptideSequence", "PeptideModifiedSequence"), names(dt))
  if (length(hit) == 0) stop("No peptide column among: ", paste(names(dt), collapse = ", "))
  hit[1]
}
plain_seq <- function(x) {
  x <- sub("^\\[[A-Za-z-]*\\]\\.", "", x)          # leading flank, e.g. "[K]."
  x <- sub("\\.\\[[A-Za-z-]*\\]$", "", x)          # trailing flank, e.g. ".[S]"
  x <- sub("_.*$", "", x)                          # modifications joined after an underscore
  x <- gsub("\\([^)]*\\)|\\[[^]]*\\]", "", x)      # bracketed modification text
  toupper(gsub("[^A-Za-z]", "", x))
}
describe <- function(f, label) {
  dt <- pick_table(readRDS(f))
  pc <- peptide_col(dt)
  cat(sprintf("%s columns: %s\n  peptide column used: %s, e.g. %s\n", label,
              paste(names(dt), collapse = ", "), pc, dt[[pc]][1]))
  dt[, feature  := paste(get(pc), PrecursorCharge, sep = "|")]
  dt[, pair_key := paste(plain_seq(get(pc)), PrecursorCharge, sep = "|")]
  obs <- dt[!is.na(Intensity) & Intensity > 0]
  cat(sprintf("  features: %d | runs: %d | rows with a positive intensity: %.1f%%\n",
              uniqueN(dt$feature), uniqueN(dt$Run), 100 * nrow(obs) / nrow(dt)))
  obs
}
cv_by <- function(obs, key) obs[, .(n = .N, cv = sd(Intensity) / mean(Intensity)), by = c(key, "Condition")][n >= 2]
q <- function(v) sprintf("%.1f%% (IQR %.1f to %.1f%%)", 100 * median(v), 100 * quantile(v, 0.25), 100 * quantile(v, 0.75))

cat("\n=== Converted data ===\n")
area <- describe("results/downstream_pxd005642_hand/cache/hand_converted.rds", "Precursor Area")
int  <- describe("results/downstream_pxd005642_llama3.1-8b_constrained_v3/cache/llm_converted.rds", "Intensity")

cat("\n=== Replicate CV within condition, each path's own features ===\n")
cva <- cv_by(area, "feature"); cvi <- cv_by(int, "feature")
cat(sprintf("Precursor Area: median %s over %d feature-condition pairs\n", q(cva$cv), nrow(cva)))
cat(sprintf("Intensity:      median %s over %d feature-condition pairs\n", q(cvi$cv), nrow(cvi)))

cat("\n=== Same peptides and charges in both paths (paired) ===\n")
pa <- cv_by(area, "pair_key"); pi <- cv_by(int, "pair_key")
pr <- merge(pa, pi, by = c("pair_key", "Condition"), suffixes = c("_area", "_int"))
cat(sprintf("Peptide-charge keys: %d in Precursor Area, %d in Intensity, %d shared\n",
            uniqueN(area$pair_key), uniqueN(int$pair_key), length(intersect(area$pair_key, int$pair_key))))
if (nrow(pr) > 0) {
  cat(sprintf("Paired condition-level CVs: %d | median Precursor Area %.1f%%, Intensity %.1f%% | Intensity higher in %.0f%%\n",
              nrow(pr), 100 * median(pr$cv_area), 100 * median(pr$cv_int), 100 * mean(pr$cv_int > pr$cv_area)))
}