# make_extdata.R — build the files the benchmark vignette reads
#
# Run from the repository root, where results/ holds the saved runs:
#   Rscript inst/scripts/make_extdata.R
#
# 1. inst/extdata/raw_trials_benchmark_20261007_163702.rds is a straight copy
#    of the saved trials, so the vignette needs no model and no API key.
# 2. inst/extdata/downstream_panels.rds is the per-protein adjusted p-values
#    behind Figure 4, cut down to what the six panels need: one row per
#    protein and contrast, both paths side by side. The full results are a
#    few hundred megabytes; this is a few hundred kilobytes.
# 3. inst/extdata/benchmark_columns.rds is the header of each deposited file,
#    which is what decides whether a mapped column exists. The files
#    themselves are not part of the installed package.

suppressPackageStartupMessages(library(data.table))

PANELS <- list(
  list(run = "downstream_puyvelde2022_reference_ProteinName_PG.ProteinNames",
       panel = "(a) Mnemonics instead of accessions",   data = "Spectronaut"),
  list(run = "downstream_solivais_llama3.1-8b_lean",
       panel = "(b) Unmodified peptide sequence",       data = "MetaMorpheus"),
  list(run = "downstream_puyvelde2022_reference_Intensity_F.NormalizedPeakArea",
       panel = "(c) Normalized peak area",              data = "Spectronaut"),
  list(run = "downstream_pxd005642_llama3.1-8b_constrained_v3",
       panel = "(d) Intensity instead of peak area",    data = "Proteome Discoverer"),
  list(run = "downstream_puyvelde2022_reference_nofilters",
       panel = "(e) No filtering or censoring",         data = "Spectronaut"),
  list(run = "downstream_pxd005642_deepseek-r1-14b_filter_aware",
       panel = "(f) Additional quality filters",        data = "Proteome Discoverer")
)

stopifnot(dir.exists("results"))
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

file.copy("results/raw_trials_benchmark_20261007_163702.rds",
          "inst/extdata/raw_trials_benchmark_20261007_163702.rds", overwrite = TRUE)

# The headers the models were shown, normalized exactly as load_test_dataset()
# normalizes them, so a predicted column can be looked up as the converter
# looks it up.
TOOLS <- c("spectronaut", "proteome_discoverer", "metamorpheus")
columns <- lapply(stats::setNames(TOOLS, TOOLS), function(tool) {
  names(MSstatsConvertLLM::load_test_dataset(tool, source = "benchmark"))
})
saveRDS(columns, "inst/extdata/benchmark_columns.rds")

# A full outer join keeps the proteins only one path tested, which the
# significant-protein counts are taken over; `in_h` and `in_l` say which path
# reported each row, so the disagreement counts can use the shared rows alone.
one_panel <- function(p) {
  f <- file.path("results", p$run, "full_results.rds")
  if (!file.exists(f)) stop("Missing ", f)
  r   <- readRDS(f)
  pos <- r$config$positive
  h   <- as.data.table(r$hand)[, .(Protein, Label, species,
                                   adjp_h = adj.pvalue, sig_h = sig, in_h = TRUE)]
  l   <- as.data.table(r$llm)[,  .(Protein, Label, species,
                                   adjp_l = adj.pvalue, sig_l = sig, in_l = TRUE)]
  m <- merge(h, l, by = c("Protein", "Label", "species"), all = TRUE)
  m[is.na(in_h), in_h := FALSE]
  m[is.na(in_l), in_l := FALSE]
  m[, `:=`(panel = p$panel, data = p$data, spiked = species == pos)]
  m[, .(panel, data, Protein, Label, spiked, in_h, adjp_h, sig_h, in_l, adjp_l, sig_l)]
}

panels <- rbindlist(lapply(PANELS, one_panel))
panels[, panel := factor(panel, levels = vapply(PANELS, `[[`, "", "panel"))]
saveRDS(panels, "inst/extdata/downstream_panels.rds", compress = "xz")

cat(sprintf("downstream_panels.rds: %d rows, %d panels, %.0f KB\n",
            nrow(panels), uniqueN(panels$panel),
            file.size("inst/extdata/downstream_panels.rds") / 1024))
