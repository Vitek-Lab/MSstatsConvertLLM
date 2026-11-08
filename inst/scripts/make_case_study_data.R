# make_case_study_data.R - build the files the case-study vignettes read
#
# Run from the repository root, after inst/scripts/make_benchmark_data.R:
#
#   Rscript inst/scripts/make_case_study_data.R
#
# Each case study is one list in inst/extdata:
#   rows           rows of the deposited export, with the column names the tool
#                  wrote; the first three are the rows the model was shown
#   specification  the conversion specification of the saved trial, repeat 1,
#                  which is the one the downstream comparison used
#   results        the Figure 4 panel of the case: one row per protein and
#                  comparison, the hand-written converter in the _h columns and
#                  the LLM-based converter in the _l columns
# The specification and the results are cut from the saved trials and
# downstream_panels.rds, so nothing is recomputed. The exports are the
# truncated deposited files in benchmark_data/; point MSSTATS_BENCHMARK_DATA
# elsewhere if they are not there.

suppressPackageStartupMessages(library(data.table))

TRIALS  <- "inst/extdata/raw_trials_benchmark_20261007_163702.rds"
PANELS  <- "inst/extdata/downstream_panels.rds"
EXPORTS <- Sys.getenv("MSSTATS_BENCHMARK_DATA", "benchmark_data")

stopifnot(file.exists(TRIALS), file.exists(PANELS))
trials <- readRDS(TRIALS)
panels <- readRDS(PANELS)

saved_specification <- function(tool, model, prompt) {
  t <- Filter(function(x) x$tool == tool && x$model == model &&
                x$prompt == prompt && x$rep == 1L, trials)
  stopifnot(length(t) == 1)
  t[[1]]$mapping
}

panel_results <- function(label) {
  r <- panels[panel == label, .(Protein, Label, spiked, in_h, adjp_h, sig_h,
                                 in_l, adjp_l, sig_l)]
  stopifnot(nrow(r) > 0)
  r
}

read_export <- function(file) {
  path <- file.path(EXPORTS, file)
  if (!file.exists(path)) stop("Missing ", path)
  fread(path, showProgress = FALSE)
}

chosen <- function(specification, field) {
  m <- specification$mappings
  m$from[m$field == field]
}

# Proteome Discoverer: the three rows the model was shown, which hold both
# candidate intensity columns.
pd <- read_export("PXD005642_PD_input.csv")
pd_case <- list(
  rows          = pd[1:3],
  specification = saved_specification("proteome_discoverer", "llama3.1-8b",
                                      "constrained_v3"),
  results       = panel_results("(d) Intensity instead of peak area"))
stopifnot(chosen(pd_case$specification, "Intensity") == "Intensity",
          all(c("Intensity", "Precursor Area") %in% names(pd_case$rows)))

# Spectronaut: the fragments of the first precursor in the export. The first
# three are the rows the model was shown; the others include a fragment with
# a neutral loss and fragments excluded from quantification.
sn <- read_export("Puyvelde2022_Spectronaut_Report-001.tsv")
sn_rows <- sn[R.FileName == R.FileName[1] & EG.PrecursorId == EG.PrecursorId[1]]
stopifnot(identical(sn_rows[1:3], sn[1:3]),
          any(sn_rows$F.FrgLossType != "noloss"),
          any(sn_rows$F.ExcludedFromQuantification))
sn_case <- list(
  rows          = sn_rows,
  specification = saved_specification("spectronaut", "llama3.1-8b",
                                      "constrained_v3"),
  results       = panel_results("(e) No filtering or censoring"))
stopifnot(NROW(sn_case$specification$filters) == 0)

# MetaMorpheus: the three rows the model was shown, where the base and the
# full sequence are the same, and the first peptide the export reports both
# with and without a modification, in one run, at one charge and in one
# protein group.
mm <- read_export("Solivais2024_QuantifiedPeaks.tsv")
FORM_KEY <- c("Base Sequence", "Precursor Charge", "File Name", "Protein Group")
forms <- mm[, .(n = uniqueN(`Full Sequence`)), by = FORM_KEY][n > 1][1]
mm_case <- list(
  rows          = mm[1:3],
  forms         = mm[forms, on = FORM_KEY],
  specification = saved_specification("metamorpheus", "llama3.1-8b", "lean"),
  results       = panel_results("(b) Unmodified peptide sequence"))
mm_case$forms[, n := NULL]
stopifnot(with(mm_case$rows, all(`Base Sequence` == `Full Sequence`)),
          uniqueN(mm_case$forms$`Full Sequence`) > 1,
          chosen(mm_case$specification, "PeptideSequence") == "Base.Sequence")

for (case in list(list(pd_case, "pd_case_study.rds"),
                  list(sn_case, "spectronaut_case_study.rds"),
                  list(mm_case, "metamorpheus_case_study.rds"))) {
  out <- file.path("inst/extdata", case[[2]])
  saveRDS(case[[1]], out, compress = "xz")
  cat(sprintf("%s: %d rows of the export, %d result rows, %.0f KB\n", out,
              nrow(case[[1]]$rows), nrow(case[[1]]$results),
              file.size(out) / 1024))
}
