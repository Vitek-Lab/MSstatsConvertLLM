# make_new_format_data.R - build the files the new-format vignette reads
#
# Run from the repository root:
#   Rscript inst/scripts/make_new_format_data.R
#
# Downloads ProteoBench's AlphaPept example, cuts it under the 1 MB an R
# package should carry, and writes the annotation its run names imply. See
# inst/extdata/README for where the data comes from and under what licence.

suppressPackageStartupMessages(library(data.table))

URL <- paste0("https://raw.githubusercontent.com/Proteobench/ProteoBench/main/",
              "test/data/quant/quant_lfq_ion_DDA_QExactive/AlphaPept_subset.csv")
OUT <- "inst/extdata"
MIN_RUNS <- 4L   # keep a protein seen in at least this many of the six runs

stopifnot(dir.exists(OUT))

src <- file.path(tempdir(), "AlphaPept_subset_full.csv")
utils::download.file(URL, src, quiet = TRUE)
full <- fread(src)
stopifnot(nrow(full) > 0,
          identical(names(full), c("sequence_naked", "shortname", "protein",
                                   "sequence", "charge", "decoy",
                                   "ms1_int_sum_apex_dn")))

# Cutting by protein rather than by row keeps each surviving protein's design
# intact. The decoy rows come along whatever their protein, because they are
# what a filter-aware prompt has to find.
seen <- full[, .(runs = uniqueN(shortname)), by = protein][runs >= MIN_RUNS, protein]
subset <- full[protein %in% seen | decoy %in% TRUE]
fwrite(subset, file.path(OUT, "AlphaPept_subset.csv"))

# Two conditions, three runs each, read off the run names. They are the three
# technical replicates of the "alpha" series, and each is given its own
# BioReplicate so that MSstats compares condition A with condition B across
# them, as ProteoBench does. The labels carry the condition, because a
# BioReplicate appearing under both conditions would make MSstats read the
# design as paired.
runs <- sort(unique(subset$shortname))
condition <- sub("^.*_Condition_([AB])_.*$", "\\1", runs)
annotation <- data.table(
  Run = runs,
  Condition = condition,
  BioReplicate = paste(condition, sub("^.*_Sample_Alpha_0?(\\d+)$", "\\1", runs), sep = "_"))
stopifnot(nrow(annotation) == 6, uniqueN(annotation$Condition) == 2,
          annotation[, .N, by = Condition]$N == c(3, 3))
fwrite(annotation, file.path(OUT, "AlphaPept_annotation.csv"))

cat(sprintf("AlphaPept_subset.csv: %d of %d rows, %d proteins, %d decoy rows, %.0f KB\n",
            nrow(subset), nrow(full), uniqueN(subset$protein),
            sum(subset$decoy %in% TRUE),
            file.size(file.path(OUT, "AlphaPept_subset.csv")) / 1024))
print(annotation)
