# figure4_downstream.R
# Figure 4 (downstream comparison) and the replicate-CV check, built from the
# saved downstream results. Nothing is recomputed: every number comes from
# results/<run>/full_results.rds and the conversion caches.
#
# Run from the repo root:  Rscript figure4_downstream.R
# Writes results/Figure4_downstream.pdf and .png, and prints the counts the
# caption needs. The replicate-CV check is in diagnose_pd_intensity.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# One panel per mapping choice, in the order the Results text discusses them.
panels <- list(
  list(run = "downstream_puyvelde2022_reference_ProteinName_PG.ProteinNames",
       choice = "(a) Mnemonics instead of accessions",
       data = "Spectronaut"),
  list(run = "downstream_solivais_llama3.1-8b_lean",
       choice = "(b) Unmodified peptide sequence",
       data = "MetaMorpheus"),
  list(run = "downstream_puyvelde2022_reference_Intensity_F.NormalizedPeakArea",
       choice = "(c) Normalized peak area",
       data = "Spectronaut"),
  list(run = "downstream_pxd005642_llama3.1-8b_constrained_v3",
       choice = "(d) Intensity instead of peak area",
       data = "Proteome Discoverer"),
  list(run = "downstream_puyvelde2022_reference_nofilters",
       choice = "(e) No filtering or censoring",
       data = "Spectronaut"),
  list(run = "downstream_pxd005642_deepseek-r1-14b_filter_aware",
       choice = "(f) Additional quality filters",
       data = "Proteome Discoverer")
)
CAP   <- 20      # cap -log10(adjusted p) for display only; stated in the caption
ALPHA <- 0.05
NEAR  <- 0.20    # a disagreeing protein is "near the threshold" if the path that
                 # missed it still has an adjusted p-value at or below this

points <- list(); counts <- list()
for (p in panels) {
  f <- file.path("results", p$run, "full_results.rds")
  if (!file.exists(f)) stop("Missing ", f, ". This run's downstream comparison has not been saved.")
  r   <- readRDS(f)
  pos <- r$config$positive
  h   <- as.data.table(r$hand)
  l   <- as.data.table(r$llm)

  m <- merge(h[, .(Protein, Label, adjp_h = adj.pvalue, sig_h = sig, species)],
             l[, .(Protein, Label, adjp_l = adj.pvalue, sig_l = sig)],
             by = c("Protein", "Label"))
  m[, truth := fifelse(species == pos, "changes (spiked)", "does not change")]
  m[, x := pmin(-log10(pmax(adjp_h, 1e-300)), CAP)]
  m[, y := pmin(-log10(pmax(adjp_l, 1e-300)), CAP)]
  m[, panel := p$choice]
  points[[p$run]] <- m[!is.na(x) & !is.na(y)]

  n_sig <- function(d) sum(d$sig, na.rm = TRUE)
  n_tp  <- function(d) sum(d$sig & d$species == pos, na.rm = TRUE)
  disc  <- m[xor(sig_h %in% TRUE, sig_l %in% TRUE)]
  near  <- disc[pmax(adjp_h, adjp_l, na.rm = TRUE) <= NEAR]
  counts[[p$run]] <- data.table(
    panel = p$choice, data = p$data, run = p$run,
    conv_sig = n_sig(h), conv_true = n_tp(h), conv_false = n_sig(h) - n_tp(h),
    llm_sig  = n_sig(l), llm_true  = n_tp(l), llm_false  = n_sig(l) - n_tp(l),
    disagree = nrow(disc), disagree_near = nrow(near))
}

d   <- rbindlist(points)
ann <- rbindlist(counts)
# Titles wrap to fit a third of the page width; the tool gets its own line
ann[, strip := sprintf("%s\n%s", vapply(strwrap(panel, 40, simplify = FALSE), paste, "", collapse = "\n"), data)]
d <- merge(d, ann[, .(panel, strip)], by = "panel")
d[, strip := factor(strip, levels = ann$strip)]

thr <- -log10(ALPHA)
g <- ggplot(d, aes(x, y, colour = truth)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = thr, linetype = "dotted", colour = "grey65") +
  geom_vline(xintercept = thr, linetype = "dotted", colour = "grey65") +
  # Points only are rasterised, so the PDF opens quickly; text and lines stay
  # vector. ragg, because ggrastr's default Cairo device needs XQuartz.
  ggrastr::rasterise(geom_point(alpha = 0.35, size = 0.6), dpi = 600, dev = "ragg") +
  facet_wrap(~ strip, ncol = 3, scales = "free") +
  scale_colour_manual(values = c("changes (spiked)" = "#B03A2E", "does not change" = "#8C9AA3"),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "-log10 adjusted p-value, hand-written converter",
       y = "-log10 adjusted p-value, LLM-based converter") +
  theme_minimal(base_size = 8.5) +
  theme(legend.position = "bottom",
        strip.text = element_text(hjust = 0, face = "bold", lineheight = 1.05),
        panel.grid.minor = element_blank())

ggsave("results/Figure4_downstream.pdf", g, width = 7.2, height = 6.8)
ggsave("results/Figure4_downstream.png", g, width = 7.2, height = 6.8, dpi = 200)

cat("\n=== Figure 4 counts (paste these back) ===\n")
print(ann[, .(run, conv_sig, conv_true, conv_false, llm_sig, llm_true, llm_false, disagree, disagree_near)])
cat(sprintf("Display cap: -log10 adjusted p clipped at %d. Near threshold: missed path's adjusted p <= %.2f.\n", CAP, NEAR))

cat("\nFigure written to results/Figure4_downstream.pdf and .png\n")