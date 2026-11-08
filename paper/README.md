# Reproducing the benchmark of the manuscript

This directory is not part of the package. It shows how to reproduce Table 2
and the downstream figure of the manuscript from the saved results in
`inst/extdata`, and how to rerun the benchmark. Run the code from the root of
the repository, with the package installed or loaded with
`devtools::load_all()`. The saved results need no language model, API key or
network connection.

```r
library(MSstatsConvertLLM)
library(data.table)
```

## What was compared

The benchmark asked four models, three local models served with Ollama and one
hosted model, to map the columns of the output of three spectral processing
tools, Spectronaut, Proteome Discoverer and MetaMorpheus, to the MSstats format.
Each combination of model and tool was run with eight prompt variants and
repeated three times, which gives 288 trials. The prompt variants are lean,
constrained, filter_aware and constrained_filter, and two revisions (v2 and v3)
of the two constrained variants.

```r
trials <- readRDS("inst/extdata/raw_trials_benchmark_20261007_163702.rds")
length(trials)

table(vapply(trials, `[[`, "", "model"), vapply(trials, `[[`, "", "tool"))
names(PROMPT_VERSIONS)
```

Each trial stores the output of `run_benchmark_trial()`, including the
conversion specification returned by the model. `score_trial()` compares the
mappings with the gold standard, which was read from the source code of the
hand-written converter. A strict match requires the default column of the
converter, and an acceptable match allows any column that the converter
accepts.

## Table 2

Table 2 of the manuscript reports, for each model and prompt, the proportion of
correct decisions over the five essential fields (ProteinName, PeptideSequence,
PrecursorCharge, Run and Intensity) among the trials whose conversion ran. A
conversion ran when each of the five fields was mapped to a column that exists
in the output. A decision is correct when the chosen column is one that the
hand-written converter accepts. The last row gives the proportion of trials
whose conversion ran.

```r
ESSENTIAL <- c("ProteinName", "PeptideSequence", "PrecursorCharge",
               "Run", "Intensity")
MODELS <- c("llama3.2-3b", "llama3.1-8b", "deepseek-r1-14b", "claude-sonnet")

# The header of each tool's output, which decides whether a chosen column exists.
columns <- readRDS("inst/extdata/benchmark_columns.rds")

scored <- rbindlist(lapply(trials, function(tr) {
  s <- score_trial(tr)[field %in% ESSENTIAL]
  s[, .(model = tr$model, tool = tr$tool, prompt = tr$prompt, rep = tr$rep,
        field, predicted, acceptable_match)]
}))

# A conversion runs when each of the five essential fields is mapped to a column
# that the tool's output has. A field left unmapped, or mapped to a column that
# is not there, leaves the converter with nothing to read.
scored[, mapped := !is.na(predicted) & predicted %in% columns[[tool[1]]], by = tool]
scored[, ran := all(mapped), by = .(model, tool, prompt, rep)]

accuracy <- dcast(scored[ran == TRUE, .(value = mean(acceptable_match)),
                         by = .(model, prompt)],
                  prompt ~ model, value.var = "value")
ran_share <- dcast(unique(scored[, .(model, tool, prompt, rep, ran)])[
                     , .(prompt = "conversion ran", value = mean(ran)), by = model],
                   prompt ~ model, value.var = "value")

table2 <- rbind(accuracy[match(names(PROMPT_VERSIONS), prompt)], ran_share)
setcolorder(table2, c("prompt", MODELS))
table2[, (MODELS) := lapply(.SD, round, 3), .SDcols = MODELS]
table2
```

The following code compares these values with the values printed in the
manuscript, and lists any value that differs.

```r
paper <- data.table(
  prompt            = c(names(PROMPT_VERSIONS), "conversion ran"),
  `llama3.2-3b`     = c(0.800, 0.667, 0.700, 0.667, 0.800, 0.800, 0.800, 0.600, 0.833),
  `llama3.1-8b`     = c(0.867, 0.800, 0.867, 0.800, 0.933, 0.933, 0.933, 0.933, 1.000),
  `deepseek-r1-14b` = c(0.867, 0.800, 1.000, 0.800, 0.933, 1.000, 0.933, 0.933, 0.958),
  `claude-sonnet`   = c(1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000))

long <- function(x) melt(x, id.vars = "prompt", variable.name = "model")
differences <- merge(long(table2), long(paper), by = c("prompt", "model"),
                     suffixes = c("_here", "_paper"))[value_here != value_paper]

if (nrow(differences) == 0) "every cell matches the paper" else differences
```

## The downstream figure

The manuscript compares the differential abundance results of the hand-written
converter and the LLM-based converter for six choices that the models made in
the conversion specification. The saved table holds, for each choice, one row
per protein and comparison, with the adjusted p-value from each converter. The
counts below are the significant proteins, the spiked proteins among them, the
proteins whose significance changed, and how many of those were near the
threshold, i.e., had an adjusted p-value of 0.2 or less with the converter that
did not call them. They are the counts reported in the manuscript. The
vignettes "2. Checking the intensity column", "3. Checking the filters" and
"4. Checking the peptide sequence" walk through three of the six choices.

```r
panels <- readRDS("inst/extdata/downstream_panels.rds")
ALPHA <- 0.05   # significance
CAP   <- 20     # display cap on -log10 adjusted p
NEAR  <- 0.20   # a disagreement is near the threshold below this adjusted p

counts <- panels[, {
  shared  <- .SD[in_h & in_l]
  changed <- shared[xor(sig_h %in% TRUE, sig_l %in% TRUE)]
  .(hand        = sum(sig_h, na.rm = TRUE),
    hand_spiked = sum(sig_h & spiked, na.rm = TRUE),
    llm         = sum(sig_l, na.rm = TRUE),
    llm_spiked  = sum(sig_l & spiked, na.rm = TRUE),
    changed     = nrow(changed),
    near        = nrow(changed[pmax(adjp_h, adjp_l, na.rm = TRUE) <= NEAR]))
}, by = .(panel, data)]

counts
```

The figure redraws the six panels. Each point is one protein in one comparison,
and points on the dashed line have the same adjusted p-value with both
converters.

```r
library(ggplot2)

strips <- unique(panels[, .(panel, data)])[order(panel)]
strips[, strip := paste(panel, data, sep = "\n")]

points <- merge(panels[in_h & in_l & !is.na(adjp_h) & !is.na(adjp_l)],
                strips, by = c("panel", "data"))
points[, `:=`(x     = pmin(-log10(pmax(adjp_h, 1e-300)), CAP),
              y     = pmin(-log10(pmax(adjp_l, 1e-300)), CAP),
              truth = fifelse(spiked, "changes (spiked)", "does not change"),
              strip = factor(strip, levels = strips$strip))]

ggplot(points, aes(x, y, colour = truth)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dotted", colour = "grey65") +
  geom_vline(xintercept = -log10(ALPHA), linetype = "dotted", colour = "grey65") +
  geom_point(alpha = 0.35, size = 0.6) +
  facet_wrap(~ strip, ncol = 3, scales = "free") +
  scale_colour_manual(values = c("changes (spiked)" = "#B03A2E",
                                 "does not change"  = "#8C9AA3"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "-log10 adjusted p-value, hand-written converter",
       y = "-log10 adjusted p-value, LLM-based converter") +
  theme_minimal(base_size = 8.5) +
  theme(legend.position = "bottom",
        strip.text = element_text(hjust = 0, face = "bold", lineheight = 1.05),
        panel.grid.minor = element_blank())
```

## Rerunning the benchmark

The results above were produced with the settings and scripts below. Rerunning
them requires Ollama with the three local models, an Anthropic API key in
`ANTHROPIC_API_KEY` for the hosted model, and, for the comparison scripts, the
full exports of the three datasets in `data/`, which are not in the repository.
The local models were run at temperature 0 with a fixed seed. Ollama reuses the
cached beginning of the previous request, so the first call after a change of
prompt can differ from the following calls, and it is discarded.

The filter prompts in the package now state that each filter describes the rows
to keep. To rerun the benchmark exactly as in the manuscript, use the prompts of
the tag benchmark-2026-10.

```r
# create_chat() passes these settings to ellmer. MSSTATS_OLLAMA_URL is read on
# every use, so it can be set after the package is attached.
Sys.setenv(MSSTATS_OLLAMA_URL = "http://localhost:11434")
ellmer::params(format = "json", temperature = 0, top_p = 1, seed = 1,
               num_ctx = 4096, num_predict = -1)
```

```r
TOOLS <- c("spectronaut", "proteome_discoverer", "metamorpheus")

for (model in names(MODEL_REGISTRY)) {
  for (tool in TOOLS) {
    for (prompt in names(PROMPT_VERSIONS)) {
      args <- list(model, tool, prompt,
                   acquisition = TOOL_ACQUISITION[[tool]],
                   allow_transforms = isTRUE(TOOL_TRANSFORMS[[tool]]),
                   source = "benchmark")
      do.call(run_benchmark_trial, args)  # discarded: first call after a change
      for (rep in 1:3) {
        trial <- do.call(run_benchmark_trial, args)
        print_scorecard(score_trial(trial), trial)
      }
    }
  }
}
```

`run_all.R` runs this loop and saves every trial. The differential abundance
results come from six runs of the three comparison scripts, and
`inst/scripts/make_benchmark_data.R` reduces their output to the tables used
above. `inst/scripts/make_case_study_data.R` then cuts the files that the
case-study vignettes read from these tables. The commands are run from the
root of the repository.

```sh
Rscript run_all.R --source benchmark

LLM_MODEL=reference REF_SWAP=ProteinName=PG.ProteinNames Rscript compare_puyvelde2022.R
LLM_MODEL=reference REF_SWAP=Intensity=F.NormalizedPeakArea Rscript compare_puyvelde2022.R
LLM_MODEL=reference REF_FILTERS=none Rscript compare_puyvelde2022.R
LLM_MODEL=llama3.1-8b LLM_PROMPT=constrained_v3 Rscript compare_pxd005642.R
LLM_MODEL=deepseek-r1-14b LLM_PROMPT=filter_aware Rscript compare_pxd005642.R
LLM_MODEL=llama3.1-8b LLM_PROMPT=lean Rscript compare_solivais.R

Rscript inst/scripts/make_benchmark_data.R
Rscript inst/scripts/make_case_study_data.R
```
