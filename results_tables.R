## Consolidates the local and hosted arms into one table, then prints the two
## views worth putting on screen. source() this, don't paste it.

library(data.table)

local_arm  <- fread("results/baseline_final_h200.csv")
hosted_arm <- rbind(
  # Spectronaut here predates the PG.Qvalue correction, so take it from the rerun
  fread("results/field_scores_20260927_161013.csv")[tool != "spectronaut"],
  fread("results/field_scores_20260927_164703.csv")
)

all_arms <- rbind(local_arm, hosted_arm, fill = TRUE)
fwrite(all_arms, "results/all_arms.csv")

model_order  <- c("llama3.2-3b", "llama3.1-8b", "deepseek-r1-14b", "claude-sonnet")
prompt_order <- c("lean", "constrained", "filter_aware", "constrained_filter",
                  "constrained_v2", "constrained_filter_v2",
                  "constrained_v3", "constrained_filter_v3")

cell <- all_arms[, .(acceptable = mean(acceptable_match),
                     exact      = mean(exact_match)),
                 by = .(model, tool, prompt)]
cell[, model  := factor(model,  levels = model_order)]
cell[, prompt := factor(prompt, levels = prompt_order)]

## View 1 — how hard each format is, and how the models rank.
cat("\n=== Acceptable accuracy by model and format (mean over 8 prompts) ===\n")
print(dcast(cell, model ~ tool, value.var = "acceptable",
            fun.aggregate = function(x) round(mean(x), 3)))

## View 2 — the paper's main claim. Local models climb across the prompt
## variants; the hosted model is flat and already at the ceiling.
cat("\n=== Acceptable accuracy by model and prompt (mean over 3 formats) ===\n")
print(dcast(cell, prompt ~ model, value.var = "acceptable",
            fun.aggregate = function(x) round(mean(x), 3)))

## View 3 — the recommended configuration, per format.
cat("\n=== constrained_v3, per format ===\n")
print(dcast(cell[prompt == "constrained_v3"], model ~ tool,
            value.var = "acceptable", fun.aggregate = function(x) round(mean(x), 3)))

cat("\nWritten: results/all_arms.csv\n")
