## Does the acquisition hint remove the cold/warm instability?
##
## The stochasticity study passed acquisition = "DDA" for Proteome Discoverer,
## which tells the model not to map FragmentIon at all. FragmentIon was the
## field that flipped in the original hand-run trials, which passed no
## acquisition argument. So the study may have pinned the variable it was
## measuring.
##
## This runs the same cell both ways. 10 reps each, ~15 minutes.

library(data.table)
library(httr2)

MODEL <- "deepseek-r1-14b"
TOOL  <- "proteome_discoverer"
N     <- 10

ollama_unload <- function(model_key) {
  request(paste0(OLLAMA_URL, "/api/generate")) |>
    req_body_json(list(model = MODEL_REGISTRY[[model_key]]$model, keep_alive = 0)) |>
    req_perform()
  Sys.sleep(2)
}

one_run <- function(acq) {
  t <- run_trial(MODEL, TOOL, "constrained", acquisition = acq)
  s <- score_trial(t)
  data.table(
    exact   = sum(s$exact_match),
    mapping = paste(s$field, ifelse(is.na(s$predicted), "null", s$predicted),
                    sep = "=", collapse = "|")
  )
}

res <- rbindlist(lapply(list(NULL, "DDA"), function(acq) {
  acq_label <- if (is.null(acq)) "none" else acq
  rbindlist(lapply(c("cold", "warm"), function(cond) {
    message("\n=== acquisition=", acq_label, " / ", cond, " ===")
    if (cond == "warm") try(one_run(acq), silent = TRUE)
    rbindlist(lapply(seq_len(N), function(i) {
      if (cond == "cold") ollama_unload(MODEL)
      r <- cbind(data.table(acquisition = acq_label, condition = cond, rep = i),
                 one_run(acq))
      message(sprintf("  %2d/%d  exact %d", i, N, r$exact))
      r
    }))
  }))
}))

fwrite(res, "results/acquisition_check.csv")

## Distinct mappings per cell. 1 means deterministic.
print(res[, .(n = .N, distinct_mappings = uniqueN(mapping),
              exact_min = min(exact), exact_max = max(exact)),
          by = .(acquisition, condition)])

## What FragmentIon actually resolved to in each cell.
print(res[, .(fragmention = paste(unique(sub(".*\\|FragmentIon=([^|]*)\\|.*", "\\1",
                                             paste0("|", mapping, "|"))),
                                  collapse = " / ")),
          by = .(acquisition, condition)])