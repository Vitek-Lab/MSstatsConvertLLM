## Stochasticity study
##
## Question: how much does a single trial result vary across repeated runs, and
## does that variation depend on the model's load state?
##
## Three conditions:
##   cold   - model unloaded from VRAM before every trial
##   warm   - model resident, identical prompt repeated (prompt cache hot)
##   reset  - model resident, unrelated prompt sent first (prompt cache busted)
##
## "reset" separates two explanations for the cold/warm split seen earlier.
## If reset matches cold, the effect is prompt caching, and every sweep cell
## depends on whichever trial ran before it. If reset matches warm, the effect
## is VRAM load state, and one warmup per model is enough.
##
## Writes one row per run, appending after each, so a dropped tunnel loses
## nothing. Re-running the script resumes where it stopped.

library(data.table)
library(httr2)

## ---------------------------------------------------------------- settings

N_REPS      <- 50
PROMPT_KEY  <- "constrained"
OUT_CSV     <- "results/stochasticity.csv"
CONDITIONS  <- c("cold", "warm", "reset")

## acquisition is left to TOOL_ACQUISITION, which derives it from the tool name.
COMBOS <- list(
  list(model = "llama3.2-3b",     tool = "spectronaut"),
  list(model = "deepseek-r1-14b", tool = "proteome_discoverer")
)

OLLAMA <- Sys.getenv("MSSTATS_OLLAMA_URL", "http://localhost:11500")

## ---------------------------------------------------------------- guard
##
## The unload and scramble calls below go straight to the tunnel, but run_trial
## goes through ellmer. On a branch without the configurable base URL, ellmer
## defaults to localhost:11434 and every trial fails against a port with nothing
## on it. Fail here instead of 150 runs later.

if (!exists("OLLAMA_URL")) {
  stop("OLLAMA_URL not found. Check out a branch that has the configurable ",
       "Ollama base URL (PR #27) before running this.")
}
if (!nzchar(Sys.getenv("MSSTATS_OLLAMA_URL"))) {
  stop("MSSTATS_OLLAMA_URL is unset. Set it, then devtools::load_all() again ",
       "so OLLAMA_URL picks it up.")
}
## OLLAMA_URL resolves when the package loads. Setting the env var after
## load_all() leaves it on the 11434 default while the tunnel is on 11500,
## and every call fails against a port with nothing on it.
if (!identical(OLLAMA_URL, Sys.getenv("MSSTATS_OLLAMA_URL"))) {
  stop("OLLAMA_URL is '", OLLAMA_URL, "' but MSSTATS_OLLAMA_URL is '",
       Sys.getenv("MSSTATS_OLLAMA_URL"), "'. Set the env var, then ",
       "devtools::load_all() again.")
}
OLLAMA <- OLLAMA_URL

## ---------------------------------------------------------------- helpers

ollama_tag <- function(model_key) {
  entry <- MODEL_REGISTRY[[model_key]]
  if (is.null(entry)) stop("no MODEL_REGISTRY entry for '", model_key, "'")
  entry$model
}

## Drop the model from VRAM. Next call pays a full load.
ollama_unload <- function(model_key) {
  request(paste0(OLLAMA, "/api/generate")) |>
    req_body_json(list(model = ollama_tag(model_key), keep_alive = 0)) |>
    req_perform()
  Sys.sleep(2)
  invisible(NULL)
}

## Keep the model resident but overwrite the cached prompt prefix. The filler
## has to differ from token zero, which is why it is not another mapping prompt:
## all the real prompts share a task-and-schema prefix that would stay cached.
FILLER <- paste(rep(
  "The harbour master kept a ledger of every vessel that entered the bay,
   noting the tide, the weather, and the cargo, in a hand that grew less
   legible with each passing winter.", 8), collapse = " ")

ollama_scramble <- function(model_key) {
  request(paste0(OLLAMA, "/api/generate")) |>
    req_body_json(list(
      model      = ollama_tag(model_key),
      prompt     = FILLER,
      stream     = FALSE,
      keep_alive = "10m",
      options    = list(num_predict = 1, temperature = 0)
    )) |>
    req_perform()
  invisible(NULL)
}

## Mirrors the call in run_all.R exactly. If the study ran a different
## configuration from the sweep, its conclusions would not transfer.
call_trial <- function(model, tool, prompt) {
  run_trial(model, tool, prompt,
            acquisition      = TOOL_ACQUISITION[[tool]],
            allow_transforms = isTRUE(TOOL_TRANSFORMS[[tool]]))
}

## One row per run. The mapping string is the point: scores can stay identical
## while the underlying prediction changes, and that would be invisible in a
## column of totals.
summarise_run <- function(trial, sc) {
  has <- function(col) col %in% names(sc)
  data.table(
    exact       = sum(sc$exact_match),
    acceptable  = sum(sc$acceptable_match),
    candidates  = if (has("candidate_hit")) sum(sc$candidate_hit, na.rm = TRUE) else NA_integer_,
    n_fields    = nrow(sc),
    elapsed_sec = trial$elapsed_sec,
    mapping     = paste(sc$field, ifelse(is.na(sc$predicted), "null", sc$predicted),
                        sep = "=", collapse = "|")
  )
}

## ---------------------------------------------------------------- resume

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

## Keyed as strings rather than a data.table filter. The column names and the
## lookup values are the same words, and comparing a column to a variable of
## the same name silently matches every row.

run_key <- function(model, tool, condition, rep) {
  paste(model, tool, condition, rep, sep = "|")
}

done_keys <- if (file.exists(OUT_CSV)) {
  d <- fread(OUT_CSV)
  run_key(d$model, d$tool, d$condition, d$rep)
} else {
  character(0)
}

## ---------------------------------------------------------------- run

for (combo in COMBOS) {
  for (condition in CONDITIONS) {

    message("\n=== ", combo$model, " / ", combo$tool, " / ", condition, " ===")

    ## warm and reset both need the model resident before the block starts.
    ## That first load is a cold call, so it is run and discarded.
    if (condition %in% c("warm", "reset")) {
      message("  warming up")
      try(call_trial(combo$model, combo$tool, PROMPT_KEY), silent = TRUE)
    }

    for (rep in seq_len(N_REPS)) {

      if (run_key(combo$model, combo$tool, condition, rep) %in% done_keys) next

      if (condition == "cold")  ollama_unload(combo$model)
      if (condition == "reset") ollama_scramble(combo$model)

      row <- tryCatch({
        trial <- call_trial(combo$model, combo$tool, PROMPT_KEY)
        cbind(summarise_run(trial, score_trial(trial)), error = NA_character_)
      }, error = function(e) {
        data.table(exact = NA_integer_, acceptable = NA_integer_,
                   candidates = NA_integer_, n_fields = NA_integer_,
                   elapsed_sec = NA_real_, mapping = NA_character_,
                   error = conditionMessage(e))
      })

      row <- cbind(
        data.table(ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                   model = combo$model, tool = combo$tool,
                   prompt = PROMPT_KEY, condition = condition, rep = rep),
        row
      )

      fwrite(row, OUT_CSV, append = file.exists(OUT_CSV))

      message(sprintf("  %2d/%d  exact %s  acceptable %s  %.1fs%s",
                      rep, N_REPS, row$exact, row$acceptable,
                      ifelse(is.na(row$elapsed_sec), 0, row$elapsed_sec),
                      ifelse(is.na(row$error), "", paste0("  ERROR: ", row$error))))
    }
  }
}

message("\nDone. Results in ", OUT_CSV)

## ---------------------------------------------------------------- summary
##
## Run summarise_stochasticity() after the loop finishes, or in a second
## R session while it is still going.

summarise_stochasticity <- function(path = OUT_CSV) {
  d <- fread(path)[is.na(error)]

  overall <- d[, .(
    n                = .N,
    distinct_mappings = uniqueN(mapping),
    exact_min        = min(exact),
    exact_max        = max(exact),
    exact_mode       = as.integer(names(sort(table(exact), decreasing = TRUE))[1]),
    median_sec       = round(median(elapsed_sec), 1)
  ), by = .(model, tool, condition)]

  ## Which fields actually move. A field that is constant across every run is
  ## not interesting; the ones that flip are the finding.
  long <- d[, {
    parts <- strsplit(unlist(strsplit(mapping, "|", fixed = TRUE)), "=", fixed = TRUE)
    .(field = vapply(parts, `[`, character(1), 1),
      value = vapply(parts, `[`, character(1), 2))
  }, by = .(model, tool, condition, rep)]

  unstable <- long[, .(distinct_values = uniqueN(value),
                       values = paste(sort(unique(value)), collapse = " / ")),
                   by = .(model, tool, condition, field)][distinct_values > 1]

  list(overall = overall[], unstable = unstable[])
}