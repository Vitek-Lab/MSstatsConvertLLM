# zzz.R — load-time setup

# OLLAMA_URL and BENCHMARK_DATA_DIR are active bindings rather than ordinary
# objects. Top-level package code runs when the package is installed, so a
# plain `Sys.getenv()` assignment baked the default into the installed
# namespace and the documented environment variables had no effect. Bound
# here, each use re-reads the environment, so setting MSSTATS_OLLAMA_URL
# after library(MSstatsConvertLLM) works.
.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  bind <- function(name, variable, default) {
    if (exists(name, envir = ns, inherits = FALSE)) return(invisible(NULL))
    makeActiveBinding(name, function() Sys.getenv(variable, default), ns)
  }
  bind("OLLAMA_URL", "MSSTATS_OLLAMA_URL", "http://localhost:11434")
  bind("BENCHMARK_DATA_DIR", "MSSTATS_BENCHMARK_DATA", "benchmark_data")
}
