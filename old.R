library(ellmer)
library(data.table)
library(MSstatsConvert)
library(jsonlite)

models <- c(
  "llama3.2:3b-instruct-q4_K_M",
  "llama3.2:3b",
  "llama3.1:8b","deepseek-r1:14b"
)

chat <- chat_ollama(
  model = models[[4]],
  params = params(
    format      = "json",  # forces JSON
    temperature = 0,       # greedy decoding
    top_p       = 1,
    seed        = 1,       # fixed RNG
    num_ctx     = 4096,    # ensure prompt fits
    num_predict = -1       # don't truncate early
  )
)




# Test data - DIA-NN
pd_raw = system.file("tinytest/raw_data/PD/pd_input.csv", 
                     package = "MSstatsConvert")
pd_raw = data.table::fread(pd_raw)

# Example data 1
example_1_data = fread(system.file(
  "tinytest/raw_data/Spectronaut/spectronaut_input.csv", package = "MSstatsConvert"))

normalize_headers <- function(df) {
  nms <- names(df)
  nms <- gsub("[^A-Za-z0-9._]+", ".", nms)
  nms <- gsub("\\.+", ".", nms)
  nms <- sub("^\\.|\\.$", "", nms)
  names(df) <- nms
  df
}

make_json_preview <- function(df, n_rows = 2, max_chars = 40) {
  df <- head(df, n_rows)
  profiles <- lapply(names(df), function(col) {
    vals <- df[[col]]
    # coerce to character, truncate
    vals <- as.character(vals)
    vals <- substr(vals, 1L, max_chars)
    vals <- vals[!is.na(vals) & vals != ""]
    vals <- unique(vals)
    vals <- head(vals, n_rows) # at most 3 examples
    list(
      column = col,
      examples = vals
    )
  })
  toJSON(profiles, auto_unbox = TRUE, pretty = TRUE)
}

# Example:


MSSTATS_FIELDS <- c(
  "ProteinName","PeptideSequence","PrecursorCharge",
  "FragmentIon","ProductCharge","Run","Intensity", "Qvalue"
)

prompt <- paste0(
  "You are a converter that maps tabular column names to the MSstats schema.\n",
  "Return ONLY valid JSON. Do all reasoning silently; never show your thoughts.\n",
  "Rules:\n",
  "1) Only use column names that appear in DATA HEADER (the new dataset).\n",
  "2) If a required field is missing, set \"from\": null and include \"candidates\".\n",
  "3) Never invent column names. Never copy anything from EXAMPLE HEADER/ROWS.\n",
  "4) Output keys must be exactly: mappings, confidence, notes, warnings.\n",
  "5) The \"mappings\" array must contain exactly 8 objects in the required order.\n",
  "6) The Qvalue column must contain the letter \"Q\". Set the Qvalue mapping’s from to null and include candidates.\n\n",
  "Violation if any \"from\" is not in DATA HEADER.\n\n",
  "Task: Map dataset columns to the MSstats schema.\n",
  "Return ONLY valid JSON. Do not invent columns.\n\n",
  "MSstats required fields in this EXACT order:\n- ",
  paste(MSSTATS_FIELDS, collapse = "\n- "), "\n\n",
  "OUTPUT REQUIREMENTS:\n",
  "- ProteinName must be a string (uniprot id/gene id/ect.), not a number \n",
  " - Top-level keys: \"mappings\", \"confidence\", \"notes\", \"warnings\" (notes/warnings can be empty arrays).\n",
  "- \"mappings\" MUST be an ARRAY of EXACTLY 8 OBJECTS (one per field above), in that exact order.\n",
  "- Each object MUST contain: \"field\", \"from\" (string or null), and \"confidence\" (0..1).\n",
  "- \"field\" is the MSstats column name/\n",
  "- \"from\" Is the name of the column that is mapped to the MSstats column name.\n",
  "- If \"from\" is null, include a \"candidates\" ARRAY sorted by descending \"score\",\n",
  "  where each candidate is an object {\"from\": <string>, \"score\": <0..1>}.\n",
  "- \"from\" must be the highest-confidence candidate available (if any).\n",
  "EXAMPLE 1  (structure-only; adjust column names to the new dataset):\n",
  "RAW DATA:\n",
  ex1_rows_text,
  "EXPECTED MAPPING:\n",
  "{\n",
  "  \"mappings\": [\n",
  "    {\"field\": \"ProteinName\",     \"from\": \"PG.ProteinGroups\", ",
  "\"confidence\": 0.8,",
    "\"candidates\": [\n",
  "       {\"from\": \"PG.ProteinGroups\", \"score\": 0.8},\n",
  "       {\"from\": \"PG.ProteinAccessions\", \"score\": 0.7}\n",
  "     ]},\n",
  "    {\"field\": \"PeptideSequence\", \"from\": \"EG.ModifiedSequence\", ",
  "\"confidence\": 0.9,",
  "\"candidates\": [\n",
  "       {\"from\": \"EG.ModifiedSequence\", \"score\": 0.9},\n",
  "       {\"from\": \"PEP.GroupingKey\", \"score\": 0.8}\n",
  "     ]},\n",
  "    {\"field\": \"PrecursorCharge\", \"from\": \"FG.Charge\", ",
  "\"confidence\": 0.92}\n",
  "    {\"field\": \"FragmentIon\", \"from\": \"F.FrgIon\", ",
  "\"confidence\": 0.9}\n",
  "    {\"field\": \"ProductCharge\", \"from\": \"F.Charge\", ",
  "\"confidence\": 0.8}\n",
  "    {\"field\": \"Run\",             \"from\": \"FileName\", ",
  "\"confidence\": 0.99},\n",
  "    {\"field\": \"Intensity\",       \"from\": \"F.PeakArea\", ",
  "\"confidence\": 0.6,",
  "\"candidates\": [\n",
  "       {\"from\": \"F.PeakArea\", \"score\": 0.6},\n",
  "       {\"from\": \"F.NormalizedPeakArea\", \"score\": 0.5}\n",
  "       {\"from\": \"F.NormalizedPeakHeight\", \"score\": 0.5}\n",
  "       {\"from\": \"F.PeakHeight\", \"score\": 0.5}\n",
  "     ]},\n",
  "    {\"field\": \"Qvalue\",       \"from\": \"PG.Qvalue\", ",
  "\"confidence\": 0.5,",
  "\"candidates\": [\n",
  "       {\"from\": \"PG.Qvalue\", \"score\": 0.5},\n",
  "       {\"from\": \"EG.Qvalue\", \"score\": 0.5}\n",
  "     ]},\n",
  "  ],\n",
  "  \"confidence\": 0.93,\n",
  "  \"notes\": [],\n",
  "  \"warnings\": []\n",
  "}\n\n",
  
  "Now map the following dataset:\n",
  "DATA SAMPLE:\n",
  new_rows_text
)


# One mapping entry per MSstats field
mapping_entry <- type_object(
  field      = type_string(),                 # required (default)
  from       = type_string(),  # may be null/omitted
  confidence = type_number(),                  # 0..1
  candidates = type_array(                     # optional alternatives
    type_object(from = type_string(), score = type_number())
  )
)

# Top-level result: array of mapping entries (+ optional metadata)
mapping_type <- type_object(
  "LLM mapping proposal for MSstats",
  mappings    = type_array(mapping_entry),                 # required
  confidence  = type_number(),             # global
  notes       = type_array(type_string()),
  warnings    = type_array(type_string()),
)


prompt <- "You are a converter that maps tabular column names to the MSstats schema.\n\nReturn ONLY valid JSON. Do all reasoning silently; never show your thoughts.\n\nHARD RULES:\n- Use ONLY columns that appear in DATA HEADER (the new dataset).\n- If a required field is missing, set \"from\": null and include \"candidates\".\n- Never invent column names. Never copy from EXAMPLE data.\n- Output keys must be exactly: mappings, confidence, notes, warnings.\n- \"mappings\" MUST contain exactly 8 objects in this exact order:\n  1) ProteinName\n  2) PeptideSequence\n  3) PrecursorCharge\n  4) FragmentIon\n  5) ProductCharge\n  6) Run\n  7) Intensity\n  8) Qvalue\n- Qvalue mapping: column NAME must contain the letter 'Q' (case-insensitive).\n  If no such column exists, set \"from\": null and include candidates that DO\n  contain 'Q'. Do NOT choose columns without 'Q' for Qvalue.\n- \"ProteinName must be a string identifier (e.g., UniProt/Gene), not numeric.\"\n\nDECISION HEURISTICS (apply in order; first match wins):\n- ProteinName: choose the first header matching any of:\n  /(Protein.*(Accessions|Accession|Groups|Group|IDs?|ID|Names?|Name))/i\n  /(PG\\.(Protein.*|.*Accession|.*Groups))/i\n  Reject columns whose values are purely numeric across examples.\n- PeptideSequence: prefer names matching:\n  /(Sequence|Peptide|StrippedSequence|ModifiedSequence)/i\n  If multiple, prefer those containing \"Modified\" or \"Stripped\".\n- PrecursorCharge: prefer headers matching:\n  /(Precursor.*Charge|FG\\.Charge|\\bCharge\\b)/i\n- FragmentIon: must look like ion labels (e.g., y7, b3). Prefer:\n  /(Fragment.*Ion|Frg(Ion|Type)|IonType|F\\.Frg(Ion|Type))/i\n  Reject pure counts like \"Matched.Ions\", \"Total.Ions\".\n  If no valid ion-label column exists, set null.\n- ProductCharge: prefer fragment-level charge:\n  /(Product.*Charge|Fragment.*Charge|F\\.Charge)/i\n  If only a generic \"Charge\" exists and it was used for PrecursorCharge,\n  do NOT reuse it; set null.\n- Run: prefer in this order if present:\n  1) /(Spectrum\\.File|FileName|R\\.FileName)/i\n  2) /(Run|R\\.Run|RawFile|Search\\.ID)/i\n- Intensity: prefer quantitative signal columns in this order:\n  /(Intensity$|Intensity[^A-Za-z]|^Intensity)/i\n  /(Precursor\\.Area|PeakArea|NormalizedPeakArea|PeakHeight)/i\n- Qvalue: choose by header name only (must contain 'Q' or 'q'):\n  /(Qvalue|Q-value|QVal|Q_?val|PG\\.Qvalue|EG\\.Qvalue)/i\n  If none, set null and list candidates that do contain 'Q'.\n\nCONFIDENCE SCORING (fixed, deterministic):\n- 0.98 exact regex match on highest-priority alias.\n- 0.92 secondary alias match.\n- 0.80 generic match (e.g., \"Charge\" for PrecursorCharge).\n- 0.00 if null.\n\nCANDIDATE LISTING:\n- For any chosen \"from\", include up to 3 descending-score candidates that\n  matched the regexes (including the chosen one).\n- If null, include up to 3 best regex matches (or empty if none).\n\nVALIDATION FILTERS:\n- ProteinName must NOT be purely numeric in examples.\n- FragmentIon must contain ion labels (letters + numbers). Exclude pure counts.\n- Qvalue must contain 'Q' IN THE HEADER NAME.\n\nOUTPUT FORMAT:\n{\n  \"mappings\": [\n    {\"field\":\"ProteinName\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"PeptideSequence\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"PrecursorCharge\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"FragmentIon\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"ProductCharge\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"Run\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"Intensity\",\"from\":..., \"confidence\":..., \"candidates\":[...]},\n    {\"field\":\"Qvalue\",\"from\":..., \"confidence\":..., \"candidates\":[...]}\n  ],\n  \"confidence\": <mean of non-null field confidences (0..1)>,\n  \"notes\": [],\n  \"warnings\": []\n}\n\nReject any output that is not valid JSON."

resp <- chat_ollama(
  model = models[[4]],
  system_prompt = prompt,
  params = params(
    format      = "json",  # forces JSON
    temperature = 0,       # greedy decoding
    top_p       = 1,
    seed        = 1,       # fixed RNG
    num_ctx     = 4096,    # ensure prompt fits
    num_predict = -1       # don't truncate early
  )
)

input = system.file("tinytest/raw_data/Metamorpheus/QuantifiedPeaks.tsv", 
                    package = "MSstatsConvert")
meta_input = make_json_preview(normalize_headers(data.table::fread(input)), n_rows=5)
pd_input = make_json_preview(normalize_headers(pd_raw), n_rows=5)
spect_input = make_json_preview(normalize_headers(example_1_data), n_rows=5)

user_prompt = paste0("DATA HEADER:\n", spect_input, "\n RETURN JSON only")

spec <- chat$chat_structured(user_prompt, type = mapping_type, convert = TRUE)
spec

chat$chat(user_prompt)

