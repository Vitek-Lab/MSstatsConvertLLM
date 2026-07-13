# ============================================================
# prompts.R — Prompt templates for schema inference
# ============================================================

# ------------------------------------------------------------------
# LEAN PROMPT: Minimal instructions, tests whether the LLM can
# infer mappings from schema knowledge + examples alone.
# This is the "can it actually reason?" condition.
# ------------------------------------------------------------------
PROMPT_LEAN <- "
You are a proteomics data converter. Your task is to map column names
from an input dataset to the MSstats analysis schema.

MSstats requires these fields (in this exact order):
1. ProteinName     — protein identifier (e.g., UniProt accession)
2. PeptideSequence — peptide or modified peptide sequence
3. PrecursorCharge — charge state of the precursor ion
4. FragmentIon     — fragment ion label (e.g., y7, b3)
5. ProductCharge   — charge of the fragment ion
6. Run             — raw file name or run identifier
7. Intensity       — quantitative signal (peak area, intensity, etc.)
8. Qvalue          — q-value or FDR for quality filtering

Rules:
- Only use column names from the DATA HEADER below.
- Each source column may be assigned to at most ONE MSstats field. No duplicates.
  If the best candidate is already used, set \"from\": null with \"candidates\".
- If no column clearly maps to a field, set \"from\": null.
- When uncertain, include a \"candidates\" array with alternatives.
- Never invent column names.

Return ONLY valid JSON with this structure:
{
  \"mappings\": [
    {\"field\": \"ProteinName\", \"from\": \"...\", \"confidence\": 0.0-1.0, \"candidates\": [...]},
    ... (8 total, one per field in order above)
  ],
  \"confidence\": <overall 0-1>,
  \"notes\": [],
  \"warnings\": []
}
"

# ------------------------------------------------------------------
# CONSTRAINED PROMPT: Regex heuristics + decision rules baked in.
# Tests whether structured guidance improves smaller models.
# ------------------------------------------------------------------
PROMPT_CONSTRAINED <- "
You are a proteomics data converter. Map column names from an input
dataset to the MSstats schema. Return ONLY valid JSON.

Do all reasoning silently. Never show your thoughts.

MSstats required fields (output in this EXACT order):
1. ProteinName
2. PeptideSequence
3. PrecursorCharge
4. FragmentIon
5. ProductCharge
6. Run
7. Intensity
8. Qvalue

HARD RULES:
- Use ONLY columns from DATA HEADER below.
- Each source column may be assigned to at most ONE field. No duplicates.
  If the best match is already taken, pick the next best or set null.
- If a field has no match, set \"from\": null with \"candidates\".
- Never invent column names.

DECISION HEURISTICS (first match wins):
- ProteinName: match /(Protein.*(Accession|Group|ID|Name))/i
  Reject columns with purely numeric values.
- PeptideSequence: match /(Sequence|Peptide|ModifiedSequence)/i
  Prefer \"Modified\" or \"Stripped\" variants.
- PrecursorCharge: match /(Precursor.*Charge|FG\\.Charge|\\bCharge\\b)/i
- FragmentIon: match /(Fragment.*Ion|Frg.*Ion|IonType)/i
  Must contain ion labels (y7, b3, etc). Reject counts.
  If no valid column, set null.
- ProductCharge: match /(Product.*Charge|Fragment.*Charge|F\\.Charge)/i
  Do NOT reuse PrecursorCharge column. If no separate match, set null.
- Run: prefer (1) /(Spectrum.File|FileName|R\\.FileName)/i then
  (2) /(Run|RawFile)/i
- Intensity: prefer /(Intensity|PeakArea|NormalizedPeakArea)/i
- Qvalue: column name MUST contain 'Q' (case-insensitive).
  Match /(Qvalue|Q.value|QVal)/i. If none, set null. Prefer columns 
  without 'Protein' or 'PG' included to avoid confusion with protein-level q-values.

CONFIDENCE SCORING:
- 0.98 exact primary match
- 0.92 secondary alias match
- 0.80 generic/weak match
- 0.00 if null

Return ONLY valid JSON:
{
  \"mappings\": [
    {\"field\": \"...\", \"from\": ..., \"confidence\": ..., \"candidates\": [...]},
    ... (exactly 8 objects)
  ],
  \"confidence\": <mean of non-null confidences>,
  \"notes\": [],
  \"warnings\": []
}
"

# ------------------------------------------------------------------
# Prompt builder: injects the data preview and acquisition context
#
# acquisition: "DIA", "DDA", or NULL (let the LLM guess)
# ------------------------------------------------------------------
ACQUISITION_CONTEXT <- list(
  DIA = paste0(
    "ACQUISITION TYPE: DIA (Data-Independent Acquisition)\n",
    "This is fragment-level quantification. All 8 MSstats fields are expected.\n",
    "FragmentIon and ProductCharge columns SHOULD exist in the data.\n"
  ),
  DDA = paste0(
    "ACQUISITION TYPE: DDA (Data-Dependent Acquisition)\n",
    "This is precursor-level quantification. FragmentIon and ProductCharge\n",
    "are NOT expected in DDA data. Set both to null with no candidates.\n",
    "Do not attempt to map any column to FragmentIon or ProductCharge.\n"
  )
)

TRANSFORM_CONTEXT <- paste0(
  "PACKED COLUMNS:\n",
  "Some columns may contain MULTIPLE values packed into a single cell,\n",
  "separated by a delimiter (e.g., semicolons). This is common in DIA\n",
  "proteomics where fragment-level quantification values are stored in\n",
  "one column per precursor row.\n\n",
  "If you detect a packed column that maps to an MSstats field (especially\n",
  "Intensity), include it in the \"transforms\" array with:\n",
  "- \"field\": the MSstats field this transform produces (e.g., \"Intensity\")\n",
  "- \"source_column\": the packed column name from DATA HEADER\n",
  "- \"action\": \"split_long\" (split the packed values into separate rows)\n",
  "- \"delimiter\": the separator character (e.g., \";\")\n",
  "- \"id_field\": the MSstats field to generate IDs for (e.g., \"FragmentIon\")\n",
  "- \"id_prefix\": prefix for generated IDs (e.g., \"Frag\" -> Frag1, Frag2, ...)\n\n",
  "When a transform handles a field, set that field's \"from\" in mappings to null\n",
  "and note in warnings that it is handled by a transform.\n",
  "Also set the id_field's \"from\" to null (it will be auto-generated).\n"
)

#' Build the user prompt injected with data preview and context
#'
#' @param data_preview_json JSON preview from [make_json_preview()].
#' @param acquisition `"DIA"`, `"DDA"`, or `NULL`.
#' @param allow_transforms Include the packed-column transform instructions.
#' @return A single prompt string.
#' @export
build_user_prompt <- function(data_preview_json, acquisition = NULL,
                              allow_transforms = FALSE) {
  parts <- character(0)
  
  if (!is.null(acquisition)) {
    acq <- toupper(acquisition)
    if (acq %in% names(ACQUISITION_CONTEXT)) {
      parts <- c(parts, ACQUISITION_CONTEXT[[acq]])
    } else {
      warning("Unknown acquisition type: ", acquisition, 
              ". Use 'DIA' or 'DDA'. Proceeding without acquisition context.")
    }
  }
  
  if (allow_transforms) {
    parts <- c(parts, TRANSFORM_CONTEXT)
  }
  
  parts <- c(parts, paste0("DATA HEADER:\n", data_preview_json))
  parts <- c(parts, "\nRETURN JSON only.")
  
  paste(parts, collapse = "\n")
}

# ------------------------------------------------------------------
# FILTER-AWARE PROMPT: Extends the lean prompt to also discover
# tool-specific filtering columns beyond the 8 MSstats fields.
# This tests whether the LLM understands proteomics QC conventions.
# ------------------------------------------------------------------
PROMPT_FILTER_AWARE <- "
You are a proteomics data converter with domain expertise.

TASK 1 — COLUMN MAPPING:
Map columns from the input dataset to the MSstats analysis schema.

MSstats requires these fields (in this exact order):
1. ProteinName     — protein identifier (e.g., UniProt accession)
2. PeptideSequence — peptide or modified peptide sequence
3. PrecursorCharge — charge state of the precursor ion
4. FragmentIon     — fragment ion label (e.g., y7, b3)
5. ProductCharge   — charge of the fragment ion
6. Run             — raw file name or run identifier
7. Intensity       — quantitative signal (peak area, intensity, etc.)
8. Qvalue          — q-value or FDR for quality filtering

TASK 2 — FILTER DISCOVERY:
Identify columns that should be used to FILTER rows before analysis.
These are tool-specific quality control columns NOT part of the 8 MSstats
fields above. Common examples in proteomics:

- Decoy/contaminant flags (e.g., columns indicating reverse sequences,
  contaminant proteins, or decoy hits that should be removed)
- Exclusion flags (e.g., features flagged as excluded from quantification)
- Fragment loss type (e.g., keeping only 'noloss' fragments)
- Score thresholds (e.g., posterior error probability, percolator scores)
- Modification flags (e.g., filtering out unwanted PTMs)

For each filter, specify:
- \"column\": the exact column name from DATA HEADER
- \"dtype\": expected data type (\"boolean\", \"string\", \"numeric\")
- \"operation\": one of \"equals\", \"not_equals\", \"less_than\",
  \"greater_than\", \"less_than_or_equals\", \"contains\", \"not_contains\"
- \"value\": the filter threshold or value
- \"description\": brief explanation of why this filter matters
- \"confidence\": 0-1 how certain you are this filter is appropriate

Only suggest filters where you are reasonably confident. Do not suggest
filtering on columns already mapped to MSstats fields (e.g., Qvalue).

RULES (apply to both tasks):
- Only use column names from the DATA HEADER below.
- Each source column may be assigned to at most ONE MSstats field.
- If no column maps to a field, set \"from\": null.
- When uncertain, include a \"candidates\" array.
- Never invent column names.

Return ONLY valid JSON:
{
  \"mappings\": [
    {\"field\": \"ProteinName\", \"from\": \"...\", \"confidence\": 0.0-1.0, \"candidates\": [...]},
    ... (8 total, one per MSstats field in order above)
  ],
  \"filters\": [
    {
      \"column\": \"...\",
      \"dtype\": \"boolean|string|numeric\",
      \"operation\": \"equals|not_equals|less_than|greater_than|...\",
      \"value\": ...,
      \"description\": \"...\",
      \"confidence\": 0.0-1.0
    }
  ],
  \"confidence\": <overall 0-1>,
  \"notes\": [],
  \"warnings\": []
}
"

# ------------------------------------------------------------------
# CONSTRAINED + FILTER PROMPT: Regex heuristics for column mapping
# combined with filter discovery. The "best effort" prompt.
# ------------------------------------------------------------------
PROMPT_CONSTRAINED_FILTER <- "
You are a proteomics data converter. Map column names from an input
dataset to the MSstats schema AND identify quality-control filters.
Return ONLY valid JSON. Do all reasoning silently.

TASK 1 — COLUMN MAPPING:

MSstats required fields (output in this EXACT order):
1. ProteinName
2. PeptideSequence
3. PrecursorCharge
4. FragmentIon
5. ProductCharge
6. Run
7. Intensity
8. Qvalue

HARD RULES:
- Use ONLY columns from DATA HEADER below.
- Each source column may be assigned to at most ONE field. No duplicates.
  If the best match is already taken, pick the next best or set null.
- If a field has no match, set \"from\": null with \"candidates\".
- Never invent column names.

DECISION HEURISTICS (first match wins):
- ProteinName: match /(Protein.*(Accession|Group|ID|Name))/i
  Reject columns with purely numeric values.
- PeptideSequence: match /(Sequence|Peptide|ModifiedSequence)/i
  Prefer \"Modified\" or \"Stripped\" variants.
- PrecursorCharge: match /(Precursor.*Charge|FG\\.Charge|\\bCharge\\b)/i
- FragmentIon: match /(Fragment.*Ion|Frg.*Ion|IonType)/i
  Must contain ion labels (y7, b3, etc). Reject counts.
  If no valid column, set null.
- ProductCharge: match /(Product.*Charge|Fragment.*Charge|F\\.Charge)/i
  Do NOT reuse PrecursorCharge column. If no separate match, set null.
- Run: prefer (1) /(Spectrum.File|FileName|R\\.FileName)/i then
  (2) /(Run|RawFile)/i
- Intensity: prefer /(Intensity|PeakArea|NormalizedPeakArea)/i
- Qvalue: column name MUST contain 'Q' (case-insensitive).
  Match /(Qvalue|Q.value|QVal)/i. If none, set null.

CONFIDENCE SCORING:
- 0.98 exact primary match
- 0.92 secondary alias match
- 0.80 generic/weak match
- 0.00 if null

TASK 2 — FILTER DISCOVERY:
Identify columns to FILTER rows before analysis. These are QC columns
NOT part of the 8 MSstats fields above. Look for:

- Decoy/contaminant flags (reverse sequences, contaminant proteins)
- Exclusion flags (features excluded from quantification: remove excluded features.)
- Fragment loss type (e.g., keeping only 'noloss' fragments)
- Score thresholds (posterior error probability, percolator scores)
- Modification flags (filtering out unwanted PTMs)

For each filter, specify:
- \"column\": exact column name from DATA HEADER
- \"dtype\": \"boolean\", \"string\", or \"numeric\"
- \"operation\": \"equals\", \"not_equals\", \"less_than\", \"greater_than\",
  \"less_than_or_equals\", \"contains\", or \"not_contains\"
- \"value\": the filter threshold or value
- \"description\": why this filter matters
- \"confidence\": 0-1

Only suggest filters where you are reasonably confident. Do not filter
on columns already mapped to MSstats fields.

Return ONLY valid JSON:
{
  \"mappings\": [
    {\"field\": \"...\", \"from\": ..., \"confidence\": ..., \"candidates\": [...]},
    ... (exactly 8 objects)
  ],
  \"filters\": [
    {\"column\": \"...\", \"dtype\": \"...\", \"operation\": \"...\",
     \"value\": ..., \"description\": \"...\", \"confidence\": ...}
  ],
  \"confidence\": <mean of non-null mapping confidences>,
  \"notes\": [],
  \"warnings\": []
}
"

# ------------------------------------------------------------------
# Available prompt versions for benchmarking
# ------------------------------------------------------------------

#' Registry of available prompt strategies
#'
#' A named list of system-prompt strings keyed by strategy
#' (`lean`, `constrained`, `filter_aware`, `constrained_filter`). Add new
#' strategies by registering them here.
#'
#' @export
PROMPT_VERSIONS <- list(
  lean                = PROMPT_LEAN,
  constrained         = PROMPT_CONSTRAINED,
  filter_aware        = PROMPT_FILTER_AWARE,
  constrained_filter  = PROMPT_CONSTRAINED_FILTER
)