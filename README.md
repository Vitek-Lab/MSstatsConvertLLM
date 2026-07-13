# Local LLMs for Automated Schema Identification in MSstatsConvert

## Overview

This project explores the use of local large language models (LLMs) to automate schema identification for mass spectrometry proteomics data as part of the [MSstatsConvert](https://github.com/Vitek-Lab/MSstatsConvert) preprocessing pipeline. Instead of relying on manual or rule-based mapping of input file columns to the standardized MSstats schema, this approach uses a locally-hosted LLM to infer column mappings automatically, reducing the manual overhead required to onboard new data formats.

This work was presented as a poster at ASMS 2026.

## Repository Structure

```
.
├── R/    # Code for schema inference pipeline (LLM prompting, parsing, integration with MSstatsConvert)
├── results/    # Output/evaluation results from testing on example datasets
└── README.md
```

## Requirements

In order to use this project you need Ollama and a local LLM installed

## Usage

*(entry point instructions to be added)*

## Next steps

TODO

## Citation

If you use this work, please cite the associated ASMS 2026 poster (citation details TBD).
