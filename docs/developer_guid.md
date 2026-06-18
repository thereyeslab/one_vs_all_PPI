# Developer Guide

## Overview

This document describes the internal architecture of pooled-PPI and is intended for developers who want to:

* Understand the codebase
* Extend existing functionality
* Add new pooling strategies
* Add new AlphaFold dialects
* Modify score extraction
* Implement new correction methods

This guide focuses on software architecture rather than end-user usage.

For usage instructions, see:

```text
workflow.md
config.md
installation.md
```

---

# Design Philosophy

The workflow was designed around several principles:

1. Modularity
2. Reproducibility
3. Config-driven execution
4. AlphaFold backend independence
5. Scalability to genome-wide screens

Each major workflow stage is implemented as an independent module.

```text
Filtering and msa cashing (ancillary modules) 
Pooling
JSON generation
Validation
Extraction
Correction
Aggregation
```

This allows stages to be developed, tested, and executed independently.

---

# Project Structure

```text
pooled-ppi/
│
├── configs/
│
├── docs/
│
├── src/
│   ├── workflow.py
│   ├── generate_pool.py (Pooling)
│   ├── write_af3_jsons.py (JSON generation)
│   ├── validate_af3_jsons.py (Validation)
│   ├── extract_af3_pair_scores.py (Extraction)
│   ├── size_correction.py (Correction)
│   ├── run_pooled_ppi_sampler.py (Not used in the current workflow)
│   ├── seq_filtering/filter_proteins.py (Filtering)
│   └── msa_cache/build_msa_cache_jsons.py (msa caching for HPC approach)
│
├── data/
│
└── environment.yml
```

The exact directory structure may evolve over time, but the architectural principles remain the same.

---

# Workflow Architecture

The workflow consists of two major phases:

```text
prepare
process
```

---

## Prepare

```text
candidate_proteins.fasta
        ↓
pooling
        ↓
AF3 json generation
        ↓
validation
```

Outputs:

```text
pools.tsv
pool_summary.tsv
pool_chain_mapping.tsv
alphafold3_jsons/
af3_json_validation_summary.tsv
```

---

## Process

```text
AF3 predictions
        ↓
score extraction
        ↓
size correction
        ↓
replicate aggregation
```

Outputs:

```text
pair_scores_raw.tsv
prediction_summary.tsv
pair_scores_size_corrected_observations.tsv
pair_scores_size_corrected_aggregated.tsv
```

---

# workflow.py

`workflow.py` is the main orchestration layer.

Responsibilities:

* Load configuration
* Build project paths
* Build/Inject module-specific configuration
* Validate inputs
* Execute workflow stages

The workflow intentionally contains very little scientific logic.

All biological and computational operations are delegated to independent modules.

---

# Configuration System

The workflow is fully configuration-driven.

Configuration is loaded from:

```text
config.yaml
```

and passed to modules as a dictionary.

The workflow dynamically injects generated paths into the configuration.

Example:

```python
config["af3_result_input"] = {
    "predictions_dir": ...
}
```

This design avoids hardcoding paths inside modules.

---

# Pooling Module

Responsibilities:

* Load FASTA sequences
* Compute sequence lengths
* Generate pools based on the decided strategy
* Track replicate coverage
* Generate pooling summaries

Outputs:

```text
pools.tsv
pool_summary.tsv
```

---

## Supported Pooling Modes

### one_vs_all

Single bait protein screened against many prey proteins.

Used for:

```text
bait vs proteome
```

screening workflows.

---

### all_vs_all

Generalized pooling strategy, not fully developed.

Designed for future support of complete interaction network reconstruction.

---

## Supported Pooling Algorithms

### probabilistic

Current default.

Characteristics:

* Length-weighted sampling
* Prioritizes unseen proteins (not pairs)
* Minimizes unnecessary pools
* Optimized for one-vs-all screening

---

### deterministic

Implements the original pooled-PPI strategy described in:

```text
Nature Structural & Molecular Biology (2025)
pooled-PPI
```

Suitable for all-vs-all workflows.

---

# AlphaFold JSON Generation

Responsibilities:

* Convert pools into AlphaFold inputs
* Create chain mappings
* Support multiple AlphaFold dialects

Outputs:

```text
alphafold3_jsons/
pool_chain_mapping.tsv
```

---

# Dialect System

The JSON builder uses a dialect abstraction layer.

Currently supported:

```text
alphafoldserver
alphafold3
```

The goal is to isolate backend-specific formatting.

This allows future support for:

```text
AlphaFold Server
AlphaFold3 local
AlphaFold3 HPC
AF3 variants
future AF models
```

without modifying upstream workflow logic.

---

# Chain Mapping

One critical design decision is the use of:

```text
pool_chain_mapping.tsv
```

AlphaFold outputs only contain chain identifiers.

The mapping file provides the translation:

```text
chain ID
      ↓
protein ID
```

and is therefore required for score extraction.

Any modifications to chain generation must maintain compatibility with this file.

---

# Validation Module

Responsibilities:

* Validate generated JSON files ( for any af3 backend)
* Validate chain mappings
* Detect malformed inputs

Checks include:

* Missing sequences
* Duplicate identifiers
* Invalid chain mappings
* Pool size limits
* AlphaFold compatibility

The validation layer exists to catch problems before expensive AlphaFold runs.

---

# Score Extraction Module

Responsibilities:

* Parse AlphaFold outputs
* Extract pairwise interaction scores
* Map chains back to proteins
* Generate observation tables

Outputs:

```text
pair_scores_raw.tsv
prediction_summary.tsv
```

---

## Backend Independence

The extractor should operate on AlphaFold outputs rather than on the execution environment.

Current support includes:

```text
AlphaFold Server outputs
AlphaFold3 HPC outputs
```

The extraction layer should remain agnostic to how predictions were generated.

---

# Size Correction Module

Motivation:

According to the paper, the af3 is biased by the size of the proteins ( the longer they are, the higher chance to detect an interaction). So they correct the outputted iptm either by using a fixed number or by fitting the data to a regression to find the formula for the size correction. 


The correction module estimates this effect and adjusts scores accordingly.

Outputs:

```text
pair_scores_size_corrected_observations.tsv
size_correction_model.tsv
```

---

## Supported Methods

### pooled_ppi_fixed

Uses fixed coefficients derived from published pooled-PPI benchmarks.

---

### fit_from_data

Fits correction coefficients directly from the current dataset.

Recommended when sufficient observations are available.

---

# Aggregation Module

Responsibilities:

* Merge replicate observations
* Compute summary statistics
* Generate final interaction scores

Current implementation:

```text
mean
```

Future strategies may include:

```text
median
weighted mean
robust estimators
```

---

# MSA Caching Support

Large-scale AlphaFold3 workflows may use MSA caching.
NOTE: The module for generating the final af3 input json files, ( incorporation of msa cach to the jsons) is not developed for now. If you want to use msa cach, you need to write the code.

The workflow intentionally separates:

```text
AF3 input generation
AF3 execution
score extraction
```

So that cached MSA pipelines can be inserted without modifying pooled-PPI internals.

Future optimization work should remain external to the core workflow whenever possible.

See:

```text
alphafold3_optimization.md
```

---

# Adding New Pooling Strategies

To add a new pooling strategy:

1. Implement the algorithm
2. Return pools in the standard schema
3. Preserve output compatibility

Required outputs:

```text
pools.tsv
pool_summary.tsv
```

Downstream modules should not require modification.

---

# Adding New AlphaFold Dialects

To add a new AlphaFold backend:

1. Create a new dialect formatter
2. Generate valid JSON
3. Preserve chain mapping compatibility

Required outputs:

```text
alphafold3_jsons/
pool_chain_mapping.tsv
```

If these outputs remain unchanged, downstream modules will continue to function.

---

# Testing Recommendations

Before large-scale runs:

Validate:

```text
1 small pool
5 pools
10 pools
```

before launching hundreds or thousands of AlphaFold jobs.

Recommended checks:

* JSON generation
* Validation
* AF3 execution
* Score extraction
* Final aggregation

Small-scale tests catch most workflow issues at minimal computational cost.

---

# Future Development Ideas

Potential future extensions include:

* Additional pooling algorithms
* Alternative correction models
* Automatic MSA caching integration
* AF3 feature reuse
* Network-level scoring
* Experimental prior integration
* Additional AlphaFold backends
* Visualization utilities

The modular architecture was designed to allow these features to be added with minimal impact on existing workflow components.

---

# Contributing

When adding new functionality:

1. Preserve existing output schemas whenever possible.
2. Avoid hardcoded paths.
3. Prefer configuration-driven behavior.
4. Keep AlphaFold-specific logic isolated.
5. Maintain compatibility with existing TSV outputs.

The workflow is designed around stable interfaces between modules. Maintaining these interfaces ensures backward compatibility and reproducibility.
