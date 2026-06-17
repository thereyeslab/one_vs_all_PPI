# Installation

## Overview

This document describes how to install and set up the pooled-PPI workflow.

The pooled-PPI codebase is independent of AlphaFold 3 and only requires a Python environment to run.

AlphaFold predictions can be generated using either:

* AlphaFold Server (recommended)
* Local AlphaFold 3 installation
* HPC AlphaFold 3 installation

The workflow itself is responsible for:

* Protein filtering
* Pool generation
* AlphaFold3 JSON generation
* Input validation
* Score extraction
* Size correction
* Replicate aggregation

---

# System Requirements

## Operating System

Tested on:

* macOS
* Linux

Windows is not officially supported.

---

## Python

Recommended:

```text
Python 3.10+
```

The workflow has been tested with:

```text
Python 3.10
```

---

# Clone the Repository

```bash
git clone <repository_url>
cd pooled-ppi
```

---

# Create the Conda Environment

The repository includes an `environment.yml` file containing all required dependencies.

Create the environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate pooled_ppi
```

---

# Environment Contents

The environment installs the core workflow dependencies:

* numpy
* pandas
* scipy
* biopython
* pyyaml
* matplotlib
* seaborn
* statsmodels
* tqdm
* numba
* af3io

AlphaFold 3 itself is **not** installed through this environment.

---

# Verify Installation

Verify that the workflow can be executed:

```bash
python src/workflow.py --help
```

You can also test that the required packages are available:

```bash
python -c "import pandas, numpy, yaml, Bio"
```

If no errors are produced, the installation is successful.

---

# Project Structure

A typical repository layout looks like:

```text
pooled-ppi/
│
├── configs/
│
├── data/
│   ├── raw/
│   ├── processed/
│   └── output/
│
├── docs/
│
├── src/
│
├── environment.yml
│
└── README.md
```

---

# AlphaFold Options

The workflow supports multiple AlphaFold execution backends.

## Option 1 — AlphaFold Server (Recommended)

For most users, AlphaFold Server is the simplest and most reliable option.

Advantages:

* No installation required
* No GPU required
* No database downloads
* Minimal setup
* Suitable for most projects

Workflow:

```text
Generate JSON files
        ↓
Upload to AlphaFold Server
        ↓
Download predictions
        ↓
Run pooled-PPI processing
```

Limitations:

* Manual upload and download
* Daily job quota
* Less automation

---

## Option 2 — Local AlphaFold 3

Requires:

* Linux
* NVIDIA GPU
* Docker or Apptainer
* AlphaFold 3 model weights
* AlphaFold 3 databases

Recommended only for advanced users.

---

## Option 3 — HPC AlphaFold 3

Requires:

* Access to an HPC cluster
* SLURM
* Apptainer/Singularity
* AlphaFold 3 model weights
* AlphaFold 3 databases

See:

```text
docs/alphafold3_hpc.md
```

for detailed setup instructions.

---

# AlphaFold Model Weights

Local and HPC AlphaFold 3 installations require model weights.

Weights must be requested directly from Google DeepMind.

The pooled-PPI repository does not distribute AlphaFold model weights.

---

# AlphaFold Databases

Local and HPC AlphaFold 3 installations require AlphaFold databases.

Depending on the AlphaFold version, the databases may require several hundred gigabytes of storage.

The pooled-PPI repository does not distribute AlphaFold databases.

---

# Recommended First Run

Before launching a large screen:

1. Run the filtering workflow on a small dataset.
2. Generate pools.
3. Generate AlphaFold3 JSON files.
4. Submit a small test job to AlphaFold Server.
5. Download the results.
6. Run the processing stage.

This validates the complete workflow before investing significant compute resources.

---

# Known Limitations

Current implementation:

* Supports one-vs-all screening
* Supports deterministic pooling
* Supports probabilistic pooling
* Supports AlphaFold Server
* Supports local/HPC AlphaFold 3

Current implementation does not support:

* Full all-vs-all workflow execution
* Automatic AlphaFold execution from within pooled-PPI
* Automatic MSA caching integration

---

# Updating the Environment

If dependencies change in future versions:

```bash
conda env update -f environment.yml --prune
```

The `--prune` option removes packages that are no longer required.

---

# Troubleshooting

## Missing Python Packages

Example:

```text
ModuleNotFoundError
```

Solution:

```bash
conda env update -f environment.yml --prune
```

---

## AlphaFold JSON Validation Errors

Run the workflow in prepare mode and inspect:

```text
af3_json_validation_summary.tsv
```

for details.

---

## AlphaFold Prediction Failures

Verify:

* JSON dialect matches the selected backend
* AlphaFold version is correct
* Pool size is within supported limits
* Input sequences are valid proteins

---

# Next Steps

After installation:

1. Read `docs/workflow.md`
2. Review `docs/config.md`
3. Prepare your input FASTA
4. Generate AlphaFold inputs
5. Run AlphaFold
6. Process the results

```
```


# Citation

If you use pooled-PPI in a publication, please cite:

1. The pooled-PPI publication.
2. AlphaFold 3.
3. This repository.

```
```
