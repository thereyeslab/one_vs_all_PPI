# Filtering Module Notes (Users and developers)

## Purpose

The filtering module prepares a clean candidate protein set before pooling.

It takes:

```text
data/raw/caulobacter_proteome.fasta
data/raw/caulobacter_uniprot_annotations.tsv
```

and produces:

```text
data/processed/all_proteins_annotated.csv
data/processed/candidate_proteins.csv
data/processed/excluded_proteins_with_reasons.csv
data/processed/candidate_proteins.fasta
```

The main goal is to remove proteins that are unlikely to be suitable for soluble pooled-PPI screening, while preserving traceability of every exclusion decision.

---

# Input Assumptions

The filtering module assumes two synchronized input files:

## 1. UniProt annotation table

Expected path:

```text
data/raw/caulobacter_uniprot_annotations.tsv
```

The table must contain at minimum:

```text
Entry
```

The `Entry` column is used as the UniProt accession ID.

The script also uses several annotation columns when available:

```python
TEXT_COLUMNS = [
    "Entry",
    "Entry Name",
    "Gene Names",
    "Protein names",
    "Organism",
    "Subcellular location [CC]",
    "Transmembrane",
    "Gene Ontology (GO)",
    "Gene Ontology (cellular component)",
    "Gene Ontology (biological process)",
    "Gene Ontology (molecular function)",
    "Keywords",
]
```

Missing optional columns are skipped silently by `combined_annotation_text()`.

---

## 2. Proteome FASTA

Expected path:

```text
data/raw/caulobacter_proteome.fasta
```

The FASTA headers are expected to follow the UniProt format:

```text
>tr|A0A0H3C2V8|A0A0H3C2V8_CAUVN ... GN=CCNA_00145 ...
```

The script uses two different identifiers from the FASTA header:

| Identifier        | Source      | Purpose                                 |
| ----------------- | ----------- | --------------------------------------- |
| UniProt accession | `record.id` | Match FASTA records to annotation table |
| Gene/locus ID     | `GN=...`    | Rewrite candidate FASTA headers         |

Example:

```text
UniProt accession: A0A0H3C2V8
Gene/locus ID:     CCNA_00145
```

---

# Identifier Design

A critical design decision is that the final `candidate_proteins.fasta` is rewritten to use gene/locus IDs as FASTA record IDs.

Example output:

```fasta
>CCNA_00145
MSEQ...
```

This is important because the rest of the pooled-PPI workflow uses these IDs consistently:

```text
candidate_proteins.fasta header = CCNA_00145
pooling ID = CCNA_00145
MSA JSON name = CCNA_00145
MSA cache file = CCNA_00145.unpaired.a3m
pooled AF3 JSON protein_id = CCNA_00145
```

Do not change this behavior unless all downstream modules are updated accordingly.

---
# Duplicate Protein Sequences and Locus IDs

## Background

In some genomes, multiple locus IDs may encode proteins with identical amino-acid sequences.

For example:

```text
CCNA_XXXXX → Protein Sequence A
CCNA_YYYYY → Protein Sequence A
```

Although the locus IDs differ, the translated protein sequence is identical.

---

## UniProt Annotation Behavior

When UniProt collapses identical proteins into a single protein entry, only one locus identifier may be retained in the annotation table.

Example:

```text
Genome:
    CCNA_00001
    CCNA_00002

Both encode identical proteins.

UniProt annotation:
    Entry → A0A1234567
    Gene Names → CCNA_00001
```

In this situation:

```text
CCNA_00002
```

may not appear anywhere in the downloaded UniProt annotation file.

---

## Current Filtering Behavior

The filtering workflow is annotation-driven.

Candidate proteins are selected using:

```python
candidate_ids = set(kept["Entry"].astype(str))
```

and then matched back to FASTA records using the UniProt accession.

As a consequence:

* Only proteins represented in the annotation table can be retained.
* If UniProt collapses multiple locus IDs into a single protein entry, only the representative locus ID will appear in the final outputs.

---

## Effect on Output Files

The following files will contain only the representative locus:

```text
candidate_proteins.csv
candidate_proteins.fasta
```

Example:

```text
Genome FASTA:
    CCNA_00001
    CCNA_00002

UniProt annotation:
    CCNA_00001

Output:
    candidate_proteins.fasta
        >CCNA_00001
```

while:

```text
CCNA_00002
```

will not appear.

---

## Biological Interpretation

This behavior is generally acceptable for pooled-PPI because:

1. The encoded protein sequence is identical.
2. AlphaFold operates on protein sequence rather than locus identity.
3. Including both loci would produce duplicate predictions.
4. Removing duplicates reduces computational cost without losing structural information.

In practice, pooled-PPI screens proteins rather than genomic loci.

---

## Important Caveat

The current workflow does not explicitly detect or report collapsed locus IDs.

If locus-level tracking is important for a project, developers should implement an additional reconciliation step between:

```text
proteome FASTA
```

and

```text
UniProt annotation table
```

to identify:

* duplicated protein sequences
* collapsed locus identifiers
* one-to-many mapping relationships

before filtering is performed.

For most protein interaction screening applications, the current behavior is acceptable and avoids redundant AlphaFold predictions.

---
# Filtering Logic

The filtering module creates an `exclusion_reason` column.

If no exclusion rule is triggered:

```text
filter_status = keep
```

If one or more exclusion rules are triggered:

```text
filter_status = exclude
```

Multiple exclusion reasons are stored as semicolon-separated values.

Example:

```text
membrane_bound;periplasmic
```

---

# Current Exclusion Rules

## Ribosomal proteins

Implemented in:

```python
is_ribosomal(row, text)
```

A protein is excluded as ribosomal if:

* Any gene name starts with:

  * `rps`
  * `rpl`
  * `rpm`
* Or annotation text contains terms such as:

  * `ribosomal`
  * `ribosome`
  * `30s ribosomal protein`
  * `50s ribosomal protein`

---

## Membrane-bound proteins

Implemented in:

```python
is_membrane_bound(row, text)
```

A protein is excluded as membrane-bound if:

* The `Transmembrane` field is non-empty
* Or annotation text contains terms such as:

  * `membrane`
  * `transmembrane`
  * `integral membrane`
  * `cell membrane`
  * `cytoplasmic membrane`
  * `inner membrane`

---

## Outer membrane or cell wall synthesis proteins

Implemented in:

```python
is_outer_membrane_or_cell_wall(row, text)
```

A protein is excluded if annotation text contains terms such as:

* `outer membrane`
* `cell wall`
* `cell envelope`
* `peptidoglycan`
* `murein`
* `lipopolysaccharide`
* `lps biosynthesis`
* `udp-n-acetylmuramate`
* `muramoyl`
* `penicillin-binding protein`
* `glycosyltransferase involved in cell wall`

---

## Periplasmic proteins

Implemented in:

```python
is_periplasmic(row, text)
```

A protein is excluded if annotation text contains:

* `periplasm`
* `periplasmic`
* `periplasmic space`

---

# Text Normalization

All annotation text is normalized using:

```python
normalize_text(value)
```

This function:

* Converts missing values to an empty string
* Converts all text to lowercase

This makes filtering rules case-insensitive and robust to missing values.

---

# Combined Annotation Text

The function:

```python
combined_annotation_text(row)
```

combines all useful annotation fields into one searchable lowercase string.

Most filtering rules operate on this combined text.

This design makes it easy to add new filters without repeatedly checking individual annotation columns.

---

# FASTA Export Logic

Candidate proteins are selected using UniProt accessions:

```python
candidate_ids = set(kept["Entry"].astype(str))
```

Each FASTA record is matched by extracting its UniProt accession:

```python
extract_uniprot_accession_from_fasta_id(record.id)
```

For retained proteins, the FASTA header is rewritten using the gene/locus ID:

```python
gene_id = extract_gene_name_from_description(record.description)

record.id = gene_id
record.name = gene_id
record.description = gene_id
```

This ensures the downstream FASTA uses biologically meaningful and stable locus identifiers.

---

# Important Developer Note

The script raises an error if a retained FASTA record does not contain `GN=`:

```python
raise ValueError(f"Could not find GN= in FASTA header: {description}")
```

This is intentional.

Without a gene/locus identifier, the downstream workflow may produce inconsistent protein IDs.

If supporting organisms or FASTA files without `GN=` fields, developers should implement a fallback identifier strategy explicitly.

Do not silently fall back to UniProt accessions unless the entire downstream workflow is designed to use UniProt IDs.

---

# Adding a New Exclusion Rule

To add a new filtering rule:

1. Implement a new Boolean function.
2. Call it inside `get_exclusion_reasons()`.
3. Append a short machine-readable reason string.

Example:

```python
def is_secreted(row: pd.Series, text: str) -> bool:
    terms = [
        "secreted",
        "signal peptide",
    ]
    return any(term in text for term in terms)
```

Then update:

```python
def get_exclusion_reasons(row: pd.Series) -> str:
    text = combined_annotation_text(row)
    reasons = []

    if is_secreted(row, text):
        reasons.append("secreted")

    return ";".join(sorted(set(reasons)))
```

Reason strings should be:

* Lowercase
* Short
* Machine-readable
* Stable across versions

Good examples:

```text
secreted
membrane_bound
low_confidence_annotation
```

Avoid:

```text
This protein seems to be secreted
```

---

# Output Files

## all_proteins_annotated.csv

Contains the full annotation table plus:

```text
exclusion_reason
filter_status
```

Useful for auditing and QC.

---

## candidate_proteins.csv

Subset of proteins with:

```text
filter_status = keep
```

---

## excluded_proteins_with_reasons.csv

Subset of proteins with:

```text
filter_status = exclude
```

This is the most important QC file for reviewing filtering decisions.

---

## candidate_proteins.fasta

FASTA file containing retained proteins.

Headers are rewritten to gene/locus IDs.

Example:

```fasta
>CCNA_00145
MSEQ...
```

This file is the primary input for the pooled-PPI workflow.

---

# Quality Control Recommendations

After running the filtering module, check:

```bash
wc -l data/processed/candidate_proteins.csv
wc -l data/processed/excluded_proteins_with_reasons.csv
grep ">" data/processed/candidate_proteins.fasta | head
```

Confirm that FASTA headers look like:

```text
>CCNA_00145
```

not:

```text
>tr|A0A0H3C2V8|A0A0H3C2V8_CAUVN
```

Also inspect exclusion reasons:

```bash
cut -d, -f exclusion_reason data/processed/excluded_proteins_with_reasons.csv | sort | uniq -c
```

---

# Known Limitations

Current implementation is organism-specific in practice because it assumes:

* UniProt-style FASTA headers
* `GN=` gene/locus identifiers
* Caulobacter-style locus IDs
* Annotation columns downloaded from UniProt

To support another organism, developers should verify:

1. FASTA headers contain usable gene/locus IDs.
2. Annotation table uses the expected column names.
3. Filtering terms are biologically appropriate.
4. Downstream pooling IDs match FASTA IDs.

---

# Recommended Refactor

The current script uses hardcoded paths:

```python
RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
```

For future development, consider converting this script into a configurable command-line tool:

```bash
python src/seq_filtering/filter_proteins.py \
  --annotation data/raw/annotations.tsv \
  --fasta data/raw/proteome.fasta \
  --output-dir data/processed
```

This would make the filtering module easier to reuse for other organisms and projects.
