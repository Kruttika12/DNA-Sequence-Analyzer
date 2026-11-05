# DNA Sequence Analyzer

## Overview

DNA Sequence Analyzer is a Python-based bioinformatics tool designed to perform basic DNA sequence analysis using multi-sequence FASTA files. It provides sequence statistics, reverse complement generation, motif detection, open reading frame (ORF) identification, protein translation, and sequence comparison.

The tool processes multiple DNA sequences, allows the user to select an individual sequence for detailed analysis, and exports sequence statistics and ORF results into CSV files.

## Features

### 1. FASTA File Processing
- Reads multi-sequence FASTA files.
- Extracts sequence identifiers and nucleotide sequences.
- Validates DNA sequences for permitted nucleotide characters (A, T, G, C, N).
- Identifies empty sequences and invalid characters.

### 2. Sequence Statistics
- Calculates total sequence length.
- Counts individual nucleotides (A, T, G, C, N).
- Calculates GC and AT content percentages.
- Determines the percentage composition of each nucleotide.

### 3. Reverse Complement Generation
- Generates the reverse complement of the selected DNA sequence.
- Handles ambiguous nucleotide N.

### 4. Motif Detection
- Searches for user-defined DNA motifs.
- Reports the number of motif occurrences.
- Identifies motif positions using 1-based indexing.

### 5. Open Reading Frame (ORF) Detection
- Identifies ORFs on forward and reverse-complement strands.
- Examines all three reading frames on each strand.
- Identifies start codons (ATG) and in-frame stop codons (TAA, TAG, TGA).
- Reports ORF strand, reading frame, start and end positions, and sequence length.

### 6. Protein Translation
- Translates identified ORF sequences into protein sequences using the standard genetic code (NCBI translation table 1).

### 7. Sequence Comparison
- Compares nucleotide composition and sequence length across all sequences in the FASTA file.
- Generates a comparative summary of sequence statistics.

## Requirements

- Python 3
- Biopython

### Python Libraries

- `Bio.Seq` (Biopython)
- `collections`
- `csv`

The `collections` and `csv` modules are included in the Python standard library.

## Installation

Install Biopython using pip:

```bash
pip install biopython
```

## Usage

Run the Python script:

```bash
python DNA_sequence_analyser.py
```

1. Enter the path to a multi-sequence FASTA file when prompted.
2. Select a sequence using its sequence ID.
3. View the sequence statistics and reverse complement.
4. Enter a DNA motif to perform motif detection.
5. View detected ORFs across both DNA strands and all six reading frames.
6. View the translated protein sequences.
7. Review the comparative statistics of all sequences.

## Output Files

The program generates two CSV files:

| File | Description |
|---|---|
| `sequence_comparison.csv` | Contains sequence IDs, lengths, nucleotide counts, GC content, and AT content for all sequences. |
| `orf_results.csv` | Contains detected ORF information, including strand, reading frame, start and end positions, length, and nucleotide sequence. |

## Limitations

- Supports DNA sequences containing A, T, G, C, and N.
- ORF detection uses ATG as the start codon and TAA, TAG, and TGA as stop codons.
- ORF detection does not establish biological gene functionality.
- Analysis is performed on one selected sequence at a time, while sequence comparison covers all sequences in the input FASTA file.

## Author

Kruttika R