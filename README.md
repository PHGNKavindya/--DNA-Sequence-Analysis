# DNA Sequence Analysis

This repository contains the Python implementations, DNA sequence data, and computational results developed for Assignment 02 of the DS5109 – Data Science for Bioinformatics module at Sabaragamuwa University of Sri Lanka.

## Project Overview

The assignment explores three important bioinformatics applications using Python and biological sequence data retrieved from the NCBI database.

### Question 01 – Human INS Gene Analysis

Retrieved and analysed the human insulin (INS) gene DNA sequence using NCBI.

**Tasks performed:**
- Retrieved the DNA sequence from NCBI.
- Downloaded and saved the sequence in FASTA format.
- Explained the purpose of accession numbers and FASTA records.
- Calculated the GC content.
- Generated the reverse-complement sequence.
- Explained the importance of reverse-complement analysis when identifying open reading frames (ORFs).

**Results:**
- NCBI Accession: `NC_060935.1`
- Sequence Length: 1,431 bp
- GC Content: 64.71%

### Question 02 – Greedy Motif Search

Implemented a greedy motif-search algorithm to identify a common sequence pattern among four human disease-associated genes.

**Genes analysed:**
- G6PD – Glucose-6-phosphate dehydrogenase
- F8 – Coagulation factor VIII
- DMD – Dystrophin
- PAH – Phenylalanine hydroxylase

**Methodology:**
- Retrieved genomic DNA records from NCBI.
- Extracted the annotated gene sequences.
- Selected the first 1,000 bases of each gene in biological 5′-to-3′ orientation.
- Applied a greedy motif-search algorithm.
- Identified the consensus motif.

**Results:**
- Motif Length: 10 bp
- Consensus Motif: `CTGGCCTCCA`
- Best Motif Mismatch Score: 4

### Question 03 – BRCA1 and BRCA2 Sequence Alignment

Developed a divide-and-conquer global sequence-alignment algorithm based on Hirschberg's algorithm.

**Methodology:**
- Retrieved BRCA1 and BRCA2 RefSeqGene records from NCBI.
- Extracted the longest annotated CDS feature from each record.
- Implemented Hirschberg's linear-space alignment algorithm.
- Evaluated matches, mismatches, gaps and alignment identity.

**Results:**
- BRCA1 CDS Length: 5,592 bp
- BRCA2 CDS Length: 10,257 bp
- Alignment Length: 10,393 columns
- Exact Matching Pairs: 4,715
- Alignment Score: -963
- Identity among non-gap columns: 86.42%

## Technologies and Tools

- Python 3.11
- Biopython
- NCBI Nucleotide Database
- NCBI Entrez E-utilities
- Greedy Motif Search
- Hirschberg's Divide-and-Conquer Algorithm

## Installation

Install the required Python libraries:

```bash
pip install biopython truststore
```

## Running the Scripts

Run the Python scripts using:

```bash
python question1_ins_analysis.py
```

```bash
python question2_greedy_motif.py
```

```bash
python question3_brca_alignment.py
```

**Note:** Configure your email address in the Entrez settings of each script before retrieving sequences from NCBI.

## Repository Structure

```text
DNA-Sequence-Analysis/
│
├── README.md
│
├── Question_01_INS_Analysis/
│   ├── question1_ins_analysis.py
│   ├── INS_sequence.fasta
│   ├── INS_reverse_complement.fasta
│   └── INS_analysis_results.txt
│
├── Question_02_Greedy_Motif_Search/
│   ├── question2_greedy_motif.py
