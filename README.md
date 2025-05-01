# Mutant_Library-Analysis
This script is designed to process deep sequencing FASTQ files. It can analyze the mutation frequencies and coverage of different variants in the pyr1-DSM-Hao library and the pyr1-Triple library. It also supports analysis of deep sequencing data generated using an improved evSeq method, enabling the identification of mutations in biosensors.

# Analysis of Biosensor Libraries from Deep Sequencing

This repository contains Python scripts for analyzing deep sequencing data of biosensor mutant libraries, including:

- **pyr1-DSM-Hao library**
- **pyr1-Triple library**
- **evSeq**
  
The scripts are designed to process FASTQ files and calculate:

- The proportion of different mutants
- The coverage of the mutant libraries
- The mutation profile of functional biosensors

The code is compatible with sequencing data generated using an **improved version of the evSeq method**.

## Scripts and related files
   ### Merge Paired-End Reads
   - `./PYR1_DSM-Hao/step1_merge_dsm.py`: merge two fastq files with forward and reverse reads from deep sequencing of DSM-Hao library
   ### Mutation Calling and Quality Control Workflow
   - `./PYR1_DSM-Hao/step2_counts_dsm.py`:
  
   This script processes merged paired-end deep sequencing FASTQ files to identify mutations in biosensor libraries. The workflow includes the following steps:

   1. **Filtering Reads Containing Ambiguous Bases**  
      - Any read containing the character `'N'` at any position is discarded.

   2. **Quality Filtering at Designed Mutation Sites**  
      - For each read, base quality is checked at all **designed mutation positions**.  
      - A read is retained only if **all mutation positions have a Phred quality score ≥ Q20**.

   3. **Strand-Aware Mutation Translation**  
      - If the read is from the **forward strand**, it is directly translated to identify amino acid mutations.  
      - If the read is from the **reverse strand**, the reverse complement is taken before translation.

   4. **Summary Table Output**  
      The script generates a summary table including:
      - ✅ Total number of high-quality reads passing all filters  
      - ❌ Number of reads discarded due to quality issues or ambiguous bases  
      - 🧬 Number of wild-type reads  
      - 🔬 A list of identified mutations and their corresponding read counts  
   
   ### Mutant Proportion and Coverage Calculation
   - `./PYR1_DSM-Hao/step3_coverage_dsm.py`:
   This script performs the following analyses:

   1. 📊 **Calculates the proportion of different mutation types**, including wild-type (WT), single, double, up to 9-site mutants.
   2. 🧩 **Computes the coverage of mutants** generated from different Golden Gate Assembly (GGA) combinations:  
      *Coverage = observed mutants in combo / expected mutants in combo*
   3. 🌐 **Estimates overall library coverage**:  
      *Total coverage = total observed mutants / total expected mutants*
   4. ⚠️ **Calculates the proportion of mutants containing deletions**.

   - `./PYR1_Triple/step1_merge_triple.py`: it has same function as `step1_merge_dsm.py`.
   - `./PYR1_Triple/step2_counts_triple.py`: it has same function as `step2_counts_dsm.py`.
   - `./PYR1_Triple/step3_coverage_triple.py`: Unlike `step3_coverage_dsm`, this script calculates the coverage of different triple-mutant combinations generated.
   - `./PYR1_evSeq/step1_merge.py`: it has same function as `step1_merge_dsm.py`.
   ### Quality Control and Read Normalization
   - `./PYR1_evSeq/step2_uniform.py`:
   This script performs quality control by removing reads with base quality scores (Q) below 20 at any of the designed mutation sites. All retained reads are then converted to the forward-strand sequence of PYR1 for consistency.
   ### Mutation Frequency Analysis
   - `./PYR1_evSeq/step3_counts.py`:
   This script analyzes mutations in deep sequencing files and their corresponding frequencies. It generates a table where each row represents a specific pair of barcodes, corresponding to individual hits.
   - For each barcode pair, the script identifies and ranks the top three mutations by frequency:
     1. Most frequent mutation
     2. Second most frequent mutation
     3. Third most frequent mutation
   - Mutations with frequencies lower than the top three are discarded from the output.


## Installation

To create the required environment using `conda`, run the following commands in your terminal:

```bash
conda env create -f environment.yml -n pyr1-analysis
conda activate pyr1-analysis
```

## Usage (Example: `PYR1_DSM-Hao` folder)

   ### 1. Merge Paired-End Reads

   Run the following command to merge paired-end FASTQ files:

   ```bash
   python3 step1_merge_dsm.py -l 409 -m 4 your_fastq_file_1.fq your_fastq_file_2.fq > your_merged_fastq_file.fq
   ```
   `-l 409`: Specifies that the theoretical merged reads length is 409 bp.
   `-m 4`: Allows up to 4 bases to be mismatched during the merging process.
   ### 2. Calculates the proportion of different mutation types
   ```bash
   python3 step2_counts_dsm.py
   ```
   ### 3. Computes the coverage of mutants
   ```bash
   python3 step3_coverage_dsm.py
