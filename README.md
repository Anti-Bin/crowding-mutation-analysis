# Variant Processing and Annotation Pipeline
This repository contains a custom variant filtering and de novo mutation identification pipeline designed to refine raw variant calls and accurately estimate mutation rates. The pipeline removes sequencing artifacts and applies stringent post-calling filters to ensure high-confidence mutation detection.

# Overview

The pipeline performs:
1. Extraction of ±3bp FASTA sequence around each variant.
2. Removal of variants in ≥3bp homopolymer regions.
3. Variant decomposition and normalization (`vt`).
4. Removal of variants present in F1 reference VCFs.
5. Querying VCFs using `bcftools`.
6. Identification of common and unique variants.
7. Counting filtered variants by quality metrics.
8. Classification of SNP, MNP, DEL, INS types.
9. Subclassification of six SNP substitution categories.
10. Annotation of variants using provided annotation tables.
11. Counting of genic, intronic, intergenic, and CDS categories.

# Requirements

| Tool | Version (tested) | Install |
|------|------------------|----------|
| bash | ≥5.1 | default |
| awk, grep, sort, wc | — | default |
| perl | ≥5.34 | default |
| bedtools | ≥2.30 | `conda install -c bioconda bedtools` |
| vt | ≥0.577 | `conda install -c bioconda vt` |
| bcftools | ≥1.13 | `conda install -c bioconda bcftools` |
| gunzip  | ≥1.10 | `sudo apt install gunzip` |







