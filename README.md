# telo_comparative_study

## Introduction

This repository contains the code used to analyze sequencing data from the Sequence Read Archive in support of the manuscript "Benchmarking Long-read Sequencing Tools for Chromosome End-specific Telomere Analysis" by Jake Reed et al. The manuscript is currently under review at bioRxiv. This pipeline was used to process all the samples in the study and generate the results presented in the manuscript. If you would like to replicate the results or use the pipeline for your own analysis, please follow the instructions below.

## Features

* Integrates telseq, telogator, and TECAT
* Supports short-read (BAM) and long-read (FASTQ) data
* Platform-specific optimizations for ONT/PacBio
* Parallel processing support

## Prerequisites

* conda/miniconda
* R ≥ 4.0
* Python ≥ 3.8
* 16GB+ RAM
* ~2x input data size for temp storage
* Telseq, Telogator2, and TECAT installed

### R Dependencies

```R
tidyverse
ggplot2
cowplot
parallel
optparse
TECAT
```

## Quick Start

```bash
# Clone repo
git clone https://github.com/jake-bioinfo/telo_comparative_study.git
cd telo_comparative_study

# Install dependencies
# The applications used in this study can be found at the following links:
# Telseq: https://github.com/zd1/telseq
# Telogator2: https://github.com/zstephens/telogator2
# TECAT: https://github.com/jake-bioinfo/tecat

# Run pipeline
./compare.sh -t ~/tmp/analysis_temp \
             -o ~/results \
             -i ~/data/sample1 \
             -s HG002 \
             -n 18 \
             -p ont
```

## Usage

```bash
./compare.sh -t <temp_directory> \
             -o <output_directory> \
             -i <input_dir> \
             -s <sample> \
             -n <threads> \
             -p <platform>
```

### Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-t` | Temp directory | Yes |
| `-o` | Output directory | Yes |
| `-i` | Input directory | Yes |
| `-s` | Sample name | Yes |
| `-n` | Thread count | No |
| `-p` | Platform (ont/pb) | No |

### Output Structure

```
output_directory/
├── sample_name/
│   ├── telseq/
│   │   └── telseq_results.txt
│   ├── telomerehunter/
│   │   └── results/
│   ├── telogator/
│   │   └── tlens_by_allele.tsv
│   └── tecat/
│       ├── frequencies/
│       ├── truncated/
│       ├── mapped/
│       └── tecat_results.csv
```

## License

[GPL-3.0](LICENSE)

## Citation

If you use this pipeline in your research, please cite the following manuscript:

```bibtex
@article{reed2021benchmarking,
  title={Benchmarking Long-read Sequencing Tools for Chromosome End-specific Telomere Analysis},
  author={Reed, Jake, Oelkuct, Mark, and Coombes, Kevin},
  journal={bioRxiv},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```