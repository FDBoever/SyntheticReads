# SyntheticReads

## Simulating synthetic long-read datasets from multi-genome inputs


SyntheticReads (`inoa-generate-reads.py`) is a Python tool for generating **synthetic long-read sequencing datasets** from multiple reference genomes.  
It enables reproducible, coverage-driven simulation of long-read data with configurable read length distributions, sequencing error profiles, and support for circular genomes.

The tool is designed for **benchmarking genome assemblers, read mappers, and metagenomic analysis pipelines** using ground-truth datasets.

---

## Features

- Simulate long reads from **multiple FASTA references**
- **Coverage-driven** read generation (no fixed read count)
- Skew-normal read length distributions
- Per-read stochastic sequencing error rates
- Support for **circular and linear genomes**
- Outputs:
  - FASTA (sequences)
  - FASTQ (with Phred scores)
  - TSV metadata (origin, coordinates, error rates)
- Optional diagnostic plots for read length and quality distributions

---

## Requirements

Python ≥ 3.8

Required Python packages:
- `numpy`
- `scipy`
- `pandas`
- `biopython`
- `plotly`
- `tabulate`

You can install dependencies with:
```bash
pip install numpy scipy pandas biopython plotly tabulate
```


## Input Files

### 1. Reference genomes

A directory containing one or more FASTA (or multi-FASTA) files.

Example:

```
genomes/
├── ecoli.fna
├── paeruginosa.fna
└── plasmid.fna
```

### 2. TSV configuration file

A tab-delimited file specifying coverage depth and circularity for each reference.

Example (`rel.ab.txt`):

```
fasta_file	coverage	is_circular
ecoli.fna	30	        True
paeruginosa.fna	10	    False
plasmid.fna	100	        True
```

- coverage defines sequencing depth per reference
- is_circular should be True or False


## Usage

Basic example

```
python inoa-generate-reads.py \
  -i genomes \
  -o output \
  --tsv_file rel.ab.txt \
  --prefix synthetic_reads
```

---

## Parameters

| Parameter                      | Flag                       | Type   | Default    | Description                                                                             |
| ------------------------------ | -------------------------- | ------ | ---------- | --------------------------------------------------------------------------------------- |
| Input directory                | `-i`, `--input_directory`  | string | *required* | Path to directory containing FASTA or multi-FASTA reference genomes                     |
| Output directory               | `-o`, `--output_directory` | string | *required* | Directory where output files will be written                                            |
| TSV configuration file         | `--tsv_file`               | string | *required* | Tab-delimited file specifying input FASTA files, coverage values, and circularity flags |
| Output prefix                  | `--prefix`                 | string | `output`   | Prefix used for naming all generated output files                                       |
| Read length peak               | `--peak`                   | float  | `13000`    | Mode (peak) of the skew-normal read length distribution (in bp)                         |
| Read length standard deviation | `--std_dev`                | float  | `5000`     | Standard deviation of the read length distribution                                      |
| Read length skewness           | `--skewness`               | float  | `3.5`      | Skewness parameter controlling asymmetry of the read length distribution                |
| Mean error rate                | `--error_mean`             | float  | `0.001`    | Mean per-read sequencing error rate                                                     |
| Error rate standard deviation  | `--error_std_dev`          | float  | `0.0005`   | Standard deviation of the per-read error rate distribution                              |

Default values approximate PacBio HiFi–like reads.

---

## Coverage-driven read generation

The total number of reads is not specified by the user.

For each reference sequence of length L and target coverage C, the number of reads N is computed as:

    N = ceil((C × L) / μ)

where μ is the mean of the sampled read length distribution.

---

## Read length model

Read lengths are sampled from a skew-normal distribution parameterized by:
- peak (mode)
- standard deviation
- skewness

This allows modeling of asymmetric long-read length profiles typical of modern sequencing platforms.

---

## Error model

- Each read is assigned an independent error rate drawn from a normal distribution
- Errors are introduced as random base substitutions
- Error rates are clipped to the range [0, 1]
- FASTQ quality scores are computed as:

    Q = -10 × log10(error_rate)

- Quality scores are capped at Q50
- Each read receives a uniform per-base quality score

---

## Circular genome handling

For circular references:
- The sequence is internally duplicated
- Reads may span the origin of replication
- Coordinates are normalized to the original sequence length

Supported use cases include bacterial chromosomes, plasmids, and organellar genomes.

---

## Output Files

For a given prefix <prefix>:

- <prefix>.*.generated_reads.fasta
- <prefix>.*.generated_reads.fastq
- <prefix>.*.generated_reads.tsv

The TSV file records:
- Read ID
- Source FASTA
- Contig ID
- Read length
- Start and end positions
- Error rate
- Q-score

---

## Example: ONT-like simulation

```
    python inoa-generate-reads.py \
      -i genomes \
      -o output \
      --tsv_file rel.ab.txt \
      --prefix ont_sim \
      --peak 10000 \
      --std_dev 30000 \
      --skewness 20 \
      --error_mean 0.01 \
      --error_std_dev 0.05
```

---


## Description

We developed SyntheticReads (SyntheticReads.py), a Python script for simulating complex synthetic long-read datasets from multiple reference genomes. By enabling bottom-up construction of long-read datasets for mock metagenomes using high-quality genomic sequences, the tool supports controlled benchmarking of genome assemblers, taxonomic profilers, and read mappers.  

The script accepts a directory of multi-FASTA files as input and supports user-defined relative abundances via a tab-delimited file, facilitating the simulation of uneven community structures reflective of real-world metagenomic samples. 

Reads are extracted according to a skew-normal distribution of lengths, parameterized by a user-defined mean, standard deviation, and skewness. This approach permits the generation of length distributions representative of contemporary long-read sequencing platforms, including PacBio HiFi and Oxford Nanopore Technologies (ONT). To simulate the latter, users can specify higher variance and greater skewness in read length distributions, along with elevated base-level error rates. Errors are introduced stochastically per read, drawn from a normal distribution with configurable mean and variance, allowing for simulation of platform-specific error profiles. 

The tool includes automated handling of circular genomic elements, which are internally duplicated to enable extraction of fragments that span the origin of replication. This ensures proper simulation of complete circular sequences, such as organellar genomes, bacterial chromosomes, and plasmids. Output includes synthetic reads in FASTA and FASTQ format, accompanied by a metadata table that tracks the origin coordinates, sequence source, and applied error characteristics for each fragment. Optional graphical summaries of read length and error rate distributions are also generated.  

By enabling the in-silico generation of mock community reads with tunable complexity, read length profiles, and sequencing error models, this tool provides a reproducible framework for generating ground truth datasets. These datasets are well-suited for the systematic benchmarking and refinement of metagenomic and genome analysis pipelines. 

# SyntheticReads