# SyntheticReads

## Simulating synthetic long-read datasets from multi-genome inputs

## Description

We developed SyntheticReads (SyntheticReads.py), a Python script for simulating complex synthetic long-read datasets from multiple reference genomes. By enabling bottom-up construction of long-read datasets for mock metagenomes using high-quality genomic sequences, the tool supports controlled benchmarking of genome assemblers, taxonomic profilers, and read mappers.  

The script accepts a directory of multi-FASTA files as input and supports user-defined relative abundances via a tab-delimited file, facilitating the simulation of uneven community structures reflective of real-world metagenomic samples. 

Reads are extracted according to a skew-normal distribution of lengths, parameterized by a user-defined mean, standard deviation, and skewness. This approach permits the generation of length distributions representative of contemporary long-read sequencing platforms, including PacBio HiFi and Oxford Nanopore Technologies (ONT). To simulate the latter, users can specify higher variance and greater skewness in read length distributions, along with elevated base-level error rates. Errors are introduced stochastically per read, drawn from a normal distribution with configurable mean and variance, allowing for simulation of platform-specific error profiles. 

The tool includes automated handling of circular genomic elements, which are internally duplicated to enable extraction of fragments that span the origin of replication. This ensures proper simulation of complete circular sequences, such as organellar genomes, bacterial chromosomes, and plasmids. Output includes synthetic reads in FASTA and FASTQ format, accompanied by a metadata table that tracks the origin coordinates, sequence source, and applied error characteristics for each fragment. Optional graphical summaries of read length and error rate distributions are also generated.  

By enabling the in-silico generation of mock community reads with tunable complexity, read length profiles, and sequencing error models, this tool provides a reproducible framework for generating ground truth datasets. These datasets are well-suited for the systematic benchmarking and refinement of metagenomic and genome analysis pipelines. 