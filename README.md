# AmpVarPro
# ONT Amplicon Sequencing Variant Calling Pipeline
# Overview
This repository contains a variant calling pipeline tailored for amplicon sequencing data generated using Oxford Nanopore Technologies (ONT). It automates the process from raw reads to annotated variants, facilitating rapid analysis of genetic variants. The pipeline integrates several bioinformatics tools, including minimap2 for mapping reads to a reference genome, samtools for processing BAM files, fastp for read quality filtering, and Clair3 for variant calling. Variants are annotated using snpEff.

# Prerequisites

- Python 3.x

- minimap2

- samtools

- fastp

- Clair3 (and its dependencies)

- snpEff

# Installation

Clone the repository:

```bash
git clone https://github.com/ajinkyakhilari/AmpVarPro.git
```

```bash
cd AmpVarPro
```

# Install dependencies:

Most dependencies can be installed via conda or bioconda channels. It is recommended to create a new conda environment.

```bash
conda create -n ont-vc python=3.8 minimap2 samtools fastp snpEff -c bioconda
```

```bash
conda activate ont-vc
```

#### Follow the installation instructions of Clair3 from its official GitHub repository.

# Usage
Prepare your environment:

Before running the pipeline, make sure all dependencies are installed and your environment is activated.

# Run the pipeline:

The pipeline is executed through a Python script that requires several arguments, including paths to input files and directories, reference genome, and quality filtering parameters.

```bash
AmpVarPro.py --input_directory /path/to/input --min_length 100 --max_length 1500 --threads 8 --reference_genome /path/to/ref.fasta --phred_quality 20 --snpeff_data /path/to/snpeff/data --model /path/to/model
```
Replace the placeholder paths and parameters with your actual data and desired settings.

### Example

Assuming you have a directory named data containing subdirectories for each barcode (e.g., barcode01, barcode02, ...), a reference genome file ref.fasta, and the necessary snpeff data and model files:

```bash
python AmpVarPro.py --input_directory ./data --min_length 100 --max_length 1500 --threads 4 --reference_genome ./ref.fasta --phred_quality 20 --snpeff_data ./snpeff/data --model ./clair/model
```

This command will process each barcode directory in parallel, perform variant calling, and output annotated variants.

# Contributing
Contributions to the pipeline are welcome. Please fork the repository, make your changes, and submit a pull request.

# License
MIT License
