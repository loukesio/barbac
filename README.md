## barbac: an R package for barocde data analysis

**barbac** is a versatile R package tailor-made for researchers and bioinformaticians working with barcoded sequences. Whether you're delving into the vast realm of genomic data or investigating the temporal dynamics of molecular barcodes, `barbac` has got you covered.

Install the package using the following commands  <img align="right" src="graphic_elements/barbac_logo.png" width=400>


### Key Features:

1. **Barcode Extraction**: Efficiently and accurately extract barcode sequences from your genomic datasets.

2. **Clustering of Barcodes**: Group similar barcode sequences with the in-built clustering algorithm, providing meaningful insights into sequence families or operational taxonomic units (OTUs).

3. **Time Series Trajectories**: Dive deep into the temporal dynamics of barcoded sequences. `barbac` not only calculates but also visually represents the trajectory of each barcode over time, making it easier to discern patterns, trends, and anomalies.

### Applications:

- **Molecular Tracking**: Monitor the abundance and variations of specific molecular barcodes over different time points.
- **Biodiversity Analysis**: Understand the diversity and prevalence of barcode sequences in different samples.
- **Evolutionary Dynamics**: Track and visualize the emergence, dominance, and decline of barcoded sequences across time.

---

Absolutely! Here's a clean and user-friendly **`README.md`** draft for your `barbac` package, based on what you've built so far.

---

## 📦 `barbac`: A Lightweight R Wrapper for Barcode-Based Bioinformatics Workflows

`barbac` is an R package that wraps commonly used command-line bioinformatics tools (FastQC, MultiQC, PEAR, minimap2, samtools) into an easy-to-use R workflow.
It is especially designed for processing paired-end reads mapped to barcode reference sequences.

---

## 🚀 Features

* 🔬 Run **FastQC** on raw reads
* 📊 Summarize results using **MultiQC**
* 🔗 Merge paired-end reads with **PEAR**
* 🎯 Map merged reads to a reference using **minimap2**
* 📎 Sort and index BAM files using **samtools**
* 📈 Summarize mapped/unmapped read counts

---

## 📁 Expected Input

A `sample_table` (data frame or CSV) with the following columns:

| sample  | R1                     | R2                     |
| ------- | ---------------------- | ---------------------- |
| sample1 | data/sample1\_R1.fastq | data/sample1\_R2.fastq |
| sample2 | data/sample2\_R1.fastq | data/sample2\_R2.fastq |

---

## ⚙️ Setup Instructions

Before using `barbac`, install the required command-line tools in a Conda environment:

```r
# Install bioinformatics tools in a clean Conda env
barbac::configure_environment()
```

This will:

* Create a Conda environment named `barbac_env`
* Install `fastqc`, `multiqc`, `pear`, `minimap2`, and `samtools` via bioconda

---

## 📌 Available Functions

### 🔬 Quality Control

```r
barbac::run_fastqc(sample_table, output_dir = "fastQC")
barbac::run_multiqc(input_dir = "fastQC", output_dir = "multiQC")
```

### 🔗 Merge Paired-End Reads

```r
barbac::run_pear_merge(sample_table, output_dir = "merged")
```

### 🧬 Map to Reference + BAM Processing

```r
barbac::run_minimap2(merged_dir = "merged", reference = "path/to/barcode.fasta")
```

### 📊 BAM Stats Summary

```r
barbac::summarise_bam_stats(bam_dir = "merged/bam")
```

---

## ✅ Optional: Check if Tools Are Available

```r
barbac::check_barbac_tools()
```

Returns a logical vector indicating if each CLI tool is available in your system path.

---

## 📦 Installation (Development Version)

```r
# You need devtools or remotes
remotes::install_github("loukesio/barbac")
```

---

## 📝 License

MIT License © \[Your Name / Institute]

---

Would you like me to also generate a `DESCRIPTION` file based on this content and current dependencies?
