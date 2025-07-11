% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_minimap2.R
\name{run_minimap2}
\alias{run_minimap2}
\title{Map merged reads to reference using minimap2 and process with samtools}
\usage{
run_minimap2(
  merged_dir = "merged",
  reference,
  output_dir = file.path(merged_dir, "bam")
)
}
\arguments{
\item{merged_dir}{Directory where merged fastq files are stored (default: "merged")}

\item{reference}{Path to reference FASTA file}

\item{output_dir}{Directory for BAM files (default: "merged/bam")}
}
\value{
Vector of commands run
}
\description{
Map merged reads to reference using minimap2 and process with samtools
}
