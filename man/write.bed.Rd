% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.bed.R
\name{write.bed}
\alias{write.bed}
\title{Custom bed files}
\usage{
write.bed(
  ldpath = NULL,
  path,
  type.len = "len",
  len = 5,
  type.gene = "gene",
  gene.region = NULL
)
}
\arguments{
\item{ldpath}{path for plink bfiles for the chromosome-wise LD reference;(e.g. 1000 Genome Project)}

\item{path}{specify the output path}

\item{type.len}{'len'/'snp'}

\item{len}{length (kb/number of SNPs)}

\item{type.gene}{'gene'/'wg': by gene TSS or whole genome}

\item{gene.region}{regulatory region of the gene (a data.frame with columns: V1(chr),V2(start),V3(end),V4(gene name))}
}
\description{
Using this function to generate blockwise LD matrix
}
