# OWAS `v1.1.0`
Update 2021.12.1
1. Add the option for clumping GWAS summary statistics, to avoid high collinearity.
2. Bug fixes.

Update 2021.6.9:
1. Add the `write.bed` function, which allows customized segment length and regulatory regions.

Update 2021.6.1:
1. Update the R package and greatly improve the computational efficiency.
2. Add the README file.

## Introduction

**OWAS** is an R implementation of a  computational approach, Openness Weighted Association Studies(OWAS), which leverages and aggregates predictions of chromosome accessibility in personal genomes to prioritize GWAS signals.



## Table of contents
* [Getting Started](#getting-started)
* [Prepare Covariance Matrix](#prepare-covariance-matrix)
* [Prepare Openness Scores](#prepare-openness-scores)
* [Run OWAS](#run-owas)
* [Output](#output)
* [A Simplified Pipeline](#a-simplified-pipeline)
* [Customized Settings](#customized-settings)


## Getting Started
OWAS is an R package which can be installed using the command:
```r
devtools::install_github('shuangsong0110/OWAS')
```

### Download the LD reference panel:
- Download the 1000 Genome Project reference panel (hg19):

`wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz`

`tar -xvzf 1000G_Phase3_plinkfiles.tgz`

- Users could also specify their own LD reference files with plink bfile format (.bim, .fam, .bed).

### Download the pre-computed openness scores:

- Download the precomputed openness scores for 12 ENCODE cell types.
`wget -O openness.tar.gz https://cloud.tsinghua.edu.cn/f/17d8b0b5de3941a2bf86/?dl=1`

We precomputed openness scores for 12 ENCODE cell types using deltaSVM (Lee et al., 2015, Nature Genetics). 

### Download the bedfile for segments:

- Download the bedfile:
`wget -O bedfile_5kb.txt https://cloud.tsinghua.edu.cn/f/7524fb35b88c468dbc02/?dl=1`


For customized opennes scores and segments length, see [Customized Settings](#customized-settings).


## Prepare Covariance Matrix
```r
library(OWAS)
pre.cov(ldpath=PATH_TO_LD_REFERENCE (required), 
	    bedfile=PATH_TO_THE_BEDFILE (required),
	    path=OUTPUT_DIR (required),
	    chr=CHROMOSOME (optional))                    
```
- PATH_TO_LD_REFERENCE (required): The LD reference plink bfile should be dividied into chromosomes and ended with the number of the chromosome. The input should include the file name but not the exact number of chromosome (e.g., ldpath='path/1000G.EUR.QC.', without .bim/.fam/.bed).

- PATH_TO_THE_BEDFILE (required): should link to the exact file. e.g., bedfile='path/bedfile_5kb.txt'

- CHROMOSOME (optional): The chromosome(s) on which the model is run, e.g., `chr=c(1, 3, 5)`. The default is `chr=1:22`.


## Prepare Openness Scores

```r
library(OWAS)
pre.sc(ctype=CELLTYPE (required, can be specified, or simply 'all'),
	   scpath=PATH_TO_OPENNESS_SCORE_FILES (required),
	   ldpath=PATH_TO_LD_REFERENCE (required), 
	   path=OUTPUT_DIR (required, and has to be same with the PATH in get.cov),
	   chr=CHROMOSOME (optional))                    
```


- CELLTYPE (required): We provided 12 precomputed openness scores in `/openness/` file, if the users hope to specify the openness scores for certain cell types, please see [Customized Settings](#customized-settings).

- PATH_TO_OPENNESS_SCORE_FILES: The PATH to the openness scores files. e.g., `scpath=PATH/openness/`

## Run OWAS
```r
library(OWAS)
run.owas(ctype=CELLTYPE (required),
	   gwas=GWAS_SUMMARY_STATISTICS (required)
	   trait=NAME_OF_THE_TRAIT (optional, default='test'),
	   ldpath=PATH_TO_LD_REFERENCE (required), 
	   path=OUTPUT_DIR (required, and has to be same with the PATH in get.cov),
	   chr=CHROMOSOME (optional),
	   plinkpath= PATH_TO_PLINK_software (required))
```


- GWAS_SUMMARY_STATISTICS (required):
Prepare the summary statistics in the following format (including the header line):
```
   chr        rsid      a1    a2       z          p
    1      rs4040617    G     A     -0.199     8.42e-01
    1      rs4075116    C     T      0.646     5.18e-01
    1      rs9442385    T     G     -0.016     9.87e-01
    ...
```

where rsid is the SNP, a1 is the reference allele, a2 is the alternative allele, z is the z scores of the a2 allele, p is the p-value of the effect. 


## Output

The `run.owas` function returns a data.frame with 8 columns:

`gene`: Gene name

`id`: Index of the segment on that gene

`snplen`: The number of SNPs on the segment

`z_g`: z statistics

`p_g`: p-value

`chr`: Chromosome

`start`: Segment start position (hg19)

`end`: Segment end position (hg19)


## A Simplified Pipeline
Clone this repository using the following git command:

$ git clone https://github.com/shuangsong0110/OWAS.git

$ cd ./OWAS

$ wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz

$ tar -xvzf 1000G_Phase3_plinkfiles.tgz

$ wget -O bedfile_5kb.txt https://cloud.tsinghua.edu.cn/f/7524fb35b88c468dbc02/?dl=1 --no-check-certificate

$ wget -O openness.tar.gz https://cloud.tsinghua.edu.cn/f/17d8b0b5de3941a2bf86/?dl=1 --no-check-certificate

$ tar -zxvf openness.tar.gz


Run with R:

```r
install.packages('OWAS_1.1.0.tar.gz')
library(OWAS)
path0 <- getwd()
pre.cov(ldpath = paste0(path0, '/1000G_EUR_Phase3_plink/1000G.EUR.QC.'), 
	bedfile = paste0(path0, '/bedfiles/bedfile_5kb.txt'), 
	path = paste0(path0, '/output/'),
	chr = 22)
             
pre.sc(ctype = 'Gm12878',
	scpath = paste0(path0, '/openness/'), 
	ldpath = paste0(path0, '/1000G_EUR_Phase3_plink/1000G.EUR.QC.'), 
	path = paste0(path0, '/output/'), 
	chr = 22)

res <- run.owas(ctype = 'Gm12878',
	gwas=fread(paste0(path0, '/example/summs.txt')),
	trait='test',
	ldpath = paste0(path0, '/1000G_EUR_Phase3_plink/1000G.EUR.QC.'), 
	path = paste0(path0, '/output/'), 
	chr = 22)
	
head(res)
             
```

## Customized Settings

### Length of segments
We prepared the bed file taking 100 KB up and down-stream from the transcription start sites (TSS) of genes as regulatory regions and divided the regions into segments of 5 KB. We also allow users to specify there own bedfiles for computation. The format of the bedfiles should be (with the header line):

```
   chr        start        end           gene        index
    1       58764864    58769863       DDX11L1         1
    1       58769864    58774863       DDX11L1         2
    1       58774864    58779863       DDX11L1         3
    ...
```

The bedfiles can also be generated with function `write.bed` in our R package. The usage is:

```r
library(OWAS)
ldpath=NULL, path,type.len='len',len=5,type.gene='gene',gene.region=NULL
write.bed(ldpath=PATH_TO_LD_REFERENCE (required), 
	    path=OUTPUT_DIR (required),
	    type.len='len'/'snp' (default = 'len'),
	    len=LENGTH OF SEGMENTS (KB)/ NUMBER OF SNPS (default = 5),
	    type.gene='gene'/'wg' (default=gene),
	    gene.region= GENE TSS (required if type.gene='gene'))                    
```
- PATH_TO_LD_REFERENCE (required): The LD reference plink bfile should be dividied into chromosomes and ended with the number of the chromosome. The input should include the file name but not the exact number of chromosome (e.g., ldpath='path/1000G.EUR.QC.', without .bim/.fam/.bed).

- OUTPUT_DIR (required): path to the output

- type.len: if type.len='len', the function will divide regions into segments with given length (len, KB). If type.len='snp', the function will divide regions into segments with 'len' SNPs in each segment.

- type.gene: if type.gene='gene', the function only focus on gene regulatory region; if type.gene='wg', the funciton will run in the whole genome.

- gene.region: A data.frame with the format below (with no header) , which can be downloaded directly with `wget -O gene.region.txt https://cloud.tsinghua.edu.cn/f/801fcbf2609a4283a26b/?dl=1`

```
      V1          V2          V3        V4    
    chr1         12314      12516     gene1   
    chr1         12872      12900     gene2    
    chr1         31592      34000     gene3    
    ...
```




### Openness scores (or other annotation weights)

We also allow users to customize the openness scores, besides our precomputed 12 cell lines. The used could write the scores in the following format (with each chromosme):
```
      rsid          score         a1    a2    
   rs4040617          2.3         G     A   
   rs4075116          0.6         C     T     
   rs9442385         -0.3         T     G    
    ...
```

The users should specify a name for the annotation (e.g., `newct`), and create a new folder in the `/openness/` file as `/openness/newct/`, and save 22 files as `result_1_newct.txt`, `result_2_newct.txt`, ... , `result_22_newct.txt`.



## Maintainer

Please contact Shuang Song (song-s19@mails.tsinghua.edu.cn) if there are any problems or questions.


## Citations

Please cite:

Song, S., Shan, N., Wang, G., Yan, X., Liu, J. S., & Hou, L. (2021). Openness weighted association studies: leveraging personal genome information to prioritize non-coding variants. Bioinformatics.

