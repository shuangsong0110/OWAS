#' @title Generate LD files
#' @description Using this function to generate blockwise LD matrix
#' @param ldpath path for plink bfiles for the chromosome-wise LD reference;(e.g. 1000 Genome Project)
#' @param bedfile path to the bedfile for OWAS segments (default: 100 KB up and down-stream of the TSS, the length of each segment is 5 KB)
#' @param path specify the output path
#' @param chr chromosome, can be a integer or vector from 1 to 22, default=1:22
#' @param shrink shrink LD matrix or not, default=T
#'
#'
#' @import  EBPRS data.table stats utils corpcor
#' @export
#'
#'

pre.cov <- function(ldpath=NULL,bedfile, path,
                    chr=1:22,shrink=T){
  ## variables
  chr0 <- chr
  path0 <- path
  #if(is.null(path)){path <- getwd()}else{dir.create(path,recursive = T)}
  if(is.null(path)){stop('please input a valid path')}
  if(is.null(bedfile)){stop('please input a valid bedfile path')}
  if(is.null(ldpath)){stop('Please specify a valid path for LD reference.')}
  ## check
  if(!is.null(ldpath)){
    for(chr in chr0){
      if(!file.exists(paste0(ldpath,chr,'.bim'))){
        stop(paste0('The specified LD file: ',ldpath,chr,'.bim does not exist!'))
      }
    }
  }
  for(chr in chr0){
  path <- paste0(path0,"/chr",chr)
  dir.create(path,recursive = T)
  setwd(path)

  bf0 <- as.data.frame(fread(bedfile))
  bf0$chr <- as.numeric(bf0$chr)
  bf <- bf0[which(bf0$chr==chr),]
  #gname <- unique(bf$gene)
  tss0 <- split(bf$start,bf$gene)
  tss1 <- split(bf$end,bf$gene)
  gname <- names(tss0)
  save(list=c("gname","tss0","tss1"),file=paste0("1_gene_ref_chr",chr,".RData"))
  #dir.create(paste0("./ctype/",ctype),recursive = TRUE)
  dir.create("./cov")
  onekg1 <- read_plink(paste0(ldpath,chr))
  onekg <- onekg1$bim
  pos <- onekg[,4]
  ind <- gene.dup <-snplen <- start <- end <- c()
  cov.all <- snp.gene.all <- list()
  k <- 1
  for(i in 1:length(gname)){
    num <- length(tss0[[i]])
    for(j in 1:num){
      snp.gene <- which((pos>=tss0[[i]][j])&(pos<=tss1[[i]][j]))
      if(length(snp.gene)==0){
        print(paste("0 snps in gene",i,",",j))
      }else{
        start[k] <- tss0[[i]][j]
        end[k] <- tss1[[i]][j]
        snp.gene.all[[k]] <- snp.gene
        snplen[k] <- length(snp.gene)
        gene.dup[k] <- gname[i]
        ind[k] <- j
        print(snplen[k])
        if(length(snp.gene)==1){
          covm <- 1
        }else{
          if(shrink){
            covm <- cov.shrink(onekg1$bed[,snp.gene])
          }else{
          covm <- cov(onekg1$bed[,snp.gene])
          }
        }
        cov.all[[k]] <- covm
        k <- k+1
        # write.table(covm,paste0("./cov/",gname[i],"_RE",j,"_cov.txt"),quote=F,
        #             col.names=F,row.names = F)
        cat(as.character(gname[i]),i,"  ",j,"\n")
      }
    }
  }
  temp <- data.frame(gene=gene.dup,ind=ind,snplen=snplen,start=start,end=end)
  temp$order <- 1:nrow(temp)
  write.table(temp,"genelist.txt",quote=F,col.names=T,row.names=F)
  save(cov.all,file=paste0('./cov/cov.RData'))
  save(snp.gene.all,file=paste0('./cov/snp.gene.RData'))
}
}
