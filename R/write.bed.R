#' @title Custom bed files
#' @description Using this function to generate blockwise LD matrix
#' @param ldpath path for plink bfiles for the chromosome-wise LD reference;(e.g. 1000 Genome Project)
#' @param path specify the output path
#' @param type.len 'len'/'snp'
#' @param len length (kb/number of SNPs)
#' @param type.gene 'gene'/'wg': by gene TSS or whole genome
#' @param gene.region regulatory region of the gene (a data.frame with columns: V1(chr),V2(start),V3(end),V4(gene name))
#'
#'
#' @import  EBPRS data.table stats utils
#' @export
#'
#'
write.bed <- function(ldpath=NULL, path,type.len='len',len=5,type.gene='gene',gene.region=NULL){
  if(is.null(ldpath)){stop('Please specify a valid path for LD reference.')}
  if(type.gene=='gene'){
    if(is.null(gene.region)){
      stop('Please specify the gene region.')
    }else{
      gene00 <- gene.region
      if(nchar(gene00$V1[1])<4){
        gene00$V1 <- paste0('chr',gene00$V1)
      }
    }
  }
  chr <- c()
  start <- c()
  end <- c()
  gene <- c()
  index <- c()
  if(type.gene=='gene'){
    if(type.len=='len'){
      if(len>100){
        stop('len should be less than 100 (kb)!')
      }
      len <- len*1000
      for(j in 1:nrow(gene00)){
        print(j)
        # if(j==1){
        # onekg1 <- fread(paste0(ldpath,gene$V1[j],'.bim'))
        # }else if(gene$V1[j]!=gene$V1[j-1]){
        #   onekg1 <- fread(paste0(ldpath,gene$V1[j],'.bim'))
        # }
        tss0 <- seq(gene00$V2[j],gene00$V3[j]-1,len)
        tss1 <- tss0+len-1
        #tss1 <- pos[seq(1,length(ind)-1,len)+1]
        #tss1 <- tss1-1
        #tss0[2:length(tss0)] <-  tss0[2:length(tss0)] -1
        tss1[length(tss1)] <- gene00$V3[j]
        chr <- c(chr,rep(gene00$V1[j],length(tss1)))
        start <- c(start,tss0)
        end <- c(end,tss1)
        gene <- c(gene,rep(gene00$V4[j],length(tss1)))
        index <- c(index,1:length(tss1))
      }
      chr <- unlist(lapply(strsplit(chr,'hr'),function(x){x[2]}))
    }else if(type.len=='snp'){
      chr.gene <- unlist(lapply(strsplit(gene00$V1,'hr'),function(x){x[2]}))
      for(j in 1:nrow(gene00)){
        if(j==1){
          onekg1 <- fread(paste0(ldpath,chr.gene[j],'.bim'))
        }else if(gene00$V1[j]!=gene00$V1[j-1]){
          print(j)
          onekg1 <- fread(paste0(ldpath,chr.gene[j],'.bim'))
        }
        ind <- which((paste0(onekg1$V1)==chr.gene[j])&
                       (onekg1$V4<gene00$V3[j])&(onekg1$V4>gene00$V2[j]))
        if(length(ind)>1){
          pos <- sort(onekg1$V4[ind])
          tss0 <- pos[seq(1,length(ind)-1,len)]
          tss1 <- tss0[-1]
          tss1 <- tss1-1
          tss1 <- c(tss1,pos[length(tss0)])
          print(length(tss0))
          print(length(tss1))
        }else if(length(ind)==1){
          tss0 <- gene00$V2[j]
          tss1 <- gene00$V3[j]
        }else{
          tss0 <- c()
          tss1 <- c()
        }
        chr <- c(chr,rep(gene00$V1[j],length(tss1)))
        start <- c(start,tss0)
        end <- c(end,tss1)
        gene <- c(gene,rep(gene00$V4[j],length(tss1)))
        index <- c(index,1:length(tss1))
      }
      chr <- unlist(lapply(strsplit(chr,'hr'),function(x){x[2]}))

    }

  }else if (type.gene=='wg'){
    if(type.len=='len'){
      if(len>100){
        stop('len should be less than 100 (kb)!')
      }
      len <- len*1000

      for(k in 1:22){
        print(k)
        onekg1 <- fread(paste0(ldpath,k,'.bim'))
        bim <- onekg1[onekg1$V1==k,]
        tss0 <- seq(bim$V4[1],bim$V4[nrow(bim)-1],len)
        tss1 <- tss0+len-1
        #tss0[2:length(tss0)] <-  tss0[2:length(tss0)] -1
        tss1[length(tss1)] <- bim$V4[nrow(bim)]
        chr <- c(chr,rep(k,length(tss1)))
        start <- c(start,tss0)
        end <- c(end,tss1)
        gene <- c(gene,rep(paste0('chr',k),length(tss1)))
        index <- c(index,1:length(tss1))
      }
    }else if(type.len=='snp'){
      for(k in 1:22){
        print(k)
        onekg1 <- fread(paste0(ldpath,k,'.bim'))
        bim <- onekg1[onekg1$V1==k,]
        ind <- which(onekg1$V1==k)
        #pos <- sort(bim$V4)
        if(length(pos)>1){
          pos <- sort(onekg1$V4[ind])
          tss0 <- pos[seq(1,length(ind)-1,len)]
          tss1 <- tss0[-1]
          tss1 <- tss1-1
          tss1 <- c(tss1,pos[length(tss0)])
        }else if(length(ind)==1){
          tss0 <- gene$V2[j]
          tss1 <- gene$V3[j]
        }else{
          tss0 <- c()
          tss1 <- c()
        }
        chr <- c(chr,rep(k,length(tss1)))
        start <- c(start,tss0)
        end <- c(end,tss1)
        gene <- c(gene,rep(paste0('chr',k),length(tss1)))
        index <- c(index,1:length(tss1))
      }
    }
  }
  dat <- data.frame(chr,start,end,gene,index)
  write.table(dat,paste0(path,'/bed_custom.txt'),quote=F,row.names = F,col.names=T)
  return(dat)
}
