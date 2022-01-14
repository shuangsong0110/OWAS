#' @title Generate openness score files
#' @description Using this function to generate blockwise openness scores
#' @param ctype celltype
#' @param scpath path to the openness scores files
#' @param ldpath path for plink bfiles for the chromosome-wise LD reference;(e.g. 1000 Genome Project)
#' @param path specify the output path
#' @param chr chromosome, can be a integer or vector from 1 to 22, default=1:22
#'
#'
#' @import  EBPRS data.table stats utils
#' @export
#'
#'
pre.sc <- function(ctype,scpath,ldpath=NULL,path,chr=1:22){
  path0 <- path
  chr0 <- chr
  for(chr in chr0){
    path <- paste0(path0,"/chr",chr)
    if(is.null(ldpath)){stop('Please specify a valid path for LD reference.')}

    if(!file.exists(paste0(path,'/cov/snp.gene.RData'))){
      stop('Please use pre.cov to generate covariance matrix first!')
    }
    if(is.null(scpath)){
      stop('Please specify the scpath (PATH to openness scores files)')
    }
    # dir.create(path,recursive = T)
    setwd(path)
    dir.create(paste0("./ctype/",ctype),recursive = T)
    load('./cov/snp.gene.RData')
    sc2 <- fread(paste0(scpath,"/",ctype,"/result_",chr,"_",ctype,".txt"))
    colnames(sc2) <- c('rsid','score')
    sc2$rsid <- paste0('r',sc2$rsid)
    # library(EBPRS)
    onekg1 <- fread(paste0(ldpath,chr,'.bim'))
    onekg1$order <- 1:nrow(onekg1)
    onekg2 <- merge(sc2,onekg1,by.x='rsid',by.y='V2',all.y=T)
    onekg2$score[which(is.na(onekg2$score))] <- 0

    # sign <- agtc(onekg2$V5,onekg2$V6,onekg2$a1,onekg2$a2)
    # onekg2$score <- onekg2$score*sign
    sc2 <- onekg2[order(onekg2$order),c('rsid','score')]
    print(length(snp.gene.all))
    sg.list <- Map(list,snp.gene.all,1:length(snp.gene.all))
    temp.dat <- data.frame(onekg1[,c(2,5,6)],sc2$score)
    try1 <- function(x,y){
      set <- unlist(x[1])
      j <-  unlist(x[2])
      pf2 <- y[set,]
      pf2$bl <- rep(j,length(set))
      pf2$bl.order <- 1:length(set)
      stats <- pf2
      return(stats)
    }
    temp.res <- lapply(sg.list,try1,y=temp.dat)
    stats <- rbindlist(temp.res)
    colnames(stats) <- c('SNP','A1','A2','sc','bl','bl.order')
    save(stats,file=paste0("./ctype/",ctype,"/stats.RData"))
  }
}
