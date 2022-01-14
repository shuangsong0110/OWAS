#' @title Main function
#' @description Run OWAS
#' @param ctype celltype
#' @param gwas GWAS summary statistics, need to include rsid, chr, a1, a2, z, p (header is necessary)
#' @param trait trait name
#' @param ldpath path for plink bfiles for the chromosome-wise LD reference;(e.g. 1000 Genome Project)
#' @param path specify the output path
#' @param chr chromosome, can be a integer or vector from 1 to 22, default=1:22
#' @param plinkpath chromosome, can be a integer or vector from 1 to 22, default=1:22
#'
#'
#' @import  EBPRS data.table stats utils
#' @export
#'
#'
run.owas <- function(ctype,gwas,trait='test',ldpath,path,plinkpath,chr=1:22,clump=T){
  path0 <- path
  chr0 <- chr
  if(is.null(ldpath)){stop('Please specify a valid path for LD reference.')}
  for(chr in chr0){
    print(chr)
   setwd(paste0(path0,'/chr',chr))
  if(clump){
    cldat <- gwas[,c('rsid','p')]
    colnames(cldat) <- c('SNP','P')
    write.table(cldat,'clump.txt',quote=F,row.names = F,col.names = T)
    system(paste0(plinkpath,' --bfile ',ldpath,chr,' --clump clump.txt  --clump-p1 1 --clump-p2 1 --clump-r2 0.99 --clump-kb 5 --out clump --noweb'))
    clumped_p <- fread(paste0("clump.clumped"),head=T)
    set <- which(gwas$rsid%in%clumped_p$SNP)
    gwas <- gwas[set,]
    system('rm clump.txt')
    system('rm clump.clumped')
  }
      name <- fread('genelist.txt')
      dir.create(paste0("./ctype/",ctype),recursive = T)
      if(!file.exists('./cov/snp.gene.RData')){
        stop('Please run pre.cov first!')
      }
      if(!file.exists('./cov/cov.RData')){
        stop('Please run pre.cov first!')
      }
      if(!file.exists(paste0('./ctype/',ctype,'/stats.RData'))){
        stop('Please run pre.sc first!')
      }
      load('./cov/snp.gene.RData')
      load('./cov/cov.RData')
      load(paste0('./ctype/',ctype,'/stats.RData'))
      n.bl <- nrow(name)
      temp.ind <- 1:length(snp.gene.all)
      ff1 <- function(x,y){return(rep(y,length(x)))}
      bl <- unlist(mapply(ff1,snp.gene.all,temp.ind))
      ff2 <- function(x){return(1:length(x))}
      bl.order <- unlist(lapply(snp.gene.all,ff2))
      stats$bl <- bl
      stats$bl.order <- bl.order
      pp <- which(gwas$chr==chr)
      gwas1 <- gwas[pp,]
      if(nrow(gwas1)==0){
        stop(paste0('No SNPs on chromosome ', chr, ' in the provided summary statistics. Please check.'))
      }
      ldbim <- fread(paste0(ldpath,chr,'.bim'))
      gwas1 <- merge(gwas1,ldbim,by.x='rsid',by.y='V2')
      sign <- agtc(gwas1$a1,gwas1$a2,gwas1$V5,gwas1$V6)
      gwas1$z <- gwas1$z*sign

      gwas1 <- gwas1[,c('rsid','chr','a1','a2','z','p')]
      sc.all <- merge(gwas1,stats,by.x="rsid",by.y='SNP')
      sc.all <- sc.all[order(sc.all$bl,sc.all$bl.order),]
      sc.all$order <- 1:nrow(sc.all)
      sc.list <- split(sc.all$bl.order,sc.all$bl)
      sc.order <- split(sc.all$order,sc.all$bl)
      ind.temp <- as.character(1:n.bl)
      sc.list1 <- sc.list[ind.temp]
      sc.order1 <- sc.order[ind.temp]
      try2 <- function(x1,or1,or2){
        if(!is.null(or1)){
          x <- as.matrix(x1)[or1,or1]
          y <- sc.all$sc[or2]
          z <- sc.all$z[or2]
          if(T){
            if(length(y)>1){
              x.sd <- sqrt(diag(x))
              x.r <- cov2cor(x)
              #x.r1 <- (eigen(x.r)$vectors)%*%diag(eigen(x.r)$values+0.1)%*%t(eigen(x.r)$vectors)
              x.r1 <- x.r*0.99+0.01*diag(nrow(x.r))
              x.r2 <- (x.r1+t(x.r1))/2
              xx <- cor2cov(x.r2, sd=x.sd)
              x <- xx
            }
          }
          temp.yxy <- t(y)%*% x %*% y
          if(temp.yxy<0){temp.yxy <- 0}
          sigma_g <- sqrt(temp.yxy)
          if(sigma_g==0){
            return(list(sigma_g=0,z_g=0,p_g=1))
          }else{
            sigma_l <- sqrt(diag(x))
            z_l <- z
            z_g <- sum(z_l*sigma_l*y)/sigma_g
            p_g <- pnorm(-abs(z_g))*2
            return(list(sigma_g=sigma_g,z_g,p_g))
          }}else{
            return(list(sigma_g=0,z_g=0,p_g=1))
          }
      }
      res <-  mapply(try2,cov.all,sc.list1,sc.order1)
      gene_stats <- as.data.frame(name)
      gene_stats$sigma_g <- unlist( res[1,])
      gene_stats$z_g <- unlist(res[2,])
      gene_stats$p_g <- unlist(res[3,])
      gene_stats$chr <- chr
      gene_stats.sig <- gene_stats[gene_stats$p_g<0.05,]
      dir.create(paste0(path0,"/ctype/",ctype,"/chr",chr,"/",trait),recursive = T)
      write.table(gene_stats,paste0(path0,"/ctype/",ctype,"/chr",chr,"/",trait,
                                    "/1_gene_stats.txt"),quote=F,col.names = T,row.names=F)
      write.table(gene_stats.sig,paste0(path0,"/ctype/",ctype,"/chr",chr,"/",trait,
                                        "/1_gene_stats_sig.txt"),quote=F,col.names = T,row.names=F)
      if(chr ==chr0[1]){
          res.sig <- gene_stats.sig
      }else{
        res.sig <- rbind(res.sig,gene_stats.sig)
      }
  }
  res.sig <- res.sig[order(res.sig$p_g),c('gene','ind','snplen','z_g','p_g','chr','start','end')]
  return(res.sig)
}
