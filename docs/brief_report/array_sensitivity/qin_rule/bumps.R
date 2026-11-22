suppressPackageStartupMessages(library(bumphunter))
suppressPackageStartupMessages(library(data.table))
d <- commandArgs(TRUE)[[1]]
beta <- fread(file.path(d,"bumps_input.tsv.gz"))
groups <- c(0,0,0,1,1,1)
design <- model.matrix(~groups)
# Enumerate all assignments, including the observed allocation. Each design has
# 36 within-group permutations, all of which fit the same two-level contrast.
assignments <- combn(1:6,3)
perms <- apply(assignments,2,function(control) c(control,setdiff(1:6,control)))
clusters <- clusterMaker(beta$chrnum,beta$pos,assumeSorted=TRUE,maxGap=500)
for (method in c("GEO","ssNoob_SWAN","ssNoob_dye")) {
  cat(method, "\n"); flush.console()
  methyl <- as.matrix(beta[,paste0(method,"_",1:6),with=FALSE])
  result <- bumphunterEngine(methyl,design,chr=beta$chrnum,pos=beta$pos,cluster=clusters,
                 cutoff=.05,maxGap=500,nullMethod="permutation",smooth=FALSE,
                 permutations=perms,B=ncol(perms))
  dt <- as.data.table(result$table)
  stopifnot(nrow(dt)>0L)
  dt[,genes:=vapply(seq_len(.N),function(i) {
    vals <- beta$UCSC_RefGene_Name[seq.int(indexStart[i],indexEnd[i])]
    paste(sort(unique(unlist(strsplit(vals[!is.na(vals)],";",fixed=TRUE)))),collapse=";")
  },FUN.VALUE=character(1L))]
  fwrite(dt,file.path(d,paste0("bumps_",method,".tsv")),sep="\t")
}
