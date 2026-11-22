suppressPackageStartupMessages(library(minfi))
args <- commandArgs(trailingOnly=TRUE)
work <- args[[1]]
idats <- sort(list.files(args[[2]], pattern="_Grn.idat.gz$", full.names=TRUE))
rg <- read.metharray(sub("_Grn.idat.gz$", "", idats))
raw <- preprocessRaw(rg)
cat("raw loaded\n")
ann <- as.data.frame(getAnnotation(raw))
cat(paste(names(ann),collapse="\n"), "\n")
ann$probe_id <- rownames(ann)
data.table::fwrite(ann, file.path(work,"epic_annotation.tsv.gz"), sep="\t")
noob <- preprocessNoob(rg, offset=15, dyeCorr=TRUE, dyeMethod="single")
set.seed(20261120)
swan <- preprocessSWAN(rg, mSet=raw)
set.seed(20261120)
combined <- preprocessSWAN(rg, mSet=noob)
names <- rownames(raw)
gg <- data.table::data.table(probe_id=names)
for (method in c("SWAN","ssNoob_SWAN","ssNoob_0")) {
  set <- switch(method, SWAN=swan, ssNoob_SWAN=combined, ssNoob_0=noob)
  beta <- getBeta(set, offset=if (method=="ssNoob_0") 0 else 100)[names,,drop=FALSE]
  for (i in seq_len(ncol(beta))) gg[[paste0(method,"_",i)]] <- beta[,i]
}
data.table::fwrite(gg,file.path(work,"swan_betas.tsv.gz"),sep="\t")
