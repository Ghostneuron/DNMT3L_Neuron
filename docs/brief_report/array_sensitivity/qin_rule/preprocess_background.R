suppressPackageStartupMessages(library(minfi))
args <- commandArgs(trailingOnly=TRUE)
work <- args[[1]]
idats <- sort(list.files(args[[2]], pattern="_Grn.idat.gz$", full.names=TRUE))
rg <- read.metharray(sub("_Grn.idat.gz$", "", idats))
raw <- preprocessRaw(rg)
cat("raw loaded\n")
noob <- preprocessNoob(rg, offset=15, dyeCorr=FALSE)
set.seed(20261120)
swan <- preprocessSWAN(rg, mSet=raw)
set.seed(20261120)
combined <- preprocessSWAN(rg, mSet=noob)
names <- rownames(raw)
gg <- data.table::data.table(probe_id=names)
for (method in c("Noob_background_SWAN")) {
  set <- combined
  beta <- getBeta(set, offset=if (method=="ssNoob_0") 0 else 100)[names,,drop=FALSE]
  for (i in seq_len(ncol(beta))) gg[[paste0(method,"_",i)]] <- beta[,i]
}
data.table::fwrite(gg,file.path(work,"background_swan_betas.tsv.gz"),sep="\t")
