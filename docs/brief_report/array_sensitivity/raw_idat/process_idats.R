args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
idat_dir <- args[1]
positions <- read.delim(args[2], colClasses = "character")
out <- args[3]
scratch <- args[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
stopifnot(identical(as.integer(positions$sample), 1:6))
positions$label <- c("Z1", "Z2", "Z3", "D1", "D2", "D3")
paths <- file.path(idat_dir, paste(positions$geo_accession, positions$chip_id, positions$position, sep = "_"))
for (i in seq_along(paths)) {
    for (color in c("Grn", "Red")) {
        idat <- illuminaio::readIDAT(paste0(paths[i], "_", color, ".idat.gz"))
        stopifnot(idat$Barcode == positions$chip_id[i], idat$Unknowns$MostlyA == positions$position[i])
        times <- idat$RunInfo[idat$RunInfo[, "BlockType"] == "Scan", "RunTime"]
        times <- as.POSIXct(times, format = "%m/%d/%Y %I:%M:%S %p", tz = "UTC")
        stopifnot(length(times) > 0L, !anyNA(times))
        if (color == "Red") {
            positions$scan_first[i] <- format(min(times), "%Y-%m-%d %H:%M:%S")
            positions$scan_last[i] <- format(max(times), "%Y-%m-%d %H:%M:%S")
            positions$scan_records[i] <- length(times)
        }
    }
}
data.table::fwrite(positions, file.path(out, "verified_idat_mapping.tsv"), sep = "\t")
cat("Verified six sample mappings in both IDAT channels.\n")

raw_ew <- ewastools::read_idats(paths, quiet = TRUE)
stopifnot(identical(as.character(raw_ew$meta$position), positions$position),
          identical(as.character(raw_ew$meta$barcode), positions$chip_id))
metrics <- ewastools::control_metrics(raw_ew)
keep_cg <- raw_ew$manifest$probe_type == "cg"
autosome_ew <- keep_cg & raw_ew$manifest$chr %in% paste0("chr",1:22)
negative <- which(raw_ew$controls$group == "NEGATIVE")
stopifnot(sum(autosome_ew) > 800000L, length(negative) > 0L)
raw_ew <- ewastools::detectionP(raw_ew)
fail_ew <- !is.finite(raw_ew$detP) | raw_ew$detP > .01

qc <- positions[c("sample","label","group","position")]
qc[["Negative control median green"]] <- apply(raw_ew$ctrlG[negative,,drop=FALSE], 2, median, na.rm=TRUE)
qc[["Negative control median red"]] <- apply(raw_ew$ctrlR[negative,,drop=FALSE], 2, median, na.rm=TRUE)
qc[["Negative control 95th percentile green"]] <- apply(raw_ew$ctrlG[negative,,drop=FALSE], 2, quantile, .95, na.rm=TRUE)
qc[["Negative control 95th percentile red"]] <- apply(raw_ew$ctrlR[negative,,drop=FALSE], 2, quantile, .95, na.rm=TRUE)
qc[["Raw median M"]] <- apply(raw_ew$M[autosome_ew,,drop=FALSE], 2, median, na.rm=TRUE)
qc[["Raw median U"]] <- apply(raw_ew$U[autosome_ew,,drop=FALSE], 2, median, na.rm=TRUE)
qc[["ewastools detection failures autosomal (%)"]] <- 100 * colMeans(fail_ew[autosome_ew,,drop=FALSE])
qc[["ewastools detection failures all CpGs (%)"]] <- 100 * colMeans(fail_ew[keep_cg,,drop=FALSE])
ew_ids <- raw_ew$manifest$probe_id
stopifnot(length(ew_ids)==nrow(fail_ew), !anyDuplicated(ew_ids))
ew_pass <- setNames(rowSums(fail_ew) == 0L, ew_ids)
long <- do.call(rbind,lapply(names(metrics),function(n) {
    value <- metrics[[n]]
    threshold <- attr(value,"threshold")
    stopifnot(length(value)==6L, length(threshold)==1L)
    data.frame(sample=positions$sample,metric=n,value=as.numeric(value),threshold=threshold,
               failed=ifelse(is.finite(value),value<threshold,NA))
}))
qc[["Control metrics below cutoff"]] <- as.integer(table(factor(long$sample[which(long$failed)],levels=positions$sample)))
data.table::fwrite(long, file.path(out,"control_metrics_long.tsv"),sep="\t",na="NA")
rm(raw_ew,metrics,fail_ew); gc()

rg <- minfi::read.metharray(paths, force=FALSE)
raw <- minfi::preprocessRaw(rg)
probe_ids <- rownames(raw)
ann <- minfi::getAnnotation(raw)
auto <- probe_ids %in% rownames(ann)[ann$chr %in% paste0("chr",1:22)]
stopifnot(sum(auto)>800000L)
minfi_p <- minfi::detectionP(rg)[probe_ids,,drop=FALSE]
fail_minfi <- !is.finite(minfi_p) | minfi_p > .01
qc[["minfi detection failures autosomal (%)"]] <- 100 * colMeans(fail_minfi[auto,,drop=FALSE])
qc[["minfi detection failures all CpGs (%)"]] <- 100 * colMeans(fail_minfi[startsWith(probe_ids,"cg"),,drop=FALSE])
data.table::fwrite(qc, file.path(out,"array_qc.tsv"),sep="\t",na="NA")
cat("Raw array QC exported. Running ssNoob.\n")

beta100 <- function(set){
    m <- minfi::getMeth(set)
    u <- minfi::getUnmeth(set)
    stopifnot(identical(rownames(m),probe_ids),identical(rownames(u),probe_ids))
    m / (m + u + 100)
}
background <- minfi::preprocessNoob(rg, offset=15, dyeCorr=FALSE, dyeMethod="single")
dye <- minfi::preprocessNoob(rg, offset=15, dyeCorr=TRUE, dyeMethod="single")
tag <- function(x, label) {
    colnames(x) <- paste0(label,"_",positions$sample)
    x
}
beta <- data.frame(probe_id = probe_ids,
    pass_ewastools = unname(ew_pass[probe_ids]),
    pass_minfi = rowSums(fail_minfi)==0L,
    tag(beta100(raw),"IDAT_raw"),tag(beta100(background),"Noob_background"),tag(beta100(dye),"ssNoob_dye"),
    check.names=FALSE)
data.table::fwrite(beta, file.path(scratch,"idat_beta.tsv.gz"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out,"R_session.txt"))
cat("Raw, background-only Noob, and ssNoob beta matrices exported.\n")
