args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
data <- read.delim(args[1L], check.names = FALSE, row.names = 1)
stopifnot(identical(names(data), paste0("sample_", 1:6)))
beta <- as.matrix(data)
beta[beta < 1e-6] <- 1e-6
beta[beta > 1 - 1e-6] <- 1 - 1e-6
m_value <- qlogis(beta) / log(2)
group <- factor(rep(c("control", "overexpression"), each = 3L))
fit <- limma::eBayes(limma::lmFit(m_value, model.matrix(~ group)), trend = TRUE)
p_value <- fit$p.value[, 2L]
results <- data.frame(probe_id = rownames(beta), m_effect = fit$coefficients[, 2L],
                      p = p_value, q = p.adjust(p_value, method = "BH"))
write.table(results, args[2L], sep = "\t", row.names = FALSE, quote = FALSE)
