# Rscript fn_out.R <lib> <out.rds>: kaks_lwl_counts() of the installed package on the test inputs (post-freeze check)
a <- commandArgs(TRUE)
.libPaths(c(a[1], .libPaths()))
f <- get("kaks_lwl_counts", envir = asNamespace("orthologr"))
x <- rbind(utils::read.delim("rand1.tsv", colClasses = "character", na.strings = character())[, c("s1", "s2", "method", "code")],
           utils::read.delim("rand2.tsv", colClasses = "character", na.strings = character())[, c("s1", "s2", "method", "code")],
           utils::read.delim("zero_all.tsv", colClasses = "character", na.strings = character())[, c("s1", "s2", "method", "code")])
res <- lapply(seq_len(nrow(x)), function(r) f(x$s1[r], x$s2[r], x$method[r], as.integer(x$code[r])))
saveRDS(res, a[2])
cat("done", length(res), "\n")
