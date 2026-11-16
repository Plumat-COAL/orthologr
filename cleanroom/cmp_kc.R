# Compare kaks_lwl_counts() (sourced from the branch) with the counts KaKs_Calculator printed, as printed text.
args <- commandArgs(TRUE)
kaks_lwl_methods <- function() c("LWL", "MLWL", "LPB", "MLPB")
source("/home/claude/w-olclean/R/lwl_counts.R")
x <- utils::read.delim(args[1], colClasses = "character", na.strings = character())
g <- function(v) formatC(v, digits = 6, format = "g")
fold <- function(v, na_when_zero) if (na_when_zero && all(v == 0)) "NA" else paste(sub("^\\s+", "", g(v)), collapse = ":")
bad <- 0
nu <- 0
t0 <- Sys.time()
for (r in seq_len(nrow(x))) {
        e <- kaks_lwl_counts(x$s1[r], x$s2[r], x$method[r], as.integer(x$code[r]))
        if (grepl("[Uu]", x$s1[r]) || grepl("[Uu]", x$s2[r])) { nu <- nu + 1; if (!is.null(e)) bad <- bad + 1; next }
        got <- c(fold(e$L, FALSE), fold(e$Si, TRUE), fold(e$Vi, TRUE))
        want <- c(x$L[r], x$S[r], x$V[r])
        if (!identical(got, want)) {
                bad <- bad + 1
                if (bad <= 10) cat("MISMATCH", x$name[r], x$method[r], "R:", got, " KC:", want, "\n")
        }
}
cat("rows", nrow(x), "of which with U (expected NULL)", nu, "mismatches", bad, "seconds", format(as.numeric(Sys.time() - t0, units = "secs")), "\n")
