# Boundary cases (1 - 2P - Q = 0 in exact arithmetic): does the sign / size of 1 - 2P - Q computed from
# kaks_lwl_counts() give KaKs_Calculator's decision?
kaks_lwl_methods <- function() c("LWL", "MLWL", "LPB", "MLPB")
source("/home/claude/w-olclean/R/lwl_counts.R")
x <- utils::read.delim(commandArgs(TRUE)[1], colClasses = "character", na.strings = character())
bad <- 0
for (r in seq_len(nrow(x))) {
        e <- kaks_lwl_counts(x$s1[r], x$s2[r], x$method[r], as.integer(x$code[r]))
        t <- as.integer(x$target[r]) + 1
        a1 <- 1 - 2 * e$P[t] - e$Q[t]
        a2 <- 1 - 2 * e$Q[t]
        if (x$kc[r] == "sat") ok <- !(a1 > 0 && a2 > 0)
        else if (x$kc[r] == "nonsat") ok <- a1 > 0 && a2 > 0
        else ok <- a1 > 0 && abs(log(1 / a1) - as.numeric(x$kc[r])) < 1e-3 * log(1 / a1)
        if (!ok) { bad <- bad + 1; if (bad < 10) cat("MISMATCH", x$name[r], x$method[r], a1, x$kc[r], "\n") }
}
cat("rows", nrow(x), "mismatches", bad, "\n")
