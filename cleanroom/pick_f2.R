kaks_lwl_methods <- function() c("LWL", "MLWL", "LPB", "MLPB")
source("/home/claude/w-olclean/R/lwl_counts.R")
f_sum <- kaks_lwl_counts
e_sum <- new.env(parent = environment(kaks_lwl_counts))
e_sum$Reduce <- function(f, x, init) sum(x)   # the only Reduce() in kaks_lwl_counts() is the running total
environment(f_sum) <- e_sum
x <- utils::read.delim("zero_all.tsv", colClasses = "character", na.strings = character())
dec <- function(e, t) { a1 <- 1 - 2 * e$P[t] - e$Q[t]; a2 <- 1 - 2 * e$Q[t]; if (a1 > 0 && a2 > 0) "nonsat" else "sat" }
out <- list(); nflip <- 0; nok <- 0
for (r in seq_len(nrow(x))) {
  t <- as.integer(x$target[r]) + 1
  kc <- if (x$kc[r] == "sat") "sat" else "nonsat"
  a <- dec(kaks_lwl_counts(x$s1[r], x$s2[r], x$method[r], as.integer(x$code[r])), t)
  s <- dec(f_sum(x$s1[r], x$s2[r], x$method[r], as.integer(x$code[r])), t)
  nok <- nok + (a == kc)
  if (a == kc && s != kc) { nflip <- nflip + 1; out[[length(out) + 1]] <- data.frame(x[r, c("name", "code", "method", "target")], kc = kc, len = nchar(x$s1[r]), s1 = x$s1[r], s2 = x$s2[r]) }
}
o <- do.call(rbind, out)
cat("rows", nrow(x), "Reduce agrees with KC", nok, "sum() flips", nflip, "\n")
print(table(o$method, o$target, o$kc))
o <- o[order(o$len), ]
pick <- do.call(rbind, lapply(split(o, paste(o$method, o$target, o$kc)), function(g) g[1, ]))
pick <- pick[order(pick$method, pick$target, pick$kc), ]
print(pick[, c("name", "code", "method", "target", "kc", "len")])
utils::write.table(pick[, c("code", "method", "target", "kc", "s1", "s2")], "f2_fixture.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
