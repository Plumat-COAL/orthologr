# Run test-lwl_counts.R with kaks_lwl_counts() patched to use sum() for the running totals: it must fail.
.libPaths(c("/home/claude/w-olclean-scratch/lib-new", .libPaths()))
library(testthat); library(orthologr); ns <- asNamespace("orthologr")
f <- get("kaks_lwl_counts", envir = ns); e <- new.env(parent = ns); e$Reduce <- function(f, x, init) sum(x)
environment(f) <- e; unlockBinding("kaks_lwl_counts", ns); assign("kaks_lwl_counts", f, envir = ns)
r <- as.data.frame(test_file("/home/claude/w-olclean/tests/testthat/test-lwl_counts.R", package = "orthologr", load_package = "installed", reporter = "silent"))
cat("WITH-SUM failed", sum(r$failed), "by test:\n"); print(r[r$failed > 0, c("test", "failed")])
