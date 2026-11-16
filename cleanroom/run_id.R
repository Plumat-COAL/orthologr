# Rscript run_id.R <lib> <tag> <ds:cfg:cores> ...   -> res/<tag>_<ds>_<cfg>_c<cores>.rds
# (driver adapted from /home/claude/w-olfix3/scripts/run.R); counts calls of kaks_lwl_counts() in res/<...>.calls
a <- commandArgs(TRUE)
lib <- a[1]; tag <- a[2]; jobs <- a[-(1:2)]
.libPaths(c(lib, "/home/claude/w-olfix3/lib-coalr", .libPaths()))
suppressMessages(library(orthologr))
S <- system.file("seqs", package = "orthologr")
Y <- "/home/claude/og-work/data/yeast_cds"
for (j in jobs) {
  p <- strsplit(j, ":")[[1]]; ds <- p[1]; cfg <- p[2]; cores <- as.integer(p[3])
  files <- switch(ds,
    d20   = c(file.path(S, "ortho_thal_cds.fasta"), file.path(S, "ortho_lyra_cds.fasta")),
    d1000 = c(file.path(S, "ortho_thal_cds_1000.fasta"), file.path(S, "ortho_lyra_cds_1000.fasta")),
    ScCg  = file.path(Y, c("Scer.fna", "Cgla.fna")),
    ScKl  = file.path(Y, c("Scer.fna", "Klac.fna")))
  args <- list(query_file = files[1], subject_file = files[2], comp_cores = cores, print_citation = FALSE)
  if (startsWith(cfg, "kc_")) args$dnds_est.method <- sub("^kc_", "", cfg)
  else if (cfg == "coal") args$engine <- "coal"
  else if (cfg == "coal_LWL") { args$engine <- "coal"; args$dnds_est.method <- "LWL" }
  out <- file.path("res", sprintf("%s_%s_%s_c%d.rds", tag, ds, cfg, cores))
  calls <- sub("rds$", "calls", out); unlink(calls)
  ns <- asNamespace("orthologr")
  if (exists("kaks_lwl_counts", envir = ns))
    suppressMessages(trace("kaks_lwl_counts", tracer = bquote(cat("call\n", file = .(calls), append = TRUE)), where = ns, print = FALSE))
  t <- system.time(r <- tryCatch(suppressWarnings(suppressMessages(do.call(dNdS, args))), error = function(e) e))
  if (exists("kaks_lwl_counts", envir = ns)) suppressMessages(untrace("kaks_lwl_counts", where = ns))
  n <- if (file.exists(calls)) length(readLines(calls)) else 0
  cat(sprintf("RUN %s %s %s c%d elapsed=%.1f rows=%s lwl_counts_calls=%d\n", tag, ds, cfg, cores, t[["elapsed"]],
              if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else nrow(r), n))
  saveRDS(list(r = r, elapsed = t[["elapsed"]]), out)
}
