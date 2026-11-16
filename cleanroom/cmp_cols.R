# raw per-column identity (no normalisation at all) and the list of attributes that differ before normalisation
.libPaths(c("/home/claude/w-olclean-scratch/lib-base", "/home/claude/w-olfix3/lib-coalr", .libPaths()))
suppressMessages(library(data.table))
ok <- TRUE
for (f in sort(Sys.glob("res/base_*.rds"))) {
  rx <- readRDS(f)$r; ry <- readRDS(sub("res/base_", "res/new_", f))$r
  same <- vapply(names(rx), function(n) identical(rx[[n]], ry[[n]]), TRUE)
  ok <- ok && all(same) && identical(names(rx), names(ry))
  an <- union(names(attributes(rx)), names(attributes(ry)))
  da <- an[!vapply(an, function(n) identical(attr(rx, n), attr(ry, n)), TRUE)]
  cat(sprintf("COLS %s: %d columns, all raw-identical=%s; attributes differing before path/date normalisation: %s\n",
              sub("res/base_", "", f), length(same), all(same), paste(da, collapse = ",")))
}
cat("ALL_COLUMNS_IDENTICAL", ok, "\n")
