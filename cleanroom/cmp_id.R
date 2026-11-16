# For every res/base_*.rds with a res/new_*.rds: are the result tables identical (all columns, bitwise), and all
# attributes after normalising library paths, temp paths, dates and timings?
.libPaths(c("/home/claude/w-olclean-scratch/lib-base", "/home/claude/w-olfix3/lib-coalr", .libPaths()))
suppressMessages(library(data.table))
norm <- function(x) { x <- gsub("/home/claude/w-olclean-scratch/lib-(base|new)", "LIB", x)
  gsub("(/tmp|/home/claude/w-olclean-scratch/tmp)/[^ '\"]*", "TMP", x) }
walk <- function(o) {
  if (is.character(o)) { a <- attributes(o); o <- norm(o); attributes(o) <- a; return(o) }
  if (is.list(o)) { a <- attributes(o); o2 <- lapply(o, walk); if (!is.null(a)) attributes(o2) <- lapply(a, walk); return(o2) }
  a <- attributes(o); if (!is.null(a)) attributes(o) <- lapply(a, walk); o }
drop_time <- function(r) { o <- attr(r, "orthologr_run"); if (!is.null(o)) { o$date <- NULL; o$seconds <- NULL; o$timing <- NULL; attr(r, "orthologr_run") <- o }
  x <- attr(r, "coal"); if (!is.null(x)) { attr(x, "seconds") <- NULL; attr(r, "coal") <- x }; r }
all_ok <- TRUE
for (f in sort(Sys.glob("res/base_*.rds"))) {
  g <- sub("res/base_", "res/new_", f); if (!file.exists(g)) next
  rx <- readRDS(f)$r; ry <- readRDS(g)$r
  cols <- identical(names(rx), names(ry)) && all(vapply(names(rx), function(n) identical(rx[[n]], ry[[n]]), TRUE))
  full <- identical(walk(drop_time(rx)), walk(drop_time(ry)))
  all_ok <- all_ok && cols && full
  cat(sprintf("CMP %s rows %d/%d columns_bitwise_identical=%s everything_identical=%s\n", sub("res/base_", "", f), NROW(rx), NROW(ry), cols, full))
  if (!full) { ax <- attributes(walk(drop_time(rx))); ay <- attributes(walk(drop_time(ry)))
    for (n in union(names(ax), names(ay))) if (!identical(ax[[n]], ay[[n]])) { cat("  attr differs:", n, "\n"); print(all.equal(ax[[n]], ay[[n]])) } }
}
cat("ALL_IDENTICAL", all_ok, "\n")
