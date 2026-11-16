# Characterise rows where the outputs of the old and new kaks_lwl_counts() differ (outputs only)
a <- readRDS("fn_base.rds"); b <- readRDS("fn_new.rds")
x <- rbind(utils::read.delim("rand1.tsv", colClasses = "character", na.strings = character())[, c("name", "s1", "s2", "method", "code")],
           utils::read.delim("rand2.tsv", colClasses = "character", na.strings = character())[, c("name", "s1", "s2", "method", "code")],
           utils::read.delim("zero_all.tsv", colClasses = "character", na.strings = character())[, c("name", "s1", "s2", "method", "code")])
d <- which(!mapply(identical, a, b))
kind <- vapply(d, function(i) {
  p <- a[[i]]; q <- b[[i]]
  if (is.null(p) || is.null(q)) return(paste("null", is.null(p), is.null(q)))
  if (!identical(names(p), names(q))) return(paste("names", paste(names(p), collapse = ",")))
  f <- names(p)[!mapply(identical, p, q)]
  num <- all(vapply(f, function(n) isTRUE(all.equal(unname(p[[n]]), unname(q[[n]]), tolerance = 1e-12)), TRUE))
  attrs <- all(vapply(f, function(n) identical(unname(p[[n]]), unname(q[[n]])), TRUE))
  paste0(paste(f, collapse = ","), if (attrs) " (attributes/names only)" else if (num) " (numeric, < 1e-12 rel)" else " (LARGER)")
}, "")
print(table(kind))
odd <- grepl("[^ACGT]", toupper(x$s1[d])) | grepl("[^ACGT]", toupper(x$s2[d]))
lower <- grepl("[a-z]", x$s1[d]) | grepl("[a-z]", x$s2[d])
u <- grepl("[Uu]", x$s1[d]) | grepl("[Uu]", x$s2[d])
cat("among differing rows: non-ACGT", sum(odd), " lowercase", sum(lower), " U", sum(u), " of", length(d), "\n")
cat("by method:\n"); print(table(x$method[d])); cat("by source:\n"); print(table(sub("_.*", "", x$name[d])))
allodd <- grepl("[^ACGT]", toupper(x$s1)) | grepl("[^ACGT]", toupper(x$s2))
cat("rows with any non-ACGT char:", sum(allodd), " of which differ:", sum(d %in% which(allodd)), "\n")
for (i in head(d[!odd], 3)) { cat("--- row", i, x$name[i], x$method[i], x$code[i], "\n"); str(a[[i]]); str(b[[i]]) }
for (i in head(d[odd], 2)) { cat("--- odd row", i, x$name[i], x$method[i], x$code[i], "\n"); str(a[[i]]); str(b[[i]]) }
