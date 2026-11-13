# How the stored tables of orthologr 0.4.3 were made. Run with orthologr 0.4.3 installed (drostlab/orthologr,
# commit 06c3098), DIAMOND, perl and KaKs_Calculator on the PATH, in an empty working directory:
#
#   Rscript make_fixture.R <output folder>
#
# It writes comeron_20.rds and yn_20.rds: the result of dNdS() on the 20-gene example of the package with the
# defaults of 0.4.3, as a plain data.frame (attributes removed), plus the versions that produced it.
out <- commandArgs(TRUE)[1]
suppressPackageStartupMessages(library(orthologr))
stopifnot(as.character(utils::packageVersion("orthologr")) == "0.4.3")
q <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
s <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
plain <- function(x) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        for (n in names(x))
                attributes(x[[n]]) <- NULL
        attributes(x) <- list(names = names(x), class = "data.frame", row.names = seq_len(nrow(x)))
        x
}
diamond <- system2("diamond", "version", stdout = TRUE)[1]
for (method in c("Comeron", "YN")) {
        res <- suppressMessages(dNdS(q, s, dnds_est.method = method, comp_cores = 1, quiet = TRUE, print_citation = FALSE))
        saveRDS(list(table = plain(res), orthologr = "0.4.3", commit = "06c3098", method = method, diamond = diamond,
                     r = R.version.string, made = format(Sys.Date())),
                file.path(out, paste0(tolower(method), "_20.rds")), version = 2)
}
