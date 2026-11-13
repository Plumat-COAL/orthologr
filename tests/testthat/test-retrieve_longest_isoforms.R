test_that("retrieve_longest_isoforms() checks its arguments", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        out <- tempfile(fileext = ".faa")
        expect_error(suppressMessages(retrieve_longest_isoforms(file.path(x$dir, "nope.faa"), x$query$gff, out)))
        expect_error(suppressMessages(retrieve_longest_isoforms(x$query$proteome, x$query$gff, out, annotation_format = "bed")))
})

for (fmt in c("gff", "gtf")) {
        test_that(paste0("retrieve_longest_isoforms() writes the longest isoform of every gene locus (", fmt, ")"), {
                skip_if_not_installed("rtracklayer")
                x <- make_test_ortholog_input()
                out <- tempfile(fileext = ".faa")
                suppressMessages(retrieve_longest_isoforms(x$query$proteome, x$query[[fmt]], out, annotation_format = fmt))
                res <- Biostrings::readAAStringSet(out)
                all <- Biostrings::readAAStringSet(x$query$proteome)
                expect_length(res, 20)                       # 25 isoforms of 20 loci
                ids <- sub(" .*", "", names(res))
                expect_false(any(grepl("_iso2$", ids)))      # the second isoforms are the shorter ones
                expect_identical(as.character(res), as.character(all[match(ids, sub(" .*", "", names(all)))]), ignore_attr = TRUE)
                # the input files are left as they were
                expect_length(Biostrings::readAAStringSet(x$query$proteome), 25)
        })
}

test_that("retrieve_longest_isoforms_all() processes every proteome of a folder", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        prot <- file.path(x$dir, "prot"); ann <- file.path(x$dir, "ann"); out <- file.path(x$dir, "longest")
        dir.create(prot); dir.create(ann)
        file.copy(c(x$query$proteome, x$subject$proteome), prot)
        file.copy(c(x$query$gff, x$subject$gff), ann)
        suppressMessages(retrieve_longest_isoforms_all(prot, ann, out, annotation_format = "gff"))
        res <- list.files(out, full.names = TRUE)
        expect_length(res, 2)
        expect_identical(unname(vapply(res, function(f) length(Biostrings::readAAStringSet(f)), 1L)), c(20L, 20L))
})
