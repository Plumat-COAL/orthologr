test_that("read.genome() reads FASTA into a data.table keyed by the ids (first word of the header)", {
        g <- read.genome(thal_cds(), format = "fasta")
        expect_s3_class(g, "data.table")
        expect_identical(names(g), c("geneids", "seqs"))
        expect_equal(nrow(g), 20)
        expect_identical(data.table::key(g), "geneids")
        expect_false(any(grepl(" ", g$geneids)))
        expect_true(all(grepl("^[ACGTN]+$", g$seqs)))        # upper case
        ref <- Biostrings::readDNAStringSet(thal_cds())
        expect_identical(unname(nchar(g$seqs)), unname(Biostrings::width(ref)[order(sub(" .*", "", names(ref)), method = "radix")]))
})

test_that("read.genome() rejects other formats, missing and empty files", {
        expect_error(read.genome(thal_cds(), format = "gbk"), "file format that is supported")
        expect_error(read.genome(file.path(tempdir(), "nope.fa"), format = "fasta"), "does not seem to exist")
        empty <- tempfile(fileext = ".fa"); file.create(empty)
        expect_error(read.genome(empty, format = "fasta"), "could not be read properly")
})
