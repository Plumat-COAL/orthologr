test_that("read.proteome() reads FASTA into a data.table keyed by the ids", {
        p <- read.proteome(thal_aa(), format = "fasta")
        expect_s3_class(p, "data.table")
        expect_identical(names(p), c("geneids", "seqs"))
        expect_equal(nrow(p), 20)
        expect_identical(data.table::key(p), "geneids")
        expect_false(any(grepl(" ", p$geneids)))
        expect_identical(p$seqs, toupper(p$seqs))
        ref <- Biostrings::readAAStringSet(thal_aa())
        expect_setequal(p$geneids, sub(" .*", "", names(ref)))
})

test_that("read.proteome() and write.proteome() round-trip", {
        p <- read.proteome(thal_aa(), format = "fasta")
        f <- tempfile(fileext = ".faa")
        write.proteome(p, f)
        expect_equal(read.proteome(f, format = "fasta"), p)
})

test_that("read.proteome() rejects other formats, missing files and files that hold no sequences", {
        expect_error(read.proteome(thal_aa(), format = "gbk"), "file format that is supported")
        expect_error(read.proteome(file.path(tempdir(), "nope.faa"), format = "fasta"), "does not seem to exist")
        empty <- tempfile(fileext = ".faa"); file.create(empty)
        expect_error(read.proteome(empty, format = "fasta"), "could not be read properly")
})
