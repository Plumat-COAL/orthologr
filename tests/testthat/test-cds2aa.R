test_that("cds2aa() adds the protein of every coding sequence", {
        dt <- data.table::data.table(geneids = c("g1", "g2"), seqs = c("atgaaatttccctaa", "ATGGGGTGA"))
        res <- cds2aa(data.table::copy(dt))
        expect_identical(res$aa, c("MKFP*", "MG*"))
        expect_identical(res$geneids, c("g1", "g2"))
})

test_that("cds2aa() leaves out empty sequences, incomplete codons (delete_corrupt_cds) and non-ACGT sequences", {
        dt <- data.table::data.table(geneids = paste0("g", 1:5),
                                     seqs = c("atgaaataa", "", "atgaaata", "atgnnntaa", NA))
        res <- suppressMessages(cds2aa(data.table::copy(dt), delete_corrupt_cds = TRUE))
        expect_identical(res$geneids, "g1")
        expect_message(cds2aa(data.table::copy(dt), delete_corrupt_cds = TRUE), "1 of 5 coding sequences has characters other than A, C, G, T.*g4")
        # the empty sequences are reported as well (they were left out silently)
        expect_message(cds2aa(data.table::copy(dt), delete_corrupt_cds = TRUE), "2 of 5 coding sequences have no sequence and are left out: g2, g5")
        expect_message(cds2aa(data.table::copy(dt), delete_corrupt_cds = TRUE), "1 of 5 coding sequences has a length that is not a multiple of 3")
        # delete_corrupt_cds = FALSE keeps the sequence with an incomplete last codon (which is not translated)
        res2 <- suppressMessages(cds2aa(data.table::copy(dt), delete_corrupt_cds = FALSE))
        expect_identical(res2$geneids, c("g1", "g3"))
        expect_identical(res2$aa, c("MK*", "MK"))
        expect_error(cds2aa(as.data.frame(dt)), "data.table")
})
