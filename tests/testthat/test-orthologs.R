test_that("orthologs() rejects unknown methods and names the valid ones", {
        expect_error(orthologs(thal_aa(), lyra_aa(), ortho_detection = "PO"), "Unsupported orthology detection method: 'PO'")
        expect_error(orthologs(thal_aa(), lyra_aa(), ortho_detection = "Orthofinder2"), "not implemented")
})

test_that("orthologs() warns about the arguments that were never used", {
        expect_warning(
                expect_error(orthologs(thal_aa(), lyra_aa(), ortho_detection = "Orthofinder2", cdd.path = "x")),
                "'cdd.path' is deprecated"
        )
})

test_that("orthologs() returns the reciprocal best hits of DIAMOND2 for proteins and for CDS", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        prot <- suppressMessages(orthologs(thal_aa(), lyra_aa(), seq_type = "protein", ortho_detection = "DIAMOND_RBH"))
        expect_equal(nrow(prot), 20)
        expect_true(all(c("query_id", "subject_id", "evalue") %in% names(prot)))
        expect_false(anyDuplicated(prot$query_id) > 0)
        expect_false(anyDuplicated(prot$subject_id) > 0)
        cds <- suppressMessages(orthologs(thal_cds(), lyra_cds(), seq_type = "cds", ortho_detection = "DIAMOND_RBH"))
        expect_equal(nrow(cds), 20)
        # the same gene pairs (the protein file ids carry a '|PACid...' suffix on the query side)
        expect_setequal(sub("\\|.*$", "", cds$subject_id), sub("\\|.*$", "", prot$subject_id))
        # it is what diamond_rec() returns
        rec <- suppressMessages(diamond_rec(thal_aa(), lyra_aa(), seq_type = "protein"))
        expect_equal(as.data.frame(prot)[ , c("query_id", "subject_id")], as.data.frame(rec)[ , c("query_id", "subject_id")],
                     ignore_attr = TRUE)
})

test_that("orthologs() best hits: one hit per query, a superset of the reciprocal best hits", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        bh  <- suppressMessages(orthologs(thal_aa(), lyra_aa(), seq_type = "protein", ortho_detection = "DIAMOND_BH"))
        rbh <- suppressMessages(orthologs(thal_aa(), lyra_aa(), seq_type = "protein", ortho_detection = "DIAMOND_RBH"))
        expect_false(anyDuplicated(bh$query_id) > 0)
        expect_true(all(paste(rbh$query_id, rbh$subject_id) %in% paste(bh$query_id, bh$subject_id)))
})

test_that("orthologs() with BLAST+ finds the same reciprocal best hits as with DIAMOND2 on the example", {
        skip_if_not(diamond_available && blast_available, "DIAMOND2 or BLAST+ is not installed")
        b <- suppressMessages(orthologs(thal_aa(), lyra_aa(), seq_type = "protein", ortho_detection = "RBH"))
        d <- suppressMessages(orthologs(thal_aa(), lyra_aa(), seq_type = "protein", ortho_detection = "DIAMOND_RBH"))
        expect_setequal(paste(b$query_id, b$subject_id), paste(d$query_id, d$subject_id))
})

test_that("orthologs() stops for several subject files (pairwise methods only)", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        expect_error(suppressMessages(orthologs(thal_aa(), c(lyra_aa(), lyra_aa()), seq_type = "protein")),
                     "only defined for pairwise comparisons")
})
