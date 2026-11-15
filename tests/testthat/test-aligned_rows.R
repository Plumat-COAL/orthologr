# pairwise_aln() writes full-length alignment rows. They used to come from pwalign::alignedPattern() /
# pwalign::alignedSubject(), which took about 40 % of the run time of a default dNdS() call; aligned_rows() builds
# the same strings from the aligned ranges. The pwalign functions stay the reference here.

reference_rows <- function(aln) {
        c(as.character(pwalign::alignedPattern(aln)), as.character(pwalign::alignedSubject(aln)))
}

rows_agree <- function(p, s, type = "global", ...) {
        aln <- pwalign::pairwiseAlignment(Biostrings::AAString(p), Biostrings::AAString(s), type = type, ...)
        identical(orthologr:::aligned_rows(aln, p, s), unname(reference_rows(aln)))
}

test_that("aligned_rows() gives the rows of alignedPattern() / alignedSubject() for overhangs on either side", {
        core <- "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"
        cases <- list(
                equal            = c(core, core),
                pattern_leading  = c(paste0("WWWWWWWWWWWWHHHHHHHH", core), core),
                pattern_trailing = c(paste0(core, "WWWWWWWWWWWWHHHHHHHH"), core),
                pattern_both     = c(paste0("CCCCCCCCCCCCCCC", core, "WWWWWWWWWWWWHHHHHHHH"), core),
                subject_leading  = c(core, paste0("WWWWWWWWWWWWHHHHHHHH", core)),
                subject_trailing = c(core, paste0(core, "WWWWWWWWWWWWHHHHHHHH")),
                subject_both     = c(core, paste0("CCCCCCCCCCCCCCC", core, "WWWWWWWWWWWWHHHHHHHH")),
                opposite_ends    = c(paste0("CCCCCCCCCCCCCCC", core), paste0(core, "WWWWWWWWWWWWHHHHHHHH")),
                opposite_ends_2  = c(paste0(core, "CCCCCCCCCCCCCCC"), paste0("WWWWWWWWWWWWHHHHHHHH", core)),
                internal_gap     = c(core, paste0(substr(core, 1, 30), substr(core, 41, nchar(core)))),
                unrelated        = c("MSTNPKPQRKTKRNTNRRPQDVKFPGG", "GAVLIMFWPSTCYNQDEKRH"),
                one_residue      = c("M", core),
                one_residue_2    = c(core, "V")
        )
        for (n in names(cases)) {
                p <- cases[[n]][1]
                s <- cases[[n]][2]
                expect_true(rows_agree(p, s, substitutionMatrix = "BLOSUM62"), info = n)
                # the scoring of 0.4.3 (no substitution matrix) and other gap costs give other ranges
                expect_true(rows_agree(p, s), info = n)
                expect_true(rows_agree(p, s, substitutionMatrix = "BLOSUM62", gapOpening = 2, gapExtension = 1), info = n)
        }
})

test_that("the rows have the same length and hold every residue of both sequences", {
        p <- "WWWWWWWWWWWWHHHHHHHHMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
        s <- "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQCCCCCCCCCCCC"
        aln <- pwalign::pairwiseAlignment(Biostrings::AAString(p), Biostrings::AAString(s), type = "global",
                                          substitutionMatrix = "BLOSUM62")
        rows <- orthologr:::aligned_rows(aln, p, s)
        expect_identical(nchar(rows[1]), nchar(rows[2]))
        expect_identical(gsub("-", "", rows[1], fixed = TRUE), p)
        expect_identical(gsub("-", "", rows[2], fixed = TRUE), s)
})

test_that("aligned_rows() agrees with pwalign on the proteins of the example files (true and mismatched pairs)", {
        q <- Biostrings::readAAStringSet(thal_aa())
        s <- Biostrings::readAAStringSet(lyra_aa())
        strip <- function(x) sub("\\*$", "", as.character(x))
        q <- strip(q)
        s <- strip(s)
        n <- min(length(q), length(s), 20)
        pairs <- rbind(cbind(seq_len(n), seq_len(n)), cbind(seq_len(n), rev(seq_len(n))))
        for (k in seq_len(nrow(pairs)))
                expect_true(rows_agree(q[[pairs[k, 1]]], s[[pairs[k, 2]]], substitutionMatrix = "BLOSUM62"),
                            info = paste(pairs[k, ], collapse = " / "))
})

test_that("pairwise_aln() writes the rows of alignedPattern() / alignedSubject()", {
        d <- new_test_dir()
        old <- options(orthologr.tmpdir = d)
        on.exit({ options(old); unlink(d, recursive = TRUE) }, add = TRUE)
        f <- file.path(d, "pair.fasta")
        writeLines(c(">query", "WWWWWWWWWWWWHHHHHHHHMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ*",
                     ">subject", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQCCCCCCCCCCCC*"), f)
        aln <- pairwise_aln(f, tool = "NW", seq_type = "protein", get_aln = TRUE, pairwise_aln_name = "t", quiet = TRUE)
        written <- Biostrings::readAAStringSet(file.path(d, "_alignment", "pairwise_aln", "t_NW_AA.aln"))
        expect_identical(toupper(unname(as.character(written))), unname(reference_rows(aln)))
})
