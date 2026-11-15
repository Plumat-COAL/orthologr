aa_aln_example  <- system.file("seqs/aa_seqs.aln", package = "orthologr")
dna_example     <- system.file("seqs/dna_seqs.fasta", package = "orthologr")

test_that("codon_aln() checks its arguments", {
        expect_error(codon_aln(file.path(tempdir(), "nope.aln"), dna_example), "file_aln")
        expect_error(codon_aln(aa_aln_example, file.path(tempdir(), "nope.fasta")), "file_nuc")
        expect_error(codon_aln(aa_aln_example, dna_example, tool = "other"), "tool that is supported")
        expect_error(codon_aln(aa_aln_example, dna_example, format = "paml"), "format that is supported")
})

test_that("codon_aln() threads the coding sequences through the protein alignment", {
        skip_if_not(perl_available, "perl is not installed")
        out <- suppressWarnings(codon_aln(aa_aln_example, dna_example, format = "fasta", codon_aln_name = "example", quiet = TRUE))
        expect_identical(out, file.path(tempdir(), "_alignment", "codon_aln", "example_pal2nal.aln"))
        codon <- Biostrings::readBStringSet(out)
        cds   <- Biostrings::readDNAStringSet(dna_example)
        expect_length(codon, 2)
        # same length, a multiple of 3, gaps in whole codons, and without gaps the codons are those of the CDS
        expect_identical(length(unique(Biostrings::width(codon))), 1L)
        expect_identical(unname(Biostrings::width(codon)[1] %% 3L), 0L)
        for (i in 1:2) {
                codons <- substring(as.character(codon[[i]]), seq(1, Biostrings::width(codon)[i], 3), seq(3, Biostrings::width(codon)[i], 3))
                expect_true(all(codons == "---" | !grepl("-", codons, fixed = TRUE)))
                ungapped <- toupper(paste(codons[codons != "---"], collapse = ""))
                expect_true(grepl(ungapped, toupper(as.character(cds[[i]])), fixed = TRUE))
        }
        # the amino acid alignment has as many columns as the codon alignment has codons
        aa <- read_clustal(aa_aln_example)
        expect_identical(unname(nchar(aa[1]) * 3L), unname(Biostrings::width(codon)[1]))
})

test_that("codon_aln() reports the gene pair when PAL2NAL cannot match protein and CDS", {
        skip_if_not(perl_available, "perl is not installed")
        d <- new_test_dir()
        bad <- file.path(d, "bad.fasta")
        writeLines(c(">a", "ATGAAATTTCCCGGG", ">b", "ATGAAATTTCCCGGG"), bad)
        aln <- file.path(d, "aa.fasta")
        writeLines(c(">a", "MWWWW", ">b", "MWWWW"), aln)
        expect_error(suppressWarnings(codon_aln(aln, bad, format = "fasta", quiet = TRUE)), "PAL2NAL .* gene pair 'a' / 'b'")
})
