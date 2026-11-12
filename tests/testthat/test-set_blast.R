test_that("set_blast() returns the proteins (and no database name) when makedb = FALSE", {
        res <- set_blast(thal_cds())
        expect_length(res, 2)
        expect_s3_class(res[[1]], "data.table")
        expect_true(all(c("geneids", "seqs", "aa") %in% names(res[[1]])))
        expect_equal(nrow(res[[1]]), 20)
        expect_identical(res[[2]], "")
        # the translation of the coding sequences is the protein file of the package
        prot <- set_blast(thal_aa(), seq_type = "protein")[[1]]
        expect_identical(sub("\\*$", "", toupper(res[[1]]$aa)), sub("\\*$", "", prot$aa[match(res[[1]]$geneids, sub("\\|.*$", "", prot$geneids))]))
})

test_that("set_blast() checks its arguments", {
        expect_error(set_blast(file.path(tempdir(), "nope.fasta")), "seems not to exist")
        expect_error(set_blast(thal_cds(), seq_type = "rna"), "'cds', 'protein', or 'dna'")
        expect_error(set_blast(thal_cds(), makedb_type = "x"), "'protein' or 'nucleotide'")
})

test_that("set_blast(makedb = TRUE) builds the database in tempdir()/_blast_db and returns its name", {
        skip_if_not(blast_available, "BLAST+ is not installed")
        db <- set_blast(thal_cds(), makedb = TRUE)[[2]]
        expect_identical(db, "blastdb_ortho_thal_cds.fasta_protein.fasta")
        dir <- file.path(tempdir(), "_blast_db")
        expect_true(file.exists(file.path(dir, db)))
        expect_true(file.exists(file.path(dir, paste0(db, ".psq"))))
        # the database can be searched
        q <- tempfile(fileext = ".faa")
        Biostrings::writeXStringSet(Biostrings::readAAStringSet(thal_aa())[1], q)
        out <- tempfile()
        res <- run_tool("blastp", args = list("-db", file.path(dir, db), "-query", q, "-outfmt", "6", "-out", out))
        expect_gt(length(readLines(out)), 0)
})
