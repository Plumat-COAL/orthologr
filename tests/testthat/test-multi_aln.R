aa_example <- system.file("seqs/multi_aln_example.fasta", package = "orthologr")

test_that("multi_aln() rejects an unknown tool and a missing file", {
        expect_error(multi_aln(aa_example, tool = "tcoffee"), "supported by this function")
        expect_error(multi_aln(file.path(tempdir(), "no_such_file.fasta"), tool = "clustalo"), "does not seem to exist")
})

for (tool in c("clustalw", "clustalo", "mafft")) {
        test_that(paste0("multi_aln(tool = '", tool, "') aligns the four example proteins and stores the alignment"), {
                available <- switch(tool, clustalw = clustalw_available, clustalo = clustalo_available, mafft = mafft_available)
                skip_if_not(available, paste(tool, "is not installed"))
                aln <- multi_aln(aa_example, tool = tool, get_aln = TRUE, quiet = TRUE)
                expect_s4_class(aln, "AAMultipleAlignment")
                expect_equal(nrow(aln), 4)
                # every row is one of the input proteins with gaps (a final '*' may be dropped by the program)
                input <- Biostrings::readAAStringSet(aa_example)
                ungapped <- gsub("-", "", unname(as.character(Biostrings::unmasked(aln))), fixed = TRUE)
                expect_setequal(sub("\\*$", "", ungapped), sub("\\*$", "", unname(as.character(input))))
                out <- file.path(tempdir(), "_alignment", "multi_aln", paste0("multi_aln_example.fasta_", tool, ".aln"))
                expect_true(file.exists(out))
                # a name for the stored file
                p <- multi_aln(aa_example, tool = tool, multi_aln_name = "mine", quiet = TRUE)
                expect_identical(basename(p), paste0("mine_", tool, ".aln"))
                expect_true(file.exists(p))
        })
}

test_that("multi_aln(path = ) starts the program from that folder", {
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        aln <- multi_aln(aa_example, tool = "clustalo", get_aln = TRUE, quiet = TRUE, path = dirname(Sys.which("clustalo")))
        expect_equal(nrow(aln), 4)
        tool_cache_reset()
        expect_error(multi_aln(aa_example, tool = "clustalo", quiet = TRUE, path = file.path(tempdir(), "no_such_bin")),
                     "Clustal-Omega installed")
})

test_that("multi_aln(params = ) replaces the default parameters", {
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        aln <- multi_aln(aa_example, tool = "clustalo", get_aln = TRUE, quiet = TRUE, params = "--outfmt clu --iter 1")
        expect_equal(nrow(aln), 4)
})

test_that("multi_aln(clean_folders = TRUE) removes the alignment folder", {
        skip_if_not(mafft_available, "MAFFT is not installed")
        aln <- multi_aln(aa_example, tool = "mafft", get_aln = TRUE, quiet = TRUE, clean_folders = TRUE)
        expect_equal(nrow(aln), 4)
        expect_false(file.exists(file.path(tempdir(), "_alignment")))
})
