write_pair <- function(a, b, file = tempfile(fileext = ".fasta")) {
        seqinr::write.fasta(sequences = list(a, b), names = c("seq_a", "seq_b"), file.out = file)
        file
}

read_aln <- function(name) {
        f <- file.path(tempdir(), "_alignment", "pairwise_aln", paste0(name, "_NW_AA.aln"))
        x <- seqinr::read.fasta(f, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
        unlist(x)
}

long_aa  <- "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"
short_aa <- "AYIAKQRQISFVKSHFSRQIEERLGLIEVQAPILSRVGDGTQDN"

test_that("pairwise_aln() writes a full-length alignment (overhangs are not cut off) ...", {
        f <- write_pair(long_aa, short_aa)
        pairwise_aln(file = f, tool = "NW", seq_type = "protein", pairwise_aln_name = "t_full", quiet = TRUE)
        aln <- read_aln("t_full")
        expect_equal(names(aln), c("seq_a", "seq_b"))
        expect_equal(nchar(aln[[1]]), nchar(aln[[2]]))
        # every residue of both sequences is in the alignment
        expect_equal(gsub("-", "", aln[[1]]), long_aa)
        expect_equal(gsub("-", "", aln[[2]]), short_aa)
        # the shorter sequence is padded with gaps where the longer one overhangs
        expect_true(startsWith(aln[[2]], "---"))
        expect_true(endsWith(aln[[2]], "---"))
        # and the same holds when the shorter sequence comes first
        f2 <- write_pair(short_aa, long_aa)
        pairwise_aln(file = f2, tool = "NW", seq_type = "protein", pairwise_aln_name = "t_full2", quiet = TRUE)
        aln2 <- read_aln("t_full2")
        expect_equal(gsub("-", "", aln2[[1]]), short_aa)
        expect_equal(gsub("-", "", aln2[[2]]), long_aa)
})

test_that("pairwise_aln() scores amino acids with BLOSUM62 by default and accepts other matrices ...", {
        f <- write_pair(long_aa, short_aa)
        aln <- pairwise_aln(file = f, tool = "NW", seq_type = "protein", get_aln = TRUE,
                            pairwise_aln_name = "t_b62", quiet = TRUE)
        ref <- pwalign::pairwiseAlignment(Biostrings::AAString(long_aa), Biostrings::AAString(short_aa),
                                          type = "global", substitutionMatrix = "BLOSUM62")
        expect_equal(pwalign::score(aln), pwalign::score(ref))
        
        old <- pairwise_aln(file = f, tool = "NW", seq_type = "protein", get_aln = TRUE,
                            pairwise_aln_name = "t_old", quiet = TRUE, aa_sub_matrix = NULL)
        ref_old <- pwalign::pairwiseAlignment(Biostrings::AAString(long_aa), Biostrings::AAString(short_aa),
                                              type = "global")
        expect_equal(pwalign::score(old), pwalign::score(ref_old))
        expect_false(isTRUE(all.equal(pwalign::score(old), pwalign::score(aln))))
        
        b45 <- pairwise_aln(file = f, tool = "NW", seq_type = "protein", get_aln = TRUE,
                            pairwise_aln_name = "t_b45", quiet = TRUE, aa_sub_matrix = "BLOSUM45")
        ref45 <- pwalign::pairwiseAlignment(Biostrings::AAString(long_aa), Biostrings::AAString(short_aa),
                                            type = "global", substitutionMatrix = "BLOSUM45")
        expect_equal(pwalign::score(b45), pwalign::score(ref45))
        
        expect_error(pairwise_aln(file = f, tool = "NW", seq_type = "protein", quiet = TRUE, aa_sub_matrix = 62),
                     "aa_sub_matrix")
        expect_error(pairwise_aln(file = f, tool = "INVALID"))
        expect_error(pairwise_aln(file = f, seq_type = "INVALID"))
})

test_that("an unknown aa_sub_matrix is an error before the search (it failed every gene pair after it)", {
        f <- tempfile(fileext = ".fasta")
        seqinr::write.fasta(sequences = list("MKVLAAGIVG", "MKVLSAGIVG"), names = c("a", "b"), file.out = f)
        expect_error(pairwise_aln(file = f, tool = "NW", seq_type = "protein", quiet = TRUE, aa_sub_matrix = "FOO"),
                     "aa_sub_matrix", class = "orthologr_error_bad_argument")
        expect_error(pairwise_aln(file = f, tool = "NW", seq_type = "protein", quiet = TRUE, aa_sub_matrix = 62),
                     class = "orthologr_error_bad_argument")
        # before DIAMOND runs: no program is needed to get this error
        expect_error(dNdS(thal_cds(), lyra_cds(), aa_sub_matrix = "FOO", aligner_path = tempfile("no_diamond_"),
                          print_citation = FALSE),
                     "BLOSUM62", class = "orthologr_error_bad_argument")
        expect_error(dNdS(thal_cds(), lyra_cds(), aa_sub_matrix = "blosum62", print_citation = FALSE),
                     class = "orthologr_error_bad_argument")
        for (m in c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"))
                expect_true(orthologr:::check_aa_sub_matrix(m), info = m)
        expect_true(orthologr:::check_aa_sub_matrix(NULL))
})

test_that("pairwise_aln() does not align a terminal stop codon symbol ...", {
        f <- write_pair(paste0(long_aa, "*"), paste0(short_aa, "*"))
        pairwise_aln(file = f, tool = "NW", seq_type = "protein", pairwise_aln_name = "t_stop", quiet = TRUE)
        aln <- read_aln("t_stop")
        expect_false(any(grepl("*", aln, fixed = TRUE)))
        expect_equal(gsub("-", "", aln[[1]]), long_aa)
        expect_equal(gsub("-", "", aln[[2]]), short_aa)
        # an internal stop symbol is kept
        f2 <- write_pair("MKTAYIAKQ*QISFVKSHF*", "MKTAYIAKQRQISFVKSHF")
        pairwise_aln(file = f2, tool = "NW", seq_type = "protein", pairwise_aln_name = "t_stop2", quiet = TRUE)
        expect_equal(gsub("-", "", read_aln("t_stop2")[[1]]), "MKTAYIAKQ*QISFVKSHF")
})

test_that("pairwise_aln() writes the score file only when store_locally = TRUE ...", {
        score_file <- file.path(tempdir(), "_pairwise_alignment_with_score", "t_nostore_NW_AA.aln")
        f <- write_pair(long_aa, short_aa)
        pairwise_aln(file = f, tool = "NW", seq_type = "protein", pairwise_aln_name = "t_nostore", quiet = TRUE)
        expect_false(file.exists(score_file))
        
        wd_new <- tempfile("orthologr_store_")
        dir.create(wd_new)
        old_wd <- setwd(wd_new)
        on.exit(setwd(old_wd), add = TRUE)
        pairwise_aln(file = f, tool = "NW", seq_type = "protein", pairwise_aln_name = "t_store", quiet = TRUE,
                     store_locally = TRUE)
        expect_true(file.exists(file.path("orthologr_alignment_files", "_alignment", "pairwise_aln", "t_store_NW_AA.aln")))
        expect_true(file.exists(file.path("orthologr_alignment_files", "_pairwise_alignment_with_score", "t_store_NW_AA.aln")))
})

test_that("options(orthologr.legacy_pairwise_aln = TRUE) restores the 0.4.3 behaviour ...", {
        old_opt <- options(orthologr.legacy_pairwise_aln = TRUE)
        on.exit(options(old_opt), add = TRUE)
        f <- write_pair(paste0(long_aa, "*"), paste0(short_aa, "*"))
        aln <- pairwise_aln(file = f, tool = "NW", seq_type = "protein", get_aln = TRUE,
                            pairwise_aln_name = "t_legacy", quiet = TRUE)
        ref_old <- pwalign::pairwiseAlignment(Biostrings::AAString(paste0(long_aa, "*")),
                                              Biostrings::AAString(paste0(short_aa, "*")), type = "global")
        expect_equal(pwalign::score(aln), pwalign::score(ref_old))
        written <- read_aln("t_legacy")
        # the 0.4.3 output: only the aligned range, as pwalign::pattern() / subject() return it
        expect_equal(toupper(written[[1]]), as.character(pwalign::pattern(ref_old)))
        expect_equal(toupper(written[[2]]), as.character(pwalign::subject(ref_old)))
        expect_true(file.exists(file.path(tempdir(), "_pairwise_alignment_with_score", "t_legacy_NW_AA.aln")))
})

test_that("pairwise_aln() aligns DNA sequences as before (aa_sub_matrix is ignored) ...", {
        f <- tempfile(fileext = ".fasta")
        seqinr::write.fasta(sequences = list("ATGGCTAAAGGTCTGGAACGT", "ATGGCTAAAGGCCTGGAA"),
                            names = c("seq_a", "seq_b"), file.out = f)
        aln <- pairwise_aln(file = f, tool = "NW", seq_type = "dna", get_aln = TRUE,
                            pairwise_aln_name = "t_dna", quiet = TRUE)
        ref <- pwalign::pairwiseAlignment(Biostrings::DNAString("ATGGCTAAAGGTCTGGAACGT"),
                                          Biostrings::DNAString("ATGGCTAAAGGCCTGGAA"), type = "global")
        expect_equal(pwalign::score(aln), pwalign::score(ref))
        out <- seqinr::read.fasta(file.path(tempdir(), "_alignment", "pairwise_aln", "t_dna_NW_DNA.aln"),
                                  as.string = TRUE, set.attributes = FALSE)
        expect_equal(toupper(gsub("-", "", out[[1]])), "ATGGCTAAAGGTCTGGAACGT")
        expect_equal(toupper(gsub("-", "", out[[2]])), "ATGGCTAAAGGCCTGGAA")
})

test_that("dNdS() takes aa_sub_matrix, and the codon alignments no longer hold the stop codon ...", {
        skip_if_not(isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE)),
                    "DIAMOND2 is not installed or not on PATH")
        query_cds   <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
        subject_cds <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
        # keep the call's working folder to look at its codon alignments
        before <- list.files(tempdir(), pattern = "^orthologr_dNdS_")
        old_opt <- options(orthologr.keep_tmp = TRUE)
        res <- suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, comp_cores = 1,
                                     quiet = TRUE, print_citation = FALSE))
        options(old_opt)
        run_dir <- file.path(tempdir(), setdiff(list.files(tempdir(), pattern = "^orthologr_dNdS_"), before))
        on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)
        expect_s3_class(res, "data.frame")
        expect_gt(nrow(res), 0)
        expect_true(all(c("query_id", "subject_id", "dN", "dS", "dNdS") %in% names(res)))
        # codon alignments of this run: full length, and no pair ends in two aligned stop codons
        codon_files <- list.files(file.path(run_dir, "_alignment", "codon_aln"), pattern = "^q[0-9]+_pal2nal\\.aln$",
                                  full.names = TRUE)
        expect_gt(length(codon_files), 0)
        last_codons <- vapply(codon_files[seq_len(min(10, length(codon_files)))], function(f) {
                x <- toupper(unlist(seqinr::read.fasta(f, as.string = TRUE, set.attributes = FALSE)))
                x <- gsub("-+$", "", x)
                paste(substr(x, nchar(x) - 2, nchar(x)), collapse = " ")
        }, character(1))
        expect_false(any(grepl("^(TAA|TAG|TGA) (TAA|TAG|TGA)$", last_codons)))
        
        res_old <- suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, comp_cores = 1,
                                         quiet = TRUE, print_citation = FALSE, aa_sub_matrix = NULL))
        expect_equal(nrow(res_old), nrow(res))
        expect_error(suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, comp_cores = 1,
                                           quiet = TRUE, print_citation = FALSE, aa_sub_matrix = 1)))
})

test_that("pairwise_aln() does not align pairs above options(orthologr.max_aln_cells); dNdS() reports them as failed (O12)", {
        f <- system.file("seqs/aa_seqs.fasta", package = "orthologr")
        rlang::local_options(orthologr.max_aln_cells = 100)
        e <- tryCatch(pairwise_aln(f, quiet = TRUE), error = function(e) e)
        expect_s3_class(e, "orthologr_error_alignment_too_large")
        expect_match(conditionMessage(e), "orthologr.max_aln_cells = Inf")
        rlang::local_options(orthologr.max_aln_cells = 1e7)
        expect_silent(pairwise_aln(f, quiet = TRUE))
        rlang::local_options(orthologr.max_aln_cells = "a lot")
        expect_error(pairwise_aln(f, quiet = TRUE), class = "orthologr_error_bad_argument")
        # Inf: no limit (also under the legacy switch)
        rlang::local_options(orthologr.max_aln_cells = NULL, orthologr.legacy = TRUE)
        expect_identical(max_aln_cells(), Inf)
        rlang::local_options(orthologr.legacy = NULL)
        expect_identical(max_aln_cells(), 5e8)
})

test_that("dNdS(): a gene pair above options(orthologr.max_aln_cells) is a failed pair with the reason (O12)", {
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        rlang::local_options(orthologr.max_aln_cells = 1e5)
        expect_warning(res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), quiet = TRUE, print_citation = FALSE)),
                       class = "orthologr_warning_failed_pairs")
        failed <- attr(res, "failed_pairs")
        expect_gt(nrow(failed), 0)
        expect_lt(nrow(failed), 20)
        expect_true(all(grepl("orthologr.max_aln_cells", failed$reason)))
        expect_true(all(res$status[paste(res$query_id, res$subject_id) %in% paste(failed$query_id, failed$subject_id)] == "failed"))
})
