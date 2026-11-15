# Options and inputs: what is documented works, what is not supported is rejected with a clear message
# (audit 2026-11-11, problems 3, 4, 5, 6 and 9; stress audit F7, F8, F9).

# ---------------------------------------------------------------- dNdS() arguments

test_that("dNdS() rejects unsupported values before any program is started, naming argument and choices", {
        q <- thal_cds(); s <- lyra_cds()
        expect_error(dNdS(q, s, aligner = "diamnd"), "aligner = 'diamnd' is not supported.*'diamond', 'blast'")
        expect_error(dNdS(q, s, ortho_detection = "PO"), "ortho_detection = 'PO' .*no longer available.*'RBH'")
        expect_error(dNdS(q, s, ortho_detection = "DELTA"), "ortho_detection = 'DELTA' .*no longer available")
        expect_error(dNdS(q, s, ortho_detection = "Orthofinder2"), "'Orthofinder2' is not available in dNdS\\(\\)")
        expect_error(dNdS(q, s, ortho_detection = "rbh"), "ortho_detection = 'rbh' is not supported.*'RBH', 'BH'")
        expect_error(dNdS(q, s, seq_type = "protein"), "seq_type = 'protein' is not possible in dNdS\\(\\)")
        expect_error(dNdS(q, s, seq_type = "dna"), "seq_type = 'dna' is not supported")
        expect_error(dNdS(q, s, format = "gbk"), "format = 'gbk' is not supported.*FASTA")
        expect_error(dNdS(q, s, format = "fastq"), "format = 'fastq' is not supported")
        expect_error(dNdS(q, s, sensitivity_mode = "sensitve"), "sensitivity_mode = 'sensitve' is not a DIAMOND2 sensitivity mode")
        expect_error(dNdS(q, s, aa_aln_type = "pairwise", aa_aln_tool = "clustalw"), "aa_aln_tool = 'clustalw' is not supported.*'NW'")
        expect_error(dNdS(q, s, aa_aln_type = "multiple", aa_aln_tool = "NW"), "aa_aln_tool = 'NW' is not supported.*'clustalw'")
        # nothing was started: no working folder is left
        expect_length(list.files(tempdir(), pattern = "^orthologr_dNdS_"), 0L)
})

test_that("the validators agree with the documentation: BH and RBH only", {
        expect_true(is.ortho_detection_method("BH"))
        expect_true(is.ortho_detection_method("RBH"))
        expect_false(is.ortho_detection_method("PO"))
        expect_false(is.ortho_detection_method("DELTA"))
        expect_identical(eval(formals(dNdS)$sensitivity_mode), "fast")
        expect_true(all(c("fast", "default") %in% diamond_sensitivity_modes()))
})

test_that("dNdS() warns about the deprecated arguments cdd.path and ortho_path", {
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        expect_warning(suppressMessages(dNdS(thal_cds(), lyra_cds(), cdd.path = "x", print_citation = FALSE)),
                       "'cdd.path' is deprecated")
})

test_that("sensitivity_mode = 'default' runs DIAMOND in its own default mode (it failed with '--default')", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        res <- suppressMessages(diamond(thal_cds(), lyra_cds(), sensitivity_mode = "default"))
        expect_equal(nrow(res), 20)
        expect_error(suppressMessages(diamond(thal_cds(), lyra_cds(), sensitivity_mode = "nonsense")),
                     "not a DIAMOND2 sensitivity mode")
        for (m in setdiff(diamond_sensitivity_modes(), "default"))
                expect_equal(nrow(suppressMessages(diamond(thal_cds(), lyra_cds(), sensitivity_mode = m))), 20, info = m)
})

test_that("comp_cores above the number of cores is capped with a warning (it was an error)", {
        n <- parallel::detectCores()
        expect_warning(res <- check_comp_cores(n + 5), "using comp_cores = ")
        expect_identical(res, n)
        expect_silent(check_comp_cores(1))
        expect_error(check_comp_cores(0), "positive")
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        # R CMD check --as-cran (_R_CHECK_LIMIT_CORES_) does not allow more than 2 processes
        skip_if(n > 2 && nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) && Sys.getenv("_R_CHECK_LIMIT_CORES_") != "false",
                "more than 2 cores and _R_CHECK_LIMIT_CORES_ is set")
        expect_warning(
                res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), comp_cores = n + 5, print_citation = FALSE)),
                "more than the .* available"
        )
        expect_equal(nrow(res), 20)
})

test_that("dNdS(clean_folders = TRUE) leaves nothing in tempdir()", {
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        before <- list.files(tempdir(), all.files = TRUE, no.. = TRUE)
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), clean_folders = TRUE, print_citation = FALSE))
        expect_equal(nrow(res), 20)
        expect_setequal(list.files(tempdir(), all.files = TRUE, no.. = TRUE), before)
})

test_that("aa_aln_params is handed to the multiple alignment program (it was ignored)", {
        skip_if_not(diamond_available && perl_available && mafft_available, "DIAMOND2, perl or MAFFT is not installed")
        ok <- suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "mafft",
                                    aa_aln_params = "--quiet --anysymbol --clustalout --op 3", print_citation = FALSE))
        expect_equal(nrow(ok), 20)
        expect_equal(nrow(attr(ok, "failed_pairs")), 0)
        # an option MAFFT does not know now reaches MAFFT: every pair fails and says so
        expect_warning(
                bad <- suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "mafft",
                                             aa_aln_params = "--no-such-option", print_citation = FALSE)),
                "gene pairs could not be processed"
        )
        expect_true(all(is.na(bad$dN)))
        expect_match(attr(bad, "failed_pairs")$reason[1], "MAFFT")
        # with the pairwise aligner there is no command line: a warning, and the result is the default one
        expect_warning(
                nw <- suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_params = "--op 3", print_citation = FALSE)),
                "aa_aln_params is not used"
        )
        expect_equal(nrow(nw), 20)
})

# ---------------------------------------------------------------- functions that failed on their own

test_that("pairwise_aln() works when called on its own", {
        aa <- system.file("seqs/aa_seqs.fasta", package = "orthologr")
        aln <- pairwise_aln(aa, tool = "NW", get_aln = TRUE, quiet = TRUE)
        expect_s4_class(aln, "PairwiseAlignmentsSingleSubject")
        pairwise_aln(aa, tool = "NW", get_aln = FALSE, quiet = TRUE, pairwise_aln_name = "standalone")
        expect_true(file.exists(file.path(tempdir(), "_alignment", "pairwise_aln", "standalone_NW_AA.aln")))
})

test_that("read.genome() reads fastq (object 'genome.dt' not found)", {
        fq <- tempfile(fileext = ".fastq")
        writeLines(c("@r1 first", "ACGTACGTAC", "+", "IIIIIIIIII", "@r2", "acgtttgtac", "+", "IIIIIIIIII"), fq)
        g <- read.genome(fq, format = "fastq")
        expect_s3_class(g, "data.table")
        expect_identical(g$geneids, c("r1", "r2"))
        expect_identical(g$seqs, c("ACGTACGTAC", "ACGTTTGTAC"))
        fa <- read.genome(thal_cds(), format = "fasta")
        expect_equal(nrow(fa), 20)
})

test_that("seq_type = 'dna' is read as it is and rejected clearly where a protein database is needed", {
        dna <- set_diamond(thal_cds(), seq_type = "dna")[[1]]
        expect_true(all(c("geneids", "seqs") %in% names(dna)))
        expect_error(set_diamond(thal_cds(), seq_type = "dna", makedb = TRUE), "seq_type = 'dna' cannot be used to build a DIAMOND database")
        expect_error(set_blast(thal_cds(), seq_type = "dna", makedb = TRUE), "seq_type = 'dna' cannot be used to build a protein database")
        expect_error(set_blast(thal_aa(), seq_type = "protein", makedb = TRUE, makedb_type = "nucleotide"),
                     "cannot be built from protein sequences")
        skip_if_not(blast_available, "BLAST+ is not installed")
        db <- set_blast(thal_cds(), seq_type = "dna", makedb = TRUE, makedb_type = "nucleotide")[[2]]
        expect_identical(db, "blastdb_ortho_thal_cds.fasta_nucleotide.fasta")
        expect_true(file.exists(file.path(tempdir(), "_blast_db", paste0(db, ".nsq"))))
})

test_that("promotor_divergence_estimation() works when all aligned pairs have the same length", {
        skip_if_not_installed("ape")
        set.seed(1)
        mk <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
        q <- vapply(1:5, function(i) mk(60), "")
        # substitutions only: every alignment has length 60
        s <- vapply(q, function(x) { substr(x, 10, 10) <- if (substr(x, 10, 10) == "A") "C" else "A"; x }, "")
        d <- new_test_dir()
        qf <- file.path(d, "query.fa"); sf <- file.path(d, "subject.fa")
        Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(q, paste0("q", 1:5))), qf)
        Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(s, paste0("s", 1:5))), sf)
        res <- suppressMessages(promotor_divergence_estimation(qf, sf))
        expect_equal(nrow(res), 5)
        expect_true(all(res$promotor_dist > 0 & res$promotor_dist < 0.1))
        expect_true(all(res$promotor_indels == 0))
        expect_length(list.files(tempdir(), pattern = "seqs_with_gaps"), 0L)
})

test_that("check_annotation() leaves the annotation file of the user as it is and keeps name and extension", {
        skip_if_not_installed("readr")
        d <- new_test_dir()
        f <- file.path(d, "out.lier.gff")
        lines <- c("##gff-version 3", "chr1\tx\tgene\t1\t20\t.\t+\t.\tID=g1",
                   paste0("chr1\tx\tregion\t1\t2\t.\t+\t.\tNote=", strrep("A", 70000)))
        writeLines(lines, f)
        md5 <- tools::md5sum(f)
        expect_warning(res <- suppressMessages(check_annotation(f)), "65000")
        expect_identical(res, f)
        clean <- suppressMessages(check_annotation(f, remove_annotation_outliers = TRUE))
        expect_identical(unname(tools::md5sum(f)), unname(md5))
        expect_identical(basename(clean), "out.lier.gff")
        expect_false(identical(normalizePath(clean), normalizePath(f)))
        expect_identical(readLines(clean), lines[1:2])
        expect_identical(list.files(d), "out.lier.gff")
        # a file without outlier lines is returned as it is
        g <- file.path(d, "ok.gff"); writeLines(lines[1:2], g)
        expect_identical(suppressMessages(check_annotation(g, remove_annotation_outliers = TRUE)), g)
})

test_that("lnc_map_counts() accepts the species column of lncRNA maps and of ortholog tables", {
        a <- lnc_map_counts(tibble::tibble(query_id = c("a", "b", "c"), species = c("s1", "s1", "s2")), species_order = c("s1", "s2"))
        b <- lnc_map_counts(tibble::tibble(query_id = c("a", "b", "c"), subject_species = c("s1", "s1", "s2")), species_order = c("s1", "s2"))
        expect_identical(a$n_orthologs, c(2L, 1L))
        expect_equal(a, b)
        expect_identical(names(a), c("subject_species", "n_orthologs"))
        expect_error(lnc_map_counts(tibble::tibble(query_id = "a"), "s1"), "'species' or 'subject_species'")
})

test_that("codon_aln(get_aln = TRUE) reads the alignment in the format it was written in", {
        skip_if_not(perl_available, "perl is not installed")
        aln <- system.file("seqs/aa_seqs.aln", package = "orthologr")
        nuc <- system.file("seqs/dna_seqs.fasta", package = "orthologr")
        fa <- suppressWarnings(codon_aln(aln, nuc, format = "fasta", get_aln = TRUE, quiet = TRUE))
        cl <- suppressWarnings(codon_aln(aln, nuc, format = "clustal", get_aln = TRUE, quiet = TRUE))
        expect_s3_class(fa, "alignment")
        expect_equal(fa$nb, 2)
        expect_equal(cl$nb, 2)
        expect_identical(toupper(unlist(fa$seq)), toupper(unlist(cl$seq)))
        # params go to PAL2NAL (before the redirection of its output)
        ng <- suppressWarnings(codon_aln(aln, nuc, format = "fasta", params = "-nogap", get_aln = TRUE, quiet = TRUE))
        expect_false(any(grepl("-", unlist(ng$seq), fixed = TRUE)))
})

# ---------------------------------------------------------------- multiple alignment programs

test_that("the availability checks return TRUE (MUSCLE returned NULL, MAFFT was asked with an invalid option)", {
        if (muscle_available) {
                expect_true(is_installed_muscle())
                expect_true(muscle_major_version() %in% c(3L, 5L))
        }
        if (mafft_available) {
                tool_cache_reset()
                expect_true(is_installed_mafft())
                expect_identical(run_tool("mafft", "--version", error = FALSE)$status, 0L)
        }
        if (clustalo_available) expect_true(is_installed_clustalo())
        if (kaks_available) expect_true(is_installed_kaks_calculator())
        skip_if_not(muscle_available || mafft_available || clustalo_available || kaks_available, "no program to check")
})

test_that("multi_aln(tool = 'muscle') works with MUSCLE 5 and keeps the order of the input", {
        skip_if_not(muscle_available, "MUSCLE is not installed")
        aa <- system.file("seqs/multi_aln_example.fasta", package = "orthologr")
        aln <- multi_aln(aa, tool = "muscle", get_aln = TRUE, quiet = TRUE)
        expect_s4_class(aln, "AAMultipleAlignment")
        expect_identical(sub("\\s.*$", "", names(Biostrings::unmasked(aln))), fasta_ids(aa))
        expect_equal(nrow(aln), 4)
})

test_that("dNdS() works with MUSCLE as aligner of the gene pairs", {
        skip_if_not(muscle_available && diamond_available && perl_available, "MUSCLE, DIAMOND2 or perl is not installed")
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "muscle", print_citation = FALSE))
        expect_equal(nrow(res), 20)
        expect_equal(nrow(attr(res, "failed_pairs")), 0)
        expect_false(anyNA(res$dN))
})

test_that("multi_aln(tool = 't_coffee') works without network in its default mode and leaves no guide tree behind", {
        skip_if_not(tcoffee_available, "T-Coffee is not installed")
        aa <- system.file("seqs/multi_aln_example.fasta", package = "orthologr")
        wd_before <- list.files(getwd())
        aln <- multi_aln(aa, tool = "t_coffee", get_aln = TRUE, quiet = TRUE)
        expect_s4_class(aln, "AAMultipleAlignment")
        expect_equal(nrow(aln), 4)
        expect_setequal(list.files(getwd()), wd_before)
        expect_identical(list.files(file.path(tempdir(), "_alignment", "multi_aln"), pattern = "t_coffee"),
                         "multi_aln_example.fasta_t_coffee.aln")
})

test_that("multi_aln(tool = 't_coffee') aligns sequences that end in '*' (they were returned unaligned)", {
        skip_if_not(tcoffee_available, "T-Coffee is not installed")
        # The example proteins end in '*', as the translated coding sequences of dNdS() do. Given with '-infile',
        # T-Coffee 12 took such a file for a finished alignment: rows without a single internal gap, padded at
        # the end, out of register.
        aa <- system.file("seqs/multi_aln_example.fasta", package = "orthologr")
        aln <- multi_aln(aa, tool = "t_coffee", get_aln = TRUE, quiet = TRUE)
        rows <- as.character(Biostrings::unmasked(aln))
        internal_gaps <- nchar(gsub("[^-]", "", gsub("^-+|-+$", "", rows)))
        expect_true(all(internal_gaps > 0))
        input <- as.character(Biostrings::readAAStringSet(aa))
        expect_identical(unname(gsub("[-*]", "", rows)), unname(gsub("\\*", "", input)))
})

test_that("dNdS() with T-Coffee as multiple aligner agrees with ClustalW on the 20-gene example", {
        skip_if_not(tcoffee_available && clustalw_available && diamond_available && perl_available,
                    "T-Coffee, ClustalW, DIAMOND2 or perl is not installed")
        run <- function(tool) suppressMessages(dNdS(thal_cds(), lyra_cds(), aa_aln_type = "multiple",
                                                    aa_aln_tool = tool, print_citation = FALSE))
        tc <- run("t_coffee")
        cw <- run("clustalw")
        expect_equal(nrow(tc), 20)
        # unaligned rows gave 10 of 20 pairs a saturated (NA) dN/dS and 7 of 20 within 0.1 of ClustalW
        expect_false(anyNA(tc$dNdS))
        m <- merge(as.data.frame(tc)[, c("query_id", "subject_id", "dNdS")],
                   as.data.frame(cw)[, c("query_id", "subject_id", "dNdS")], by = c("query_id", "subject_id"))
        expect_equal(nrow(m), 20)
        expect_gte(sum(abs(m$dNdS.x - m$dNdS.y) <= 0.1, na.rm = TRUE), 18)
        expect_gt(stats::cor(m$dNdS.x, m$dNdS.y, use = "complete.obs"), 0.95)
})

# ---------------------------------------------------------------- coding sequences that are left out

test_that("is.dnaSequence() accepts A, C, G, T only", {
        expect_true(is.dnaSequence("ACGTacgt"))
        expect_false(is.dnaSequence("ACGNT"))
        expect_false(is.dnaSequence("ACGRYT"))
        # the old pattern "[^(A|C|G|T|a|c|g|t)]" let these pass
        expect_false(is.dnaSequence("ACG|T"))
        expect_false(is.dnaSequence("ACG(T)"))
})

test_that("coding sequences with N or IUPAC codes are left out with a message that names them", {
        d <- new_test_dir()
        f <- file.path(d, "cds.fasta")
        writeLines(c(">g1", "ATGAAATTTTAA", ">g2", "ATGAANTTTTAA", ">g3", "ATGAARTTTTAA", ">g4", "ATGCCCGGGTAA"), f)
        expect_message(res <- set_diamond(f)[[1]], "2 of 4 coding sequences of 'cds.fasta' have characters other than A, C, G, T.*g2, g3")
        expect_identical(res$geneids, c("g1", "g4"))
        expect_identical(toupper(res$aa), c("MKF*", "MPG*"))
        # nothing to report for clean input
        expect_message(set_diamond(thal_cds()), NA)
})

test_that("the file with corrupted coding sequences is not written to the working directory", {
        d <- new_test_dir()
        f <- file.path(d, "cds.fasta")
        writeLines(c(">g1", "ATGAAATTTTAA", ">g2", "ATGAAATTTTA"), f)
        wd_before <- list.files(getwd(), all.files = TRUE)
        res <- suppressMessages(read.cds(f, format = "fasta"))
        expect_identical(res$geneids, "g1")
        expect_setequal(list.files(getwd(), all.files = TRUE), wd_before)
        out <- file.path(tempdir(), "cds.fasta_corrupted_cds_seqs.fasta")
        expect_true(file.exists(out))
        expect_identical(names(Biostrings::readDNAStringSet(out)), "g2")
        # ... or to the folder the user names
        mine <- file.path(d, "inspect")
        old <- options(orthologr.corrupted_cds_dir = mine)
        on.exit(options(old), add = TRUE)
        expect_message(read.cds(f, format = "fasta", delete_corrupt_cds = FALSE), "inspect")
        expect_true(file.exists(file.path(mine, "cds.fasta_corrupted_cds_seqs.fasta")))
})

# ---------------------------------------------------------------- OrthoFinder (see also test-orthofinder2.R)

test_that("orthologs(ortho_detection = 'Orthofinder2') is an error instead of NULL", {
        expect_error(orthologs(thal_aa(), lyra_aa(), ortho_detection = "Orthofinder2"), "not implemented in orthologs\\(\\)")
})
