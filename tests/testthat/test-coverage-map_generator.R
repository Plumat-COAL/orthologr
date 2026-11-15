# Coverage of map_generator_dnds() (R/map.generator.R), map_generator_lnc() (R/map.generator.lnc.R) and
# orthologs_lnc() (R/orthologs_lnc.R). The searches they start (dNdS(); BLAST through the suggested package
# metablastr) are replaced by stand-ins, so that the loops, filters, output files and checks of these functions
# are tested without DIAMOND, BLAST or metablastr. One test runs the real BLAST search when metablastr and
# BLAST+ are installed.

# A hit table as orthologs_lnc() returns it (metablastr's columns that map_generator_lnc() uses).
cov_lnc_hits <- function(prefix = "") {
        tibble::tibble(
                query_id      = paste0(prefix, c("q1", "q2", "q3", "q4")),
                subject_id    = paste0(prefix, c("s1", "s2", "s3", "s4")),
                perc_identity = c(50, 20, 90, 10),
                num_ident_matches = c(100L, 40L, 27L, 3L),
                alig_length   = c(200L, 200L, 30L, 30L),
                qcovhsp       = c(80, 80, 10, 10),
                evalue        = c(1e-50, 1e-20, 1e-5, 1e-3),
                bit_score     = c(300, 120, 40, 20)
        )
}

cov_subjects_folder <- function(files = c("Alyrata.fa", "Bnapus.v2.fa")) {
        dir <- new_test_dir("cov_subjects_")
        for (f in files)
                writeLines(c(">x", "ACGT"), file.path(dir, f))
        dir
}

cov_run_lnc <- function(..., hits = cov_lnc_hits(), env = parent.frame()) {
        calls <- new.env()
        calls$args <- list()
        testthat::local_mocked_bindings(orthologs_lnc = function(...) {
                calls$args[[length(calls$args) + 1]] <- list(...)
                hits
        }, .package = "orthologr", .env = env)
        out <- new_test_dir("cov_lnc_out_")
        unlink(out, recursive = TRUE)   # created by the function
        res <- suppressMessages(map_generator_lnc(query_file = thal_cds(), output_folder = out, ...))
        list(res = res, out = out, calls = calls$args)
}

test_that("map_generator_lnc() checks its input", {
        expect_error(map_generator_lnc(file.path(tempdir(), "no-query.fa"), cov_subjects_folder(), tempdir()),
                     "valid path to the 'query_file'")
        expect_error(map_generator_lnc(thal_cds(), cov_subjects_folder(), tempdir(), logical_connective = "XOR"),
                     "logical_connective")
        expect_error(suppressMessages(map_generator_lnc(thal_cds(), new_test_dir("cov_empty_"), tempdir())),
                     "seems to be empty")
})

test_that("map_generator_lnc() runs one search per subject file and writes one table each", {
        subj <- cov_subjects_folder()
        r <- cov_run_lnc(subjects_folder = subj, progress_bar = FALSE, task = "dc-megablast", eval = "1E-3",
                         path = "/opt/blast/bin", max.target.seqs = 5)
        expect_length(r$calls, 2)
        expect_identical(r$calls[[1]]$subject_file, file.path(subj, "Alyrata.fa"))
        expect_identical(r$calls[[2]]$subject_file, file.path(subj, "Bnapus.v2.fa"))
        expect_identical(r$calls[[1]]$task, "dc-megablast")
        expect_identical(r$calls[[1]]$eval, "1E-3")
        expect_identical(r$calls[[1]]$path, "/opt/blast/bin")
        expect_identical(r$calls[[1]]$max.target.seqs, 5)
        expect_setequal(list.files(r$out), c(
                "lncRNA_orthologs_q=ortho_thal_cds.fasta_s=Alyrata.fa_orthodetec=RBH_eval=1E-3.csv",
                "lncRNA_orthologs_q=ortho_thal_cds.fasta_s=Bnapus.v2.fa_orthodetec=RBH_eval=1E-3.csv"))
        # default filter: query coverage >= 30 AND identity >= 30
        expect_identical(r$res$query_id, c("q1", "q1"))
        # the species is the file name up to the first dot
        expect_identical(sort(r$res$species), c("Alyrata", "Bnapus"))
        expect_identical(names(r$res)[1:3], c("species", "query_id", "subject_id"))
})

test_that("map_generator_lnc() filters by coverage, identity and alignment length (AND / OR)", {
        subj <- cov_subjects_folder("Alyrata.fa")
        q <- function(...) cov_run_lnc(subjects_folder = subj, progress_bar = FALSE, ...)$res$query_id
        expect_identical(q(logical_connective = "AND"), "q1")
        expect_identical(q(logical_connective = "OR"), c("q1", "q2", "q3"))
        expect_identical(q(logical_connective = "AND", min_alig_length = 100), "q1")
        expect_identical(q(logical_connective = "OR", min_alig_length = 100), c("q1", "q2"))
        expect_identical(q(logical_connective = "OR", min_alig_length = 100, min_qry_coverage_hsp = 5), c("q1", "q2", "q3"))
        expect_identical(q(min_qry_coverage_hsp = 0, min_qry_perc_identity = 0), c("q1", "q2", "q3", "q4"))
})

test_that("map_generator_lnc() shows a progress bar also for a single subject file, and keeps ids as text", {
        subj <- cov_subjects_folder("Alyrata.fa")
        hits <- cov_lnc_hits()
        hits$query_id <- c("001", "002", "003", "004")
        hits$subject_id <- c("1e5", "2", "3", "4")
        out <- capture.output(r <- cov_run_lnc(subjects_folder = subj, progress_bar = TRUE, hits = hits,
                                               min_qry_coverage_hsp = 0, min_qry_perc_identity = 0, sep = "\t"))
        expect_true(any(grepl("100%", out)))
        expect_identical(r$res$query_id, c("001", "002", "003", "004"))
        expect_identical(r$res$subject_id, c("1e5", "2", "3", "4"))
})

test_that("map_generator_lnc() stops when a table holds a query or a subject twice", {
        subj <- cov_subjects_folder("Alyrata.fa")
        dup_q <- cov_lnc_hits()
        dup_q$query_id[2] <- "q1"
        expect_error(cov_run_lnc(subjects_folder = subj, progress_bar = FALSE, hits = dup_q, min_qry_perc_identity = 0),
                     "non-unique query_ids")
        dup_s <- cov_lnc_hits()
        dup_s$subject_id[2] <- "s1"
        expect_error(cov_run_lnc(subjects_folder = subj, progress_bar = FALSE, hits = dup_s, min_qry_perc_identity = 0),
                     "non-unique subject_ids")
})

# ---------------------------------------------------------------------------------------------------------------

test_that("orthologs_lnc() checks the method and the number of subject files before searching", {
        expect_error(suppressMessages(orthologs_lnc(thal_cds(), lyra_cds(), ortho_detection = "PO")),
                     "orthology detection method that is supported")
        expect_error(suppressMessages(orthologs_lnc(thal_cds(), c(lyra_cds(), thal_cds()), ortho_detection = "BH")),
                     "only defined for pairwise")
        expect_error(suppressMessages(orthologs_lnc(thal_cds(), c(lyra_cds(), thal_cds()), ortho_detection = "RBH")),
                     "only defined for pairwise")
})

test_that("orthologs_lnc() without metablastr says how to install it", {
        skip_if(requireNamespace("metablastr", quietly = TRUE), "metablastr is installed")
        expect_error(suppressMessages(orthologs_lnc(thal_cds(), lyra_cds(), ortho_detection = "BH")),
                     "install_github\\('drostlab/metablastr'\\)")
        expect_error(suppressMessages(orthologs_lnc(thal_cds(), lyra_cds(), ortho_detection = "RBH")),
                     "install_github\\('drostlab/metablastr'\\)")
})

test_that("orthologs_lnc() hands its settings to metablastr's best (reciprocal) hit search", {
        skip_if_not_installed("metablastr")
        seen <- new.env()
        fake <- function(which) function(...) {
                seen[[which]] <- list(...)
                cov_lnc_hits(which)
        }
        testthat::local_mocked_bindings(blast_best_hit = fake("bh"), blast_best_reciprocal_hit = fake("rbh"),
                                        .package = "metablastr")
        res <- suppressMessages(orthologs_lnc(thal_cds(), lyra_cds(), ortho_detection = "BH", task = "megablast",
                                              eval = "1E-10", max.target.seqs = 3, comp_cores = 1,
                                              path = file.path(tempdir(), "blast")))
        expect_identical(res$query_id[1], "bhq1")
        expect_identical(seen$bh$search_type, "nucleotide_to_nucleotide")
        expect_identical(seen$bh$task, "megablast")
        expect_identical(seen$bh$evalue, "1E-10")
        expect_identical(seen$bh$max.target.seqs, 3)
        expect_identical(seen$bh$blast.path, normalizePath(file.path(tempdir(), "blast"), mustWork = FALSE))
        res <- suppressMessages(orthologs_lnc(thal_cds(), lyra_cds()))
        expect_identical(res$query_id[1], "rbhq1")
        expect_null(seen$rbh$blast.path)
        expect_identical(seen$rbh$task, "blastn")
})

test_that("map_generator_lnc() with metablastr and BLAST+ (a small real search)", {
        skip_on_cran()
        skip_if_not_installed("metablastr")
        skip_if_not(nzchar(Sys.which("blastn")) && nzchar(Sys.which("makeblastdb")), "BLAST+ is not installed")
        dir <- new_test_dir("cov_lnc_real_")
        q <- Biostrings::readDNAStringSet(thal_cds())[1:6]
        s <- Biostrings::readDNAStringSet(lyra_cds())[1:6]
        names(q) <- sub(" .*", "", names(q))
        names(s) <- sub(" .*", "", names(s))
        qf <- file.path(dir, "thal.fa")
        Biostrings::writeXStringSet(q, qf)
        subj <- file.path(dir, "subjects")
        dir.create(subj)
        Biostrings::writeXStringSet(s, file.path(subj, "lyra.fa"))
        out <- file.path(dir, "out")
        res <- suppressWarnings(suppressMessages(map_generator_lnc(qf, subj, out, progress_bar = FALSE)))
        expect_s3_class(res, "data.frame")
        expect_true(all(is.element(c("species", "query_id", "subject_id", "qcovhsp", "perc_identity"), names(res))))
        expect_true(all(res$species == "lyra"))
        expect_true(nrow(res) >= 1)
        expect_true(all(res$qcovhsp >= 30 & res$perc_identity >= 30))
        expect_false(anyDuplicated(res$query_id) > 0)
})

# ---------------------------------------------------------------------------------------------------------------

cov_dnds_table <- function() {
        tibble::tibble(query_id = c("a", "b", "c"), subject_id = c("x", "y", "z"), dN = c(0.1, 0.2, 0.3),
                       dS = c(1, 2, 3), dNdS = c(0.1, 0.1, 0.1), perc_identity = c(90, 20, 90),
                       qcovhsp = c(90, 90, 20))
}

test_that("map_generator_dnds() checks its thresholds and the subject folder", {
        expect_error(map_generator_dnds(thal_cds(), new_test_dir("cov_empty_"), tempdir()), "seems to be empty")
        subj <- cov_subjects_folder("Alyrata.fa")
        expect_error(map_generator_dnds(thal_cds(), subj, tempdir(), min_qry_coverage_hsp = 0), "min_qry_coverage_hsp")
        expect_error(map_generator_dnds(thal_cds(), subj, tempdir(), min_qry_perc_identity = 101), "min_qry_perc_identity")
})

test_that("map_generator_dnds() runs dNdS() per subject file and writes the filtered maps", {
        subj <- cov_subjects_folder()
        seen <- new.env()
        seen$calls <- list()
        testthat::local_mocked_bindings(dNdS = function(...) {
                seen$calls[[length(seen$calls) + 1]] <- list(...)
                cov_dnds_table()
        }, .package = "orthologr")
        out <- file.path(new_test_dir("cov_dnds_out_"), "maps")
        capture.output(suppressMessages(map_generator_dnds(thal_cds(), subj, out, eval = "1E-4", comp_cores = 1,
                                                           progress_bar = TRUE, sep = ",", aligner = "diamond")))
        expect_length(seen$calls, 2)
        expect_identical(seen$calls[[1]]$eval, "1E-4")
        expect_identical(seen$calls[[1]]$dnds_est.method, "Comeron")
        expect_identical(seen$calls[[1]]$aligner, "diamond")
        expect_false(seen$calls[[1]]$print_citation)
        expect_identical(seen$calls[[2]]$subject_file, file.path(subj, "Bnapus.v2.fa"))
        files <- list.files(out)
        expect_setequal(files, c("map_q=ortho_thal_cds.fasta_s=Alyrata.fa_orthodetec=RBH_eval=1E-4.csv",
                                 "map_q=ortho_thal_cds.fasta_s=Bnapus.v2.fa_orthodetec=RBH_eval=1E-4.csv"))
        m <- utils::read.csv(file.path(out, files[1]))
        # coverage >= 50 and identity >= 30
        expect_identical(m$query_id, "a")
        # without a progress bar
        out2 <- file.path(new_test_dir("cov_dnds_out_"), "maps")
        expect_silent(suppressMessages(map_generator_dnds(thal_cds(), subj, out2, progress_bar = FALSE)))
        expect_length(list.files(out2), 2)
})
