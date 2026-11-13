diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))
kaks_available <- isTRUE(tryCatch(suppressWarnings(is_installed_kaks_calculator()), error = function(e) FALSE))
perl_available <- nzchar(Sys.which("perl"))

query_cds   <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")

write_fasta <- function(seqs, ids, file = tempfile(fileext = ".fasta")) {
        writeLines(as.vector(rbind(paste0(">", ids), seqs)), file)
        file
}

cds_a <- "ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACC"
cds_b <- "ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTGAAGGCGTGGTTGACC"
aa_a  <- "MKTAYIAKQRLEGDLKAWLT"
aa_b  <- "MKTAYIAKQRLDGDLKAWLT"

# ---------------------------------------------------------------- exit status of external programs

test_that("check_external() names the program and what it ran for ...", {
        expect_true(check_external(0L, tool = "PAL2NAL"))
        expect_error(check_external(1L, tool = "PAL2NAL", context = "gene pair 'A' / 'B' (q3)"),
                     "PAL2NAL failed for gene pair 'A' / 'B' \\(q3\\): exit status 1")
        expect_error(check_external(127L, tool = "KaKs_Calculator", context = "gene pair 'A' / 'B'"),
                     "not found")
        log <- tempfile()
        writeLines(c("reading", "#---  ERROR: inconsistency between the following pep and nuc seqs  ---#", ">A"), log)
        expect_error(check_external(2L, tool = "PAL2NAL", log = log), "PAL2NAL said: #---  ERROR: inconsistency")
        # status 0 but no output
        out <- tempfile()
        expect_error(check_external(0L, tool = "PAL2NAL", context = "gene pair 'A' / 'B'", output = out),
                     "PAL2NAL wrote no output for gene pair 'A' / 'B'")
        file.create(out)
        expect_error(check_external(0L, tool = "PAL2NAL", output = out), "wrote no output")
        writeLines("x", out)
        expect_true(check_external(0L, tool = "PAL2NAL", output = out))
})

test_that("codon_aln() stops with a message naming the pair when PAL2NAL cannot build the alignment ...", {
        skip_if_not(perl_available, "perl is not installed")
        # the protein of gene-B is longer than its coding sequence: PAL2NAL prints an ERROR, writes nothing,
        # and exits with status 0
        aln <- write_fasta(c(aa_a, paste0(aa_b, "WWWWWWWWWW")), c("gene-A", "gene-B"))
        nuc <- write_fasta(c(cds_a, cds_b), c("gene-A", "gene-B"))
        expect_error(
                codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", codon_aln_name = "bad_pair", quiet = TRUE),
                "PAL2NAL (failed|wrote no output) for gene pair 'gene-A' / 'gene-B' \\(bad_pair\\)"
        )
        expect_error(
                codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", codon_aln_name = "bad_pair", quiet = TRUE),
                "inconsistency"
        )
        # a good pair still works, also when an output of that name exists already
        aln_ok <- write_fasta(c(aa_a, aa_b), c("gene-A", "gene-B"))
        for (k in 1:2) {
                codon_aln(file_aln = aln_ok, file_nuc = nuc, format = "fasta", codon_aln_name = "good_pair", quiet = TRUE)
                out <- file.path(tempdir(), "_alignment", "codon_aln", "good_pair_pal2nal.aln")
                expect_true(file.exists(out))
                expect_equal(fasta_ids(out), c("gene-A", "gene-B"))
        }
        # a stale output of a previous call is not left behind when the call fails
        file.copy(out, file.path(tempdir(), "_alignment", "codon_aln", "bad_pair_pal2nal.aln"), overwrite = TRUE)
        expect_error(codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", codon_aln_name = "bad_pair", quiet = TRUE))
        stale <- file.path(tempdir(), "_alignment", "codon_aln", "bad_pair_pal2nal.aln")
        expect_true(!file.exists(stale) || file.size(stale) == 0)
})

test_that("substitutionrate() reports a KaKs_Calculator that cannot be run ...", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        f <- write_fasta(c(cds_a, cds_b), c("gene-A", "gene-B"))
        expect_error(
                suppressWarnings(substitutionrate(f, est.method = "YN", quiet = TRUE, subst_name = "q3",
                                                  kaks_calc_path = file.path(tempdir(), "no_such_folder"))),
                "KaKs_Calculator failed for gene pair 'gene-A' / 'gene-B' \\(q3\\): exit status"
        )
})

# ---------------------------------------------------------------- one pair does not end the run

test_that("compute_dnds() returns an NA row with a reason for a pair that fails, and finishes the others ...", {
        skip_if_not(perl_available, "perl is not installed")
        tbl <- data.table::data.table(
                query_id    = c("q-ok1", "q-bad", "q-ok2"),
                subject_id  = c("s-ok1", "s-bad", "s-ok2"),
                query_cds   = c(cds_a, cds_a, cds_b),
                subject_cds = c(cds_b, cds_b, cds_a),
                # the protein of 'q-bad' has 12 residues that its CDS does not encode
                query_aa    = c(aa_a, paste0(substr(aa_a, 1, 10), "WWWWWWWWWWWW", substr(aa_a, 11, 20)), aa_b),
                subject_aa  = c(aa_b, aa_b, aa_a)
        )
        for (cores in c(1, 2)) {
                skip_if(cores > parallel::detectCores(), "one core only")
                expect_warning(
                        res <- compute_dnds(tbl, aa_aln_type = "pairwise", aa_aln_tool = "NW", codon_aln_tool = "pal2nal",
                                            dnds_est.method = "Comeron", quiet = TRUE, comp_cores = cores),
                        "1 of 3 gene pairs could not be processed"
                )
                expect_equal(nrow(res), 3)
                expect_setequal(res$query_id, c("q-ok1", "q-bad", "q-ok2"))
                bad <- res[res$query_id == "q-bad", ]
                expect_equal(bad$subject_id, "s-bad")
                expect_true(is.na(bad$dN) && is.na(bad$dS) && is.na(bad$dNdS))
                expect_true(is.numeric(res$dN) && is.numeric(res$dS) && is.numeric(res$dNdS))
                expect_false(anyNA(res$dN[res$query_id != "q-bad"]))
                failed <- attr(res, "failed_pairs")
                expect_s3_class(failed, "data.frame")
                expect_equal(failed$query_id, "q-bad")
                expect_equal(failed$subject_id, "s-bad")
                expect_match(failed$reason, "PAL2NAL")
                expect_match(failed$reason, "'q-bad' / 's-bad'")
        }
        # nothing failed: no warning, an empty table of failures
        expect_silent(res <- compute_dnds(tbl[c(1, 3), ], aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                          dnds_est.method = "Comeron", quiet = TRUE, comp_cores = 1))
        expect_equal(nrow(attr(res, "failed_pairs")), 0)
        expect_false(anyNA(res$dNdS))
})

test_that("an NA row has the columns of the method's result ...", {
        expect_equal(names(dnds_na_row("a", "b", "Comeron", "why")), c("query_id", "subject_id", "dN", "dS", "dNdS", ".reason"))
        expect_equal(names(dnds_na_row("a", "b", "YN", "why")), c("query_id", "subject_id", "dN", "dS", "dNdS", "method", ".reason"))
        expect_true(is.numeric(dnds_na_row("a", "b", "Li", "why")$dN))
})

# ---------------------------------------------------------------- one working folder per call

test_that("dNdS() works in a folder of its own and removes it ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        before <- list.files(tempdir(), pattern = "^orthologr_dNdS_")
        res <- suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, comp_cores = 1,
                                     quiet = TRUE, print_citation = FALSE))
        expect_gt(nrow(res), 0)
        expect_equal(list.files(tempdir(), pattern = "^orthologr_dNdS_"), before)
        expect_null(getOption("orthologr.tmpdir"))
        expect_equal(nrow(attr(res, "failed_pairs")), 0)
        # the folder is also removed, and the option reset, when the call fails
        expect_error(suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds,
                                           ortho_detection = "INVALID")))
        expect_equal(list.files(tempdir(), pattern = "^orthologr_dNdS_"), before)
        expect_null(getOption("orthologr.tmpdir"))
        # options(orthologr.keep_tmp = TRUE) keeps it
        old <- options(orthologr.keep_tmp = TRUE)
        on.exit(options(old), add = TRUE)
        res2 <- suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, comp_cores = 1,
                                      quiet = TRUE, print_citation = FALSE))
        kept <- setdiff(list.files(tempdir(), pattern = "^orthologr_dNdS_"), before)
        expect_length(kept, 1)
        expect_true(file.exists(file.path(tempdir(), kept, "_alignment", "codon_aln", "q1_pal2nal.aln")))
        unlink(file.path(tempdir(), kept), recursive = TRUE)
        expect_equal(as.data.frame(res2), as.data.frame(res), ignore_attr = TRUE)
})

test_that("orthologr_tmpdir() is tempdir() unless a call set its own folder ...", {
        expect_equal(orthologr_tmpdir(), tempdir())
        end <- orthologr_tmpdir_start("orthologr_test_")
        d <- orthologr_tmpdir()
        expect_true(dir.exists(d))
        expect_false(identical(d, tempdir()))
        expect_equal(dirname(d), tempdir())
        end()
        expect_false(dir.exists(d))
        expect_equal(orthologr_tmpdir(), tempdir())
})

# ---------------------------------------------------------------- gene ids with a hyphen (KaKs_Calculator methods)

test_that("gene ids with hyphens survive the KaKs_Calculator result table ...", {
        # what KaKs_Calculator writes for a pair (the label used to be '<query>-<subject>' and was split at '-')
        kaks_tbl <- data.frame(Sequence = "AT1G01010-1-Al-scaffold-0001_1", Method = "YN", Ka = 0.0340338,
                               Ks = 0.484553, Ka.Ks = 0.0702375, stringsAsFactors = FALSE)
        res <- kaks_result_table(kaks_tbl, c("AT1G01010-1", "Al-scaffold-0001_1"))
        expect_equal(res$query_id, "AT1G01010-1")
        expect_equal(res$subject_id, "Al-scaffold-0001_1")
        expect_equal(res$dN, 0.0340338)
        expect_equal(as.character(res$method), "YN")
        # NA values of KaKs_Calculator stay NA, and the columns stay numeric
        kaks_tbl$Ks <- "NA"
        kaks_tbl$Ka.Ks <- "NA"
        res <- kaks_result_table(kaks_tbl, c("a-1", "b-2"))
        expect_true(is.numeric(res$dS) && is.na(res$dS) && is.na(res$dNdS))
        expect_error(kaks_result_table(kaks_tbl[0, ], c("a", "b"), "gene pair 'a' / 'b'"), "no result for gene pair 'a' / 'b'")
})

test_that("substitutionrate() keeps gene ids with hyphens for KaKs_Calculator methods ...", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        f <- write_fasta(c(cds_a, cds_b), c("AT1G01010-1", "Al-scaffold-0001_1"))
        res <- substitutionrate(f, est.method = "YN", quiet = TRUE, subst_name = "hyphen")
        expect_equal(res$query_id, "AT1G01010-1")
        expect_equal(res$subject_id, "Al-scaffold-0001_1")
        expect_true(is.numeric(res$dN) && is.finite(res$dN))
        ng <- substitutionrate(f, est.method = "NG", quiet = TRUE, subst_name = "hyphen")
        expect_equal(as.character(ng$method), "NG")
        expect_equal(ng$query_id, "AT1G01010-1")
})

test_that("dNdS() returns every pair when the gene ids contain hyphens (KaKs_Calculator method) ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        # the example files with '.' in the gene ids replaced by '-' (AT1G01010.1 -> AT1G01010-1)
        hyphenate <- function(file) {
                x <- readLines(file)
                h <- startsWith(x, ">")
                x[h] <- gsub(".", "-", sub("\\s.*$", "", x[h]), fixed = TRUE)
                out <- tempfile(fileext = ".fasta")
                writeLines(x, out)
                out
        }
        q <- hyphenate(query_cds)
        res <- suppressMessages(dNdS(query_file = q, subject_file = subject_cds, dnds_est.method = "YN",
                                     comp_cores = 1, quiet = TRUE, print_citation = FALSE))
        ref <- suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, dnds_est.method = "YN",
                                     comp_cores = 1, quiet = TRUE, print_citation = FALSE))
        expect_equal(nrow(res), nrow(ref))
        expect_gt(nrow(res), 0)
        expect_true(all(grepl("-", res$query_id, fixed = TRUE)))
        expect_equal(sort(gsub("-", ".", res$query_id, fixed = TRUE)), sort(ref$query_id))
        m <- match(gsub("-", ".", res$query_id, fixed = TRUE), ref$query_id)
        expect_equal(res$subject_id, ref$subject_id[m])
        expect_equal(res$dN, ref$dN[m])
        expect_equal(res$dS, ref$dS[m])
})

# ---------------------------------------------------------------- query and subject files of the same name

test_that("dNdS() gives the same result when both files have the same name in different folders ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        da <- tempfile("same_a_")
        db <- tempfile("same_b_")
        dir.create(da)
        dir.create(db)
        file.copy(query_cds, file.path(da, "cds.fa"))
        file.copy(subject_cds, file.path(db, "cds.fa"))
        same <- suppressMessages(dNdS(query_file = file.path(da, "cds.fa"), subject_file = file.path(db, "cds.fa"),
                                      comp_cores = 1, quiet = TRUE, print_citation = FALSE))
        ref <- suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds,
                                     comp_cores = 1, quiet = TRUE, print_citation = FALSE))
        expect_gt(nrow(ref), 0)
        expect_equal(nrow(same), nrow(ref))
        expect_equal(same$query_id, ref$query_id)
        expect_equal(same$subject_id, ref$subject_id)
        expect_equal(same$dN, ref$dN)
        expect_equal(same$dS, ref$dS)
        expect_equal(same$dNdS, ref$dNdS)
        expect_equal(nrow(attr(same, "failed_pairs")), 0)
})

# ---------------------------------------------------------------- Li's method: numbers, compared as numbers

test_that("Li's method returns numeric columns and compares numbers, not strings ...", {
        # a tiny value is not saturation ("9.95e-05" > "9.9" is TRUE for strings)
        v <- li_values(ka = 9.95e-05, ks = 9.99e-03)
        expect_equal(v$dN, 9.95e-05)
        expect_equal(v$dS, 9.99e-03)
        expect_equal(v$dNdS, 9.95e-05 / 9.99e-03)
        # seqinr's markers: about 10 = saturation, negative = not computable
        expect_true(all(is.na(unlist(li_values(ka = 9.999999, ks = 0.5)[c("dN", "dNdS")]))))
        expect_equal(li_values(ka = 9.999999, ks = 0.5)$dS, 0.5)
        expect_true(all(is.na(unlist(li_values(ka = 0.1, ks = 10)[c("dS", "dNdS")]))))
        expect_true(all(is.na(unlist(li_values(ka = -1, ks = 0.2)[c("dN", "dNdS")]))))
        expect_true(is.na(li_values(ka = 0.1, ks = 0)$dNdS))
        expect_true(is.na(li_values(ka = 0, ks = 0)$dNdS))
        
        f <- write_fasta(c(cds_a, cds_b), c("gene-A-1", "gene-B.2-x"))
        res <- substitutionrate(f, est.method = "Li", quiet = TRUE)
        expect_true(is.numeric(res$dN) && is.numeric(res$dS) && is.numeric(res$dNdS))
        expect_equal(res$query_id, "gene-A-1")
        expect_equal(res$subject_id, "gene-B.2-x")
        ref <- seqinr::kaks(seqinr::read.alignment(f, format = "fasta"))
        expect_equal(res$dN, as.numeric(ref$ka))
        expect_equal(res$dS, as.numeric(ref$ks))
        expect_equal(res$dNdS, as.numeric(ref$ka) / as.numeric(ref$ks))
        
        # two unrelated sequences: saturated, NA instead of seqinr's 9.999999
        set.seed(1)
        rnd <- function() paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
        f2 <- write_fasta(c(rnd(), rnd()), c("x", "y"))
        sat <- substitutionrate(f2, est.method = "Li", quiet = TRUE)
        expect_true(is.numeric(sat$dN))
        expect_true(is.na(sat$dNdS))
})
