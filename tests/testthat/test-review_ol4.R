# REVIEW-OL-4: functions that wrote over their input (M2), filter_dNdS()'s threshold (M3), the switches of dNdS() and
# dNdS_coal() checked first (M4), read_dNdS() without its record (L1), IUPAC codes in the check of protein files
# (L2), and tests of rules that had none (L4).

ol4_md5 <- function(f) unname(tools::md5sum(f))

# ---- M2: an output that is the input -----------------------------------------------------------------------------

test_that("translate_cds_to_protein() refuses to write over its input file (M2)", {
        d <- new_test_dir("tr_")
        f <- file.path(d, "cds.fa")
        file.copy(thal_cds(), f)
        before <- ol4_md5(f)
        expect_error(suppressMessages(translate_cds_to_protein(f, f)), "would overwrite the input",
                     class = "orthologr_error_bad_argument")
        old <- setwd(d)
        on.exit(setwd(old), add = TRUE)
        expect_error(suppressMessages(translate_cds_to_protein("cds.fa", file.path(".", "cds.fa"))),
                     class = "orthologr_error_bad_argument")
        expect_identical(ol4_md5(f), before)
        # another output file is written as before
        out <- file.path(d, "aa.fa")
        suppressMessages(translate_cds_to_protein(f, out))
        expect_true(file.exists(out))
        expect_error(translate_cds_to_protein(file.path(d, "missing.fa"), out), class = "orthologr_error_bad_input")
})

test_that("translate_cds_to_protein_all() refuses output_folder = input_folder (M2)", {
        d <- new_test_dir("tra_")
        file.copy(thal_cds(), file.path(d, "a.fa"))
        before <- ol4_md5(file.path(d, "a.fa"))
        expect_error(suppressMessages(translate_cds_to_protein_all(d, d)), "would overwrite",
                     class = "orthologr_error_bad_argument")
        expect_error(suppressMessages(translate_cds_to_protein_all(d, file.path(d, "."))),
                     class = "orthologr_error_bad_argument")
        expect_identical(ol4_md5(file.path(d, "a.fa")), before)
        out <- file.path(new_test_dir("tra_out_"), "aa")
        suppressMessages(translate_cds_to_protein_all(d, out))
        expect_true(file.exists(file.path(out, "a.fa")))
})

test_that("retrieve_longest_isoforms() and _all() refuse to write over their input files (M2)", {
        x <- make_test_ortholog_input()
        p <- x$query$proteome
        before <- ol4_md5(c(p, x$query$gff))
        expect_error(suppressMessages(retrieve_longest_isoforms(p, x$query$gff, p)), "would overwrite the input",
                     class = "orthologr_error_bad_argument")
        expect_error(suppressMessages(retrieve_longest_isoforms(p, x$query$gff, x$query$gff)),
                     class = "orthologr_error_bad_argument")
        expect_identical(ol4_md5(c(p, x$query$gff)), before)
        # _all(): the proteome folder as output folder writes '<organism>.faa' over the proteome 'Athaliana.faa'
        prot <- new_test_dir("prot_")
        ann <- new_test_dir("ann_")
        file.copy(p, prot)
        file.copy(x$query$gff, ann)
        before <- ol4_md5(file.path(prot, basename(p)))
        expect_error(suppressMessages(retrieve_longest_isoforms_all(prot, ann, prot)), class = "orthologr_error_bad_argument")
        expect_identical(ol4_md5(file.path(prot, basename(p))), before)
        # another output folder is written as before
        out <- new_test_dir("iso_out_")
        suppressMessages(retrieve_longest_isoforms_all(prot, ann, out))
        expect_true(file.exists(file.path(out, "Athaliana.faa")))
        expect_error(retrieve_longest_isoforms(p, x$query$gff, tempfile(), annotation_format = "bed"),
                     class = "orthologr_error_bad_argument")
})

# ---- M3: filter_dNdS() -------------------------------------------------------------------------------------------

test_that("filter_dNdS() needs one number as threshold; a pair at the threshold is kept (M3, Q2)", {
        t <- data.frame(query_id = c("a", "b", "c"), subject_id = c("x", "y", "z"), dN = c(0.2, 0.1, 1),
                        dS = c(0.1, 0.1, 0.1), dNdS = c(2, 1, 10))
        expect_identical(nrow(suppressMessages(filter_dNdS(t, 2))), 2L)
        for (bad in list("2", NULL, NA_real_, c(1, 2), TRUE))
                expect_error(suppressMessages(filter_dNdS(t, bad)), "dnds.threshold must be one number",
                             class = "orthologr_error_bad_argument")
        expect_identical(nrow(suppressMessages(filter_dNdS(t, 2L))), 2L)
})

# ---- M4: the switches of dNdS() and dNdS_coal() ------------------------------------------------------------------

test_that("dNdS() and dNdS_coal() check their switches before any work (M4)", {
        for (a in list(list(clean_folders = NA), list(quiet = NA), list(store_locally = "yes"),
                       list(delete_corrupt_cds = NA), list(print_citation = NA), list(quiet = c(TRUE, FALSE)))) {
                t0 <- proc.time()[["elapsed"]]
                e <- tryCatch(do.call(dNdS, c(list(thal_cds(), lyra_cds()), a)), error = function(e) e)
                expect_s3_class(e, "orthologr_error_bad_argument")
                expect_match(conditionMessage(e), paste0(names(a), " must be TRUE or FALSE"), fixed = TRUE)
                expect_lt(proc.time()[["elapsed"]] - t0, 5)
        }
        for (a in list(list(quiet = NA), list(delete_corrupt_cds = "no"), list(print_citation = NA))) {
                e <- tryCatch(do.call(dNdS_coal, c(list(thal_cds(), lyra_cds()), a)), error = function(e) e)
                expect_s3_class(e, "orthologr_error_bad_argument")
                expect_match(conditionMessage(e), paste0(names(a), " must be TRUE or FALSE"), fixed = TRUE)
        }
        # several values for ortho_detection or aligner: a classed error, not R's "the condition has length > 1"
        expect_error(dNdS(thal_cds(), lyra_cds(), ortho_detection = c("RBH", "BH")), "ortho_detection",
                     class = "orthologr_error_bad_argument")
        expect_error(dNdS(thal_cds(), lyra_cds(), ortho_detection = c("PO", "DELTA")), class = "orthologr_error_bad_argument")
        expect_error(dNdS_coal(thal_cds(), lyra_cds(), ortho_detection = c("RBH", "BH")), class = "orthologr_error_bad_argument")
        expect_error(dNdS_coal(thal_cds(), lyra_cds(), aligner = c("diamond", "blast")), class = "orthologr_error_bad_argument")
})

# ---- L1: read_dNdS() without its record --------------------------------------------------------------------------

test_that("read_dNdS() without the record keeps IDs and other text columns as text (L1)", {
        f <- tempfile(fileext = ".tsv")
        writeLines(c("orthogroup\tquery_id\tsubject_id\tmethod\tdN\tdS\tdNdS\tevalue\tstatus",
                     "1\t007\t1e3\tYN\t0.1\t0.2\t0.5\t1e-30\tok",
                     "2\tNaN\t42\tYN\tNA\tNaN\tNA\t2.5e-10\tsaturated"), f)
        expect_warning(read_dNdS(f), "without the record")
        x <- suppressWarnings(read_dNdS(f))
        expect_identical(x$query_id, c("007", "NaN"))
        expect_identical(x$subject_id, c("1e3", "42"))
        expect_identical(x$orthogroup, c("1", "2"))
        expect_identical(x$method, c("YN", "YN"))
        expect_identical(x$status, c("ok", "saturated"))
        expect_identical(x$dN, c(0.1, NA))
        expect_true(is.nan(x$dS[2]))
        expect_identical(x$evalue, c(1e-30, 2.5e-10))
})

# ---- L2: IUPAC codes count as nucleotides ------------------------------------------------------------------------

test_that("coding sequences with IUPAC ambiguity codes are refused as proteins; proteins are read (L2)", {
        f <- tempfile(fileext = ".fasta")
        writeLines(c(">c1", paste(rep("atggctaaaggcRYKMSW", 10), collapse = "")), f)
        expect_error(read.proteome(f, format = "fasta"), "nucleotide codes", class = "orthologr_error_bad_input")
        rlang::local_options(orthologr.check_protein = FALSE)
        expect_identical(read.proteome(f, format = "fasta")$geneids, "c1")
        rlang::local_options(orthologr.check_protein = NULL)
        expect_identical(nrow(read.proteome(thal_aa(), format = "fasta")), 20L)
        expect_identical(nrow(read.proteome(lyra_aa(), format = "fasta")), 20L)
})

# ---- L4: rules without a test ------------------------------------------------------------------------------------

test_that("an in-paralog is attached to one seed only: the highest score, of equal ones the first gene (L4)", {
        genes <- data.table::data.table(species = c("A", "A", "A", "B", "B"))
        rbh <- data.table::data.table(q = c(1L, 2L), s = c(4L, 5L))
        hits <- function(s23) data.table::data.table(
                q = c(1L, 4L, 2L, 5L, 3L, 1L, 3L, 2L, 3L),
                s = c(4L, 1L, 5L, 2L, 4L, 3L, 1L, 3L, 2L),
                bit_score = c(100, 100, 100, 100, 50, 300, 300, s23, s23))
        # gene 3 (no reciprocal best hit) hits both seeds of its species above their best between-species score
        expect_identical(as.data.frame(orthogroups_inparalogs(hits(250), genes, rbh)), data.frame(q = 1L, s = 3L))
        expect_identical(as.data.frame(orthogroups_inparalogs(hits(400), genes, rbh)), data.frame(q = 2L, s = 3L))
        expect_identical(as.data.frame(orthogroups_inparalogs(hits(300), genes, rbh)), data.frame(q = 1L, s = 3L))
        # below the between-species score of the seed: not attached
        expect_identical(nrow(orthogroups_inparalogs(hits(90)[-(6:7)], genes, rbh)), 0L)
})

test_that("the note on saturation comes above 20 % saturated pairs or a median dS above 1.5, from 10 pairs on (L4)", {
        st <- function(n_sat, n) c(rep("saturated", n_sat), rep("ok", n - n_sat))
        expect_false(dnds_saturation_hint(st(2, 10), 1))
        expect_true(dnds_saturation_hint(st(3, 10), 1))
        expect_false(dnds_saturation_hint(st(9, 9), 1))
        expect_false(dnds_saturation_hint(st(0, 10), 1.5))
        expect_true(dnds_saturation_hint(st(0, 10), 1.51))
        expect_false(dnds_saturation_hint(st(0, 10), NA))
        expect_match(saturation_hint_text(3, 10, 1.2), "median of the dS values that were estimated is 1.2", fixed = TRUE)
})

test_that("best_hits(rule = \"bitscore\") refuses a hit table with missing bit scores (L4)", {
        h <- data.frame(query_id = "a", subject_id = c("x", "y"), bit_score = c(NA, 5), evalue = 1e-5, alig_length = 10L,
                        perc_identity = 90)
        expect_error(best_hits(h, rule = "bitscore"), "missing bit scores", class = "orthologr_error_bad_input")
        h$bit_score <- c(4, 5)
        expect_identical(best_hits(h, rule = "bitscore")$subject_id, "y")
})
