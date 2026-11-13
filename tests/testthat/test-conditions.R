# Classed conditions (?orthologr_options, section "Errors and warnings"): the common kinds of failure can be told
# apart by class, without matching message text.

expect_orthologr_error <- function(expr, kind) {
        e <- tryCatch(force(expr), error = function(e) e)
        expect_s3_class(e, paste0("orthologr_error_", kind))
        expect_s3_class(e, "orthologr_error")
        expect_null(conditionCall(e))
        invisible(e)
}

test_that("orthologr_abort() / orthologr_warn() give classed conditions with the message of stop() / warning()", {
        e <- tryCatch(orthologr_abort("a", 1, c("b", "c"), class = c("alignment_failed", "tool_failed")), error = function(e) e)
        expect_identical(conditionMessage(e), tryCatch(stop("a", 1, c("b", "c"), call. = FALSE), error = conditionMessage))
        expect_identical(class(e)[1:3], c("orthologr_error_alignment_failed", "orthologr_error_tool_failed", "orthologr_error"))
        expect_warning(orthologr_warn("x", class = "isoforms"), class = "orthologr_warning_isoforms")
        expect_warning(orthologr_warn("x"), class = "orthologr_warning")
        expect_true(all(c("tool_missing", "no_hits", "bad_input", "alignment_failed") %in% orthologr_condition_classes()$error))
})

test_that("a missing program is orthologr_error_tool_missing", {
        expect_orthologr_error(tool_available("orthologr_no_such_program_xyz", "--version", "version"), "tool_missing")
        e <- expect_orthologr_error(check_external(127L, tool = "PAL2NAL", context = "gene pair 'a' / 'b'"), "tool_missing")
        expect_s3_class(e, "orthologr_error_alignment_failed")
        old <- Sys.getenv("PATH")
        on.exit({ Sys.setenv(PATH = old); tool_cache_reset() }, add = TRUE)
        Sys.setenv(PATH = new_test_dir("cond_empty_"))
        tool_cache_reset()
        expect_orthologr_error(dnds_require_programs(aa_aln_tool = "NW", dnds_est.method = "Comeron"), "tool_missing")
})

test_that("a missing R package is orthologr_error_package_missing", {
        expect_orthologr_error(need_package("orthologrNoSuchPackage", "this test"), "package_missing")
})

test_that("a failed program is orthologr_error_tool_failed; an aligner also _alignment_failed", {
        e <- expect_orthologr_error(check_external(1L, tool = "MAFFT"), "alignment_failed")
        expect_s3_class(e, "orthologr_error_tool_failed")
        e <- expect_orthologr_error(check_external(2L, tool = "DIAMOND"), "tool_failed")
        expect_false(inherits(e, "orthologr_error_alignment_failed"))
        out <- tempfile()
        expect_orthologr_error(check_external(0L, tool = "PAL2NAL", output = out), "alignment_failed")
})

test_that("bad input files and bad arguments have their classes", {
        prot <- tempfile(fileext = ".fasta")
        writeLines(c(">p1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAACIBVM"), prot)
        expect_orthologr_error(read.cds(prot, format = "fasta"), "bad_input")
        expect_orthologr_error(check_duplicate_ids(c("a", "b", "a"), file = "x.fasta"), "bad_input")
        expect_orthologr_error(check_sensitivity_mode("nonsense"), "bad_argument")
        expect_orthologr_error(check_comp_cores(0), "bad_argument")
        expect_orthologr_error(check_evalue("1e-5; echo hi"), "bad_argument")
})

test_that("no hits is orthologr_error_no_hits; failed pairs and isoforms are classed warnings", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        set.seed(1)
        rnd <- function(n) paste(sample(c("A", "C", "G", "T"), 3 * n, replace = TRUE), collapse = "")
        cds <- function(n) { x <- rnd(n); gsub("TAA|TAG|TGA", "GCC", x) }
        d <- new_test_dir("cond_nohits_")
        q <- file.path(d, "q.fasta"); s <- file.path(d, "s.fasta")
        writeLines(c(">q1", paste0("ATG", cds(150), "TAA"), ">q2", paste0("ATG", cds(120), "TAA")), q)
        writeLines(c(">s1", paste0("ATG", cds(150), "TAA"), ">s2", paste0("ATG", cds(120), "TAA")), s)
        expect_orthologr_error(suppressMessages(dNdS(q, s, print_citation = FALSE)), "no_hits")
})

test_that("a gene pair that fails gives orthologr_warning_failed_pairs", {
        skip_if_not(perl_available, "perl is not installed")
        cds_a <- "ATGGCTAAAGGTCTTGCAACTGGTCTTGCAACTGGTAAA"
        cds_b <- "ATGGCCAAAGGCCTGGCTACCGGTTTGGCAACAGGCAAA"
        tbl <- data.table::data.table(query_id = c("q1", "q2"), subject_id = c("s1", "s2"),
                                      query_cds = c(cds_a, cds_a), subject_cds = c(cds_b, cds_b),
                                      query_aa = c("MAKGLATGLATGK", "MAKGLWWWWWWATGK"), subject_aa = c("MAKGLATGLATGK", "MAKGLATGLATGK"))
        expect_warning(compute_dnds(tbl, aa_aln_type = "pairwise", aa_aln_tool = "NW", dnds_est.method = "Comeron",
                                    quiet = TRUE, comp_cores = 1),
                       class = "orthologr_warning_failed_pairs")
})
