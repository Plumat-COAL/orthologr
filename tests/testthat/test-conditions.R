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
        Sys.setenv(PATH = with_coreutils(new_test_dir("cond_empty_")))
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

test_that("sensitivity_hint(): a message for distant species with sensitivity_mode = 'fast' only", {
        hits <- data.frame(perc_identity = rep(c(50, 60), 15))
        expect_message(sensitivity_hint(hits, "diamond", "fast"), "sensitivity_mode = \"default\"")
        expect_silent(sensitivity_hint(hits, "diamond", "default"))
        expect_silent(sensitivity_hint(hits, "blast", "fast"))
        expect_silent(sensitivity_hint(data.frame(perc_identity = rep(95, 30)), "diamond", "fast"))
        expect_silent(sensitivity_hint(hits[1:5, , drop = FALSE], "diamond", "fast"))
})

# One row per documented kind of error (?orthologr_options, "Errors and warnings"): the argument and input-file
# checks of the main functions. All of them fail before any external program is needed.
test_that("argument and input-file errors of the main functions have the documented class", {
        q <- thal_cds()
        s <- lyra_cds()
        nope <- file.path(new_test_dir("cond_nope_"), "no_such_file.fasta")
        probes <- list(
                list("dNdS codon_aln_tool",         quote(dNdS(q, s, codon_aln_tool = "x")),                  "bad_argument"),
                list("dNdS dnds_est.method",        quote(dNdS(q, s, dnds_est.method = "x")),                 "bad_argument"),
                list("dNdS two methods",            quote(dNdS(q, s, dnds_est.method = c("YN", "NG"))),       "bad_argument"),
                list("dNdS best_hit_rule",          quote(dNdS(q, s, best_hit_rule = "x")),                   "bad_argument"),
                list("dNdS comp_cores",             quote(dNdS(q, s, comp_cores = 0)),                        "bad_argument"),
                list("dNdS aa_sub_matrix",          quote(dNdS(q, s, aa_sub_matrix = "FOO")),                 "bad_argument"),
                list("dNdS aligner",                quote(dNdS(q, s, aligner = "x")),                         "bad_argument"),
                list("dNdS missing query file",     quote(dNdS(nope, s)),                                     "bad_input"),
                list("dNdS missing subject file",   quote(dNdS(q, nope)),                                     "bad_input"),
                list("read.cds missing file",       quote(read.cds(nope, format = "fasta")),                  "bad_input"),
                list("read.cds format",             quote(read.cds(q, format = "x")),                         "bad_argument"),
                list("read.proteome missing file",  quote(read.proteome(nope, format = "fasta")),             "bad_input"),
                list("multi_aln tool",              quote(multi_aln(nope, tool = "x")),                       "bad_argument"),
                list("pairwise_aln tool",           quote(pairwise_aln(nope, tool = "x")),                    "bad_argument"),
                list("filter_dNdS dS_zero",         quote(filter_dNdS(data.frame(), dS_zero = "x")),          "bad_argument"),
                list("orthogroups one file",        quote(orthogroups(c(a = q))),                             "bad_argument"),
                list("orthogroups missing file",    quote(orthogroups(c(a = q, b = nope))),                   "bad_input"),
                list("read_dNdS missing file",      quote(read_dNdS(nope)),                                   "bad_input")
        )
        # the COAL engine (needs no program for these checks either)
        probes <- c(probes, list(
                list("dNdS engine",                 quote(dNdS(q, s, engine = "x")),                          "bad_argument"),
                list("dNdS kaks_tool",              quote(dNdS(q, s, kaks_tool = "x")),                       "bad_argument"),
                list("dNdS needle with native",     quote(dNdS(q, s, aa_aln_scoring = "needle")),             "bad_argument"),
                list("dNdS_coal aligner",           quote(dNdS_coal(q, s, aligner = "x")),                    "bad_argument"),
                list("dNdS_coal missing file",      quote(dNdS_coal(nope, s)),                                "bad_input"),
                list("dNdS_clusters missing file",  quote(dNdS_clusters(nope)),                               "bad_input"),
                list("dNdS_clusters aa_aln_tool",   quote(dNdS_clusters(q, aa_aln_tool = "x")),               "bad_argument"),
                list("COAL backend option",         quote(local({ old <- options(orthologr.coal_backend = "x")
                                                                  on.exit(options(old))
                                                                  dNdS(q, s, codon_aln_tool = "coal") })),    "bad_argument"),
                list("COAL not there",              quote(local({ old <- options(orthologr.coal_backend = "none")
                                                                  on.exit(options(old))
                                                                  dNdS(q, s, dnds_est.method = "YN", kaks_tool = "coal") })), "package_missing")
        ))
        if (requireNamespace("coalr", quietly = TRUE))
                probes <- c(probes, list(
                        list("dNdS coal method GY", quote(dNdS(q, s, engine = "coal", dnds_est.method = "GY")), "bad_argument")))
        for (p in probes) {
                e <- tryCatch(suppressMessages(eval(p[[2]])), error = function(e) e)
                expect_s3_class(e, paste0("orthologr_error_", p[[3]]))
                if (!inherits(e, paste0("orthologr_error_", p[[3]])))
                        fail(paste0(p[[1]], ": ", paste(class(e), collapse = ", "), ": ", conditionMessage(e)))
        }
        expect_true(all(vapply(probes, `[[`, "", 3) %in% orthologr_condition_classes()$error))
})

test_that("DIAMOND: an old version warns, a program that crashes is tool_failed, a failed search is tool_failed (O14)", {
        skip_on_os("windows")
        bin <- new_test_dir("old_diamond_")
        writeLines(c("#!/bin/sh", "echo 'diamond version 0.9.24'"), file.path(bin, "diamond"))
        Sys.chmod(file.path(bin, "diamond"), "755")
        tool_cache_reset()
        on.exit(tool_cache_reset(), add = TRUE)
        expect_warning(v <- diamond_version(path = bin), class = "orthologr_warning_tool_version")
        expect_match(v, "0.9.24")
        expect_silent(diamond_version(path = bin))           # once per session and program
        crash <- new_test_dir("crash_diamond_")
        writeLines(c("#!/bin/sh", "echo 'Segmentation fault' >&2", "exit 139"), file.path(crash, "diamond"))
        Sys.chmod(file.path(crash, "diamond"), "755")
        e <- expect_orthologr_error(diamond_version(path = crash), "tool_failed")
        expect_match(conditionMessage(e), "could not be started.*status 139")
        expect_false(inherits(e, "orthologr_error_tool_missing"))
        expect_orthologr_error(diamond_version(path = new_test_dir("no_diamond_")), "tool_missing")
        # a version 2.1 or newer: no warning
        new <- new_test_dir("new_diamond_")
        writeLines(c("#!/bin/sh", "echo 'diamond version 2.1.8'"), file.path(new, "diamond"))
        Sys.chmod(file.path(new, "diamond"), "755")
        expect_silent(diamond_version(path = new))
})

test_that("comp_cores above the number of cores is a classed warning (O14)", {
        cores <- parallel::detectCores()
        skip_if(is.na(cores))
        expect_warning(n <- check_comp_cores(cores + 5), class = "orthologr_warning_comp_cores")
        expect_identical(as.numeric(n), as.numeric(cores))
})
