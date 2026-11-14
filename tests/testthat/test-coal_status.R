# The column 'status', the summary message and the record of the run with the COAL engine: the same words and
# the same structure as with the default engine (test-run_summary.R).

coal_status_possible <- function() requireNamespace("coalr", quietly = TRUE) && diamond_available

coal_messages <- function(expr) {
        msgs <- character(0)
        out <- utils::capture.output(
                value <- withCallingHandlers(expr, message = function(m) {
                        msgs <<- c(msgs, conditionMessage(m))
                        invokeRestart("muffleMessage")
                })
        )
        list(value = value, messages = msgs, stdout = out)
}

test_that("dNdS(engine = \"coal\") ends with the column 'status', the summary and the record of the run", {
        skip_if_not(coal_status_possible(), "coalr or DIAMOND is not installed")
        run <- coal_messages(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "YN", engine = "coal", print_citation = FALSE))
        res <- run$value
        expect_identical(names(res)[ncol(res)], "status")
        expect_true(all(is.element(res$status, orthologr:::dnds_status_levels())))
        expect_identical(res$status == "ok", is.finite(res$dNdS))
        expect_identical(run$stdout, character(0))
        expect_lte(length(run$messages), 3)
        expect_match(run$messages[1], "^dNdS_coal\\(\\): reciprocal best hits of 'ortho_thal_cds.fasta'")
        expect_match(run$messages[length(run$messages)], paste0("dNdS_coal(): ", nrow(res), " gene pairs - "), fixed = TRUE)
        rec <- attr(res, "orthologr_run")
        expect_identical(rec[["function"]], "dNdS_coal")
        expect_identical(rec$engine, "coal")
        expect_identical(rec$tools$coalr, as.character(utils::packageVersion("coalr")))
        expect_identical(rec$arguments$dnds_est.method, "YN")
        expect_identical(sum(unlist(rec$status_counts)), nrow(res))
        expect_identical(rec$input_files$md5, unname(tools::md5sum(c(thal_cds(), lyra_cds()))))
        # COAL's own table is still there
        expect_true(is.data.frame(attr(res, "coal")))
})

test_that("both engines give the same status for the same gene pairs (Comeron, 20 genes)", {
        skip_if_not(coal_status_possible() && perl_available, "coalr, DIAMOND or perl is not installed")
        skip_if_not(orthologr:::coal_has_native_methods(), "the installed COAL does not compute Comeron's method")
        native <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        coal <- suppressMessages(dNdS(thal_cds(), lyra_cds(), engine = "coal", print_citation = FALSE))
        expect_identical(names(coal), names(native))
        m <- merge(as.data.frame(native)[ , c("query_id", "subject_id", "status")],
                   as.data.frame(coal)[ , c("query_id", "subject_id", "status")], by = c("query_id", "subject_id"))
        expect_identical(nrow(m), nrow(native))
        expect_identical(m$status.x, m$status.y)
})

test_that("a saturated pair has status 'saturated' with the COAL engine (1,000 genes, Comeron)", {
        skip_if_not(coal_status_possible(), "coalr or DIAMOND is not installed")
        skip_if_not(orthologr:::coal_has_native_methods(), "the installed COAL does not compute Comeron's method")
        q <- system.file("seqs/ortho_thal_cds_1000.fasta", package = "orthologr")
        s <- system.file("seqs/ortho_lyra_cds_1000.fasta", package = "orthologr")
        run <- coal_messages(dNdS(q, s, engine = "coal", comp_cores = 2, print_citation = FALSE))
        res <- run$value
        expect_identical(nrow(res), 995L)
        expect_identical(as.vector(table(factor(res$status, levels = c("ok", "saturated")))), c(994L, 1L))
        expect_true(all(is.na(res$dNdS[res$status == "saturated"])))
        expect_match(run$messages[length(run$messages)], "995 gene pairs - 994 with dN/dS, 1 saturated", fixed = TRUE)
})

test_that("options(orthologr.status_column = FALSE) leaves the column out with the COAL engine", {
        skip_if_not(coal_status_possible(), "coalr or DIAMOND is not installed")
        old <- options(orthologr.status_column = FALSE)
        on.exit(options(old), add = TRUE)
        res <- suppressMessages(dNdS_coal(thal_cds(), lyra_cds(), dnds_est.method = "YN", print_citation = FALSE))
        expect_false(is.element("status", names(res)))
        expect_identical(sum(unlist(attr(res, "orthologr_run")$status_counts)), nrow(res))
})

test_that("dNdS_clusters() uses the same words in 'status' and carries the record of the run", {
        skip_if_not(coal_status_possible(), "coalr or DIAMOND is not installed")
        run <- coal_messages(dNdS_clusters(c(thal_cds(), lyra_cds()), aa_aln_tool = "pairwise"))
        res <- run$value
        expect_identical(names(res), c("pairs", "summary", "failed"))
        expect_true(all(is.element(res$pairs$status, orthologr:::dnds_status_levels())))
        expect_identical(res$pairs$status == "ok", is.finite(res$pairs$dNdS))
        expect_true(any(grepl(paste0("dNdS_clusters(): ", nrow(res$pairs), " gene pairs in "), run$messages, fixed = TRUE)))
        expect_identical(run$stdout, character(0))
        rec <- attr(res, "orthologr_run")
        expect_identical(rec[["function"]], "dNdS_clusters")
        expect_identical(rec$engine, "coal")
        expect_identical(sum(unlist(rec$status_counts)), nrow(res$pairs))
        expect_identical(nrow(rec$input_files), 2L)
})

test_that("dNdS_orthogroups(engine = \"coal\") has the status column and the record of the run", {
        skip_if_not(coal_status_possible(), "coalr or DIAMOND is not installed")
        files <- c(Athaliana = thal_cds(), Alyrata = lyra_cds())
        res <- suppressMessages(suppressWarnings(dNdS_orthogroups(files, engine = "coal", aa_aln_tool = "pairwise")))
        expect_identical(names(res$pairs)[ncol(res$pairs)], "status")
        expect_true(all(is.element(res$pairs$status, orthologr:::dnds_status_levels())))
        expect_identical(res$pairs$status == "ok", is.finite(res$pairs$dNdS))
        expect_identical(attr(res, "orthologr_run")$engine, "coal")
        expect_true(nzchar(attr(res, "orthologr_run")$tools$coalr))
})

test_that("dNdS_clusters(): a cluster whose pairs are only saturated is 'saturated', not 'failed' (O13)", {
        sm <- data.frame(cluster_id = c("a", "b", "c", "d", "e"), status = c("failed", "failed", "ok", "failed", "skipped"),
                         n_pairs = c(1L, 3L, 2L, 2L, 0L), n_ok = c(0L, 0L, 1L, 0L, 0L), n_not_converged = 0L,
                         n_saturated = c(1L, 2L, 1L, 1L, 0L), n_no_data = c(0L, 1L, 0L, 0L, 0L), n_failed = c(0L, 0L, 0L, 1L, 0L),
                         note = c("", NA, "", "", "too large"), stringsAsFactors = FALSE)
        out <- dnds_cluster_status(sm)
        expect_identical(out$status, c("saturated", "saturated", "ok", "failed", "skipped"))
        expect_match(out$note[1], "no pair has an estimate \\(1 pair: 1 saturated\\)")
        expect_match(out$note[2], "2 saturated, 1 without codons to compare")
        expect_match(out$note[4], "1 failed")
        expect_identical(out$note[c(3, 5)], c("", "too large"))
        sm$n_saturated[2] <- 0L
        expect_identical(dnds_cluster_status(sm)$status[2], "no_data")
})
