# orthologr 0.6.0 performance: PAL2NAL and KaKs_Calculator with one process for a chunk of gene pairs
# (pal2nal_batch(), kaks_calc_batch()) give exactly what one process per pair gives
# (options(orthologr.batch_tools = FALSE), the code of 0.6.0 before the batches).

# The gene pairs of the 20-gene example in the form dNdS() hands to compute_dnds() (see test-perf.R).
tools_pair_table <- function(n = 20) {
        q <- suppressMessages(read.cds(thal_cds(), format = "fasta", delete_corrupt_cds = FALSE))
        s <- suppressMessages(read.cds(lyra_cds(), format = "fasta", delete_corrupt_cds = FALSE))
        n <- min(n, nrow(q), nrow(s))
        data.table::data.table(query_id = q$geneids[seq_len(n)], subject_id = s$geneids[seq_len(n)],
                               query_cds = q$seqs[seq_len(n)], subject_cds = s$seqs[seq_len(n)],
                               query_aa = toupper(translate_cds(q$seqs[seq_len(n)])),
                               subject_aa = toupper(translate_cds(s$seqs[seq_len(n)])))
}

# compute_dnds() of the table with options(orthologr.batch_tools = batch); the name of the temporary file a failed
# PAL2NAL call was to write (a new one in every call) is replaced
tools_dnds <- function(tbl, batch, ...) {
        old <- options(orthologr.batch_tools = batch)
        on.exit(options(old), add = TRUE)
        r <- suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                           quiet = TRUE, ...))
        f <- attr(r, "failed_pairs")
        f$reason <- gsub("aln_[0-9a-f]+\\.tmp", "aln_<tmp>.tmp", f$reason)
        f$reason <- gsub("[^ ']*_calculation/[^ ']*", "<calc file>", f$reason)
        attr(r, "failed_pairs") <- f
        r
}

# How often a function of the package is called while `expr` is evaluated.
count_calls <- function(fun, expr) {
        n <- 0L
        ns <- asNamespace("orthologr")
        f <- get(fun, envir = ns)
        counter <- function(...) {
                n <<- n + 1L
                f(...)
        }
        testthat::local_mocked_bindings(!!fun := counter, .package = "orthologr")
        force(expr)
        n
}

test_that("pal2nal_batch.pl gives the output, exit status and nothing else of one 'perl pal2nal.pl' per job", {
        skip_on_cran()
        skip_if_not(perl_available, "perl is not installed")
        skip_if_not(pal2nal_batch_possible(), "no fork()")
        d <- new_test_dir()
        on.exit(unlink(d, recursive = TRUE), add = TRUE)
        p2n <- system.file("pal2nal", "pal2nal.v14", "pal2nal.pl", package = "orthologr")
        aln <- system.file("pal2nal", "pal2nal.v14", "test.aln", package = "orthologr")
        nuc <- system.file("pal2nal", "pal2nal.v14", "test.nuc", package = "orthologr")
        jobs <- list(c(aln, nuc, "-output", "fasta"), c(aln, file.path(d, "missing.fasta"), "-output", "fasta"),
                     c(aln, nuc, "-output", "clustal"), c(aln, nuc, "-output", "nonsense"), c(aln, nuc, "-output", "paml"))
        out <- file.path(d, paste0("out", seq_along(jobs)))
        writeLines(vapply(seq_along(jobs), function(k) paste(c(out[k], jobs[[k]]), collapse = "\t"), ""), file.path(d, "jobs"))
        res <- system2("perl", shQuote(c(system.file("pal2nal", "pal2nal_batch.pl", package = "orthologr"), p2n,
                                         file.path(d, "jobs"), file.path(d, "status"))))
        expect_identical(res, 0L)
        st <- as.integer(readLines(file.path(d, "status")))
        for (k in seq_along(jobs)) {
                ref <- file.path(d, paste0("ref", k))
                s <- suppressWarnings(system2("perl", shQuote(c(p2n, jobs[[k]])), stdout = ref, stderr = FALSE))
                expect_identical(st[k] %/% 256L, as.integer(s), info = k)
                expect_identical(unname(tools::md5sum(out[k])), unname(tools::md5sum(ref)), info = k)
        }
        expect_identical(st[c(1, 3, 5)], c(0L, 0L, 0L))
        expect_true(st[2] != 0)
})

test_that("pal2nal_batch() writes the codon alignments of codon_aln() and leaves the pairs it cannot do", {
        skip_on_cran()
        skip_if_not(perl_available, "perl is not installed")
        skip_if_not(pal2nal_batch_possible(), "no fork()")
        d <- new_test_dir()
        old <- options(orthologr.tmpdir = d)
        on.exit({ options(old); unlink(d, recursive = TRUE) }, add = TRUE)
        aln <- system.file("pal2nal", "pal2nal.v14", "test.aln", package = "orthologr")
        nuc <- system.file("pal2nal", "pal2nal.v14", "test.nuc", package = "orthologr")
        bad <- file.path(d, "bad.fasta")
        writeLines(c(">a", "ATGATG", ">b", "ATGATG"), bad)
        done <- pal2nal_batch(c(aln, aln, aln), c(nuc, bad, file.path(d, "missing")), c("p1", "p2", "p3"))
        expect_identical(done, c(TRUE, FALSE, FALSE))
        out <- file.path(d, "_alignment", "codon_aln", "p1_pal2nal.aln")
        got <- readLines(out)
        codon_aln(aln, nuc, format = "fasta", tool = "pal2nal", codon_aln_name = "ref", quiet = TRUE)
        expect_identical(got, readLines(file.path(d, "_alignment", "codon_aln", "ref_pal2nal.aln")))
        expect_false(file.exists(file.path(d, "_alignment", "codon_aln", "p2_pal2nal.aln")))
        # nothing is left behind but the alignments
        expect_setequal(list.files(file.path(d, "_alignment", "codon_aln")), c("p1_pal2nal.aln", "ref_pal2nal.aln"))
})

test_that("compute_dnds(): PAL2NAL and KaKs_Calculator per chunk give the table of one process per pair", {
        skip_on_cran()
        skip_if_not(perl_available, "perl is not installed")
        tbl <- tools_pair_table()
        # an unrelated pair, a pair whose protein does not fit its CDS (PAL2NAL fails) and a short pair
        tbl$subject_aa[5] <- tbl$subject_aa[6]
        tbl$subject_cds[5] <- tbl$subject_cds[6]
        tbl$query_aa[7] <- paste0(substr(tbl$query_aa[7], 1, 10), "WWWWWWWWWWWW", substr(tbl$query_aa[7], 11, 400))
        methods <- c("Comeron", if (kaks_available) c("NG", "LWL", "LPB", "MLWL", "MLPB", "YN", "MYN", "GY", "GMYN"))
        for (m in methods) {
                new <- tools_dnds(tbl, TRUE, codon_aln_tool = "pal2nal", dnds_est.method = m)
                old <- tools_dnds(tbl, FALSE, codon_aln_tool = "pal2nal", dnds_est.method = m)
                expect_identical(new, old, info = m)
                expect_identical(nrow(attr(new, "failed_pairs")), 1L, info = m)
                expect_match(attr(new, "failed_pairs")$reason, "PAL2NAL", info = m)
        }
        if (kaks_available && requireNamespace("coalr", quietly = TRUE)) {
                new <- tools_dnds(tbl, TRUE, codon_aln_tool = "coal", dnds_est.method = "YN")
                expect_identical(new, tools_dnds(tbl, FALSE, codon_aln_tool = "coal", dnds_est.method = "YN"))
        }
})

test_that("compute_dnds(): one process per chunk, and no per-pair call where the batch computed the pair", {
        skip_on_cran()
        skip_if_not(perl_available && kaks_available, "perl or KaKs_Calculator is not installed")
        skip_if_not(pal2nal_batch_possible(), "no fork()")
        tbl <- tools_pair_table(12)
        old <- options(orthologr.batch_tools = TRUE)
        on.exit(options(old), add = TRUE)
        run <- function() suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise",
                                                         aa_aln_tool = "NW", codon_aln_tool = "pal2nal",
                                                         dnds_est.method = "YN", quiet = TRUE))
        expect_identical(count_calls("substitutionrate", run()), 0L)
        expect_identical(count_calls("codon_aln", run()), 0L)
        expect_identical(count_calls("kaks_calc_batch", run()), 1L)
        expect_identical(count_calls("pal2nal_batch", run()), 1L)
        local({
                testthat::local_mocked_bindings(dnds_chunk_size = function() 5L, .package = "orthologr")
                expect_identical(count_calls("kaks_calc_batch", run()), 3L)
        })
        # the settings of 0.4.3: one pair at a time, as there
        suppressMessages(legacy <- orthologr_legacy_options())
        on.exit(options(legacy), add = TRUE)
        expect_identical(count_calls("kaks_calc_batch", run()), 0L)
        expect_identical(count_calls("substitutionrate", run()), 12L)
        expect_identical(count_calls("pal2nal_batch", run()), 0L)
})

test_that("kaks_calc_batch(): a pair the batch cannot match is computed on its own, with the same values and errors", {
        skip_on_cran()
        skip_if_not(kaks_available && requireNamespace("coalr", quietly = TRUE), "KaKs_Calculator or coalr is not installed")
        skip_if(.Platform$OS.type == "windows", "shell scripts")
        tbl <- tools_pair_table(12)
        ref <- tools_dnds(tbl, FALSE, codon_aln_tool = "coal", dnds_est.method = "NG")
        real <- Sys.which("KaKs_Calculator")
        d <- new_test_dir()
        on.exit(unlink(d, recursive = TRUE), add = TRUE)
        fake <- function(body) {
                dir <- tempfile("kc_", tmpdir = d)
                dir.create(dir)
                writeLines(c("#!/bin/sh", body), file.path(dir, "KaKs_Calculator"))
                Sys.chmod(file.path(dir, "KaKs_Calculator"), "755")
                dir
        }
        # a record of more than one pair loses its rows (KaKs_Calculator leaves out a record it rejects); a batch that
        # fails; a batch that writes rows in another order
        many <- 'n=$(grep -c "^q" "$2"); '
        wrappers <- list(
                drop = fake(paste0(many, shQuote(real), ' "$@" > /dev/null; if [ "$n" -gt 1 ]; then grep -v "^q3\t" "$4" > "$4.x"; mv "$4.x" "$4"; fi')),
                fail = fake(paste0(many, 'if [ "$n" -gt 1 ]; then exit 3; fi; exec ', shQuote(real), ' "$@"')),
                order = fake(paste0(many, shQuote(real), ' "$@" > /dev/null; if [ "$n" -gt 1 ]; then (head -1 "$4"; tail -n +2 "$4" | sort -r) > "$4.x"; mv "$4.x" "$4"; fi')))
        for (w in names(wrappers)) {
                n <- 0L
                got <- local({
                        f <- substitutionrate
                        testthat::local_mocked_bindings(substitutionrate = function(...) { n <<- n + 1L; f(...) },
                                                        .package = "orthologr")
                        tools_dnds(tbl, TRUE, codon_aln_tool = "coal", dnds_est.method = "NG", kaks_calc_path = wrappers[[w]])
                })
                expect_identical(got, ref, info = w)
                # every pair of the chunk was computed on its own
                expect_identical(n, 12L, info = w)
        }
        # KaKs_Calculator is not there: every pair fails with the error of the per-pair call
        missing <- file.path(d, "no_such_folder")
        expect_identical(tools_dnds(tbl, TRUE, codon_aln_tool = "coal", dnds_est.method = "NG", kaks_calc_path = missing),
                         tools_dnds(tbl, FALSE, codon_aln_tool = "coal", dnds_est.method = "NG", kaks_calc_path = missing))
})

test_that("compute_dnds(): the batched table does not depend on the chunk size or the number of processes", {
        skip_on_cran()
        skip_if_not(perl_available && kaks_available, "perl or KaKs_Calculator is not installed")
        tbl <- tools_pair_table()
        tbl$subject_cds[5] <- tbl$subject_cds[6]
        for (cfg in list(list(codon_aln_tool = "pal2nal", dnds_est.method = "YN"),
                         list(codon_aln_tool = "pal2nal", dnds_est.method = "MLPB"))) {
                ref <- do.call(tools_dnds, c(list(tbl, FALSE), cfg))
                for (size in c(1L, 3L, 7L)) {
                        local({
                                testthat::local_mocked_bindings(dnds_chunk_size = function() size, .package = "orthologr")
                                expect_identical(do.call(tools_dnds, c(list(tbl, TRUE), cfg)), ref, info = size)
                        })
                }
                if (parallel::detectCores() >= 2 && .Platform$OS.type != "windows")
                        expect_identical(do.call(tools_dnds, c(list(tbl, TRUE), cfg, list(comp_cores = 2))), ref)
        }
})

test_that("options(orthologr.batch_tools) is not part of the record of a run", {
        old <- options(orthologr.batch_tools = FALSE)
        on.exit(options(old), add = TRUE)
        expect_false(is.element("orthologr.batch_tools", names(orthologr_options_in_force())))
})
