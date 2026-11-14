# The end of a dN/dS run: the column 'status', the summary message, progress messages, the record of the run
# (attr "orthologr_run"), write_dNdS() / read_dNdS(), and message hygiene (nothing on stdout).

collect_messages <- function(expr) {
        msgs <- character(0)
        out <- utils::capture.output(
                value <- withCallingHandlers(expr, message = function(m) {
                        msgs <<- c(msgs, conditionMessage(m))
                        invokeRestart("muffleMessage")
                })
        )
        list(value = value, messages = msgs, stdout = out)
}

native_possible <- function() diamond_available && perl_available

# ---------------------------------------------------------------- status

test_that("dnds_status() reads a pair by its values", {
        st <- orthologr:::dnds_status(
                dN     = c(0.01, 0.02, 0.03, NA,   0.05, NA, NA),
                dS     = c(0.10, 0,    NA,   0.30, 0.50, NA, NA),
                dNdS   = c(0.10, NA,   NA,   NA,   0.10, NA, NA),
                failed = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
        )
        expect_identical(st, c("ok", "dS_zero", "saturated", "saturated", "ok", "saturated", "failed"))
        expect_true(all(is.element(st, orthologr:::dnds_status_levels())))
})

test_that("dnds_status() uses what the engine says for pairs without data and failed pairs", {
        st <- orthologr:::dnds_status(dN = c(NA, NA, 0.1, NA), dS = c(NA, NA, 0.2, 0.3), dNdS = c(NA, NA, 0.5, NA),
                                      engine_status = c("no_data", "error", "warning", "saturated"))
        expect_identical(st, c("no_data", "failed", "ok", "saturated"))
})

test_that("the summary names the numbers by status, the median dS and the sequences left out", {
        status <- c(rep("ok", 90), rep("saturated", 5), rep("dS_zero", 3), rep("failed", 2))
        dropped <- data.frame(file = "a.fa", id = c("x", "y"), reason = c("not_acgt", "not_acgt"))
        txt <- orthologr:::dnds_summary_text(status, dS = c(rep(0.2, 90), rep(NA, 10)), dropped = dropped)
        expect_match(txt, "dNdS(): 100 gene pairs", fixed = TRUE)
        expect_match(txt, "90 with dN/dS", fixed = TRUE)
        expect_match(txt, "5 saturated", fixed = TRUE)
        expect_match(txt, "3 with dS = 0", fixed = TRUE)
        expect_match(txt, "2 failed", fixed = TRUE)
        expect_match(txt, "Median dS: 0.2", fixed = TRUE)
        expect_match(txt, "2 input sequences were left out before the search", fixed = TRUE)
        expect_false(grepl("too distant", txt, fixed = TRUE))
})

test_that("the summary says so when the species are too distant for dN/dS", {
        # many saturated pairs
        status <- c(rep("ok", 30), rep("saturated", 70))
        txt <- orthologr:::dnds_summary_text(status, dS = c(rep(1.9, 30), rep(NA, 70)))
        expect_match(txt, "70 of 100 pairs (70 %) are saturated", fixed = TRUE)
        expect_match(txt, "too distant", fixed = TRUE)
        # no saturated pair, but a high median dS
        txt <- orthologr:::dnds_summary_text(rep("ok", 50), dS = rep(2.4, 50))
        expect_match(txt, "too distant", fixed = TRUE)
        # a handful of pairs is not judged
        txt <- orthologr:::dnds_summary_text(rep("saturated", 5), dS = rep(NA_real_, 5))
        expect_false(grepl("too distant", txt, fixed = TRUE))
})

# ---------------------------------------------------------------- dNdS(): status column, messages, record

test_that("dNdS() ends with the column 'status', one summary message and the record of the run", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        orthologr:::tool_cache_reset()
        run <- collect_messages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        res <- run$value

        # status: last column, known words, in agreement with the values
        expect_identical(names(res)[ncol(res)], "status")
        expect_true(all(is.element(res$status, orthologr:::dnds_status_levels())))
        expect_identical(res$status == "ok", is.finite(res$dNdS))

        # nothing on stdout; a handful of messages, none of them empty
        expect_identical(run$stdout, character(0))
        expect_lte(length(run$messages), 3)
        expect_true(all(nzchar(trimws(run$messages))))
        expect_match(run$messages[1], "^dNdS\\(\\): reciprocal best hits of 'ortho_thal_cds.fasta' against 'ortho_lyra_cds.fasta'")
        expect_match(run$messages[length(run$messages)], paste0("dNdS(): ", nrow(res), " gene pairs - ", sum(res$status == "ok"), " with dN/dS"),
                     fixed = TRUE)

        # the record of the run
        rec <- attr(res, "orthologr_run")
        expect_true(all(is.element(c("orthologr", "function", "arguments", "options", "engine", "tools", "date",
                                     "input_files", "status_counts"), names(rec))))
        expect_identical(rec$orthologr, as.character(utils::packageVersion("orthologr")))
        expect_identical(rec[["function"]], "dNdS")
        expect_identical(rec$engine, "orthologr")
        expect_identical(rec$arguments$dnds_est.method, "Comeron")
        expect_identical(rec$arguments$best_hit_rule, "bitscore")
        expect_identical(rec$arguments$query_file, thal_cds())
        expect_identical(rec$options$orthologr.best_hit_rule, "bitscore")
        expect_identical(rec$options$orthologr.legacy_pairwise_aln, FALSE)
        expect_identical(rec$options$orthologr.status_column, TRUE)
        expect_identical(rec$options$orthologr.ambiguous_cds, "drop")
        expect_match(rec$tools$diamond, "diamond")
        expect_true(nzchar(rec$tools$pwalign))
        expect_identical(rec$input_files$md5, unname(tools::md5sum(c(thal_cds(), lyra_cds()))))
        expect_identical(rec$input_files$size, as.numeric(file.size(c(thal_cds(), lyra_cds()))))
        expect_identical(sum(unlist(rec$status_counts)), nrow(res))
        expect_identical(rec$status_counts$ok, sum(res$status == "ok"))
})

test_that("sequences trimmed to a multiple of 3 are counted in the summary, not in a block of messages (quiet = TRUE)", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        cds <- Biostrings::readDNAStringSet(thal_cds())
        cds[1:2] <- Biostrings::DNAStringSet(paste0(as.character(cds[1:2]), c("A", "AC")))
        f <- file.path(new_test_dir("trimmed_"), "thal_trim.fasta")
        Biostrings::writeXStringSet(cds, f)
        run <- collect_messages(dNdS(f, lyra_cds(), print_citation = FALSE))
        res <- run$value
        # one line at the start and the summary (was: 10 messages, the old block about corrupted CDS among them)
        expect_length(run$messages, 2)
        expect_match(run$messages[2], "2 coding sequences have a length that is not a multiple of 3 and were trimmed", fixed = TRUE)
        tr <- attr(res, "trimmed_sequences")
        expect_identical(sort(tr$id), sort(fasta_header_id(names(cds)[1:2])))
        expect_identical(tr$removed_nt[match(fasta_header_id(names(cds)[1:2]), tr$id)], c(1L, 2L))
        expect_true(all(file.exists(tr$corrupted_file)))
        # quiet = FALSE: the long explanation as before, the same numbers
        loud <- collect_messages(dNdS(f, lyra_cds(), print_citation = FALSE, quiet = FALSE))
        expect_true(any(grepl("cannot be properly divided in base triplets", loud$messages)))
        expect_equal(as.data.frame(loud$value)[, c("dN", "dS", "dNdS")], as.data.frame(res)[, c("dN", "dS", "dNdS")])
})

test_that("options(orthologr.status_column = FALSE) gives the table without the column and the same values", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        with_col <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        old <- options(orthologr.status_column = FALSE)
        on.exit(options(old), add = TRUE)
        without <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_false(is.element("status", names(without)))
        expect_identical(names(without), setdiff(names(with_col), "status"))
        expect_identical(as.data.frame(without)[ , names(without)], as.data.frame(with_col)[ , names(without)],
                         ignore_attr = TRUE)
        # the counts are still in the record
        expect_identical(attr(without, "orthologr_run")$status_counts, attr(with_col, "orthologr_run")$status_counts)
        expect_identical(attr(without, "orthologr_run")$options$orthologr.status_column, FALSE)
})

test_that("the citation note is a message and is shown once per session", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        state <- get(".orthologr_state", envir = asNamespace("orthologr"))
        old <- state$citation_shown
        on.exit(assign("citation_shown", old, envir = state), add = TRUE)
        assign("citation_shown", NULL, envir = state)
        first <- collect_messages(dNdS(thal_cds(), lyra_cds()))
        second <- collect_messages(dNdS(thal_cds(), lyra_cds()))
        expect_identical(first$stdout, character(0))
        expect_identical(sum(grepl("Please cite: Drost et al.", first$messages, fixed = TRUE)), 1L)
        expect_true(all(endsWith(first$messages, "\n")))
        expect_identical(sum(grepl("Please cite", second$messages, fixed = TRUE)), 0L)
})

test_that("quiet = FALSE reports progress, serial and parallel, without changing the result", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        quiet_res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        # a progress message every 5 gene pairs instead of every tenth of the run (20 pairs here)
        testthat::local_mocked_bindings(progress_step = function(n) 5L, .package = "orthologr")
        run <- collect_messages(dNdS(thal_cds(), lyra_cds(), quiet = FALSE, print_citation = FALSE))
        expect_identical(run$stdout, character(0))
        progress <- grep("^dN/dS: [0-9]+ of [0-9]+ gene pairs \\(", run$messages, value = TRUE)
        expect_identical(sub(" \\(.*", "", sub("\n$", "", progress)),
                         paste0("dN/dS: ", c(5, 10, 15), " of ", nrow(quiet_res), " gene pairs"))
        expect_true(any(grepl("Running diamond version", run$messages, fixed = TRUE)))
        expect_identical(as.data.frame(run$value), as.data.frame(quiet_res), ignore_attr = TRUE)

        skip_if(parallel::detectCores() < 2, "fewer than 2 cores")
        skip_if(.Platform$OS.type == "windows", "forked processes are needed to see the changed step in the workers")
        run2 <- collect_messages(dNdS(thal_cds(), lyra_cds(), quiet = FALSE, comp_cores = 2, print_citation = FALSE))
        progress2 <- grep("^dN/dS: [0-9]+ of [0-9]+ gene pairs \\(", run2$messages, value = TRUE)
        expect_identical(length(progress2), 3L)
        expect_identical(as.data.frame(run2$value), as.data.frame(quiet_res), ignore_attr = TRUE)
})

test_that("a failed gene pair has status 'failed' and is counted in the summary", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        # the estimation step fails for one pair
        real <- get("substitutionrate", envir = asNamespace("orthologr"))
        calls <- 0
        testthat::local_mocked_bindings(substitutionrate = function(...) {
                calls <<- calls + 1
                if (calls == 3)
                        stop("made to fail")
                real(...)
        }, .package = "orthologr")
        run <- collect_messages(suppressWarnings(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE)))
        res <- run$value
        expect_identical(sum(res$status == "failed"), 1L)
        fp <- attr(res, "failed_pairs")
        expect_identical(res$query_id[res$status == "failed"], fp$query_id)
        expect_match(run$messages[length(run$messages)], "1 failed (reasons in attr(<result>, \"failed_pairs\"))", fixed = TRUE)
        expect_identical(attr(res, "orthologr_run")$status_counts$failed, 1L)
})

# ---------------------------------------------------------------- write_dNdS() / read_dNdS()

test_that("write_dNdS() and read_dNdS() keep the table, its column types and its attributes", {
        skip_if_not_installed("jsonlite")
        tbl <- tibble::tibble(
                query_id = c("007", "1e5", "TRUE", "a|b-c"), subject_id = c("x 1", "NA.1", "0010", "y"),
                dN = c(0.0123456789012345, 0, NA, 1e-300), dS = c(0.1, 0, 0.3, NA), dNdS = c(0.123456789012345, NA, NA, NA),
                alig_length = c(100L, 200L, NA, 4L), status = c("ok", "dS_zero", "saturated", "saturated")
        )
        attr(tbl, "failed_pairs") <- data.frame(query_id = "q", subject_id = "s", reason = "PAL2NAL failed:\nline 2", stringsAsFactors = FALSE)
        attr(tbl, "dropped_sequences") <- data.frame(file = "a.fa", id = c("g1", "g2"), reason = "not_acgt", stringsAsFactors = FALSE)
        attr(tbl, "orthologr_run") <- orthologr:::run_record(
                fun = "dNdS", args = list(query_file = "a.fa", eval = "1E-5", comp_cores = 1, aligner_path = NULL),
                engine = "orthologr", tools = list(diamond = "diamond version 2.1.0"), files = thal_cds(), status = tbl$status)
        f <- file.path(new_test_dir(), "dnds table.tsv")
        files <- write_dNdS(tbl, f)
        expect_true(all(file.exists(files)))
        expect_identical(unname(files[["run"]]), paste0(f, ".run.json"))

        back <- read_dNdS(f)
        expect_s3_class(back, "tbl_df")
        expect_equal(as.data.frame(back), as.data.frame(tbl), ignore_attr = TRUE, tolerance = 1e-6)
        expect_identical(vapply(back, function(v) class(v)[1], ""), vapply(tbl, function(v) class(v)[1], ""))
        expect_identical(attr(back, "failed_pairs"), attr(tbl, "failed_pairs"))
        expect_identical(attr(back, "dropped_sequences"), attr(tbl, "dropped_sequences"))
        rec <- attr(tbl, "orthologr_run")
        got <- attr(back, "orthologr_run")
        for (n in c("orthologr", "function", "engine", "date"))
                expect_identical(got[[n]], rec[[n]], info = n)
        expect_identical(got$options$orthologr.best_hit_rule, rec$options$orthologr.best_hit_rule)
        expect_identical(got$arguments$eval, "1E-5")
        expect_identical(got$tools$diamond, "diamond version 2.1.0")
        expect_identical(got$input_files$md5, rec$input_files$md5)
        expect_identical(as.integer(unlist(got$status_counts)), as.integer(unlist(rec$status_counts)))

        # IDs with quotes and other characters that are special in a table (O10: 'd"q' came back as 'd""q')
        odd <- tibble::tibble(query_id = c('d"q', '"x"', "q'uote", "x,y", "h#sh", "a|b", "\u00fcn\u00efcode", '"'),
                              subject_id = c("s", "t", "u", "v", "w", "x", "y", "z"), dN = 0.1, dS = 0.2, dNdS = 0.5,
                              status = "ok")
        g <- file.path(new_test_dir(), "odd.tsv")
        write_dNdS(odd, g)
        expect_identical(read_dNdS(g)$query_id, odd$query_id)

        # the table file alone can be read by anything; read_dNdS() says that the record is missing
        plain <- utils::read.delim(f, stringsAsFactors = FALSE, colClasses = "character")
        expect_identical(plain$query_id, tbl$query_id)
        unlink(files[["run"]])
        expect_warning(alone <- read_dNdS(f), "run.json' was not found")
        expect_identical(nrow(alone), nrow(tbl))
        expect_null(attr(alone, "orthologr_run"))
})

test_that("read_dNdS() returns non-ASCII IDs equal to those written, in a UTF-8 and in a C locale (V03)", {
        skip_if_not_installed("jsonlite")
        skip_on_os("windows")
        # IDs as dNdS() returns them (the bytes of the FASTA file, unmarked), marked as UTF-8, and in latin1
        raw_id <- rawToChar(as.raw(c(0xc3, 0xbc, 0x6e, 0xc3, 0xaf, 0x63, 0x6f, 0x64, 0x65)))
        utf8_id <- raw_id
        Encoding(utf8_id) <- "UTF-8"
        latin_id <- iconv(rawToChar(as.raw(c(0xc3, 0xa9))), from = "UTF-8", to = "latin1")
        tbl <- data.frame(query_id = c(raw_id, utf8_id, "plain"), subject_id = c("s1", latin_id, raw_id),
                          dN = 0.1, dS = 0.2, dNdS = 0.5, status = "ok", stringsAsFactors = FALSE)
        attr(tbl, "failed_pairs") <- data.frame(query_id = c(raw_id, "x"), subject_id = c(latin_id, utf8_id),
                                                reason = c("r", raw_id), stringsAsFactors = FALSE)
        attr(tbl, "dropped_sequences") <- data.frame(file = "f.fa", id = raw_id, reason = "r", stringsAsFactors = FALSE)
        bytes <- function(v) lapply(v, charToRaw)
        old <- Sys.getlocale("LC_CTYPE")
        on.exit(Sys.setlocale("LC_CTYPE", old), add = TRUE)
        for (loc in c("C.UTF-8", "C")) {
                # (the comparisons are made in that locale too; equal bytes and marks are equal in every locale)
                if (!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc))))
                        next
                f <- file.path(new_test_dir(), paste0("enc_", loc, ".tsv"))
                write_dNdS(tbl, f)
                back <- read_dNdS(f)
                expect_true(all(back$query_id == tbl$query_id), info = loc)
                expect_false(anyNA(match(tbl$query_id, back$query_id)), info = loc)
                for (col in c("query_id", "subject_id")) {
                        expect_identical(back[[col]], tbl[[col]], info = paste(loc, col))
                        expect_identical(bytes(back[[col]]), bytes(tbl[[col]]), info = paste(loc, col))
                        expect_identical(Encoding(back[[col]]), Encoding(tbl[[col]]), info = paste(loc, col))
                }
                expect_identical(attr(back, "failed_pairs"), attr(tbl, "failed_pairs"), info = loc)
                expect_identical(bytes(attr(back, "failed_pairs")$subject_id), bytes(attr(tbl, "failed_pairs")$subject_id), info = loc)
                expect_identical(Encoding(attr(back, "failed_pairs")$query_id), Encoding(attr(tbl, "failed_pairs")$query_id), info = loc)
                expect_identical(attr(back, "dropped_sequences"), attr(tbl, "dropped_sequences"), info = loc)
                expect_identical(Encoding(attr(back, "dropped_sequences")$id), "unknown", info = loc)
        }
})

test_that("write_dNdS() checks its input", {
        skip_if_not_installed("jsonlite")
        expect_error(write_dNdS(1:3, tempfile()), "must be a dN/dS table")
        expect_error(write_dNdS(data.frame(a = 1), NA_character_), "'file' must be the path")
        expect_error(read_dNdS(file.path(tempdir(), "no_such_file.tsv")), "was not found")
})

test_that("a dNdS() result survives write_dNdS() / read_dNdS()", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        skip_if_not_installed("jsonlite")
        res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        f <- tempfile(fileext = ".tsv")
        on.exit(unlink(c(f, paste0(f, ".run.json"))), add = TRUE)
        write_dNdS(res, f)
        back <- read_dNdS(f)
        # (tolerance: R's reading of numbers of extreme exponent, e.g. 1e-300, is not exact on macOS arm64)
        expect_equal(as.data.frame(back), as.data.frame(res), ignore_attr = TRUE, tolerance = 1e-6)
        expect_identical(attr(back, "orthologr_run")$tools, attr(res, "orthologr_run")$tools)
        expect_identical(attr(back, "orthologr_run")$input_files, attr(res, "orthologr_run")$input_files)
        expect_identical(attr(back, "orthologr_run")$arguments$dnds_est.method, "Comeron")
})

# ---------------------------------------------------------------- other functions

test_that("divergence_stratigraphy() carries the record of its run", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        ds <- suppressMessages(divergence_stratigraphy(thal_cds(), lyra_cds(), eval = "1E-5", ortho_detection = "RBH",
                                                       comp_cores = 1, quiet = TRUE, clean_folders = TRUE))
        rec <- attr(ds, "orthologr_run")
        expect_identical(rec[["function"]], "divergence_stratigraphy")
        expect_identical(rec$arguments$ortho_detection, "RBH")
        expect_identical(rec$gene_pairs_in_map, nrow(ds))
        expect_identical(rec$input_files$md5, unname(tools::md5sum(c(thal_cds(), lyra_cds()))))
        expect_identical(names(attr(ds, "dNdS_tbl"))[ncol(attr(ds, "dNdS_tbl"))], "status")
})

test_that("dNdS_orthogroups() (engine = \"orthologr\") has the status column, the summary and the record", {
        skip_if_not(native_possible(), "DIAMOND or perl is not installed")
        files <- c(Athaliana = thal_cds(), Alyrata = lyra_cds())
        run <- collect_messages(dNdS_orthogroups(files, aa_aln_tool = "pairwise"))
        res <- run$value
        expect_identical(names(res$pairs)[ncol(res$pairs)], "status")
        expect_true(all(is.element(res$pairs$status, orthologr:::dnds_status_levels())))
        expect_identical(res$pairs$status == "ok", is.finite(res$pairs$dNdS))
        expect_true(any(grepl(paste0("dNdS_orthogroups(): ", nrow(res$pairs), " gene pairs in "), run$messages, fixed = TRUE)))
        expect_identical(run$stdout, character(0))
        rec <- attr(res, "orthologr_run")
        expect_identical(rec[["function"]], "dNdS_orthogroups")
        expect_identical(rec$engine, "orthologr")
        expect_identical(sum(unlist(rec$status_counts)), nrow(res$pairs))
        expect_identical(basename(rec$input_files$file), basename(unname(files)))
})

test_that("orthogroups() records the package version and checks for splice variants", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        og <- suppressMessages(orthogroups(c(A = thal_aa(), B = lyra_aa()), seq_type = "protein"))
        expect_identical(og$params$orthologr_version, as.character(utils::packageVersion("orthologr")))

        # a proteome with a second transcript for 8 of 20 genes (IDs with a transcript suffix)
        aa <- Biostrings::readAAStringSet(thal_aa())
        ids <- sub("\\..*$", "", sub(" .*", "", names(aa)))
        iso <- aa[1:8]
        names(aa) <- paste0(ids, ".1")
        names(iso) <- paste0(ids[1:8], ".2")
        f <- file.path(new_test_dir(), "with_isoforms.faa")
        Biostrings::writeXStringSet(c(aa, iso), f)
        expect_warning(suppressMessages(orthogroups(c(A = f, B = lyra_aa()), seq_type = "protein")),
                       "the file of species 'A', 'with_isoforms.faa' seems to hold several transcripts")
        old <- options(orthologr.isoform_warning = FALSE)
        on.exit(options(old), add = TRUE)
        expect_no_warning(suppressMessages(orthogroups(c(A = f, B = lyra_aa()), seq_type = "protein")))
})

test_that("a protein file given as coding sequences is an error that says so", {
        expect_true(orthologr:::fasta_looks_like_protein(thal_aa()))
        expect_false(orthologr:::fasta_looks_like_protein(thal_cds()))
        expect_error(read.cds(thal_aa(), format = "fasta"), "does not look like coding DNA")
        skip_if_not(diamond_available, "DIAMOND is not installed")
        expect_error(suppressMessages(dNdS(thal_aa(), lyra_cds())), "Is this a protein file")
        # coding sequences with a few N are not taken for proteins
        f <- file.path(new_test_dir(), "with_n.fa")
        writeLines(c(">a", "ATGNNNAAACCCGGGTTTNNNACGTRYACGTACGTACGTAA"), f)
        expect_false(orthologr:::fasta_looks_like_protein(f))
})

test_that("the isoform check needs 1 % of the genes named in the headers to have several sequences", {
        d <- new_test_dir()
        # 300 genes, 2 of them with two coding sequences (as in a bacterial genome): no warning
        genes <- sprintf("b%04d", 1:300)
        one <- c(genes, genes[1:2])
        f1 <- file.path(d, "few.fna")
        writeLines(as.vector(rbind(sprintf(">lcl|cds_%d [locus_tag=%s] [protein_id=P%d]", seq_along(one), one, seq_along(one)), "ATGAAATAA")), f1)
        expect_no_warning(orthologr:::isoform_warning(f1, "query_file"))
        # 30 of 300 genes with two sequences: a warning
        many <- c(genes, genes[1:30])
        f2 <- file.path(d, "many.fna")
        writeLines(as.vector(rbind(sprintf(">lcl|cds_%d [locus_tag=%s] [protein_id=P%d]", seq_along(many), many, seq_along(many)), "ATGAAATAA")), f2)
        expect_warning(orthologr:::isoform_warning(f2, "query_file"), "30 genes with more than one sequence")
})

test_that("dNdS(dnds_est.method = 'kaks_calc') is an error before the search; 'ALL' is summarised as pairs x methods (O07)", {
        expect_error(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "kaks_calc"), class = "orthologr_error_bad_argument")
        expect_error(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "kaks_calc"), "only substitutionrate")
        res <- data.frame(query_id = rep(c("a", "b"), each = 3), subject_id = rep(c("x", "y"), each = 3),
                          method = rep(c("NG", "YN", "MA"), 2), dN = 0.1, dS = 0.5, dNdS = 0.2)
        msgs <- character()
        withCallingHandlers(dnds_finish(res, failed_pairs = NULL, dropped = NULL, print_citation = FALSE,
                                        record = list(fun = "dNdS", args = list(), engine = "orthologr", tools = list(), files = character(0))),
                            message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
        expect_true(any(grepl("dNdS(): 6 rows (2 gene pairs x 3 methods) - 6 with dN/dS", msgs, fixed = TRUE)))
})
