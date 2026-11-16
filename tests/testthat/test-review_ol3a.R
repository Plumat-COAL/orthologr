# REVIEW-OL-3 fixes: nucleotide files given as proteins, time limits of DIAMOND's steps, ungrouped best hits, the
# record of a run (dnds_equal(), write_dNdS() / read_dNdS() with every table), orthologr_sitrep(), and smaller ones.

ol3_fake_program <- function(name, lines, dir = new_test_dir("fake_tool_")) {
        writeLines(c("#!/bin/sh", lines), file.path(dir, name))
        Sys.chmod(file.path(dir, name), "755")
        dir
}

ol3_fasta <- function(seqs, file = tempfile(fileext = ".fasta")) {
        writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), file)
        file
}

# ---- H1: nucleotide sequences given as proteins; time limits of DIAMOND's steps --------------------------------

test_that("read.proteome() rejects nucleotide sequences (90 % or more A, C, G, T, U, N) with a bad_input error", {
        e <- tryCatch(read.proteome(thal_cds(), format = "fasta"), error = function(e) e)
        expect_s3_class(e, "orthologr_error_bad_input")
        expect_match(conditionMessage(e), "does not look like protein sequences", fixed = TRUE)
        expect_match(conditionMessage(e), "seq_type = \"cds\"", fixed = TRUE)
        # proteins are read as before
        expect_gt(nrow(read.proteome(thal_aa(), format = "fasta")), 0)
        # a protein rich in A, G, C, T and N stays a protein; RNA and gzip are recognised
        prot <- ol3_fasta(c(p1 = paste(rep("MAGTCNAGTCLKEWQ", 20), collapse = "")))
        expect_identical(read.proteome(prot, format = "fasta")$geneids, "p1")
        rna <- ol3_fasta(c(r1 = paste(rep("AUGGCUAAAGAUGAACCC", 10), collapse = "")))
        expect_error(read.proteome(rna, format = "fasta"), class = "orthologr_error_bad_input")
        gz <- tempfile(fileext = ".fasta.gz")
        con <- gzfile(gz, "w"); writeLines(readLines(thal_cds()), con); close(con)
        expect_error(read.proteome(gz, format = "fasta"), class = "orthologr_error_bad_input")
        # very short files are not judged
        tiny <- ol3_fasta(c(t1 = "ACGT"))
        expect_identical(read.proteome(tiny, format = "fasta")$geneids, "t1")
        expect_true(is.na(fasta_nucleotide_share(tiny)))
        expect_gt(fasta_nucleotide_share(thal_cds()), 0.99)
        expect_lt(fasta_nucleotide_share(thal_aa()), 0.5)
})

test_that("a CDS file searched as proteins stops before DIAMOND is started (protein_table(), diamond_best(), deepclust())", {
        expect_error(protein_table(thal_cds(), seq_type = "protein"), class = "orthologr_error_bad_input")
        # the stand-in would hang: it must not be reached
        bin <- ol3_fake_program("diamond", c("case \"$1\" in version) echo 'diamond version 2.2.8';; *) exec sleep 30;; esac"))
        skip_on_os("windows")
        tool_cache_reset()
        on.exit(tool_cache_reset(), add = TRUE)
        t0 <- proc.time()[["elapsed"]]
        expect_error(suppressMessages(diamond_best(thal_cds(), lyra_cds(), seq_type = "protein", path = bin)),
                     class = "orthologr_error_bad_input")
        expect_error(suppressMessages(deepclust(thal_cds(), path = bin)), class = "orthologr_error_bad_input")
        expect_lt(proc.time()[["elapsed"]] - t0, 20)
})

test_that("diamond makedb has a limit of its own (1 h); deepclust and realign have a finite one (7 days)", {
        rlang::local_options(orthologr.tool_timeout = NULL)
        expect_identical(tool_timeout("diamond makedb", "diamond"), 3600)
        expect_identical(tool_timeout_key("diamond makedb", "diamond"), "diamond makedb")
        expect_identical(tool_timeout("diamond deepclust", "diamond"), 604800)
        expect_identical(tool_timeout("diamond realign", "diamond"), 604800)
        expect_identical(tool_timeout_key("diamond deepclust", "diamond"), "diamond deepclust")
        expect_identical(tool_timeout("diamond blastp", "diamond"), 604800)
        expect_identical(tool_timeout_key("diamond blastp", "diamond"), "DIAMOND")
        # a limit for DIAMOND applies to its steps unless they are named themselves
        rlang::local_options(orthologr.tool_timeout = c(DIAMOND = 20))
        expect_identical(tool_timeout("diamond makedb", "diamond"), 20)
        expect_identical(tool_timeout("diamond deepclust", "diamond"), 20)
        rlang::local_options(orthologr.tool_timeout = c(DIAMOND = 20, "diamond makedb" = 5))
        expect_identical(tool_timeout("diamond makedb", "diamond"), 5)
        expect_identical(tool_timeout("diamond blastp", "diamond"), 20)
        rlang::local_options(orthologr.tool_timeout = c("Diamond Makedb" = 7))
        expect_identical(tool_timeout("diamond makedb", "diamond"), 7)
        expect_identical(tool_timeout("diamond blastp", "diamond"), 604800)
})

test_that("a DIAMOND that hangs in makedb or deepclust is stopped at its limit", {
        skip_on_cran()
        skip_on_os("windows")
        bin <- ol3_fake_program("diamond", "exec sleep 30")
        fa <- ol3_fasta(c(p1 = "MKLVAGTWQERTYIPASDFGHKLCVNM"))
        rlang::local_options(orthologr.tool_timeout = c("diamond makedb" = 1, "diamond deepclust" = 1))
        t0 <- proc.time()[["elapsed"]]
        e <- tryCatch(diamond_makedb(fa, tempfile("db_"), path = bin), error = function(e) e)
        expect_s3_class(e, "orthologr_error_tool_timeout")
        expect_match(conditionMessage(e), "diamond makedb did not finish within 1 seconds", fixed = TRUE)
        expect_match(conditionMessage(e), "c(\"diamond makedb\" = 4)", fixed = TRUE)
        e <- tryCatch(run_tool("diamond", args = list("deepclust"), path = bin, tool = "diamond deepclust"), error = function(e) e)
        expect_s3_class(e, "orthologr_error_tool_timeout")
        expect_lt(proc.time()[["elapsed"]] - t0, 20)
})

test_that("time limits: fractions are rounded up, too large numbers are no limit, unnamed elements are refused", {
        expect_identical(tool_timeout_value(0.4), 1)
        expect_identical(tool_timeout_value(2.5), 3)
        expect_identical(tool_timeout_value(.Machine$integer.max + 1), Inf)
        expect_identical(tool_timeout_value(1e12), Inf)
        expect_identical(tool_timeout_value(.Machine$integer.max), as.numeric(.Machine$integer.max))
        rlang::local_options(orthologr.tool_timeout = 0.2)
        expect_identical(tool_timeout("KaKs_Calculator"), 1)
        rlang::local_options(orthologr.tool_timeout = 1e10)
        expect_identical(tool_timeout("KaKs_Calculator"), Inf)
        rlang::local_options(orthologr.tool_timeout = c(KaKs_Calculator = 10, 5))
        expect_error(tool_timeout("KaKs_Calculator"), "needs the name of its program", class = "orthologr_error_bad_argument")
        x <- c(10, 5); names(x) <- c("KaKs_Calculator", NA)
        rlang::local_options(orthologr.tool_timeout = x)
        expect_error(tool_timeout("PAL2NAL"), class = "orthologr_error_bad_argument")
})

test_that("timeout_cap() shortens a limit for the time of an expression only", {
        rlang::local_options(orthologr.tool_timeout = NULL)
        expect_identical(timeout_cap("KaKs_Calculator", 300, tool_timeout("KaKs_Calculator")), 300)
        expect_identical(timeout_cap("KaKs_Calculator", 300, tool_timeout("PAL2NAL")), 3600)
        expect_identical(timeout_cap("KaKs_Calculator", 1e6, tool_timeout("KaKs_Calculator")), 3600)
        expect_identical(tool_timeout("KaKs_Calculator"), 3600)
        expect_identical(dnds_retry("KaKs_Calculator", FALSE, tool_timeout("KaKs_Calculator")), 3600)
        expect_identical(dnds_retry("KaKs_Calculator", TRUE, tool_timeout("KaKs_Calculator")), 300)
        rlang::local_options(orthologr.tool_timeout = c(KaKs_Calculator = 20))
        expect_identical(dnds_retry("KaKs_Calculator", TRUE, tool_timeout("KaKs_Calculator")), 20)
        # an error inside restores the limit
        try(timeout_cap("KaKs_Calculator", 2, stop("x")), silent = TRUE)
        expect_identical(tool_timeout("KaKs_Calculator"), 20)
})

test_that("after a chunk of KaKs_Calculator timed out, each pair is retried with the shorter limit", {
        skip_on_cran()
        skip_on_os("windows")
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        rlang::local_options(orthologr.tool_timeout = NULL)
        q <- c("ATGGCTAAAGATGAACCCGGG", "ATGAAACCCGGGTTTAAA")
        s <- c("ATGGCCAAAGATGAGCCCGGA", "ATGAAGCCCGGGTTCAAA")
        tbl <- data.table::data.table(query_id = c("a", "b"), subject_id = c("x", "y"), query_cds = q, subject_cds = s,
                                      query_aa = translate_cds(q), subject_aa = translate_cds(s))
        bin <- ol3_fake_program("KaKs_Calculator", "exec sleep 30")
        # the chunk "times out" (as after the 1 h limit); the retry of each pair gets dnds_retry_limit() = 1 s here
        local_mocked_bindings(
                kaks_calc_batch = function(files, ...) {
                        .orthologr_state$timeouts <- tool_timeouts() + 1L
                        vector("list", length(files))
                },
                dnds_retry_limit = function() 1)
        t0 <- proc.time()[["elapsed"]]
        res <- suppressWarnings(compute_dnds(data.table::copy(tbl), aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                             codon_aln_tool = "coal", dnds_est.method = "YN", kaks_calc_path = bin,
                                             quiet = TRUE))
        expect_lt(proc.time()[["elapsed"]] - t0, 20)
        fp <- attr(res, "failed_pairs")
        expect_identical(fp$query_id, c("a", "b"))
        expect_true(all(grepl("KaKs_Calculator did not finish within 1 seconds", fp$reason, fixed = TRUE)))
        # the limit is back
        expect_identical(tool_timeout("KaKs_Calculator"), 3600)
})

test_that("a run that times out is counted (tool_timeouts())", {
        skip_on_cran()
        skip_on_os("windows")
        bin <- ol3_fake_program("KaKs_Calculator", "exec sleep 30")
        rlang::local_options(orthologr.tool_timeout = 1)
        n <- tool_timeouts()
        expect_error(run_tool("KaKs_Calculator", path = bin), class = "orthologr_error_tool_timeout")
        expect_identical(tool_timeouts(), n + 1L)
})

# ---- M5: best hits are not grouped ------------------------------------------------------------------------------

test_that("best_hits(), diamond_best(), diamond_rec() and orthologs() return tibbles that are not grouped", {
        h <- data.frame(query_id = c("a", "a", "b"), subject_id = c("x", "y", "z"), evalue = c(1e-5, 1e-9, 1e-3),
                        alig_length = c(10L, 10L, 5L), perc_identity = c(90, 80, 70), bit_score = c(50, 60, 30),
                        stringsAsFactors = FALSE)
        for (rule in c("bitscore", "evalue")) {
                b <- best_hits(h, rule = rule)
                expect_false(inherits(b, "grouped_df"), info = rule)
                expect_s3_class(b, "tbl_df")
                expect_identical(b$subject_id, c("y", "z"), info = rule)
        }
        expect_false(inherits(best_hits(h[0, ]), "grouped_df"))
        h$evalue[1] <- NA
        expect_false(inherits(best_hits(h, rule = "evalue"), "grouped_df"))
        skip_if_not(diamond_available, "DIAMOND is not installed")
        for (od in c("DIAMOND_RBH", "DIAMOND_BH")) {
                o <- suppressMessages(orthologs(thal_aa(), lyra_aa(), ortho_detection = od))
                expect_false(inherits(o, "grouped_df"), info = od)
                expect_gt(nrow(o), 0)
        }
        expect_false(inherits(suppressMessages(diamond_rec(thal_cds(), lyra_cds())), "grouped_df"))
})

# ---- M6: two runs compare equal with dnds_equal(); version checks are no commands --------------------------------

test_that("dnds_equal(): two runs of the same analysis are equal; version checks are not in the commands", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        tool_cache_reset()
        r1 <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        r2 <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        c1 <- attr(r1, "orthologr_run")$commands
        expect_gt(length(c1), 0)
        expect_false(any(grepl("'version'$", c1)))
        expect_identical(length(c1), length(attr(r2, "orthologr_run")$commands))
        expect_true(dnds_equal(r1, r2))
        # with nothing ignored it is identical()
        expect_identical(dnds_equal(r1, r2, ignore = character(0)), identical(r1, r2))
        # a different value is a difference
        r3 <- r2
        r3$dN[1] <- r3$dN[1] + 1e-12
        expect_false(dnds_equal(r1, r3))
        r4 <- r2
        attr(r4, "orthologr_run")$options$orthologr.dS_zero <- "keep"
        expect_false(dnds_equal(r1, r4))
        # lists and tables in attributes
        expect_true(dnds_equal(list(pairs = r1, n = 1), list(pairs = r2, n = 1)))
        d1 <- data.frame(a = 1); attr(d1, "dNdS_tbl") <- r1
        d2 <- data.frame(a = 1); attr(d2, "dNdS_tbl") <- r2
        expect_true(dnds_equal(d1, d2))
        expect_error(dnds_equal(r1, r2, ignore = 1), class = "orthologr_error_bad_argument")
})

# ---- M7: write_dNdS() / read_dNdS() keep every table of the run ----------------------------------------------------

test_that("write_dNdS() / read_dNdS() round-trip every table attribute with its types and class", {
        skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is not installed")
        tbl <- tibble::tibble(query_id = c("a", "b"), subject_id = c("x", "y"), dN = c(0.1, 1 / 3), dS = c(0.2, NA),
                              dNdS = c(0.5, NA), status = c("ok", "saturated"))
        attr(tbl, "failed_pairs") <- data.frame(query_id = "b", subject_id = "y", reason = "r é", stringsAsFactors = FALSE)
        attr(tbl, "dropped_sequences") <- data.frame(file = "f.fa", id = "c", reason = "short", stringsAsFactors = FALSE)
        attr(tbl, "trimmed_sequences") <- data.frame(file = "f.fa", id = "a", removed_nt = 2L, corrupted_file = "/tmp/x.fasta",
                                                     stringsAsFactors = FALSE)
        attr(tbl, "stop_codon_sequences") <- data.frame(file = "f.fa", id = c("a", "b"), stop_codons = c(1L, 3L),
                                                        stringsAsFactors = FALSE)
        attr(tbl, "ambiguous_sequences") <- data.frame(file = "f.fa", id = "b", ambiguous_codons = 4L, flag = TRUE,
                                                       stringsAsFactors = FALSE)
        attr(tbl, "coal") <- tibble::tibble(seq1 = c("a", "b"), seq2 = c("x", "y"), dN = c(0.1000000000000000055, 2 / 3),
                                            dS = c(pi, NA), note = c("", "saturated"))
        attr(tbl, "orthologr_run") <- list(orthologr = "0.6.0", engine = "orthologr")
        file <- tempfile(fileext = ".tsv")
        write_dNdS(tbl, file)
        back <- read_dNdS(file)
        for (a in c("failed_pairs", "dropped_sequences", "trimmed_sequences", "stop_codon_sequences",
                    "ambiguous_sequences", "coal"))
                expect_identical(attr(back, a), attr(tbl, a), info = a)
        expect_identical(as.data.frame(back), as.data.frame(tbl))
        # empty tables stay empty tables with their columns
        attr(tbl, "failed_pairs") <- data.frame(query_id = character(0), subject_id = character(0), reason = character(0),
                                                stringsAsFactors = FALSE)
        write_dNdS(tbl, file)
        expect_identical(attr(read_dNdS(file), "failed_pairs"), attr(tbl, "failed_pairs"))
        # a sidecar of 0.6.0 (failed_pairs and dropped_sequences, no types) is still read
        side <- jsonlite::fromJSON(paste0(file, ".run.json"), simplifyVector = TRUE, simplifyDataFrame = FALSE)
        side$tables <- NULL
        side$trimmed_sequences <- NULL
        writeLines(jsonlite::toJSON(side, auto_unbox = TRUE, null = "null", na = "null"), paste0(file, ".run.json"))
        old <- read_dNdS(file)
        expect_null(attr(old, "trimmed_sequences"))
        expect_identical(attr(old, "dropped_sequences")$id, "c")
})

test_that("a real run with trimmed, stop-codon and ambiguous sequences reads back with all its tables", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is not installed")
        q <- Biostrings::readDNAStringSet(thal_cds())
        s <- as.character(q)
        s[1] <- paste0(s[1], "A")                                      # not a multiple of 3: trimmed
        substr(s[2], 31, 33) <- "TAA"                                  # an in-frame stop codon
        substr(s[3], 40, 40) <- "N"                                    # an ambiguous base
        qf <- tempfile(fileext = ".fasta")
        writeLines(as.vector(rbind(paste0(">", names(q)), s)), qf)
        r <- suppressWarnings(suppressMessages(dNdS(qf, lyra_cds(), print_citation = FALSE)))
        for (a in c("trimmed_sequences", "stop_codon_sequences", "ambiguous_sequences"))
                expect_false(is.null(attr(r, a)), info = a)
        file <- tempfile(fileext = ".tsv")
        write_dNdS(r, file)
        back <- read_dNdS(file)
        for (a in c("failed_pairs", "trimmed_sequences", "stop_codon_sequences", "ambiguous_sequences"))
                expect_identical(attr(back, a), attr(r, a), info = a)
        cols <- function(x) stats::setNames(lapply(names(x), function(n) x[[n]]), names(x))
        expect_identical(cols(back), cols(r))
})

# ---- M8: orthologr_sitrep() ---------------------------------------------------------------------------------------

test_that("orthologr_sitrep() reports what is found and hints at what is missing, without remembering a miss", {
        skip_on_os("windows")
        old <- Sys.getenv(c("PATH", "R_USER_CACHE_DIR", "COAL_CPP"), unset = NA)
        on.exit({
                for (v in names(old)) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
                tool_cache_reset()
        }, add = TRUE)
        bin <- with_coreutils(ol3_fake_program("diamond", "echo 'diamond version 2.2.8'"))
        Sys.setenv(PATH = bin, R_USER_CACHE_DIR = new_test_dir("cache_"))
        Sys.unsetenv("COAL_CPP")
        rlang::local_options(orthologr.coal_path = NULL)
        tool_cache_reset()
        rep <- orthologr_sitrep()
        expect_s3_class(rep, "tbl_df")
        expect_identical(names(rep), c("component", "needed_for", "found", "path", "version", "hint"))
        expect_type(rep$found, "logical")
        row <- function(x, name) x[x$component == name, ]
        d <- row(rep, "DIAMOND")
        expect_true(d$found)
        expect_identical(d$version, "diamond version 2.2.8")
        expect_identical(d$path, normalizePath(file.path(bin, "diamond")))
        expect_true(is.na(d$hint))
        expect_false(row(rep, "DIAMOND (copy of diamond_install())")$found)
        for (name in c("BLAST+ (blastp, makeblastdb)", "Perl", "KaKs_Calculator", "Clustal Omega", "MAFFT", "MUSCLE",
                       "T-Coffee", "ClustalW", "coal-cpp", "OrthoFinder")) {
                r <- row(rep, name)
                expect_identical(nrow(r), 1L, info = name)
                expect_false(r$found, info = name)
                expect_true(is.na(r$path) && is.na(r$version), info = name)
                expect_true(nzchar(r$hint), info = name)
        }
        expect_match(row(rep, "KaKs_Calculator")$hint, "kakscalculator2", fixed = TRUE)
        expect_true(all(c("coalr (COAL in R)", "coalr: GY method", "coalr: Comeron's method", "coalr: Li's method")
                        %in% rep$component))
        # installed later in the session: found at the next call
        ol3_fake_program("KaKs_Calculator", c("echo 'Program: KaKs_Calculator Toolbox'", "echo 'Version: 2.0, June. 2009'"),
                         dir = bin)
        ol3_fake_program("coal-cpp", "echo 'coal 0.6.0 (C++)'", dir = bin)
        rep2 <- orthologr_sitrep()
        expect_true(row(rep2, "KaKs_Calculator")$found)
        expect_match(row(rep2, "KaKs_Calculator")$version, "KaKs_Calculator")
        expect_true(row(rep2, "coal-cpp")$found)
        expect_identical(row(rep2, "coal-cpp")$version, "coal 0.6.0 (C++)")
        # a program that hangs on its version check is given up after the limit of the check
        skip_on_cran()
        ol3_fake_program("mafft", "exec sleep 120", dir = bin)
        local_mocked_bindings(sitrep_probe = local({
                real <- sitrep_probe
                function(program, version_args, pattern, path = NULL) {
                        if (program == "mafft") return(timeout_cap("mafft", 1, real(program, version_args, pattern, path)))
                        real(program, version_args, pattern, path)
                }
        }))
        t0 <- proc.time()[["elapsed"]]
        rep3 <- orthologr_sitrep()
        expect_false(row(rep3, "MAFFT")$found)
        expect_lt(proc.time()[["elapsed"]] - t0, 60)
})

test_that("coalr's methods are reported from coal_info()", {
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        rep <- orthologr_sitrep()
        info <- coalr::coal_info()
        r <- rep[rep$component == "coalr (COAL in R)", ]
        expect_true(r$found)
        expect_identical(r$version, as.character(utils::packageVersion("coalr")))
        expect_identical(rep$found[rep$component == "coalr: GY method"], is.element("gy94_ml", info$features))
        expect_identical(rep$found[rep$component == "coalr: Li's method"], is.element("Li", info$methods))
})

# ---- lows -------------------------------------------------------------------------------------------------------

test_that("the note on saturation names the species pairs of a run with several species", {
        status <- c(rep("saturated", 30), rep("ok", 70))
        dS <- c(rep(NA, 30), rep(c(0.1, 2.5, 3.2, 0.2, 1.9, 0.3, 0.2), 10))
        sp <- c(rep("A / B", 30), rep(c("A / C", "B / D", "C / D", "A / C", "B / C", "A / C", "A / C"), 10))
        two <- dnds_summary_text(status, dS = dS)
        expect_match(two, "the two species are probably too distant", fixed = TRUE)
        many <- dnds_summary_text(status, dS = dS, species_pair = sp)
        expect_false(grepl("the two species", many, fixed = TRUE))
        expect_match(many, "some of the species compared are probably too distant", fixed = TRUE)
        expect_match(many, "highest median dS: A / B (-, 100 % saturated); C / D (3.2, 0 % saturated); B / D (2.5, 0 % saturated)",
                     fixed = TRUE)
        one <- dnds_summary_text(status, dS = dS, species_pair = rep("A / B", 100))
        expect_match(one, "the two species (A / B) are probably too distant", fixed = TRUE)
        # a species_pair of the wrong length is not used
        expect_identical(dnds_summary_text(status, dS = dS, species_pair = "A / B"), two)
})

test_that("read.cds(): the file of corrupted sequences is named after the whole path of the input file", {
        d1 <- new_test_dir("in1_"); d2 <- new_test_dir("in2_"); out <- new_test_dir("corrupted_")
        rlang::local_options(orthologr.corrupted_cds_dir = out)
        f1 <- ol3_fasta(c(a = "ATGAAACCCGGGTA", b = "ATGAAATAA"), file.path(d1, "cds.fasta"))
        f2 <- ol3_fasta(c(c = "ATGAAACCCGG", d = "ATGAAATAA"), file.path(d2, "cds.fasta"))
        suppressMessages(read.cds(f1, format = "fasta"))
        suppressMessages(read.cds(f2, format = "fasta"))
        files <- list.files(out)
        expect_length(files, 2)
        expect_true(all(grepl("^cds\\.fasta_[0-9a-f]{8}_corrupted_cds_seqs\\.fasta$", files)))
        expect_identical(corrupted_cds_file(f1, out), corrupted_cds_file(f1, out))
        expect_false(identical(corrupted_cds_file(f1, out), corrupted_cds_file(f2, out)))
        expect_identical(names(Biostrings::readDNAStringSet(corrupted_cds_file(f2, out))), "c")
        # read again: the same file, no temporary files left
        suppressMessages(read.cds(f1, format = "fasta"))
        expect_length(list.files(out), 2)
})

test_that("dNdS() refuses engine = NULL and kaks_tool of length 2 with a classed error", {
        expect_error(dNdS(thal_cds(), lyra_cds(), engine = NULL), "engine", class = "orthologr_error_bad_argument")
        expect_error(dNdS(thal_cds(), lyra_cds(), engine = c("orthologr", "coal")), class = "orthologr_error_bad_argument")
        expect_error(dNdS(thal_cds(), lyra_cds(), kaks_tool = c("coal", "kaks_calculator")), "kaks_tool",
                     class = "orthologr_error_bad_argument")
        expect_error(dNdS(thal_cds(), lyra_cds(), kaks_tool = NA_character_), class = "orthologr_error_bad_argument")
})

test_that("dnds_bind() refuses a record without rows", {
        good <- dnds_values(dN = 0.1, dS = 0.2, dNdS = 0.5, method = NULL, reason = NA_character_)
        good$query_id <- "a"; good$subject_id <- "x"
        empty <- dnds_values(dN = numeric(0), dS = numeric(0), dNdS = numeric(0), method = NULL, reason = character(0))
        empty$query_id <- "b"; empty$subject_id <- "y"
        expect_identical(nrow(dnds_bind(list(good))[[1]]), 1L)
        expect_error(dnds_bind(list(good, empty)), "gene pair 'b' / 'y' hold no row", fixed = TRUE)
})

test_that("plain errors of these functions are classed", {
        expect_error(split_params("-a 'b"), class = "orthologr_error_bad_argument")
        h <- data.frame(query_id = "a", subject_id = "x", evalue = 1, alig_length = 1L)
        expect_error(best_hits(h, rule = "bitscore"), "needs the columns", class = "orthologr_error_bad_input")
        skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is not installed")
        expect_error(write_dNdS(data.frame(a = 1), NA), class = "orthologr_error_bad_argument")
        expect_error(write_dNdS(1:3, tempfile()), class = "orthologr_error_bad_argument")
})

test_that("diamond_install(): a staged copy has a name of its own; a cached copy of the release asked for is no warning", {
        skip_on_os("windows")
        old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA)
        on.exit(if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR") else Sys.setenv(R_USER_CACHE_DIR = old), add = TRUE)
        Sys.setenv(R_USER_CACHE_DIR = new_test_dir("cache_"))
        dir <- diamond_cache_dir()
        dir.create(dir, recursive = TRUE)
        writeLines(c("#!/bin/sh", "echo 'diamond version 2.2.8'"), file.path(dir, "diamond"))
        Sys.chmod(file.path(dir, "diamond"), "755")
        expect_identical(diamond_cached_version(file.path(dir, "diamond")), "2.2.8")
        expect_message(diamond_install(version = "v2.2.8"), "already installed")
        expect_message(diamond_install(version = "2.2.8"), "already installed")
        expect_silent(diamond_install(version = "v2.2.8", quiet = TRUE))
        expect_silent(diamond_install(sha256 = strrep("a", 64), quiet = TRUE))
        expect_warning(diamond_install(version = "v2.1.0"), "DIAMOND v2.2.8 is already installed.*v2.1.0 was asked for",
                       class = "orthologr_warning_tool_version")
        # stale staged copies of both name forms are removed, fresh ones and other files are kept
        old1 <- file.path(dir, ".new-999-diamond")
        old2 <- file.path(dir, ".new-999-1a2b3c-diamond")
        fresh <- file.path(dir, ".new-998-4d5e6f-diamond")
        other <- file.path(dir, ".new-x-diamond")
        file.create(old1, old2, fresh, other)
        Sys.setFileTime(c(old1, old2, other), Sys.time() - 7200)
        removed <- diamond_remove_stale_staged(dir, "diamond")
        expect_setequal(basename(removed), basename(c(old1, old2)))
        expect_true(all(file.exists(c(fresh, other))))
})

# ---- after merging fix/review-ol-3b: AUDIT-OLB lows L1-L3 -----------------------------------------------------------

test_that("dNdS_orthogroups(): an orthogroups object reports the same input tables as the files (L1)", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        toy <- make_toy_species()
        a <- readLines(toy$files[["A"]])
        seq_lines <- which(!startsWith(a, ">"))
        substr(a[seq_lines[1]], 13, 15) <- "TAA"                       # an in-frame stop codon
        substr(a[seq_lines[2]], 20, 20) <- "N"                         # an ambiguous base
        writeLines(a, toy$files[["A"]])
        from_files <- suppressWarnings(suppressMessages(dNdS_orthogroups(toy$files, comp_cores = 1)))
        og <- suppressMessages(orthogroups(toy$files))
        from_object <- suppressWarnings(suppressMessages(dNdS_orthogroups(og, comp_cores = 1)))
        for (at in c("ambiguous_sequences", "stop_codon_sequences", "dropped_sequences", "trimmed_sequences"))
                expect_identical(attr(from_object, at), attr(from_files, at), info = at)
        expect_identical(nrow(attr(from_object, "ambiguous_sequences")), 1L)
        expect_identical(nrow(attr(from_object, "stop_codon_sequences")), 1L)
        expect_identical(from_object$pairs, from_files$pairs)
})

test_that("a coal-cpp without GY: the message also names 'pip install -U coal-align' (L2)", {
        local_mocked_bindings(coal_backend = function(need = NULL) list(type = "binary", path = "/x/coal-cpp"),
                              coal_has_gy = function() FALSE)
        e <- tryCatch(coal_check_methods("GY", "engine = \"coal\""), error = function(e) e)
        expect_s3_class(e, "orthologr_error_package_missing")
        expect_match(conditionMessage(e), "pip install -U coal-align", fixed = TRUE)
        expect_match(conditionMessage(e), coalr_install_hint(), fixed = TRUE)
})

test_that("substitutionrate() names est.method, kaks_tool = \"coal\" lists every method it takes (L3)", {
        f <- system.file("seqs/pal2nal.aln", package = "orthologr")
        expect_error(substitutionrate(f, est.method = "XX"), "^est\\.method = \"XX\" is not a method of substitutionrate\\(\\)",
                     class = "orthologr_error_bad_argument")
        e <- tryCatch(substitutionrate(f, est.method = "GMYN", kaks_tool = "coal"), error = function(e) e)
        expect_s3_class(e, "orthologr_error_bad_argument")
        for (m in kaks_tool_coal_methods())
                expect_match(conditionMessage(e), paste0("\"", m, "\""), fixed = TRUE, info = m)
        expect_true(all(c("Comeron", "Li") %in% kaks_tool_coal_methods()))
        expect_match(conditionMessage(e), "takes est.method = ", fixed = TRUE)
        e <- tryCatch(dNdS(thal_cds(), lyra_cds(), dnds_est.method = "GMYN", kaks_tool = "coal"), error = function(e) e)
        expect_match(conditionMessage(e), "takes dnds_est.method = \"Comeron\", \"Li\"", fixed = TRUE)
})

test_that("number_text(): every double reads back as the same number, also very small e-values", {
        x <- c(1.37e-215, 2.17e-180, 1.2e-215, 4.47e-161, 0.1, 1 / 3, pi, 1e-300, 5e-324, .Machine$double.xmax, 0, -2.5, NA,
               Inf, -Inf, NaN)
        txt <- number_text(x)
        back <- suppressWarnings(as.numeric(txt))
        fin <- is.finite(x)
        expect_identical(back[fin], x[fin])
        expect_true(all(is.na(txt[is.na(x)])))
        # the hexadecimal fallback is read exactly by R
        expect_identical(as.numeric(sprintf("%a", x[fin])), x[fin])
        expect_identical(number_text(0.25), "0.25")
})
