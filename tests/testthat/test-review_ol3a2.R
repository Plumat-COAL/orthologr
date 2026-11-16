# AUDIT-OLA fixes: the time limit of diamond makedb grows with its input (M1), deepclust(direct_input = TRUE)
# refuses nucleotide files (M2), options(orthologr.check_protein) and letters only proteins have (L1),
# write_dNdS() / read_dNdS() keep the types of the record (L2), engine = NULL is named (L3), and a broken OrthoFinder
# is not taken for an installed one.

ola_fake_program <- function(name, lines, dir = new_test_dir("fake_tool_")) {
        writeLines(c("#!/bin/sh", lines), file.path(dir, name))
        Sys.chmod(file.path(dir, name), "755")
        dir
}

ola_fasta <- function(seqs, file = tempfile(fileext = ".fasta")) {
        writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), file)
        file
}

# ---- M1: diamond makedb ------------------------------------------------------------------------------------------

test_that("the default limit of diamond makedb is 1 hour, or the time to read the input at 2 MB/s (M1)", {
        expect_identical(makedb_time_limit(bytes = 1e6), 3600)
        expect_identical(makedb_time_limit(bytes = 7.2e9), 3600)
        expect_identical(makedb_time_limit(bytes = 1e11), 50000)
        f <- ola_fasta(c(p1 = "MKLVAGTWQERTYIPASDFGHKLCVNM"))
        expect_identical(makedb_time_limit(f), 3600)
        rlang::local_options(orthologr.tool_timeout = NULL)
        # the default of the run replaces that of the program; options still decide
        expect_identical(tool_timeout("diamond makedb", "diamond", default = 50000), 50000)
        expect_identical(tool_timeout("diamond blastp", "diamond", default = 50000), 50000)
        expect_identical(tool_timeout("diamond makedb", "diamond"), 3600)
        rlang::local_options(orthologr.tool_timeout = c("diamond makedb" = 5))
        expect_identical(tool_timeout("diamond makedb", "diamond", default = 50000), 5)
        rlang::local_options(orthologr.tool_timeout = c(DIAMOND = 20))
        expect_identical(tool_timeout("diamond makedb", "diamond", default = 50000), 20)
        rlang::local_options(orthologr.tool_timeout = 7)
        expect_identical(tool_timeout("diamond makedb", "diamond", default = 50000), 7)
})

test_that("diamond_makedb() gives run_tool() the limit for the size of its input (M1)", {
        seen <- NULL
        local_mocked_bindings(run_tool = function(program, args = list(), ..., limit = NULL) {
                seen <<- limit
                invisible(list(status = 0L, stdout = NULL, stderr = character(0)))
        })
        f <- ola_fasta(c(p1 = "MKLVAGTWQERTYIPASDFGHKLCVNM"))
        diamond_makedb(f, tempfile("db_"))
        expect_identical(seen, 3600)
        local_mocked_bindings(makedb_time_limit = function(in_file, bytes = NULL) 12345)
        diamond_makedb(f, tempfile("db_"))
        expect_identical(seen, 12345)
})

# ---- M2: deepclust(direct_input = TRUE) ------------------------------------------------------------------------

test_that("deepclust(direct_input = TRUE) refuses nucleotide files before DIAMOND is started (M2)", {
        skip_on_os("windows")
        # the stand-in would hang in makedb: it must not be reached
        bin <- ola_fake_program("diamond", c("case \"$1\" in version) echo 'diamond version 2.2.8';; *) exec sleep 30;; esac"))
        tool_cache_reset()
        on.exit(tool_cache_reset(), add = TRUE)
        t0 <- proc.time()[["elapsed"]]
        e <- tryCatch(suppressMessages(deepclust(thal_cds(), path = bin, direct_input = TRUE)), error = function(e) e)
        expect_s3_class(e, "orthologr_error_bad_input")
        expect_match(conditionMessage(e), "ortho_thal_cds.fasta' does not look like protein sequences", fixed = TRUE)
        # the second of several files, and a compressed file
        expect_error(suppressMessages(deepclust(c(thal_aa(), lyra_cds()), path = bin, direct_input = TRUE)),
                     "ortho_lyra_cds.fasta", class = "orthologr_error_bad_input")
        gz <- tempfile(fileext = ".fasta.gz")
        con <- gzfile(gz, "w"); writeLines(readLines(thal_cds()), con); close(con)
        expect_error(suppressMessages(deepclust(gz, path = bin, direct_input = TRUE)), class = "orthologr_error_bad_input")
        expect_lt(proc.time()[["elapsed"]] - t0, 20)
})

test_that("with options(orthologr.check_protein = FALSE) the files of direct_input go to DIAMOND unchecked (M2, L1)", {
        skip_on_os("windows")
        bin <- ola_fake_program("diamond", c("case \"$1\" in version) echo 'diamond version 2.2.8';; *) echo 'refused' >&2; exit 3;; esac"))
        tool_cache_reset()
        on.exit(tool_cache_reset(), add = TRUE)
        # refused, and the message names the option
        expect_error(suppressMessages(deepclust(thal_cds(), path = bin, direct_input = TRUE)),
                     "options(orthologr.check_protein = FALSE)", fixed = TRUE, class = "orthologr_error_bad_input")
        rlang::local_options(orthologr.check_protein = FALSE)
        e <- tryCatch(suppressMessages(deepclust(thal_cds(), path = bin, direct_input = TRUE)), error = function(e) e)
        expect_s3_class(e, "orthologr_error_tool_failed")
        expect_false(inherits(e, "orthologr_error_bad_input"))
})

# ---- L1: the check of protein files ------------------------------------------------------------------------------

test_that("a low-complexity protein can be read with options(orthologr.check_protein = FALSE), which the error names (L1)", {
        silk <- ola_fasta(c(silk = strrep("GAGAGSGAGAGT", 10)))
        e <- tryCatch(read.proteome(silk, format = "fasta"), error = function(e) e)
        expect_s3_class(e, "orthologr_error_bad_input")
        expect_match(conditionMessage(e), "options(orthologr.check_protein = FALSE)", fixed = TRUE)
        rlang::local_options(orthologr.check_protein = FALSE)
        expect_identical(read.proteome(silk, format = "fasta")$geneids, "silk")
        expect_identical(nrow(read.proteome(thal_cds(), format = "fasta")), 20L)
        rlang::local_options(orthologr.check_protein = "no")
        expect_error(read.proteome(silk, format = "fasta"), "TRUE or FALSE", class = "orthologr_error_bad_argument")
})

test_that("letters only proteins have (E, F, I, L, P, Q) make a file of proteins; nucleotide files are still refused (L1)", {
        # 95 % A, C, G, T, U, N, and Q / P: an asparagine-rich protein (as in Dictyostelium), not nucleotides
        polyn <- ola_fasta(c(d1 = paste0("M", strrep("NNNNTNNNNA", 3), "NNNNGNNNQ"),
                             d2 = paste0("M", strrep("NNNNNNNNNA", 3), "NNNNPNNNN")))
        expect_identical(read.proteome(polyn, format = "fasta")$geneids, c("d1", "d2"))
        # one stray protein letter in a nucleotide file is below 1 %: refused, as the coding sequences and RNA
        cds <- readLines(thal_cds())
        seq_line <- which(!startsWith(cds, ">"))[1]
        substr(cds[seq_line], 5, 5) <- "E"
        stray <- tempfile(fileext = ".fasta")
        writeLines(cds, stray)
        for (f in c(thal_cds(), lyra_cds(), stray))
                expect_error(read.proteome(f, format = "fasta"), class = "orthologr_error_bad_input")
        rna <- ola_fasta(c(r1 = paste(rep("AUGGCUAAAGAUGAACCC", 10), collapse = "")))
        expect_error(read.proteome(rna, format = "fasta"), class = "orthologr_error_bad_input")
        # coding sequences with a few IUPAC ambiguity codes (none of them is E, F, I, L, P or Q) as well
        amb <- ola_fasta(c(a1 = paste(rep("ATGGCTAAAGATGAACCR", 10), collapse = "")))
        expect_error(read.proteome(amb, format = "fasta"), class = "orthologr_error_bad_input")
})

# ---- L2: write_dNdS() / read_dNdS() ------------------------------------------------------------------------------

test_that("write_dNdS() / read_dNdS(): the record keeps its types, NaN stays NaN, the table coal keeps 'seconds' (L2)", {
        skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is not installed")
        tbl <- tibble::tibble(query_id = c("a", "b", "c"), subject_id = c("x", "y", "z"), dN = c(0.1, NaN, 1 / 3),
                              dS = c(0.2, NA, Inf), dNdS = c(0.5, NA, 0), status = c("ok", "failed", "saturated"))
        coal <- tibble::tibble(seq1 = c("a", "b"), seq2 = c("x", "y"), dN = c(0.1, NaN), dS = c(pi, NA))
        attr(coal, "seconds") <- c(align_codon = 0.00969, kaks = 1 / 3)
        attr(tbl, "coal") <- coal
        attr(tbl, "orthologr_run") <- list(
                orthologr = "0.6.0", "function" = "dNdS_coal",
                arguments = list(comp_cores = 1, eval = "1E-5", aligner_path = NULL, dnds_est.method = "Comeron"),
                called_as = list("function" = "dNdS", arguments = list(comp_cores = 2, files = character(0),
                                                                       empty = list(), named_empty = stats::setNames(list(), character(0)))),
                options = list(orthologr.legacy = FALSE, orthologr.comeron_min_sites = 3, orthologr.max_aln_cells = 5e8,
                               orthologr.tool_timeout = c(KaKs_Calculator = Inf, DIAMOND = 600), orthologr.flag = NA),
                engine = "coal",
                input_files = data.frame(file = c("q.fa", "s.fa"), size = c(100, 2e10), md5 = c("abc", NA),
                                         stringsAsFactors = FALSE),
                status_counts = list(ok = 1L, failed = 1L, saturated = 1L))
        file <- tempfile(fileext = ".tsv")
        write_dNdS(tbl, file)
        back <- read_dNdS(file)
        expect_identical(attr(back, "orthologr_run"), attr(tbl, "orthologr_run"))
        expect_identical(attr(back, "coal"), attr(tbl, "coal"))
        expect_identical(attr(attr(back, "coal"), "seconds"), c(align_codon = 0.00969, kaks = 1 / 3))
        expect_true(is.nan(back$dN[2]))
        expect_true(is.na(back$dS[2]) && !is.nan(back$dS[2]))
        expect_identical(back$dS[3], Inf)
        expect_true(dnds_equal(tbl, back, ignore = character(0)))
        # the readable record is still in the file
        side <- jsonlite::fromJSON(paste0(file, ".run.json"), simplifyVector = TRUE, simplifyDataFrame = FALSE)
        expect_identical(side$run$orthologr, "0.6.0")
        expect_identical(side$run[["function"]], "dNdS_coal")
        # a file without the typed record (written before it) is read from the readable one
        side$run_r <- NULL
        writeLines(jsonlite::toJSON(side, auto_unbox = TRUE, null = "null", na = "null", digits = NA),
                   paste0(file, ".run.json"))
        old <- read_dNdS(file)
        expect_identical(attr(old, "orthologr_run")$engine, "coal")
        expect_identical(attr(old, "orthologr_run")$input_files$size, c(100, 2e10))
})

test_that("a run of dNdS(engine = \"coal\") reads back equal; called_as has no internal mark (L2)", {
        skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is not installed")
        skip_if_not(requireNamespace("coalr", quietly = TRUE), "coalr is not installed")
        skip_if_not(has_program("diamond"), "DIAMOND is not installed")
        r <- suppressMessages(dNdS(thal_cds(), lyra_cds(), engine = "coal", print_citation = FALSE))
        expect_null(attr(attr(r, "orthologr_run")$called_as$arguments, "command_mark"))
        file <- tempfile(fileext = ".tsv")
        write_dNdS(r, file)
        back <- read_dNdS(file)
        expect_true(dnds_equal(r, back, ignore = character(0)))
        expect_false(is.null(attr(attr(back, "coal"), "seconds")))
})

# ---- L3: engine = NULL -------------------------------------------------------------------------------------------

test_that("the message for a value that is not one text shows the value as R code (L3)", {
        expect_error(dNdS(thal_cds(), lyra_cds(), engine = NULL), "engine = NULL is not supported", fixed = TRUE)
        expect_error(dNdS(thal_cds(), lyra_cds(), engine = c("orthologr", "coal")),
                     "engine = c(\"orthologr\", \"coal\") is not supported", fixed = TRUE)
        expect_error(dNdS(thal_cds(), lyra_cds(), kaks_tool = NA_character_), "kaks_tool = NA_character_ is not supported",
                     fixed = TRUE)
        expect_error(dNdS(thal_cds(), lyra_cds(), aligner = "diamnd"), "aligner = 'diamnd' is not supported", fixed = TRUE)
})

# ---- OrthoFinder that does not run -------------------------------------------------------------------------------

test_that("an OrthoFinder that stops with a traceback is not installed, also not for orthologr_sitrep() after the check", {
        skip_on_os("windows")
        old <- Sys.getenv(c("PATH", "R_USER_CACHE_DIR", "COAL_CPP"), unset = NA)
        on.exit({
                for (v in names(old)) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
                tool_cache_reset()
        }, add = TRUE)
        broken <- with_coreutils(ola_fake_program("orthofinder", c(
                "echo 'Traceback (most recent call last):' >&2",
                "echo '  File \"/opt/conda/bin/OrthoFinder/orthofinder.py\", line 12, in <module>' >&2",
                "echo \"ModuleNotFoundError: No module named 'scipy'\" >&2",
                "exit 1")))
        tool_cache_reset()
        e <- tryCatch(is_installed_orthofinder(path = broken), error = function(e) e)
        expect_s3_class(e, "orthologr_error_tool_failed")
        Sys.setenv(PATH = broken, R_USER_CACHE_DIR = new_test_dir("cache_"))
        Sys.unsetenv("COAL_CPP")
        rlang::local_options(orthologr.coal_path = NULL)
        expect_error(is_installed_orthofinder(), class = "orthologr_error_tool_failed")
        rep <- orthologr_sitrep()
        of <- rep[rep$component == "OrthoFinder", ]
        expect_false(of$found)
        expect_true(is.na(of$version))
        expect_match(of$hint, "could not be started", fixed = TRUE)
        # the version line with an exit status other than 0 is no installed OrthoFinder either
        tool_cache_reset()
        ola_fake_program("orthofinder", c("echo 'OrthoFinder version 2.5.5 Copyright (C) 2014 David Emms'", "exit 1"),
                         dir = broken)
        expect_error(is_installed_orthofinder(), class = "orthologr_error_tool_failed")
        expect_false(orthologr_sitrep()$found[rep$component == "OrthoFinder"])
        # one that runs is found, by both
        tool_cache_reset()
        ola_fake_program("orthofinder", c("echo 'OrthoFinder version 2.5.5 Copyright (C) 2014 David Emms'", "exit 0"),
                         dir = broken)
        expect_true(is_installed_orthofinder())
        rep2 <- orthologr_sitrep()
        expect_true(rep2$found[rep2$component == "OrthoFinder"])
        expect_identical(rep2$version[rep2$component == "OrthoFinder"], "OrthoFinder version 2.5.5 Copyright (C) 2014 David Emms")
})
