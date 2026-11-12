# External programs are run through one helper: status and output checks, quoting, no stale results, no change of
# the working directory or the PATH (audit 2026-11-11, problem 1; stress audit F12).

rscript_dir <- R.home("bin")

test_that("run_tool() passes arguments with spaces as one argument and returns status, stdout and stderr", {
        res <- run_tool("Rscript", args = list("-e", "cat(commandArgs(TRUE), sep = '|'); message('to stderr')", "a b", "c"),
                        path = rscript_dir)
        expect_identical(res$status, 0L)
        expect_identical(res$stdout, "a b|c")
        expect_true(any(grepl("to stderr", res$stderr)))
})

test_that("run_tool() passes I() arguments unquoted (free-text parameter strings)", {
        res <- run_tool("Rscript", args = list("-e", "cat(length(commandArgs(TRUE)))", I("x y z")), path = rscript_dir)
        expect_identical(res$stdout, "3")
})

test_that("run_tool() stops on a non-zero exit status and quotes what the program said", {
        expect_error(
                run_tool("Rscript", args = list("-e", "message('Error: bad input'); quit(status = 3)"),
                         path = rscript_dir, tool = "MyTool", context = "gene pair 'a' / 'b'"),
                "MyTool failed for gene pair 'a' / 'b': exit status 3.*bad input"
        )
})

test_that("run_tool() stops when the program is not there", {
        expect_error(run_tool("orthologr_no_such_program_xyz", args = "-h", tool = "NoSuch"),
                     "NoSuch failed: exit status 127")
        expect_error(run_tool("orthologr_no_such_program_xyz", args = "-h", path = tempdir(), tool = "NoSuch"),
                     "NoSuch failed")
})

test_that("run_tool() stops when the expected output file is missing or empty although the status is 0", {
        out <- tempfile()
        expect_error(
                run_tool("Rscript", args = list("-e", "invisible(1)"), path = rscript_dir, tool = "MyTool", output = out),
                "MyTool wrote no output"
        )
        expect_error(run_tool("Rscript", args = list("-e", "invisible(1)"), path = rscript_dir, tool = "MyTool",
                              stdout = out, output = out), "MyTool wrote no output")
        expect_silent(run_tool("Rscript", args = list("-e", "cat('x')"), path = rscript_dir, stdout = out, output = out))
        expect_identical(readLines(out, warn = FALSE), "x")
})

test_that("run_tool(error = FALSE) returns the failure instead of stopping", {
        res <- run_tool("Rscript", args = list("-e", "quit(status = 2)"), path = rscript_dir, error = FALSE)
        expect_identical(res$status, 2L)
})

test_that("run_tool() leaves the working directory and the PATH as they were", {
        wd <- getwd()
        p <- Sys.getenv("PATH")
        run_tool("Rscript", args = list("-e", "invisible(1)"), path = rscript_dir)
        try(run_tool("Rscript", args = list("-e", "quit(status = 1)"), path = rscript_dir), silent = TRUE)
        expect_identical(getwd(), wd)
        expect_identical(Sys.getenv("PATH"), p)
})

test_that("tool_available() starts a program once per session and path, and reports a missing program", {
        skip_on_os("windows")
        d <- new_test_dir()
        log <- file.path(d, "calls.log")
        script <- file.path(d, "faketool")
        writeLines(c("#!/bin/sh", paste0("echo call >> ", shQuote(log)), "echo 'FakeTool version 1.2'"), script)
        Sys.chmod(script, "755")
        tool_cache_reset()
        expect_identical(tool_available("faketool", "--version", "FakeTool", path = d), "FakeTool version 1.2")
        expect_identical(tool_available("faketool", "--version", "FakeTool", path = d), "FakeTool version 1.2")
        expect_length(readLines(log), 1L)
        expect_error(tool_available("faketool", "--version", "FakeTool", path = file.path(d, "nowhere"),
                                    not_found = "FakeTool is not installed"), "FakeTool is not installed")
        # a failed lookup is not remembered
        expect_error(tool_available("faketool", "--version", "FakeTool", path = file.path(d, "nowhere")), "not found")
        tool_cache_reset()
})

test_that("the is_installed_*() checks return TRUE and do not start the program again", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        expect_true(is_installed_diamond())
        t <- system.time(for (i in 1:20) is_installed_diamond())[["elapsed"]]
        expect_lt(t, 0.5)
})

test_that("safe_file_name() and blast_db_arg() handle spaces", {
        expect_identical(safe_file_name("my file (1).fasta"), "my_file__1_.fasta")
        expect_identical(blast_db_arg("/a/b"), "/a/b")
        expect_identical(blast_db_arg("/a b/c"), "\"/a b/c\"")
})

# ---------------------------------------------------------------- stale results

test_that("blast() stops when BLAST fails instead of returning the table of the previous call", {
        skip_if_not(blast_available, "BLAST+ is not installed")
        ok <- suppressMessages(blast(thal_cds(), lyra_cds()))
        expect_equal(nrow(ok), 20)
        # the example of the old man page: BLAST rejects the duplicated flag
        expect_error(
                suppressMessages(blast(thal_cds(), lyra_cds(), blast_params = "-max_target_seqs 1")),
                "blastp failed for query 'ortho_thal_cds.fasta': exit status.*max_target_seqs"
        )
        expect_false(file.exists(file.path(tempdir(), "_blast_db", "blastresult_ortho_thal_cds.fasta.csv")))
})

test_that("blast() keeps its hit table where it was documented and copies it to save.output", {
        skip_if_not(blast_available, "BLAST+ is not installed")
        out <- new_test_dir()
        suppressMessages(blast(thal_cds(), lyra_cds(), save.output = out))
        expect_true(file.exists(file.path(tempdir(), "_blast_db", "blastresult_ortho_thal_cds.fasta.csv")))
        expect_true(file.exists(file.path(out, "blastresult_ortho_thal_cds.fasta.csv")))
        expect_equal(nrow(utils::read.delim(file.path(out, "blastresult_ortho_thal_cds.fasta.csv"), header = FALSE)), 20)
})

test_that("diamond() stops when DIAMOND fails, says what DIAMOND said, and leaves no table behind", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        ok <- suppressMessages(diamond(thal_cds(), lyra_cds()))
        expect_equal(nrow(ok), 20)
        expect_error(
                suppressMessages(diamond(thal_cds(), lyra_cds(), diamond_params = "--no-such-option")),
                "exited with non-zero status.*diamond said"
        )
        expect_false(file.exists(file.path(tempdir(), "_blast_db", "blastresult_ortho_thal_cds.fasta.csv")))
})

test_that("substitutionrate() does not return the rows of an earlier KaKs_Calculator run", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        aln <- system.file("seqs/pal2nal.aln", package = "orthologr")
        first <- substitutionrate(aln, est.method = "YN", quiet = TRUE)
        expect_equal(nrow(first), 1)
        # "kaks_calc" without kaks_calc.params is not a KaKs_Calculator model: an error, not the rows of 'first'
        expect_error(substitutionrate(aln, est.method = "kaks_calc", quiet = TRUE), "KaKs_Calculator")
        ng <- substitutionrate(aln, est.method = "NG", quiet = TRUE)
        expect_identical(ng$method, "NG")
        # nothing is left in the calculation folder
        expect_length(list.files(file.path(tempdir(), "_calculation")), 0L)
})

test_that("multi_aln() stops when the aligner fails and leaves no alignment of an earlier call behind", {
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        aa <- system.file("seqs/multi_aln_example.fasta", package = "orthologr")
        expect_s4_class(multi_aln(aa, tool = "clustalo", get_aln = TRUE, quiet = TRUE), "AAMultipleAlignment")
        out <- file.path(tempdir(), "_alignment", "multi_aln", "multi_aln_example.fasta_clustalo.aln")
        expect_true(file.exists(out))
        expect_error(multi_aln(aa, tool = "clustalo", get_aln = TRUE, quiet = TRUE, params = "--no-such-option"),
                     "Clustal Omega failed for 'multi_aln_example.fasta'")
        expect_false(file.exists(out))
        expect_length(list.files(dirname(out), pattern = "\\.tmp$"), 0L)
})

test_that("multi_aln() returns the path of the alignment invisibly when get_aln = FALSE", {
        skip_if_not(mafft_available, "MAFFT is not installed")
        aa <- system.file("seqs/multi_aln_example.fasta", package = "orthologr")
        res <- withVisible(multi_aln(aa, tool = "mafft", quiet = TRUE, multi_aln_name = "named"))
        expect_false(res$visible)
        expect_identical(res$value, file.path(tempdir(), "_alignment", "multi_aln", "named_mafft.aln"))
        expect_gt(file.size(res$value), 0)
})

test_that("select_orthologs() does not reuse the table of an earlier call with another dnds_tbl", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        all <- suppressMessages(select_orthologs(x$dnds, x$query$gtf, x$subject$gtf, format = c("gtf", "gtf")))
        expect_equal(nrow(all), 20)
        three <- suppressMessages(select_orthologs(x$dnds[6:8, ], x$query$gtf, x$subject$gtf, format = c("gtf", "gtf")))
        expect_equal(nrow(three), 3)
        expect_setequal(three$query_id, x$dnds$query_id[6:8])
        expect_length(list.files(tempdir(), pattern = "final_join_dnds_annotation_tbl"), 0L)
})

# ---------------------------------------------------------------- file names with spaces

test_that("file and folder names with spaces work in the search wrappers", {
        d <- new_test_dir("orthologr test ")
        q <- file.path(d, "thal cds.fasta")
        s <- file.path(d, "lyra cds.fasta")
        file.copy(thal_cds(), q)
        file.copy(lyra_cds(), s)
        if (diamond_available) {
                expect_identical(suppressMessages(diamond(q, s)), suppressMessages(diamond(thal_cds(), lyra_cds())))
                expect_equal(nrow(suppressMessages(diamond_rec(q, s))), 20)
                db <- set_diamond(q, makedb = TRUE)[[2]]
                expect_identical(db, "diamonddb_thal_cds.fasta_protein.fasta")
                expect_true(file.exists(file.path(tempdir(), "_blast_db", paste0(db, ".dmnd"))))
                expect_equal(nrow(suppressMessages(deepclust(c(q, s), seq_type = "cds"))), 40)
        }
        if (blast_available) {
                expect_identical(suppressMessages(blast(q, s)), suppressMessages(blast(thal_cds(), lyra_cds())))
                expect_equal(nrow(suppressMessages(blast_rec(q, s))), 20)
                expect_identical(set_blast(q, makedb = TRUE)[[2]], "blastdb_thal_cds.fasta_protein.fasta")
        }
        skip_if_not(diamond_available || blast_available, "neither DIAMOND2 nor BLAST+ is installed")
})

test_that("dNdS() works with input files and a working directory whose names hold spaces", {
        skip_if_not(diamond_available && perl_available, "DIAMOND2 or perl is not installed")
        d <- new_test_dir("orthologr test ")
        q <- file.path(d, "thal cds.fasta")
        s <- file.path(d, "lyra cds.fasta")
        file.copy(thal_cds(), q)
        file.copy(lyra_cds(), s)
        wd <- getwd()
        p <- Sys.getenv("PATH")
        res <- suppressMessages(dNdS(q, s, print_citation = FALSE))
        ref <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE))
        expect_equal(as.data.frame(res), as.data.frame(ref), ignore_attr = TRUE)
        expect_identical(getwd(), wd)
        expect_identical(Sys.getenv("PATH"), p)
})

test_that("multi_aln() and codon_aln() accept file names with spaces", {
        skip_if_not(clustalo_available && perl_available, "Clustal Omega or perl is not installed")
        d <- new_test_dir("orthologr test ")
        aa <- file.path(d, "aa seqs.fasta")
        file.copy(system.file("seqs/multi_aln_example.fasta", package = "orthologr"), aa)
        expect_s4_class(multi_aln(aa, tool = "clustalo", get_aln = TRUE, quiet = TRUE), "AAMultipleAlignment")
        aln <- file.path(d, "aa seqs.aln")
        nuc <- file.path(d, "dna seqs.fasta")
        file.copy(system.file("seqs/aa_seqs.aln", package = "orthologr"), aln)
        file.copy(system.file("seqs/dna_seqs.fasta", package = "orthologr"), nuc)
        out <- suppressWarnings(codon_aln(aln, nuc, format = "fasta", codon_aln_name = "spaces", quiet = TRUE))
        expect_identical(out, file.path(tempdir(), "_alignment", "codon_aln", "spaces_pal2nal.aln"))
        expect_gt(file.size(out), 0)
})

# ---------------------------------------------------------------- per-call folders

test_that("search wrappers leave no per-call folder behind and do not change the working directory", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        wd <- getwd()
        before <- list.files(tempdir())
        suppressMessages(diamond_rec(thal_cds(), lyra_cds()))
        suppressMessages(deepclust(thal_aa()))
        expect_identical(getwd(), wd)
        expect_setequal(setdiff(list.files(tempdir()), c(before, "_blast_db")), character(0))
        expect_length(list.files(tempdir(), pattern = "^(diamond|blast|deepclust)_"), 0L)
})
