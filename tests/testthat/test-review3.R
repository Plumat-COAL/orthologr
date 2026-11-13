# Regression tests for the findings of the third review (R01 - R23); one block per finding.

# ---------------------------------------------------------------- R07: parameter strings are argument lists, not shell text

test_that("R07: split_params() splits like a shell and gives no character a meaning", {
        expect_identical(split_params(NULL), character(0))
        expect_identical(split_params(""), character(0))
        expect_identical(split_params("  -seg no   -word_size 3 "), c("-seg", "no", "-word_size", "3"))
        expect_identical(split_params("--force; touch INJ"), c("--force;", "touch", "INJ"))
        expect_identical(split_params("-a 'b c' \"d e\" f\\ g"), c("-a", "b c", "d e", "f g"))
        expect_identical(split_params("-a \"x\\\"y\" ''"), c("-a", "x\"y", ""))
        expect_identical(split_params("$(touch X) `touch Y` | cat > f && ls *"),
                         c("$(touch", "X)", "`touch", "Y`", "|", "cat", ">", "f", "&&", "ls", "*"))
        expect_identical(split_params("-a\t-b\n-c"), c("-a", "-b", "-c"))
        expect_identical(split_params(c("-a 1", "-b 2")), c("-a", "1", "-b", "2"))
        expect_error(split_params("-a 'open"), "Unbalanced quote")
        # every token is quoted for the shell
        skip_on_os("windows")
        expect_identical(quote_args(list("-i", I("--id 10; touch X"))), c("'-i'", "'--id'", "'10;'", "'touch'", "'X'"))
        expect_identical(quote_args(list(I(""))), character(0))
        res <- run_tool("echo", args = list(I("a; echo INJECTED $(echo sub) `echo bq` > /dev/null | cat")), error = FALSE)
        expect_identical(res$stdout, "a; echo INJECTED $(echo sub) `echo bq` > /dev/null | cat")
})

test_that("R07: an injection attempt in any parameter string runs no command", {
        skip_on_cran()
        skip_on_os("windows")
        d <- new_test_dir("r07_")
        n <- 0
        attempt <- function(expr_fun) {
                n <<- n + 1
                marker <- file.path(d, paste0("INJ", n))
                # (the injected command makes a folder: a program that takes the path for the name of an output
                # file of its own - T-Coffee writes its log there - is not mistaken for an executed command)
                for (inj in c(paste0("; mkdir ", marker), paste0("$(mkdir ", marker, ")"), paste0("`mkdir ", marker, "`"),
                              paste0("| mkdir ", marker), paste0("&& mkdir ", marker), paste0("\n mkdir ", marker))) {
                        try(suppressWarnings(suppressMessages(utils::capture.output(expr_fun(inj)))), silent = TRUE)
                        expect_false(dir.exists(marker), info = paste(n, inj))
                        unlink(marker)
                }
        }
        aln <- system.file("seqs/aa_seqs.fasta", package = "orthologr")
        if (clustalo_available) attempt(function(p) multi_aln(aln, tool = "clustalo", params = paste("--force", p), quiet = TRUE))
        if (mafft_available) attempt(function(p) multi_aln(aln, tool = "mafft", params = paste("--quiet", p), quiet = TRUE))
        if (muscle_available) attempt(function(p) multi_aln(aln, tool = "muscle", params = p, quiet = TRUE))
        if (clustalw_available) attempt(function(p) multi_aln(aln, tool = "clustalw", params = paste("-quiet", p), quiet = TRUE))
        if (tcoffee_available) attempt(function(p) multi_aln(aln, tool = "t_coffee", params = paste("-quiet", p), quiet = TRUE))
        if (perl_available) attempt(function(p) codon_aln(
                file_aln = system.file("seqs/aa_seqs.aln", package = "orthologr"),
                file_nuc = system.file("seqs/dna_seqs.fasta", package = "orthologr"),
                format = "clustal", codon_aln_name = "r07", params = p, quiet = TRUE))
        if (kaks_available) attempt(function(p) substitutionrate(
                system.file("seqs/pal2nal.aln", package = "orthologr"), est.method = "YN", kaks_calc.params = paste("-c 1", p), quiet = TRUE))
        if (diamond_available) {
                attempt(function(p) diamond(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = paste("--id 10", p)))
                attempt(function(p) diamond_best(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = p))
                attempt(function(p) diamond_rec(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = p))
                if (perl_available) {
                        attempt(function(p) dNdS(thal_cds(), lyra_cds(), aligner_params = p, print_citation = FALSE))
                        if (clustalo_available) attempt(function(p) dNdS(
                                thal_cds(), lyra_cds(), aa_aln_type = "multiple", aa_aln_tool = "clustalo",
                                aa_aln_params = paste("--force", p), print_citation = FALSE))
                }
        }
        if (blast_available) {
                attempt(function(p) blast(thal_aa(), lyra_aa(), seq_type = "protein", blast_params = paste("-seg no", p)))
                attempt(function(p) blast_best(thal_aa(), lyra_aa(), seq_type = "protein", blast_params = p))
                attempt(function(p) blast_rec(thal_aa(), lyra_aa(), seq_type = "protein", blast_params = p))
        }
        expect_length(list.dirs(d, recursive = FALSE), 0)
        expect_true(n > 0 || !perl_available)
})

test_that("R07: ordinary parameter strings still reach the program", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed")
        a <- suppressMessages(diamond(thal_aa(), lyra_aa(), seq_type = "protein"))
        b <- suppressMessages(diamond(thal_aa(), lyra_aa(), seq_type = "protein", diamond_params = "--id 99.5  --max-hsps 1"))
        expect_true(nrow(b) < nrow(a))
        expect_true(all(b$perc_identity >= 99.5))
})
