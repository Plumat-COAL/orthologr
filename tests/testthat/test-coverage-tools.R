# Coverage of code paths no other test reaches: the run path of orthofinder2() (with a stand-in for OrthoFinder,
# a POSIX shell script in fixtures/coverage; skipped on Windows), gestimator()'s verbose output, its output file and
# its removal of gapped columns, the "which methods are there" helpers, and argument checks and error paths of
# several file readers. No external program is needed.

# A copy of the fake OrthoFinder in a folder of its own (a fresh path per test: orthologr remembers per path
# whether a program is there). Returns the folder (of_path).
cov_fake_orthofinder <- function() {
        testthat::skip_on_os("windows")
        testthat::skip_if_not(nzchar(Sys.which("sh")), "no POSIX shell")
        dir <- new_test_dir("fake_of_")
        exe <- file.path(dir, "orthofinder")
        file.copy(testthat::test_path("fixtures", "coverage", "fake-orthofinder.sh"), exe)
        Sys.chmod(exe, "0755")
        dir
}

cov_proteome_folder <- function() {
        d <- new_test_dir("cov_of_proteomes_")
        writeLines(c(">a1", "MKTAYIAKQR"), file.path(d, "A.faa"))
        writeLines(c(">b1", "MKTAYIAKQR"), file.path(d, "B.faa"))
        d
}

test_that("orthofinder2() runs OrthoFinder and returns its results folder", {
        of <- cov_fake_orthofinder()
        old <- Sys.getenv("FAKE_OF_MODE", unset = NA)
        on.exit(if (is.na(old)) Sys.unsetenv("FAKE_OF_MODE") else Sys.setenv(FAKE_OF_MODE = old), add = TRUE)
        Sys.setenv(FAKE_OF_MODE = "ok")
        expect_true(is_installed_orthofinder(path = of))

        d <- cov_proteome_folder()
        res <- suppressWarnings(suppressMessages(orthofinder2(d, of_path = of, comp_cores = 1)))
        expected <- file.path(normalizePath(d), "OrthoFinder", paste0("Results_", basename(d)))
        expect_identical(res, expected)
        args <- readLines(file.path(res, "args.txt"))
        expect_match(args, paste("-t 1 -a 1 -S diamond -n", basename(d), "-og"), fixed = TRUE)
        # a second run would overwrite the results: refused
        expect_error(suppressWarnings(suppressMessages(orthofinder2(d, of_path = of))), "exists already")

        # without -og
        d2 <- cov_proteome_folder()
        res2 <- suppressWarnings(suppressMessages(orthofinder2(d2, of_path = of, orthogroup_only = FALSE)))
        expect_false(grepl("-og", readLines(file.path(res2, "args.txt")), fixed = TRUE))

        # OrthoFinder fails, or ends without its results folder
        Sys.setenv(FAKE_OF_MODE = "fail")
        expect_error(suppressWarnings(suppressMessages(orthofinder2(cov_proteome_folder(), of_path = of))),
                     "failed on purpose", class = "orthologr_error_tool_failed")
        Sys.setenv(FAKE_OF_MODE = "noresults")
        expect_error(suppressWarnings(suppressMessages(orthofinder2(cov_proteome_folder(), of_path = of))),
                     "results folder .* was not found")
})

# ---------------------------------------------------------------------------------------------------------------

cov_gest_file <- function(seqs) {
        f <- tempfile("cov_gest_", fileext = ".fasta")
        writeLines(as.vector(rbind(paste0(">", names(seqs)), seqs)), f)
        f
}

cov_gest_seqs <- c(
        s1 = "ATGAAAACCGCGTATATTGCG---CAGCGTCTGGAAGGTGATCTG",
        s2 = "ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTG",
        s3 = "ATGAAGACAGCCTACATTGCGAAACAACGTCTGGATGGCGACCTG"
)

test_that("gestimator() writes the verbose table to a file, one row per pair", {
        f <- cov_gest_file(cov_gest_seqs)
        out <- tempfile("cov_gest_out_")
        gestimator(f, file_out = out, verbose = TRUE)
        x <- utils::read.delim(out, stringsAsFactors = FALSE)
        expect_identical(names(x), c("seq1", "seq2", "Ka", "Ks", "omega", "p0", "p2S", "p2V", "p4", "q0", "q2S", "q2V",
                                     "q4", "Aa", "As", "Ba", "Bs", "L0", "L2S", "L2V", "L4"))
        expect_identical(paste(x$seq1, x$seq2), c("s1 s2", "s1 s3", "s2 s3"))
        expect_true(all(x$L0 > 0 & x$L4 > 0))
        # the counts are those the plain table is computed from
        plain <- tempfile("cov_gest_plain_")
        gestimator(f, file_out = plain)
        y <- utils::read.delim(plain, header = FALSE, stringsAsFactors = FALSE)
        expect_equal(y[[3]], x$Ka)
        expect_equal(y[[4]], x$Ks)
})

test_that("gestimator() prints to the console without an output file (plain and verbose)", {
        f <- cov_gest_file(cov_gest_seqs[1:2])
        plain <- capture.output(gestimator(f))
        expect_length(plain, 1)
        expect_identical(strsplit(plain, "\t")[[1]][1:2], c("s1", "s2"))
        verbose <- capture.output(gestimator(f, verbose = TRUE))
        expect_length(verbose, 2)
        expect_match(verbose[1], "^seq1\tseq2\tKa\tKs\tomega")
        expect_length(strsplit(verbose[2], "\t")[[1]], 21)
})

test_that("gestimator(remove_all_gaps = TRUE) drops every column with a gap before counting", {
        f <- cov_gest_file(cov_gest_seqs)
        with_gaps <- tempfile()
        without <- tempfile()
        gestimator(f, file_out = with_gaps, verbose = TRUE)
        gestimator(f, file_out = without, verbose = TRUE, remove_all_gaps = TRUE)
        a <- utils::read.delim(with_gaps)
        b <- utils::read.delim(without)
        # the gapped codon of s1 is gone from all three pairs: fewer sites for the pair s2 / s3 as well
        expect_lt(b$L0[3] + b$L2S[3] + b$L2V[3] + b$L4[3], a$L0[3] + a$L2S[3] + a$L2V[3] + a$L4[3])
        # the ungapped alignment gives the values the pair s2 / s3 has without the gapped column
        cut <- function(s) paste0(substr(s, 1, 21), substr(s, 25, nchar(s)))
        f2 <- cov_gest_file(vapply(cov_gest_seqs[2:3], cut, ""))
        c2 <- tempfile()
        gestimator(f2, file_out = c2, verbose = TRUE)
        expect_equal(utils::read.delim(c2)$Ka, b$Ka[3])
        # an alignment without gaps is left as it is
        f3 <- cov_gest_file(cov_gest_seqs[2:3])
        o3a <- tempfile(); o3b <- tempfile()
        gestimator(f3, file_out = o3a)
        gestimator(f3, file_out = o3b, remove_all_gaps = TRUE)
        expect_identical(readLines(o3a), readLines(o3b))
})

test_that("gestimator(): a pair without synonymous sites (Met and Trp codons only)", {
        f <- cov_gest_file(c(a = "ATGTGGATGTGG", b = "ATGTGGATGTGG"))
        out <- tempfile()
        gestimator(f, file_out = out, verbose = TRUE)
        x <- utils::read.delim(out)
        expect_equal(x$L2V + x$L4, 0)
        expect_equal(x$Ka, 0)
        expect_equal(x$Bs, 0)
})

test_that("gestimator() checks maxHits and the input file", {
        f <- cov_gest_file(cov_gest_seqs[1:2])
        expect_error(gestimator(f, maxHits = 4), "maxHits must be 1, 2 or 3")
        expect_error(gestimator(file.path(tempdir(), "no-such-alignment.fasta")), "Input File not found")
})

# ---------------------------------------------------------------------------------------------------------------

test_that("the is.*() helpers list what is provided when asked without an argument", {
        helpers <- list(is.ortho_detection_method, is.multiple_aln_tool, is.pairwise_aln_tool, is.codon_aln_tool,
                        is.dnds_est_method, is.alignment_search_tool)
        for (h in helpers) {
                expect_message(r <- h(), "The following methods are provided: ")
                expect_false(r)
        }
        expect_true(is.alignment_search_tool("ggsearch"))
        expect_false(is.alignment_search_tool("blastp"))
        expect_true(is.codon_aln_tool("coal"))
        expect_true(is.dnds_est_method("YN"))
        expect_output(test("x"), "Test x passed.")
})

test_that("read_clustal() refuses files that are not CLUSTAL alignments", {
        f <- tempfile(fileext = ".aln")
        writeLines(c(">a", "ACGT"), f)
        expect_error(read_clustal(f), "is not a file in CLUSTAL format")
        writeLines(character(0), f)
        expect_error(read_clustal(f), "is not a file in CLUSTAL format")
        writeLines(c("CLUSTAL W", "", "lonelyname"), f)
        expect_error(read_clustal(f), "is not a file in CLUSTAL format")
        writeLines(c("CLUSTAL W", "", "a  ACGT", "b  AC"), f)
        expect_error(read_clustal(f), "do not have the same length")
        writeLines(c("CLUSTAL W", "", "a  ACGT", "b  AC-T", "   ** *", "", "a  AA", "b  AA"), f)
        expect_identical(read_clustal(f), c(a = "ACGTAA", b = "AC-TAA"))
})

test_that("readers of ortholog tables check their input", {
        expect_error(read.dnds.tbl(file.path(tempdir(), "no-dnds.csv")), "does not seem to exist",
                     class = "orthologr_error_bad_input")
        empty <- tempfile(fileext = ".csv")
        writeLines(paste(c("query_id", "subject_id", "dN", "dS", "dNdS", "perc_identity", "alig_length", "mismatches",
                           "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"), collapse = ";"),
                   empty)
        expect_error(read.dnds.tbl(empty), "seems to be empty")
        expect_error(import_ortholog_tables_all(file.path(tempdir(), "no-ortholog-tables")), "does not seem to exist")
})

test_that("generate_ortholog_tables() checks its files and output_type first", {
        dnds <- tempfile(fileext = ".csv")
        writeLines("x", dnds)
        anno <- tempfile(fileext = ".gtf")
        writeLines("x", anno)
        nope <- file.path(tempdir(), "no-such-file")
        expect_error(generate_ortholog_tables(nope, anno, anno, tempdir(), format = c("gtf", "gtf")), "dNdS_file")
        expect_error(generate_ortholog_tables(dnds, nope, anno, tempdir(), format = c("gtf", "gtf")), "annotation_file_query")
        expect_error(generate_ortholog_tables(dnds, anno, nope, tempdir(), format = c("gtf", "gtf")), "annotation_file_subject")
        expect_error(generate_ortholog_tables(dnds, anno, anno, tempdir(), output_type = "isoform", format = c("gtf", "gtf")),
                     "output_type")
})

test_that("generate_ortholog_tables(output_type = \"splice_variant\") returns the splice-variant table", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        dnds_file <- file.path(x$dir, "Athaliana_vs_Alyrata.csv")
        utils::write.table(x$dnds, dnds_file, sep = ";", quote = FALSE, row.names = FALSE)
        out <- file.path(x$dir, "tables_sv")
        dir.create(out)
        sv <- suppressMessages(suppressWarnings(generate_ortholog_tables(dnds_file, x$query$gtf, x$subject$gtf,
                                                                         output_folder = out, output_type = "splice_variant",
                                                                         format = c("gtf", "gtf"))))
        expect_identical(nrow(sv), 25L)
})

test_that("retrieve_longest_isoforms_all() checks its folders and pairs proteome and annotation files by name", {
        p <- new_test_dir("cov_prot_")
        a <- new_test_dir("cov_anno_")
        expect_error(retrieve_longest_isoforms_all(file.path(p, "nope"), a, tempdir()), "seems not to exist")
        expect_error(retrieve_longest_isoforms_all(p, file.path(a, "nope"), tempdir()), "seems not to exist")
        writeLines(c(">x", "MK"), file.path(p, "Athaliana.faa"))
        writeLines("x", file.path(a, "Alyrata.gff"))
        out <- new_test_dir("cov_iso_out_")
        msgs <- character(0)
        expect_error(withCallingHandlers(retrieve_longest_isoforms_all(p, a, out),
                                         message = function(m) {
                                                 msgs <<- c(msgs, conditionMessage(m))
                                                 invokeRestart("muffleMessage")
                                         }),
                     "seem not to match")
        expect_true(any(grepl("exists already", msgs)))
})

test_that("retrieve_core_orthologs() checks the species and the uniqueness of the core set", {
        tbl <- tibble::tibble(query_id = c("g1", "g1", "g2"), subject_species = c("s1", "s2", "s1"),
                              q_len = 100, alig_length = 100L)
        expect_error(suppressMessages(retrieve_core_orthologs(tbl, c("s1", "s2", "s3"))), "does not match")
        dup <- tibble::tibble(query_id = c("g1", "g1", "g1"), subject_species = c("s1", "s2", "s2"),
                              q_len = 100, alig_length = 100L)
        expect_error(suppressMessages(retrieve_core_orthologs(dup, c("s1", "s2"))), "non-unique query_id")
})

test_that("blast() checks the BLAST mode and the remote database before it starts BLAST", {
        expect_error(blast(thal_aa(), lyra_aa(), blast_algorithm = "blastn"), "Only 'blastp'")
        expect_error(blast(thal_aa(), lyra_aa(), remote = TRUE), "specify the 'db' argument")
        expect_error(blast(thal_aa(), lyra_aa(), db = "swissprot"), "supported by remote BLAST")
})
