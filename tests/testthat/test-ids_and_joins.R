# Sequence IDs are carried from the input files to the result, whatever they look like; sequences that are left
# out are reported; searches without hits are reported as such. (Second audit: N02, N04, N13, N16, N20, N21, N22,
# N24, N25, N28, N31, N32.)

q_all <- Biostrings::readDNAStringSet(thal_cds())
s_all <- Biostrings::readDNAStringSet(lyra_cds())

write_renamed <- function(x, ids, file = tempfile(fileext = ".fa")) {
        names(x) <- ids
        Biostrings::writeXStringSet(x, file)
        file
}

run_dnds <- function(...) suppressMessages(dNdS(..., print_citation = FALSE))

# ---------------------------------------------------------------- N04: IDs that look like numbers

test_that("read_gestimator_table() reads IDs as text ...", {
        f <- tempfile()
        writeLines(c("1001\t2001\t0.1\t0.2\t0.5", "007\t1e5\t0\t0.3\t0", "NA\tTRUE\t999\t999\t999", "1.10\tx|y\t0.01\t0\t999"), f)
        x <- read_gestimator_table(f)
        expect_identical(x$query_id, c("1001", "007", "NA", "1.10"))
        expect_identical(x$subject_id, c("2001", "1e5", "TRUE", "x|y"))
        expect_type(x$dN, "double")
        expect_equal(x$dNdS, c(0.5, 0, 999, 999))
})

test_that("dNdS() works with purely numeric IDs and with IDs a type guesser rewrites (N04) ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        ref <- run_dnds(write_renamed(q_all[1:6], paste0("g", 1:6)), write_renamed(s_all, paste0("s", 1:20)))
        expect_equal(nrow(ref), 6)
        
        num <- run_dnds(write_renamed(q_all[1:6], as.character(1001:1006)), write_renamed(s_all, as.character(2001:2020)))
        expect_type(num$query_id, "character")
        expect_identical(num$query_id, as.character(1001:1006))
        expect_identical(num$subject_id, sub("^s", "20", sub("^s([0-9])$", "s0\\1", ref$subject_id)))
        expect_equal(num$dN, ref$dN)
        expect_equal(num$dS, ref$dS)
        
        odd_ids <- c("1e5", "TRUE", "NA", "007", "1.10", "geneX")
        odd <- run_dnds(write_renamed(q_all[1:6], odd_ids), write_renamed(s_all, paste0("s", 1:20)))
        expect_setequal(odd$query_id, odd_ids)
        expect_equal(odd$dN[match(odd_ids, odd$query_id)], ref$dN)
        expect_equal(odd$subject_id[match(odd_ids, odd$query_id)], ref$subject_id)
})

# ---------------------------------------------------------------- N02: long IDs and MAFFT

test_that("multi_aln(tool = 'mafft') keeps IDs longer than 15 characters (N02) ...", {
        skip_if_not(mafft_available, "MAFFT is not installed")
        aa <- Biostrings::readAAStringSet(thal_aa())[1:2]
        ids <- c("Athaliana_TAIR10_AT1G01010_transcript_1", "ENSMUST00000012345.7")
        f <- tempfile(fileext = ".fa")
        names(aa) <- ids
        Biostrings::writeXStringSet(aa, f)
        out <- multi_aln(f, tool = "mafft", get_aln = FALSE, multi_aln_name = "longids", quiet = TRUE)
        expect_identical(names(read_clustal(out)), ids)
})

test_that("dNdS() returns the input IDs for every protein aligner, also when they are long (N02) ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        qid <- paste0("Athaliana_TAIR10_", sub("[.|].*", "", names(q_all)[1:8]))     # 26 characters
        sid <- paste0("Alyrata_v2_gene_", sub("[|].*", "", names(s_all)))            # 22 characters
        qf <- write_renamed(q_all[1:8], qid)
        sf <- write_renamed(s_all, sid)
        ref <- run_dnds(qf, sf)
        expect_setequal(ref$query_id, qid)
        expect_true(all(ref$subject_id %in% sid))
        for (tool in c("mafft", "clustalo")) {
                if (!has_program(tool))
                        next
                x <- run_dnds(qf, sf, aa_aln_type = "multiple", aa_aln_tool = tool)
                expect_equal(nrow(x), 8, info = tool)
                expect_setequal(x$query_id, qid)
                expect_identical(x$subject_id[match(qid, x$query_id)], ref$subject_id[match(qid, ref$query_id)], info = tool)
                expect_false(anyNA(x$dN), info = tool)
        }
})

test_that("the IDs of a gene pair come from the pair table, not from the programs ...", {
        skip_if_not(perl_available, "perl is not installed")
        tbl <- data.table::data.table(
                query_id    = "a-very-long-query-identifier|with:pipes.1",
                subject_id  = "000123",
                query_cds   = tolower("ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACC"),
                subject_cds = tolower("ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTGAAGGCGTGGTTGACC"),
                query_aa    = "MKTAYIAKQRLEGDLKAWLT",
                subject_aa  = "MKTAYIAKQRLDGDLKAWLT"
        )
        for (m in c("Comeron", "Li")) {
                x <- suppressMessages(compute_dnds(data.table::copy(tbl), dnds_est.method = m, quiet = TRUE))
                expect_identical(x$query_id, tbl$query_id, info = m)
                expect_identical(x$subject_id, tbl$subject_id, info = m)
                expect_true(is.numeric(x$dN) && !is.na(x$dN), info = m)
        }
})

# ---------------------------------------------------------------- N28: duplicated IDs

test_that("duplicated sequence IDs are an error that names them (N28) ...", {
        f <- write_renamed(q_all[1:6], c("dup", "dup", "g3", "g4", "g4", "g6"))
        expect_error(suppressMessages(cds2aa(read.cds(f, format = "fasta"))), "2 sequence IDs occur more than once.*'dup', 'g4'")
        skip_if_not(diamond_available, "DIAMOND is not installed")
        expect_error(suppressMessages(set_diamond(f)), "more than once among the coding sequences of '.*': 'dup', 'g4'")
        expect_error(run_dnds(f, lyra_cds()), "must be unique within a file")
        aa <- Biostrings::readAAStringSet(thal_aa())[1:3]
        names(aa) <- c("p1", "p1", "p2")
        fa <- tempfile(fileext = ".faa")
        Biostrings::writeXStringSet(aa, fa)
        expect_error(suppressMessages(set_diamond(fa, seq_type = "protein")), "1 sequence ID occurs more than once among the protein sequences")
})

# ---------------------------------------------------------------- N25: sequences that are left out

test_that("sequences that are left out are reported with number and reason and returned (N25) ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        x <- as.character(q_all[1:6])
        substr(x[2], 200, 200) <- "N"                    # one ambiguous base
        substr(x[3], 1, 3) <- "NNN"                      # 5'-incomplete CDS padded with N
        substr(x[4], 10, 10) <- "R"                      # IUPAC code
        x[5] <- substr(x[5], 2, nchar(x[5]))             # starts in phase 1: length not a multiple of 3
        f <- tempfile(fileext = ".fa")
        writeLines(as.vector(rbind(paste0(">", paste0("g", 1:6)), x)), f)
        
        msg <- testthat::capture_messages(res <- dNdS(f, lyra_cds(), delete_corrupt_cds = TRUE, print_citation = FALSE))
        expect_match(paste(msg, collapse = ""), "3 of 5 coding sequences of '.*' have characters other than A, C, G, T")
        expect_match(paste(msg, collapse = ""), "2 of them with N as the only such character")
        expect_match(paste(msg, collapse = ""), "4 input sequences were left out before the search")
        d <- attr(res, "dropped_sequences")
        expect_s3_class(d, "data.frame")
        expect_identical(names(d), c("file", "id", "reason"))
        expect_setequal(d$id, c("g2", "g3", "g4", "g5"))
        expect_match(d$reason[d$id == "g5"], "not a multiple of 3")
        expect_match(d$reason[d$id == "g2"], "other than A, C, G, T")
        expect_setequal(res$query_id, c("g1", "g6"))
        
        # nothing left out: an empty table
        ok <- run_dnds(thal_cds(), lyra_cds())
        expect_equal(nrow(attr(ok, "dropped_sequences")), 0)
        
        # delete_corrupt_cds = FALSE: the phase-1 sequence is trimmed at its 3' end (and kept); the message says
        # what that means for 5'-incomplete sequences
        msg2 <- testthat::capture_messages(res2 <- dNdS(f, lyra_cds(), delete_corrupt_cds = FALSE, print_citation = FALSE))
        expect_match(paste(msg2, collapse = ""), "trimmed at their 3' end")
        expect_setequal(attr(res2, "dropped_sequences")$id, c("g2", "g3", "g4"))
})

# ---------------------------------------------------------------- N24: several transcripts per gene

test_that("isoform_map() finds genes in Ensembl and NCBI headers and in ID patterns, conservatively (N24) ...", {
        wr <- function(headers) {
                f <- tempfile(fileext = ".fa")
                writeLines(as.vector(rbind(paste0(">", headers), "ATGAAATAA")), f)
                f
        }
        ens <- wr(c("ENST01.2 cds chromosome:GRCh38:21:1:9:1 gene:ENSG01.5 gene_biotype:protein_coding",
                    "ENST02.1 cds chromosome:GRCh38:21:1:9:1 gene:ENSG01.5 gene_biotype:protein_coding",
                    "ENST03.1 cds chromosome:GRCh38:21:1:9:1 gene:ENSG02.1 gene_biotype:protein_coding"))
        m <- isoform_map(ens)
        expect_identical(unname(m[c("ENST01.2", "ENST02.1", "ENST03.1")]), c("ENSG01", "ENSG01", "ENSG02"))
        expect_identical(attr(m, "source"), "header")
        expect_warning(isoform_warning(ens, "query_file"), "3 sequences for 2 genes, 1 genes with more than one sequence")
        
        ncbi <- wr(c("lcl|NC_1_cds_NP_1.1_1 [gene=thrL] [locus_tag=b0001] [protein_id=NP_1.1]",
                     "lcl|NC_1_cds_NP_2.1_2 [locus_tag=b0002] [protein_id=NP_2.1]",
                     "lcl|NC_1_cds_NP_3.1_3 [gene=thrB] [locus_tag=b0003] [protein_id=NP_3.1]"))
        expect_identical(as.vector(isoform_map(ncbi)), c("b0001", "b0002", "b0003"))
        expect_silent(isoform_warning(ncbi))
        
        # ID patterns: 6 of 30 genes with two transcripts are reported, 4 are not (fewer than 5)
        ids <- paste0("AT1G", sprintf("%05d", 1:30), ".1")
        expect_warning(isoform_warning(wr(c(ids, sub("1$", "2", ids[1:6])))), "judged by the sequence IDs")
        expect_silent(isoform_warning(wr(c(ids, sub("1$", "2", ids[1:4])))))
        # versioned IDs and the example files: nothing
        expect_null(isoform_map(wr(paste0("ENST", sprintf("%011d", 1:40), ".", 1:4))))
        expect_silent(isoform_warning(thal_cds()))
        expect_silent(isoform_warning(lyra_cds()))
        expect_silent(isoform_warning(system.file("seqs/ortho_thal_cds_1000.fasta", package = "orthologr")))
        # switched off
        old <- options(orthologr.isoform_warning = FALSE)
        expect_silent(isoform_warning(ens))
        options(old)
})

test_that("dNdS() warns about several transcripts per gene and says how many rows per gene it returns (N24) ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        ids <- sub("[|].*", "", names(q_all))
        x <- c(as.character(q_all), as.character(Biostrings::subseq(q_all[1:6], 1, 3 * (Biostrings::width(q_all[1:6]) %/% 6))))
        f <- tempfile(fileext = ".fa")
        writeLines(as.vector(rbind(paste0(">", c(ids, sub("\\.[0-9]+$", ".9", ids[1:6]))), x)), f)
        msg <- NULL
        expect_warning(
                msg <- testthat::capture_messages(res <- dNdS(f, lyra_cds(), ortho_detection = "BH", print_citation = FALSE)),
                "query_file '.*' seems to hold several transcripts \\(splice variants\\) per gene: 26 sequences for 20 genes, 6 genes"
        )
        expect_match(paste(msg, collapse = ""), "The result holds 26 rows for 20 query genes: 6 genes have more than one row")
        expect_equal(nrow(res), 26)
})

# ---------------------------------------------------------------- N16: a search without hits

test_that("a search without hits returns an empty hit table and dNdS() says 'no hits' (N16) ...", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        set.seed(1)
        rnd <- vapply(1:3, function(i) paste(sample(c("GCT", "GAA", "AAA", "CTG", "GGT", "TCT"), 150, TRUE), collapse = ""), "")
        rf <- tempfile(fileext = ".fa")
        writeLines(as.vector(rbind(paste0(">r", 1:3), rnd)), rf)
        wf <- tempfile(fileext = ".fa")
        writeLines(c(">z1", paste(rep("TGGTGTCATATG", 40), collapse = "")), wf)
        
        expect_message(h <- diamond(wf, rf), "DIAMOND found no hits")
        expect_s3_class(h, "tbl_df")
        expect_equal(nrow(h), 0)
        expect_identical(names(h)[c(1, 2, 19, 20)], c("query_id", "subject_id", "evalue", "bit_score"))
        expect_type(h$query_id, "character")
        expect_type(h$alig_length, "integer")
        expect_equal(nrow(suppressMessages(diamond_best(wf, rf))), 0)
        expect_equal(nrow(suppressMessages(diamond_rec(wf, rf))), 0)
        expect_error(run_dnds(wf, rf), "No orthologs could be found: the DIAMOND search of .* returned no reciprocal best hits at e-value 1E-5")
})

# ---------------------------------------------------------------- N20, N21: KaKs_Calculator arguments

test_that("kaks_calc.params is added to the method, it does not replace it (N20) ...", {
        expect_identical(unlist(lapply(kaks_calc_args("YN", NULL), as.character)), c("-m", "YN"))
        expect_identical(unlist(lapply(kaks_calc_args("YN", "-c 1"), as.character)), c("-m", "YN", "-c 1"))
        expect_identical(unlist(lapply(kaks_calc_args("kaks_calc", "-m NG -m YN"), as.character)), "-m NG -m YN")
        expect_error(kaks_calc_args("kaks_calc", NULL), "needs the methods in 'kaks_calc.params'")
        
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        aln <- system.file("seqs/pal2nal.aln", package = "orthologr")
        a <- substitutionrate(aln, est.method = "YN", quiet = TRUE)
        b <- substitutionrate(aln, est.method = "YN", kaks_calc.params = "-c 1", quiet = TRUE)
        expect_identical(as.character(b$method), "YN")
        expect_equal(b$dN, a$dN)
        expect_equal(b$dS, a$dS)
        two <- substitutionrate(aln, est.method = "kaks_calc", kaks_calc.params = "-m NG -m YN", quiet = TRUE)
        expect_setequal(as.character(two$method), c("NG", "YN"))
})

test_that("kaks_calc_path is used when KaKs_Calculator is not on the PATH (N21) ...", {
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        skip_on_os("windows")
        d <- new_test_dir("kaks only here ")
        file.copy(Sys.which("KaKs_Calculator"), file.path(d, "KaKs_Calculator"))
        aln <- system.file("seqs/pal2nal.aln", package = "orthologr")
        ref <- substitutionrate(aln, est.method = "NG", quiet = TRUE)
        old <- Sys.getenv("PATH")
        on.exit({ Sys.setenv(PATH = old); tool_cache_reset() }, add = TRUE)
        Sys.setenv(PATH = "/usr/bin:/bin")
        tool_cache_reset()
        skip_if(nzchar(Sys.which("KaKs_Calculator")), "KaKs_Calculator is in /usr/bin")
        expect_error(substitutionrate(aln, est.method = "NG", quiet = TRUE), "KaKs_Calculator")
        x <- substitutionrate(aln, est.method = "NG", kaks_calc_path = d, quiet = TRUE)
        expect_equal(x$dN, ref$dN)
})

# ---------------------------------------------------------------- N22: the message a new user sees first

test_that("the 'DIAMOND not found' message names diamond_install() and real arguments (N22) ...", {
        d <- new_test_dir("no_diamond_")
        tool_cache_reset()
        on.exit(tool_cache_reset(), add = TRUE)
        err <- tryCatch(diamond_version(path = d), error = function(e) conditionMessage(e))
        expect_match(err, "diamond_install()", fixed = TRUE)
        expect_match(err, "'aligner_path'", fixed = TRUE)
        expect_match(err, "'path'", fixed = TRUE)
        expect_false(grepl("diamond_exec_path", err, fixed = TRUE))
        for (a in c("path", "aligner_path"))
                expect_true(a %in% c(names(formals(diamond)), names(formals(dNdS))))
})

# ---------------------------------------------------------------- N13: progress bar with one subject file

test_that("map_generator_dnds() works with one subject file and the default progress bar (N13) ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        subj <- new_test_dir("subj_")
        file.copy(lyra_cds(), file.path(subj, "Alyrata.fa"))
        out <- file.path(new_test_dir("maps_"), "out")
        utils::capture.output(
                res <- suppressMessages(map_generator_dnds(thal_cds(), subj, out, progress_bar = TRUE))
        )
        expect_length(list.files(out, pattern = "^map_q=.*Alyrata.*csv$"), 1)
})

# ---------------------------------------------------------------- N31, N32

test_that("the ID of a sequence ends at the first white space, also a tab (N31) ...", {
        expect_identical(fasta_header_id(c("g1\tdescription here", "g2 more text", "g3", "")), c("g1", "g2", "g3", ""))
        f <- tempfile(fileext = ".fa")
        writeLines(c(">g1\tdescription here", "ATGGCTAAATTGTAA", ">g2 x", "ATGGCTAAATTGTAA"), f)
        expect_identical(suppressMessages(read.cds(f, "fasta"))$geneids, c("g1", "g2"))
        expect_identical(read.genome(f, "fasta")$geneids, c("g1", "g2"))
        fa <- tempfile(fileext = ".faa")
        writeLines(c(">p1\tdescription", "MAKL"), fa)
        expect_identical(read.proteome(fa, "fasta")$geneids, "p1")
})

test_that("orthologs(seq_type = 'cds') searches the input files and leaves no '_translated' copies (N32) ...", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        unlink(list.files(tempdir(), pattern = "_translated\\.fasta$", full.names = TRUE))
        a <- suppressMessages(orthologs(thal_cds(), lyra_cds(), seq_type = "cds", ortho_detection = "DIAMOND_RBH"))
        expect_length(list.files(tempdir(), pattern = "_translated\\.fasta$"), 0)
        b <- suppressMessages(diamond_rec(thal_cds(), lyra_cds(), seq_type = "cds"))
        expect_equal(as.data.frame(a), as.data.frame(b), ignore_attr = TRUE)
        expect_equal(nrow(a), 20)
})
