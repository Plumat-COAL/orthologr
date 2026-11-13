diamond_available <- isTRUE(tryCatch(is_installed_diamond(), error = function(e) FALSE))
kaks_available <- isTRUE(tryCatch(suppressWarnings(is_installed_kaks_calculator()), error = function(e) FALSE))
perl_available <- nzchar(Sys.which("perl"))
coalr_available <- requireNamespace("coalr", quietly = TRUE)
coal_cpp <- coal_binary()

query_cds   <- system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
subject_cds <- system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")

write_fasta <- function(seqs, ids, file = tempfile(fileext = ".fasta")) {
        writeLines(as.vector(rbind(paste0(">", ids), seqs)), file)
        file
}

cds_a <- "ATGAAAACCGCGTATATTGCGAAACAGCGTCTGGAAGGTGATCTGAAAGCGTGGCTGACCTAA"
cds_b <- "ATGAAGACCGCCTATATTGCGAAACAACGTGATCTGAAGGCGTGGTTGACCTAG"
aa_a  <- "MKTAYIAKQRLEGDLKAWLT"
aa_b  <- "MKTAYIAKQR--DLKAWLT-"
aa_a2 <- "MKTAYIAKQRLEGDLKAWLT"

run_dNdS <- function(...) {
        suppressMessages(dNdS(query_file = query_cds, subject_file = subject_cds, comp_cores = 1, quiet = TRUE,
                              print_citation = FALSE, ...))
}

# every test runs once per COAL backend that is installed: the R package 'coalr' (in process) and the coal-cpp program
backends <- c(if (coalr_available) "coalr", if (!is.null(coal_cpp)) "binary")

with_backend <- function(backend, code) {
        old <- options(orthologr.coal_backend = backend)
        on.exit(options(old), add = TRUE)
        force(code)
}

test_that("COAL is optional: a helpful message when it is not installed ...", {
        old <- options(orthologr.coal_backend = "none")
        on.exit(options(old), add = TRUE)
        expect_false(is_installed_coal())
        expect_error(coal_backend(), "install_github\\(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\"\\)")
        aln <- write_fasta(c(aa_a, aa_b), c("gene-A", "gene-B"))
        nuc <- write_fasta(c(cds_a, cds_b), c("gene-A", "gene-B"))
        expect_error(codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", tool = "coal", quiet = TRUE),
                     "COAL is not available")
        # dNdS() says so before the search starts
        expect_error(run_dNdS(codon_aln_tool = "coal"), "COAL is not available")
        expect_error(run_dNdS(dnds_est.method = "YN", kaks_tool = "coal"), "COAL is not available")
        expect_error(run_dNdS(dnds_est.method = "YN", engine = "coal"), "COAL is not available")
})

test_that("arguments of the COAL options are checked ...", {
        expect_true(is.codon_aln_tool("coal"))
        expect_true(is.dnds_est_method("MLPB"))
        expect_error(run_dNdS(engine = "INVALID"), "engine")
        expect_error(run_dNdS(kaks_tool = "INVALID"), "kaks_tool")
        expect_error(run_dNdS(dnds_est.method = "GY", kaks_tool = "coal"), "needs kaks_tool = \"kaks_calculator\"")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        # Comeron (the default method) and Li need a COAL that has them; other methods are not COAL's at all
        if (!coal_has_native_methods())
                expect_error(run_dNdS(engine = "coal"), "needs a COAL that computes Comeron's and Li's method")
        expect_error(run_dNdS(engine = "coal", dnds_est.method = "MA"), "\"MA\" is only available with engine = \"orthologr\"")
        expect_error(run_dNdS(engine = "coal", dnds_est.method = "YN", aa_aln_scoring = "blast"), "aa_aln_scoring")
        expect_error(run_dNdS(aa_aln_scoring = "needle"), "needs engine = \"coal\"")
        expect_error(run_dNdS(engine = "coal", dnds_est.method = "YN", ortho_detection = "Orthofinder2"),
                     "RBH")
        expect_error(dNdS_clusters(query_cds, dnds_est.method = "GY"), "dnds_est.method")
        expect_error(dNdS_clusters(query_cds, aa_aln_scoring = "x"), "aa_aln_scoring")
        expect_error(dNdS_clusters(query_cds, aa_aln_tool = "clustalw"), "aa_aln_tool")
        expect_error(dNdS_clusters("no_such_file.fasta"), "not found")
})

test_that("codon_aln(tool = \"coal\") gives PAL2NAL's codon alignment, without Perl ...", {
        skip_if(length(backends) == 0, "COAL is not installed")
        aln <- write_fasta(c(aa_a, aa_b), c("gene-A", "gene-B"))
        nuc <- write_fasta(c(cds_a, cds_b), c("gene-A", "gene-B"))
        for (backend in backends) with_backend(backend, {
                codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", tool = "coal",
                          codon_aln_name = paste0("t_", backend), quiet = TRUE)
                out <- file.path(tempdir(), "_alignment", "codon_aln", paste0("t_", backend, "_coal.aln"))
                expect_true(file.exists(out))
                got <- toupper(read_fasta_strings(out))
                expect_equal(names(got), c("gene-A", "gene-B"))
                # the stop codons are not part of the alignment; the gaps of the protein alignment are codon gaps
                expect_equal(unname(got[1]), substr(cds_a, 1, 60))
                expect_equal(unname(got[2]), "ATGAAGACCGCCTATATTGCGAAACAACGT------GATCTGAAGGCGTGGTTGACC---")
                if (perl_available) {
                        codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", tool = "pal2nal",
                                  codon_aln_name = paste0("t_", backend), quiet = TRUE)
                        ref <- toupper(read_fasta_strings(file.path(tempdir(), "_alignment", "codon_aln",
                                                                    paste0("t_", backend, "_pal2nal.aln"))))
                        expect_equal(got, ref)
                }
                # get_aln returns the alignment
                res <- codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", tool = "coal",
                                 codon_aln_name = paste0("t2_", backend), get_aln = TRUE, quiet = TRUE)
                expect_equal(res$nb, 2)
                expect_equal(toupper(res$seq[[1]]), substr(cds_a, 1, 60))
                # a protein that does not fit its coding sequence: an error naming the pair and COAL
                bad <- write_fasta(c(aa_a, "MKTAYIAKQRWWWWWWWWWWWWWWWWLEGDLKAWLT"), c("gene-A", "gene-B"))
                bad_nuc <- write_fasta(c(cds_a, cds_a), c("gene-A", "gene-B"))
                expect_error(codon_aln(file_aln = bad, file_nuc = bad_nuc, format = "fasta", tool = "coal",
                                       codon_aln_name = "bad", quiet = TRUE),
                             "COAL.*gene pair 'gene-A' / 'gene-B' \\(bad\\)")
                expect_error(codon_aln(file_aln = aln, file_nuc = nuc, format = "fasta", tool = "coal",
                                       params = "-nogap", quiet = TRUE), "params")
        })
})

test_that("codon_aln(tool = \"coal\") writes clustal format as PAL2NAL does ...", {
        skip_if(length(backends) == 0, "COAL is not installed")
        aln <- write_fasta(c(aa_a, aa_b), c("gene-A", "gene-B"))
        # a coding sequence longer than one block of 60 columns
        long_a <- paste0(substr(cds_a, 1, 60), substr(cds_a, 4, 60))
        long_b <- paste0(substr(cds_b, 1, 51), substr(cds_b, 4, 51))
        aln2 <- write_fasta(c(paste0(aa_a, substr(aa_a, 2, 20)), paste0(aa_b, substr(aa_b, 2, 20))), c("gene-A", "gene-B"))
        nuc2 <- write_fasta(c(long_a, long_b), c("gene-A", "gene-B"))
        for (backend in backends) with_backend(backend, {
                codon_aln(file_aln = aln2, file_nuc = nuc2, format = "clustal", tool = "coal",
                          codon_aln_name = paste0("t_clustal_", backend), quiet = TRUE)
                got <- readLines(file.path(tempdir(), "_alignment", "codon_aln", paste0("t_clustal_", backend, "_coal.aln")))
                expect_match(got[1], "^CLUSTAL")
                rows <- got[grepl("^gene-", got)]
                expect_equal(length(rows), 4)
                seq_a <- paste(sub("^gene-A\\s+", "", rows[grepl("^gene-A", rows)]), collapse = "")
                expect_equal(toupper(seq_a), long_a)
                if (perl_available) {
                        codon_aln(file_aln = aln2, file_nuc = nuc2, format = "clustal", tool = "pal2nal",
                                  codon_aln_name = paste0("t_clustal_", backend), quiet = TRUE)
                        ref <- readLines(file.path(tempdir(), "_alignment", "codon_aln",
                                                   paste0("t_clustal_", backend, "_pal2nal.aln")))
                        expect_equal(got[nzchar(got)], ref[nzchar(ref)])
                }
        })
})

test_that("substitutionrate(kaks_tool = \"coal\") computes the KaKs_Calculator methods ...", {
        skip_if(length(backends) == 0, "COAL is not installed")
        f <- write_fasta(c(substr(cds_a, 1, 60), "ATGAAGACCGCCTATATTGCGAAACAACGTCTGGATGGTGATCTGAAGGCGTGGTTGACC"),
                         c("AT1G01010-1", "Al-scaffold-0001_1"))
        for (backend in backends) with_backend(backend, {
                for (m in c("YN", "NG", "LWL", "MLWL", "LPB", "MLPB", "MYN")) {
                        res <- substitutionrate(f, est.method = m, quiet = TRUE, subst_name = "t", kaks_tool = "coal")
                        expect_equal(names(res), c("query_id", "subject_id", "dN", "dS", "dNdS", "method"))
                        expect_equal(res$query_id, "AT1G01010-1")
                        expect_equal(res$subject_id, "Al-scaffold-0001_1")
                        expect_equal(as.character(res$method), m)
                        expect_true(is.numeric(res$dN) && is.finite(res$dN) && is.finite(res$dS))
                        expect_equal(res$dNdS, res$dN / res$dS)
                        if (kaks_available && m != "MLPB") {
                                ref <- substitutionrate(f, est.method = m, quiet = TRUE, subst_name = "t")
                                # KaKs_Calculator prints 6 significant digits (and "-nan" for MYN on a pair as short
                                # as this one)
                                if (is.finite(ref$dN) && is.finite(ref$dS)) {
                                        expect_equal(res$dN, ref$dN, tolerance = 1e-5)
                                        expect_equal(res$dS, ref$dS, tolerance = 1e-5)
                                }
                        }
                }
                # Comeron and Li are not affected by kaks_tool
                expect_equal(substitutionrate(f, est.method = "Comeron", quiet = TRUE, subst_name = "t", kaks_tool = "coal")$dN,
                             substitutionrate(f, est.method = "Comeron", quiet = TRUE, subst_name = "t")$dN)
                expect_error(substitutionrate(f, est.method = "GY", quiet = TRUE, kaks_tool = "coal"), "kaks_calculator")
                # sequences that are not a codon alignment: an error naming the pair
                bad <- write_fasta(c("ATGAAAACCG", "ATGAAGACCGCC"), c("x", "y"))
                expect_error(substitutionrate(bad, est.method = "YN", quiet = TRUE, kaks_tool = "coal"),
                             "COAL.*gene pair 'x' / 'y'")
        })
})

test_that("dNdS(codon_aln_tool = \"coal\") gives the values of PAL2NAL ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if_not(perl_available, "perl is not installed")
        skip_if(length(backends) == 0, "COAL is not installed")
        ref <- run_dNdS()
        for (backend in backends) with_backend(backend, {
                res <- run_dNdS(codon_aln_tool = "coal")
                expect_equal(as.data.frame(res), as.data.frame(ref), ignore_attr = TRUE)
                expect_equal(nrow(attr(res, "failed_pairs")), 0)
        })
})

test_that("dNdS(kaks_tool = \"coal\") gives KaKs_Calculator's values without KaKs_Calculator or Perl ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if(length(backends) == 0, "COAL is not installed")
        for (backend in backends) with_backend(backend, {
                res <- run_dNdS(dnds_est.method = "YN", codon_aln_tool = "coal", kaks_tool = "coal")
                expect_gt(nrow(res), 0)
                expect_true(all(c("query_id", "subject_id", "dN", "dS", "dNdS", "method") %in% names(res)))
                expect_true(is.numeric(res$dN) && is.numeric(res$dS) && is.numeric(res$dNdS))
                expect_true(all(res$method == "YN"))
                if (kaks_available && perl_available) {
                        ref <- run_dNdS(dnds_est.method = "YN")
                        expect_equal(res$query_id, ref$query_id)
                        expect_equal(res$subject_id, ref$subject_id)
                        expect_equal(names(res), names(ref))
                        expect_equal(unname(vapply(res, function(x) class(x)[1], "")), unname(vapply(ref, function(x) class(x)[1], "")))
                        both <- !is.na(ref$dN) & !is.na(ref$dS)
                        expect_equal(res$dN[both], ref$dN[both], tolerance = 1e-5)
                        expect_equal(res$dS[both], ref$dS[both], tolerance = 1e-5)
                        # KaKs_Calculator prints values below 1e-6 as NA; COAL returns them
                        expect_true(all(res$dN[is.na(ref$dN)] < 1e-6))
                }
        })
})

test_that("dNdS(engine = \"coal\") returns the table of dNdS(), computed in one batch ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        res <- run_dNdS(dnds_est.method = "YN", engine = "coal")
        direct <- suppressMessages(dNdS_coal(query_file = query_cds, subject_file = subject_cds, dnds_est.method = "YN",
                                             comp_cores = 1, print_citation = FALSE))
        expect_equal(as.data.frame(res), as.data.frame(direct), ignore_attr = TRUE)
        expect_s3_class(res, "tbl_df")
        expect_equal(nrow(attr(res, "failed_pairs")), 0)
        expect_true(all(c("SE_dN", "codons_compared", "status") %in% names(attr(res, "coal"))))
        # the same gene pairs, columns and column types as the pair-by-pair path with the same methods
        ref <- run_dNdS(dnds_est.method = "YN", codon_aln_tool = "coal", kaks_tool = "coal")
        expect_equal(names(res), names(ref))
        expect_equal(unname(vapply(res, function(x) class(x)[1], "")), unname(vapply(ref, function(x) class(x)[1], "")))
        expect_equal(res$query_id, ref$query_id)
        expect_equal(res$subject_id, ref$subject_id)
        expect_equal(as.data.frame(res)[ , -(3:5)], as.data.frame(ref)[ , -(3:5)], ignore_attr = TRUE)
        # the values: the same alignment scoring, but co-optimal alignments can differ in where a gap is placed
        expect_gt(mean(abs(res$dN - ref$dN) < 1e-9), 0.5)
        expect_gt(stats::cor(res$dN, ref$dN), 0.999)
        expect_gt(stats::cor(res$dS, ref$dS), 0.99)
        # two threads and several methods
        two <- run_dNdS(dnds_est.method = c("YN", "NG"), engine = "coal")
        expect_equal(nrow(two), 2 * nrow(res))
        expect_equal(two$dN[two$method == "YN"], res$dN)
        # best hits instead of reciprocal best hits, and needle's scoring through dNdS_coal()
        bh <- suppressMessages(dNdS_coal(query_file = query_cds, subject_file = subject_cds, ortho_detection = "BH",
                                         aa_aln_scoring = "needle", print_citation = FALSE))
        expect_gte(nrow(bh), nrow(res))
        expect_true(all(is.finite(bh$dS)))
})

test_that("dNdS(engine = \"coal\") keeps gene ids with hyphens and files of the same name apart ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        ref <- run_dNdS(dnds_est.method = "YN", engine = "coal")
        x <- readLines(query_cds)
        h <- startsWith(x, ">")
        x[h] <- gsub(".", "-", sub("\\s.*$", "", x[h]), fixed = TRUE)
        da <- tempfile("same_a_")
        db <- tempfile("same_b_")
        dir.create(da)
        dir.create(db)
        writeLines(x, file.path(da, "cds.fa"))
        file.copy(subject_cds, file.path(db, "cds.fa"))
        res <- suppressMessages(dNdS(query_file = file.path(da, "cds.fa"), subject_file = file.path(db, "cds.fa"),
                                     dnds_est.method = "YN", engine = "coal", print_citation = FALSE))
        expect_equal(nrow(res), nrow(ref))
        expect_equal(sort(gsub("-", ".", res$query_id, fixed = TRUE)), sort(ref$query_id))
        m <- match(gsub("-", ".", res$query_id, fixed = TRUE), ref$query_id)
        expect_equal(res$dN, ref$dN[m])
        expect_equal(res$subject_id, ref$subject_id[m])
})

test_that("dNdS_clusters() computes dN/dS for the pairs of deepclust() clusters ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        cl <- suppressMessages(deepclust(c(query_cds, subject_cds), seq_type = "cds"))
        res <- suppressMessages(dNdS_clusters(c(query_cds, subject_cds), clusters = cl, aa_aln_tool = "pairwise"))
        expect_named(res, c("pairs", "summary", "failed"))
        expect_s3_class(res$pairs, "tbl_df")
        expect_true(all(c("cluster_id", "query_id", "subject_id", "method", "dN", "dS", "dNdS", "status") %in% names(res$pairs)))
        expect_equal(sort(unique(res$summary$cluster_id)), sort(unique(cl$representative_id)))
        sizes <- table(cl$representative_id)
        expect_equal(res$summary$n_members, as.integer(sizes[res$summary$cluster_id]))
        expect_equal(sum(res$summary$n_pairs), nrow(res$pairs))
        expect_equal(sum(choose(res$summary$n_used, 2)), nrow(res$pairs))
        expect_true(all(res$pairs$status %in% c("ok", "warning", "saturated", "no_data")))
        # every pair is a pair of members of its cluster
        member_of <- stats::setNames(cl$representative_id, cl$member_id)
        expect_equal(unname(member_of[res$pairs$query_id]), res$pairs$cluster_id)
        expect_equal(unname(member_of[res$pairs$subject_id]), res$pairs$cluster_id)
        # the summary is the summary of the pairs
        one <- res$summary[res$summary$n_pairs > 0, ][1, ]
        p <- res$pairs[res$pairs$cluster_id == one$cluster_id, ]
        expect_equal(one$dN_median, stats::median(p$dN, na.rm = TRUE))
        expect_equal(one$n_dN, sum(!is.na(p$dN)))
        expect_equal(one$dNdS_mean, mean(p$dNdS[!is.na(p$dNdS)]))
        # clusters = NULL runs deepclust() itself; the realign table (seven columns) is accepted too
        auto <- suppressMessages(dNdS_clusters(c(query_cds, subject_cds), aa_aln_tool = "pairwise"))
        expect_equal(auto$pairs, res$pairs)
        re <- suppressMessages(deepclust_realign(c(query_cds, subject_cds), clusters = cl, seq_type = "cds"))
        from_re <- suppressMessages(dNdS_clusters(c(query_cds, subject_cds), clusters = re, aa_aln_tool = "pairwise",
                                                  dnds_est.method = c("YN", "NG")))
        expect_equal(nrow(from_re$pairs), 2 * nrow(res$pairs))
        # the size cap
        capped <- suppressMessages(dNdS_clusters(c(query_cds, subject_cds), clusters = cl, aa_aln_tool = "pairwise",
                                                 max_cluster_size = 2, large_clusters = "skip"))
        expect_true(all(capped$summary$n_pairs[capped$summary$n_used > 2] == 0))
})

test_that("dNdS_clusters() aligns larger clusters with Clustal Omega ...", {
        skip_if_not(diamond_available, "DIAMOND2 is not installed or not on PATH")
        skip_if_not(coalr_available, "the R package coalr is not installed")
        skip_if(!nzchar(Sys.which("clustalo")), "clustalo is not on the PATH")
        # three clusters built by hand from the example genes: each query gene with its ortholog and one more gene
        hits <- suppressMessages(diamond_rec(query_cds, subject_cds))
        q_ids <- hits$query_id[1:6]
        cl <- data.frame(representative_id = rep(q_ids[c(1, 3, 5)], each = 3),
                         member_id = c(q_ids[1], hits$subject_id[1], q_ids[2],
                                       q_ids[3], hits$subject_id[3], q_ids[4],
                                       q_ids[5], hits$subject_id[5], q_ids[6]),
                         stringsAsFactors = FALSE)
        res <- suppressMessages(dNdS_clusters(c(query_cds, subject_cds), clusters = cl, aa_aln_tool = "clustalo"))
        expect_equal(nrow(res$pairs), 9)
        expect_true(all(res$summary$aligner == "clustalo"))
        expect_equal(res$summary$n_pairs, c(3L, 3L, 3L))
        expect_true(all(res$pairs$status != "error"))
        expect_error(dNdS_clusters(c(query_cds, subject_cds), clusters = cl, aa_aln_tool = "clustalo",
                                   aa_aln_path = file.path(tempdir(), "no_such_clustalo")),
                     "was not found")
})
