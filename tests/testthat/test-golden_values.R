# Golden values: results that are pinned, so that a change of numbers, of ortholog pairs or of IDs cannot go
# unnoticed. (Second audit, N39: no test asserted a single numeric dN/dS value, a known ortholog pair, or an ID
# with '|', digits only or more than 15 characters.)
#
# The codon alignments in fixtures/codon_aln are the PAL2NAL alignments of four gene pairs of the example data
# (inst/seqs/ortho_thal_cds.fasta vs ortho_lyra_cds.fasta), as dNdS() builds them with its defaults
# (see fixtures/README). The expected values were computed with the implementation of this branch and compared
# with an independent program, KaKs_Calculator 2.0 (method MLPB for Comeron and Li, which are methods of the same
# family: Li 1993 / Pamilo & Bianchi 1993, with Comeron's and Tzeng's modifications).

golden_file <- function(f) test_path("fixtures", "codon_aln", f)

golden <- utils::read.delim(test_path("fixtures", "codon_aln", "golden_values.tsv"), stringsAsFactors = FALSE)


# ---------------------------------------------------------------- dN / dS values

test_that("Comeron: dN, dS and dN/dS of fixed codon alignments are the pinned values ...", {
        for (i in seq_len(nrow(golden))) {
                x <- substitutionrate(golden_file(golden$file[i]), est.method = "Comeron", quiet = TRUE)
                expect_identical(x$query_id, golden$query_id[i])
                expect_identical(x$subject_id, golden$subject_id[i])
                # gestimator prints four significant digits
                expect_equal(x$dN, golden$comeron_dN[i], tolerance = 1e-12, info = golden$file[i])
                expect_equal(x$dS, golden$comeron_dS[i], tolerance = 1e-12, info = golden$file[i])
                expect_equal(x$dNdS, golden$comeron_dNdS[i], tolerance = 1e-12, info = golden$file[i])
        }
})

test_that("Li: dN, dS and dN/dS of fixed codon alignments are the pinned values ...", {
        for (i in seq_len(nrow(golden))) {
                x <- substitutionrate(golden_file(golden$file[i]), est.method = "Li", quiet = TRUE)
                expect_equal(x$dN, golden$li_dN[i], tolerance = 1e-6, info = golden$file[i])
                expect_equal(x$dS, golden$li_dS[i], tolerance = 1e-6, info = golden$file[i])
                expect_equal(x$dNdS, golden$li_dNdS[i], tolerance = 1e-6, info = golden$file[i])
        }
})

test_that("the pinned Comeron and Li values agree with KaKs_Calculator (MLPB) within the stated tolerance ...", {
        # Stated tolerance: 15% of the value or 0.005, whichever is larger (the methods differ in how they weight
        # transitions and transversions and in the treatment of multiple substitutions per codon).
        close_to <- function(a, b) abs(a - b) <= pmax(0.15 * abs(b), 0.005)
        expect_true(all(close_to(golden$comeron_dN, golden$mlpb_dN)))
        expect_true(all(close_to(golden$comeron_dS, golden$mlpb_dS)))
        expect_true(all(close_to(golden$li_dN, golden$mlpb_dN)))
        expect_true(all(close_to(golden$li_dS, golden$mlpb_dS)))
        # the reference values themselves are pinned from a run of KaKs_Calculator 2.0; when the program is
        # installed they are computed again
        skip_if_not(kaks_available, "KaKs_Calculator is not installed")
        for (i in seq_len(nrow(golden))) {
                k <- substitutionrate(golden_file(golden$file[i]), est.method = "kaks_calc", kaks_calc.params = "-m MLPB", quiet = TRUE)
                expect_equal(k$dN, golden$mlpb_dN[i], tolerance = 1e-4, info = golden$file[i])
                expect_equal(k$dS, golden$mlpb_dS[i], tolerance = 1e-4, info = golden$file[i])
        }
})

test_that("hand-made alignments: identical sequences, one synonymous and one nonsynonymous change ...", {
        f <- tempfile(fileext = ".fa")
        base <- "ATGGCTAAACTGGGTTTTGAACGTCCGATTCAGTCTGTTCATTGGTATAACGACTGC"
        wr <- function(b) { writeLines(c(">a", base, ">b", b), f); f }
        same <- substitutionrate(wr(base), est.method = "Comeron", quiet = TRUE)
        expect_equal(c(same$dN, same$dS), c(0, 0))
        expect_true(is.na(same$dNdS))                           # dS = 0: not defined
        syn <- base; substr(syn, 6, 6) <- "C"                   # GCT -> GCC (Ala -> Ala)
        x <- substitutionrate(wr(syn), est.method = "Comeron", quiet = TRUE)
        expect_equal(x$dN, 0)
        expect_gt(x$dS, 0)
        expect_equal(x$dNdS, 0)
        non <- base; substr(non, 5, 5) <- "A"                   # GCT -> GAT (Ala -> Asp)
        y <- substitutionrate(wr(non), est.method = "Comeron", quiet = TRUE)
        expect_gt(y$dN, 0)
        expect_equal(y$dS, 0)
        expect_true(is.na(y$dNdS))
        li <- substitutionrate(wr(syn), est.method = "Li", quiet = TRUE)
        expect_equal(li$dN, 0)
        expect_gt(li$dS, 0)
})

# ---------------------------------------------------------------- known ortholog pairs, end to end

known_pairs <- utils::read.delim(test_path("fixtures", "codon_aln", "known_pairs.tsv"), stringsAsFactors = FALSE)

test_that("dNdS() on the example data returns the known ortholog pairs with the pinned values ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        # The values are those of the protein alignments of orthologr <= 0.4.3. Branches that change the default
        # alignment scoring (BLOSUM62, options(orthologr.legacy_pairwise_aln = )) reproduce them with this option;
        # where the option does not exist it has no effect.
        opt <- options(orthologr.legacy_pairwise_aln = TRUE)
        on.exit(options(opt), add = TRUE)
        for (rule in c("bitscore", "evalue")) {
                res <- suppressMessages(dNdS(thal_cds(), lyra_cds(), print_citation = FALSE, best_hit_rule = rule))
                expect_equal(nrow(res), 20, info = rule)
                expect_identical(res$query_id, known_pairs$query_id, info = rule)
                expect_identical(res$subject_id, known_pairs$subject_id, info = rule)
                expect_equal(res$dN, known_pairs$dN, tolerance = 1e-12, info = rule)
                expect_equal(res$dS, known_pairs$dS, tolerance = 1e-12, info = rule)
                expect_equal(res$dNdS, known_pairs$dNdS, tolerance = 1e-12, info = rule)
        }
        # best hits alone give the same pairs here (every gene has its ortholog in the other file)
        bh <- suppressMessages(diamond_best(thal_cds(), lyra_cds()))
        expect_identical(bh$subject_id, known_pairs$subject_id)
})

test_that("reciprocal best hits are symmetric: query and subject swapped give the same pairs ...", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        for (rule in c("bitscore", "evalue")) {
                ab <- suppressMessages(diamond_rec(thal_cds(), lyra_cds(), best_hit_rule = rule))
                ba <- suppressMessages(diamond_rec(lyra_cds(), thal_cds(), best_hit_rule = rule))
                expect_setequal(paste(ab$query_id, ab$subject_id), paste(ba$subject_id, ba$query_id))
                expect_setequal(paste(ab$query_id, ab$subject_id), paste(known_pairs$query_id, known_pairs$subject_id))
                # proteins give the pairs of the coding sequences
                aa <- suppressMessages(diamond_rec(thal_aa(), lyra_aa(), seq_type = "protein", best_hit_rule = rule))
                expect_equal(nrow(aa), 20)
        }
        skip_if_not(blast_available, "BLAST+ is not installed")
        ab <- suppressMessages(blast_rec(thal_cds(), lyra_cds()))
        ba <- suppressMessages(blast_rec(lyra_cds(), thal_cds()))
        expect_setequal(paste(ab$query_id, ab$subject_id), paste(ba$subject_id, ba$query_id))
        expect_setequal(paste(ab$query_id, ab$subject_id), paste(known_pairs$query_id, known_pairs$subject_id))
})

test_that("a tie-break fixture read from a hit table file gives the documented best hits ...", {
        h <- read_hit_table(test_path("fixtures", "hits", "tie_break_hits.tsv"),
                            c("query_id", "subject_id", "perc_identity", "num_ident_matches", "alig_length", "mismatches", "gap_openings",
                              "n_gaps", "pos_match", "ppos", "q_start", "q_end", "q_len", "qcov", "qcovhsp", "s_start", "s_end", "s_len",
                              "evalue", "bit_score", "score_raw"), tool = "DIAMOND")
        expected <- utils::read.delim(test_path("fixtures", "hits", "tie_break_expected.tsv"), stringsAsFactors = FALSE, colClasses = "character")
        # (the table is sorted by query_id and subject_id when it is read: 'aaa' is the first row of q5 although
        # 'bbb' comes first in the file; the IDs '007' and '1e5' stay text)
        expect_true(all(c("007", "1e5") %in% c(h$query_id, h$subject_id)))
        new <- best_hits(h, rule = "bitscore")
        old <- best_hits(h, rule = "evalue")
        expect_identical(new$subject_id[match(expected$query_id, new$query_id)], expected$bitscore)
        expect_identical(old$subject_id[match(expected$query_id, old$query_id)], expected$evalue)
        expect_identical(old, best_hits_do(h))
})

# ---------------------------------------------------------------- IDs end to end

test_that("pipe, numeric, long and hyphen IDs come through dNdS() unchanged, with the values of plain IDs ...", {
        skip_if_not(diamond_available && perl_available, "DIAMOND or perl is not installed")
        q <- Biostrings::readDNAStringSet(thal_cds())[1:8]
        s <- Biostrings::readDNAStringSet(lyra_cds())
        qid <- c("sp|P12345|NAC1_ARATH", "gi|15223276|ref|NP_171609.1|", "1001", "007", "AT1G01040-1",
                 "lcl|NC_003070.9_cds_NP_171609.1_1", paste0("Athaliana_TAIR10_chromosome1_gene_", strrep("x", 30), "_1"), "g8.1")
        sid <- c(paste0("ENSMUST", sprintf("%011d", 1:5)), as.character(2001:2005), paste0("Al-", 1:5), paste0("tr|A0A0", 1:5, "|X_ARALY"))
        names(q) <- qid; names(s) <- sid
        qf <- tempfile(fileext = ".fa"); sf <- tempfile(fileext = ".fa")
        Biostrings::writeXStringSet(q, qf); Biostrings::writeXStringSet(s, sf)
        ref <- known_pairs[1:8, ]
        check <- function(res, label, values = TRUE, allow_na = FALSE) {
                expect_equal(nrow(res), 8, info = label)
                expect_setequal(res$query_id, qid)
                res <- res[match(qid, res$query_id), ]
                expect_true(all(res$subject_id %in% sid), info = label)
                # the same partners as with the original IDs (by position in the subject file)
                orig_sid <- sub(" .*", "", names(Biostrings::readDNAStringSet(lyra_cds())))
                expect_identical(match(res$subject_id, sid), match(ref$subject_id, orig_sid), info = label)
                if (values) {
                        expect_equal(res$dN, ref$dN, tolerance = 1e-12, info = label)
                        expect_equal(res$dS, ref$dS, tolerance = 1e-12, info = label)
                }
                # (KaKs_Calculator prints NA for a dN of 0)
                if (!allow_na)
                        expect_false(anyNA(res$dN), info = label)
                expect_false(anyNA(res$dS), info = label)
        }
        opt <- options(orthologr.legacy_pairwise_aln = TRUE)       # see above: the alignments of the pinned values
        check(suppressMessages(dNdS(qf, sf, print_citation = FALSE)), "Comeron / NW")
        options(opt)
        check(suppressMessages(dNdS(qf, sf, dnds_est.method = "Li", print_citation = FALSE)), "Li / NW", values = FALSE)
        if (kaks_available)
                check(suppressMessages(dNdS(qf, sf, dnds_est.method = "YN", print_citation = FALSE)), "YN / NW", values = FALSE, allow_na = TRUE)
        for (tool in c("mafft", "clustalo")) {
                if (has_program(tool))
                        check(suppressMessages(dNdS(qf, sf, aa_aln_type = "multiple", aa_aln_tool = tool, print_citation = FALSE)),
                              tool, values = FALSE)
        }
        ds <- suppressWarnings(suppressMessages(divergence_stratigraphy(qf, sf, quiet = TRUE, subject.id = TRUE, n_quantile = 4)))
        expect_setequal(ds$query_id, qid)
})
