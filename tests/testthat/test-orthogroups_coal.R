# dNdS_orthogroups(engine = "coal"): the gene pairs of the groups computed by COAL (R package 'coalr')

test_that("dNdS_orthogroups(engine = 'coal') checks its arguments", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        og <- suppressMessages(orthogroups(c(thal = thal_cds(), lyra = lyra_cds())))
        old <- options(orthologr.coal_backend = "none")
        expect_error(dNdS_orthogroups(og, engine = "coal"), "COAL is not available")
        options(old)
        skip_if_not_installed("coalr")
        expect_error(dNdS_orthogroups(og, engine = "coal", dnds_est.method = "GY"), "only available with engine = \"orthologr\"")
        if (!coal_has_native_methods())
                expect_error(dNdS_orthogroups(og, engine = "coal", dnds_est.method = "Comeron"), "needs a COAL that computes")
        expect_error(dNdS_orthogroups(og, aa_aln_tool = "pairwise", aa_aln_scoring = "needle"), "needs engine = \"coal\"")
        expect_error(dNdS_orthogroups(og, engine = "coal", aa_aln_tool = "muscle"), "aa_aln_tool = 'muscle' is not supported")
})

test_that("dNdS_orthogroups(engine = 'coal') returns COAL's values for the between-species pairs of every group", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        skip_if_not_installed("coalr")
        toy <- make_toy_species()
        og <- suppressMessages(orthogroups(toy$files, inparalogs = TRUE))
        fam <- function(k) unique(og$groups$orthogroup[is.element(og$groups$gene, toy$truth$gene[toy$truth$family == k])])
        ids <- c(fam(1), fam(11), fam(12), fam(13))
        
        # pairwise: the pairs handed to coalr::coal_dnds_pairs()
        res <- suppressMessages(dNdS_orthogroups(og, groups = ids, engine = "coal", aa_aln_tool = "pairwise"))
        expect_identical(names(res), c("pairs", "summary", "failed"))
        expect_identical(names(res$pairs), c("orthogroup", "query_species", "query_id", "subject_species", "subject_id",
                                             "method", "dN", "dS", "dNdS", "note", "status"))
        p <- res$pairs
        s <- res$summary[match(ids, res$summary$orthogroup), ]
        expect_identical(s$n_pairs, c(6L, 9L, 3L, 1L))
        expect_true(all(p$method == "YN"))
        expect_true(all(match(p$query_species, og$species) < match(p$subject_species, og$species)))
        expect_false(anyNA(p$dN))
        expect_true(all(is.na(p$note)))
        cds <- unlist(unname(lapply(toy$files, orthologr:::read_fasta_vector)))
        ref <- coalr::coal_dnds_pairs(query_cds = cds[unique(p$query_id)], subject_cds = cds[unique(p$subject_id)],
                                      pairs = data.frame(query_id = p$query_id, subject_id = p$subject_id), methods = "YN")
        at <- match(paste(p$query_id, p$subject_id), paste(ref$seq1, ref$seq2))
        expect_false(anyNA(at))
        expect_equal(p$dN, ref$dN[at])
        expect_equal(p$dS, ref$dS[at])
        expect_equal(p$dNdS, ref$dNdS[at])
        expect_equal(s$dN_median[1], stats::median(p$dN[p$orthogroup == ids[1]]))
        
        # several methods: one row per pair and method, one summary row per group and method
        two <- suppressMessages(dNdS_orthogroups(og, groups = ids[3:4], engine = "coal", aa_aln_tool = "pairwise",
                                                 dnds_est.method = c("NG", "YN"), comp_cores = 1))
        expect_identical(nrow(two$pairs), 8L)
        expect_identical(sort(unique(two$pairs$method)), c("NG", "YN"))
        expect_identical(nrow(two$summary), 4L)
        expect_equal(two$pairs$dN[two$pairs$method == "YN"], p$dN[is.element(p$orthogroup, ids[3:4])])
        
        # one multiple alignment per group (the machinery of dNdS_clusters()): the same pairs; the sequences differ
        # by substitutions only, so the alignments and the values are the same
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        msa <- suppressMessages(dNdS_orthogroups(og, groups = ids, engine = "coal", aa_aln_tool = "clustalo"))
        expect_identical(msa$pairs[ , 1:6], p[ , 1:6])
        expect_equal(msa$pairs$dN, p$dN, tolerance = 1e-8)
        expect_equal(msa$pairs$dS, p$dS, tolerance = 1e-8)
        expect_identical(msa$summary$n_pairs, res$summary$n_pairs)
        # dNdS_clusters() on the same groups gives these values (it also computes the pair of the two copies of C)
        cl <- og$groups[is.element(og$groups$orthogroup, ids), ]
        clusters <- data.frame(representative_id = cl$gene[match(cl$orthogroup, cl$orthogroup)], member_id = cl$gene)
        dc <- suppressMessages(dNdS_clusters(unname(toy$files), clusters = clusters, aa_aln_tool = "clustalo"))
        key <- function(a, b) paste(pmin(a, b), pmax(a, b))
        at <- match(key(msa$pairs$query_id, msa$pairs$subject_id), key(dc$pairs$query_id, dc$pairs$subject_id))
        expect_false(anyNA(at))
        expect_equal(msa$pairs$dN, dc$pairs$dN[at])
        expect_equal(msa$pairs$dNdS, dc$pairs$dNdS[at])
        expect_identical(nrow(dc$pairs), nrow(msa$pairs) + 1L)
})
