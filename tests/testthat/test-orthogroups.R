# orthogroups(): DIAMOND-based orthogroups of 2..N species (RBH graph, DeepClust route), writer, dN/dS per group

quiet_og <- function(...) suppressMessages(orthogroups(...))

truth_strings <- function(truth, drop = character(0)) {
        t <- truth[!is.element(truth$gene, drop), ]
        g <- split(t$gene, t$family)
        g <- g[vapply(split(t$species, t$family), function(s) length(unique(s)) >= 2, logical(1))]
        sort(unname(vapply(g, function(v) paste(sort(v), collapse = " "), character(1))))
}

check_invariants <- function(og) {
        # every gene in at most one group, every gene either in a group or unassigned
        expect_false(any(duplicated(og$groups$gene)))
        expect_setequal(c(og$groups$gene, og$unassigned$gene), og$genes$gene)
        expect_length(intersect(og$groups$gene, og$unassigned$gene), 0)
        # a group has genes of at least two species
        n_sp <- tapply(og$groups$species, og$groups$orthogroup, function(s) length(unique(s)))
        expect_true(all(n_sp >= 2))
        # wide and long table hold the same
        long <- unlist(lapply(og$species, function(k) unlist(strsplit(og$wide[[k]], ", ", fixed = TRUE))))
        expect_setequal(long, og$groups$gene)
        expect_identical(og$wide$orthogroup, sprintf("OG%07d", seq_len(nrow(og$wide)) - 1L))
        expect_identical(names(og$wide), c("orthogroup", og$species))
        # single-copy groups: exactly one gene of every species
        sc <- og$groups[is.element(og$groups$orthogroup, og$single_copy$orthogroup), ]
        if (nrow(sc) > 0) {
                tab <- table(sc$orthogroup, factor(sc$species, levels = og$species))
                expect_true(all(tab == 1))
        }
        counts <- table(og$groups$orthogroup, factor(og$groups$species, levels = og$species))
        expect_setequal(og$single_copy$orthogroup, rownames(counts)[apply(counts == 1, 1, all)])
        # pairs: both directions, different species, same group, genes of the group
        p <- og$pairs
        expect_true(all(p$query_species != p$subject_species))
        expect_setequal(paste(p$query_id, p$subject_id), paste(p$subject_id, p$query_id))
        expect_false(any(duplicated(paste(p$query_id, p$subject_id, sep = "\r"))))
        grp <- stats::setNames(og$groups$orthogroup, og$groups$gene)
        expect_identical(unname(grp[p$query_id]), p$orthogroup)
        expect_identical(unname(grp[p$subject_id]), p$orthogroup)
        # number of pairs: all pairs of genes of different species within a group
        expected <- sum(apply(counts, 1, function(v) sum(v)^2 - sum(v^2)))
        expect_identical(nrow(p), as.integer(expected))
        expect_true(all(is.element(p$relation, c("1:1", "co-ortholog"))))
        one <- counts[cbind(p$orthogroup, p$query_species)] == 1 & counts[cbind(p$orthogroup, p$subject_species)] == 1
        expect_identical(p$relation == "1:1", as.vector(one))
}

test_that("orthogroups() rejects wrong input before any program is started", {
        f <- c(a = thal_aa(), b = lyra_aa())
        expect_error(orthogroups(thal_aa()), "at least two species")
        expect_error(orthogroups(c(a = thal_aa(), a = lyra_aa())), "species names must be unique")
        expect_error(orthogroups(c(thal_aa(), file.path(tempdir(), "no_such_file.fa"))), "not found")
        expect_error(orthogroups(f, method = "mcl"), "method = 'mcl' is not supported")
        expect_error(orthogroups(f, seq_type = "dna"), "seq_type = 'dna' is not supported")
        expect_error(orthogroups(f, best_hit_rule = "length"), "best_hit_rule")
        expect_error(orthogroups(f, sensitivity_mode = "turbo"), "not a DIAMOND2 sensitivity mode")
        expect_error(orthogroups(f, eval = "1e-5; echo hi"), "not an e-value threshold")
        expect_error(orthogroups(f, min_qcov = 150), "between 0 and 100")
        expect_error(orthogroups(f, inparalogs = NA), "TRUE or FALSE")
        expect_error(orthogroups(f, approx_id = 50), "Unused argument")
        expect_error(orthogroups(f, refine = TRUE), "belong to method = \"deepclust\"")
        expect_error(orthogroups(f, method = "deepclust", inparalogs = TRUE), "refine = TRUE")
        d <- new_test_dir()
        expect_error(orthogroups(d), "holds no FASTA files")
        expect_error(write_orthogroups(list(), d), "result of orthogroups")
        # species names: given names, else file names without extension (and without .gz)
        inp <- orthologr:::orthogroups_input(list(x = thal_aa(), lyra_aa()))
        expect_identical(inp$species, c("x", "ortho_lyra_aa"))
        file.copy(thal_aa(), file.path(d, "Sp_one.v2.faa"))
        file.copy(lyra_aa(), file.path(d, "sp_two.fasta"))
        writeLines("x", file.path(d, "notes.txt"))
        expect_identical(orthologr:::orthogroups_input(d)$species, c("Sp_one.v2", "sp_two"))
})

test_that("connected components and the RBH graph are right on a hand-made hit table (no DIAMOND needed)", {
        cc <- orthologr:::connected_components
        expect_identical(cc(5, integer(0), integer(0)), 1:5)
        expect_identical(cc(6, c(2L, 5L, 3L), c(3L, 6L, 4L)), c(1L, 2L, 2L, 2L, 5L, 5L))
        # a long chain given in an unfavourable order
        n <- 200L
        expect_identical(cc(n, (n - 1L):1L, n:2L), rep(1L, n))
        set.seed(1)
        from <- sample.int(300, 150, replace = TRUE); to <- sample.int(300, 150, replace = TRUE)
        lab <- cc(300, from, to)
        expect_true(all(lab[from] == lab[to]))
        expect_true(all(lab <= seq_len(300)) && all(lab[lab] == lab))
        # as many components as a union-find gives
        parent <- seq_len(300)
        find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
        for (k in seq_along(from)) { a <- find(from[k]); b <- find(to[k]); if (a != b) parent[max(a, b)] <- min(a, b) }
        expect_identical(length(unique(lab)), length(unique(vapply(seq_len(300), find, numeric(1)))))
        
        # genes: species A = a1, a2, a3; B = b1, b2; C = c1
        genes <- data.table::data.table(species = c("A", "A", "A", "B", "B", "C"), id = c("a1", "a2", "a3", "b1", "b2", "c1"),
                                        gene = c("a1", "a2", "a3", "b1", "b2", "c1"), length = 100L)
        h <- function(q, s, bits, e = 1e-50) data.table::data.table(q = q, s = s, perc_identity = 90, alig_length = 100L,
                                                                    qcov = 100, scov = 100, evalue = e, bit_score = bits)
        hits <- rbind(
                h(1L, 4L, 200), h(4L, 1L, 200),      # a1 <-> b1 reciprocal best
                h(1L, 5L, 150), h(5L, 1L, 150),      # b2's best hit in A is a1, but a1 prefers b1
                h(2L, 4L, 120), h(4L, 2L, 120),      # a2's best hit in B is b1, not reciprocal
                h(1L, 2L, 300), h(2L, 1L, 290),      # a1 / a2: within-species, better than any between-species hit
                h(1L, 6L, 180), h(6L, 1L, 180),      # a1 <-> c1
                h(4L, 6L, 90),  h(6L, 4L, 95),       # b1 <-> c1
                h(3L, 5L, 80),                       # a3 -> b2 one direction only
                h(4L, 5L, 160), h(5L, 4L, 160)       # b1 / b2 within-species, but below b1's best between-species hit
        )
        best <- orthologr:::orthogroups_best_hits(hits, genes, "bitscore")
        expect_setequal(paste(best$q, best$s), c("1 4", "1 6", "2 4", "3 5", "4 1", "4 6", "5 1", "6 1", "6 4"))
        g <- orthologr:::orthogroups_graph(hits, genes, "bitscore")
        expect_identical(paste(g$rbh$q, g$rbh$s), c("1 4", "1 6", "4 6"))
        expect_identical(g$label, c(1L, 2L, 3L, 1L, 5L, 1L))
        # in-paralogs: a2 joins a1 (290 > 200 and > 120); b2 does not join b1 (160 < 200)
        g2 <- orthologr:::orthogroups_graph(hits, genes, "bitscore", inparalogs = TRUE)
        expect_identical(g2$label, c(1L, 1L, 3L, 1L, 5L, 1L))
        res <- orthologr:::orthogroups_assemble(genes, g2$label, hits, g2$rbh, species = c("C", "A", "B"))
        expect_identical(res$wide$orthogroup, "OG0000000")
        expect_identical(names(res$wide), c("orthogroup", "C", "A", "B"))
        expect_identical(res$wide$A, "a1, a2")
        expect_identical(nrow(res$single_copy), 0L)
        expect_setequal(res$unassigned$gene, c("a3", "b2"))
        p <- res$pairs
        expect_identical(p$relation[p$query_id == "c1" & p$subject_id == "b1"], "1:1")
        expect_identical(p$relation[p$query_id == "a2" & p$subject_id == "b1"], "co-ortholog")
        expect_identical(p$bit_score[p$query_id == "c1" & p$subject_id == "b1"], 95)
        expect_identical(p$bit_score[p$query_id == "b1" & p$subject_id == "c1"], 90)
        expect_true(is.na(p$bit_score[p$query_id == "a2" & p$subject_id == "c1"]))
        expect_identical(p$rbh[p$query_id == "a2" & p$subject_id == "b1"], FALSE)
        expect_identical(p$rbh[p$query_id == "b1" & p$subject_id == "a1"], TRUE)
        # ties: equal bit scores are decided by e-value, alignment length, identity, then the ID - as in diamond_rec()
        tie <- rbind(h(6L, 1L, 100), h(6L, 2L, 100), h(1L, 6L, 100), h(2L, 6L, 100))
        b <- orthologr:::orthogroups_best_hits(tie, genes, "bitscore")
        expect_identical(b$s[b$q == 6L], 1L)
        tie$evalue[2] <- 1e-60
        b <- orthologr:::orthogroups_best_hits(tie, genes, "bitscore")
        expect_identical(b$s[b$q == 6L], 2L)
})

test_that("orthogroups(method = 'rbh') finds the known groups of the simulated species", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        toy <- make_toy_species()
        og <- quiet_og(toy$files)
        expect_s3_class(og, "orthogroups")
        expect_identical(og$species, c("A", "B", "C", "D"))
        check_invariants(og)
        
        # the duplicated gene of species C: one copy is the reciprocal best hit, the other one is unassigned
        dup <- toy$truth$gene[toy$truth$family == 11 & toy$truth$species == "C"]
        out <- intersect(dup, og$unassigned$gene)
        expect_length(out, 1)
        expect_identical(og_group_strings(og), truth_strings(toy$truth, drop = out))
        # the genes without ortholog are unassigned
        expect_setequal(og$unassigned$gene, c(out, toy$truth$gene[is.element(toy$truth$family, 14:15)]))
        # 10 families in all species, plus family 11 with the copy that was kept
        expect_identical(nrow(og$single_copy), 11L)
        expect_identical(nrow(og$wide), 13L)
        expect_true(all(og$pairs$rbh))
        expect_true(all(og$pairs$relation == "1:1"))
        expect_false(anyNA(og$pairs$bit_score))
        expect_true(all(og$pairs$evalue <= 1e-5))
        # IDs of every style come back unchanged
        expect_setequal(og$genes$gene, toy$truth$gene)
        expect_identical(og$genes$gene, og$genes$id)
        
        # with in-paralogs both copies of C are in the group: co-orthologs
        ogi <- quiet_og(toy$files, inparalogs = TRUE)
        check_invariants(ogi)
        expect_identical(og_group_strings(ogi), truth_strings(toy$truth))
        expect_identical(nrow(ogi$single_copy), 10L)
        fam11 <- ogi$pairs[is.element(ogi$pairs$query_id, toy$truth$gene[toy$truth$family == 11]), ]
        expect_true(all(fam11$relation[fam11$query_species == "C" | fam11$subject_species == "C"] == "co-ortholog"))
        expect_true(all(fam11$relation[fam11$query_species != "C" & fam11$subject_species != "C"] == "1:1"))
        expect_identical(sum(!fam11$rbh), 6L)
        
        # print and summary
        expect_output(print(og), "13 groups with 49 of 52 genes; 11 single-copy")
        s <- summary(og)
        expect_identical(s$per_species$unassigned, c(1L, 1L, 1L, 0L))
        expect_identical(s$group_sizes$groups[s$group_sizes$genes == 4], 11L)
        expect_identical(s$species_per_group$groups, c(1L, 1L, 11L))
        expect_output(print(s), "Groups by number of species")
        df <- as.data.frame(og)
        expect_identical(names(df), c("orthogroup", "species", "gene"))
        expect_identical(nrow(df), 49L)
        expect_identical(sum(is.na(as.data.frame(og, unassigned = TRUE)$orthogroup)), 3L)
})

test_that("orthogroups() does not depend on comp_cores, on the order of the files, or on a folder as input", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        toy <- make_toy_species()
        og1 <- quiet_og(toy$files, inparalogs = TRUE)
        og2 <- suppressWarnings(quiet_og(toy$files, inparalogs = TRUE, comp_cores = 2))
        for (tbl in c("groups", "wide", "pairs", "single_copy", "unassigned", "genes"))
                expect_identical(og1[[tbl]], og2[[tbl]])
        og3 <- quiet_og(rev(toy$files), inparalogs = TRUE)
        expect_identical(og3$species, rev(og1$species))
        # the same groups with the same numbers; only the order of the species columns and rows follows the input
        expect_identical(og3$wide[ , names(og1$wide)], og1$wide)
        key <- function(p) p[order(p$orthogroup, p$query_id, p$subject_id), ]
        expect_identical(key(og3$pairs), key(og1$pairs))
        expect_identical(og3$groups[order(og3$groups$orthogroup, og3$groups$gene), ],
                         og1$groups[order(og1$groups$orthogroup, og1$groups$gene), ])
        # a folder: species names from the file names
        og4 <- quiet_og(toy$dir, inparalogs = TRUE)
        expect_identical(og4$wide, og1$wide)
        expect_identical(og4$pairs, og1$pairs)
})

test_that("orthogroups() with two species returns the reciprocal best hits of diamond_rec()", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        q <- system.file("seqs/ortho_thal_cds_1000.fasta", package = "orthologr")
        s <- system.file("seqs/ortho_lyra_cds_1000.fasta", package = "orthologr")
        for (rule in c("bitscore", "evalue")) {
                rec <- suppressMessages(diamond_rec(q, s, best_hit_rule = rule))
                og <- quiet_og(c(thal = q, lyra = s), best_hit_rule = rule)
                p <- og$pairs[og$pairs$query_species == "thal" & og$pairs$rbh, ]
                expect_setequal(paste(p$query_id, p$subject_id), paste(rec$query_id, rec$subject_id))
                expect_identical(nrow(p), nrow(rec))
                # the statistics are those of the search of thal against lyra: bit score, identity, alignment length
                # and coverage exactly, the e-value up to the three digits DIAMOND prints
                m <- merge(as.data.frame(rec), as.data.frame(p), by = c("query_id", "subject_id"), suffixes = c(".rec", ".og"))
                expect_identical(m$bit_score.rec, m$bit_score.og)
                expect_identical(m$perc_identity.rec, m$perc_identity.og)
                expect_identical(m$alig_length.rec, m$alig_length.og)
                expect_identical(m$qcov.rec, m$qcov.og)
                expect_true(all(abs(m$evalue.rec - m$evalue.og) <= 0.02 * m$evalue.rec))
        }
        # the 20-gene example with protein input and the species in the other order
        rec <- suppressMessages(diamond_rec(lyra_aa(), thal_aa(), seq_type = "protein"))
        og <- quiet_og(c(lyra = lyra_aa(), thal = thal_aa()), seq_type = "protein")
        p <- og$pairs[og$pairs$query_species == "lyra", ]
        expect_setequal(paste(p$query_id, p$subject_id), paste(rec$query_id, rec$subject_id))
        expect_identical(nrow(og$single_copy), nrow(rec))
})

test_that("orthogroups(): hit filters, IDs shared by species, three species from protein files", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        toy <- make_toy_species()
        # filters that no hit passes: no groups, all genes unassigned, tables with their columns
        none <- quiet_og(toy$files, min_identity = 100)
        expect_identical(nrow(none$wide), 0L)
        expect_identical(nrow(none$unassigned), 52L)
        expect_identical(names(none$pairs), names(quiet_og(toy$files)$pairs))
        expect_identical(names(none$wide), c("orthogroup", "A", "B", "C", "D"))
        expect_output(print(none), "0 groups")
        expect_identical(nrow(quiet_og(toy$files, min_qcov = 50, min_identity = 50)$wide), 13L)
        
        # the same IDs in two species: all genes are named <species>|<ID>, and it is said
        d <- new_test_dir()
        f1 <- file.path(d, "one.fa"); f2 <- file.path(d, "two.fa")
        aa <- Biostrings::readAAStringSet(thal_aa())
        names(aa) <- paste0("g", seq_along(aa))
        Biostrings::writeXStringSet(aa, f1)
        Biostrings::writeXStringSet(rev(aa), f2)
        suppressMessages(expect_message(dup <- orthogroups(c(one = f1, two = f2), seq_type = "protein"), "more than one species"))
        check_invariants(dup)
        expect_true(dup$params$ids_prefixed)
        expect_identical(nrow(dup$single_copy), 20L)
        expect_true(all(dup$single_copy$one == sub("^two", "one", dup$single_copy$two)))
        expect_identical(dup$genes$gene, paste0(dup$genes$species, "|", dup$genes$id))
        
        # four species of the example proteomes (protein input)
        files <- c(athal = system.file("seqs/example_athal_aa.faa", package = "orthologr"),
                   alyra = system.file("seqs/example_alyra_aa.faa", package = "orthologr"),
                   brapa = system.file("seqs/example_brapa_aa.faa", package = "orthologr"),
                   thalo = system.file("seqs/example_thalo_aa.faa", package = "orthologr"))
        og <- quiet_og(files, seq_type = "protein")
        check_invariants(og)
        expect_identical(og$species, names(files))
        expect_gt(nrow(og$wide), 0)
})

test_that("orthogroups(method = 'deepclust') clusters the simulated species; refine = TRUE splits by reciprocal best hits", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        toy <- make_toy_species()
        dc <- quiet_og(toy$files, method = "deepclust")
        check_invariants(dc)
        # every family is in one cluster (the duplicated gene with both copies)
        expect_identical(og_group_strings(dc), truth_strings(toy$truth))
        expect_true(all(is.na(dc$pairs$bit_score)))
        expect_true(all(is.na(dc$pairs$rbh)))
        # realign: statistics of the members with their representative
        re <- quiet_og(toy$files, method = "deepclust", realign = TRUE)
        expect_identical(re$wide, dc$wide)
        expect_gt(sum(!is.na(re$pairs$bit_score)), 0)
        # refined: as method = "rbh"
        rbh <- quiet_og(toy$files)
        ref <- quiet_og(toy$files, method = "deepclust", refine = TRUE)
        check_invariants(ref)
        expect_identical(og_group_strings(ref), og_group_strings(rbh))
        expect_identical(nrow(ref$single_copy), 11L)
        refi <- quiet_og(toy$files, method = "deepclust", refine = TRUE, inparalogs = TRUE, realign = TRUE)
        expect_identical(og_group_strings(refi), truth_strings(toy$truth))
        # arguments of deepclust() are handed on
        strict <- quiet_og(toy$files, method = "deepclust", approx_id = 99.9)
        expect_lt(nrow(strict$groups), nrow(dc$groups))
})

test_that("write_orthogroups() writes OrthoFinder's layout and the importers read it back", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        toy <- make_toy_species()
        og <- quiet_og(toy$files, inparalogs = TRUE)
        d <- file.path(new_test_dir(), "out")
        files <- write_orthogroups(og, d, ortholog_tables = TRUE)
        expect_true(all(file.exists(files)))
        expect_identical(basename(files[c("orthogroups", "single_copy", "unassigned", "pairs")]),
                         c("Orthogroups.tsv", "Orthogroups_SingleCopyOrthologues.txt", "Orthogroups_UnassignedGenes.tsv", "pairs.tsv"))
        
        # Orthogroups.tsv: header 'Orthogroup' and the species, genes separated by ", "
        lines <- readLines(files[["orthogroups"]])
        expect_identical(lines[1], "Orthogroup\tA\tB\tC\tD")
        expect_identical(length(lines), nrow(og$wide) + 1L)
        back <- utils::read.delim(files[["orthogroups"]], check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
        expect_identical(back$C, og$wide$C)
        expect_true(any(grepl(", ", back$C, fixed = TRUE)))
        # the importer of the package: the single-copy groups (IDs such as '007' stay text)
        core <- orthofinder2_retrieve_core_orthologs(files[["orthogroups"]], files[["single_copy"]])
        expect_identical(core$orthogroup, og$single_copy$orthogroup)
        for (k in c("A", "C", "D"))
                expect_identical(core[[k]], og$single_copy[[k]])
        expect_identical(sprintf("%03d", as.integer(core$B)), og$single_copy$B)
        # the layout orthofinder2(use_existing_output = TRUE) looks for
        found <- orthologr:::orthofinder_core_files(d)
        expect_identical(normalizePath(found$orthogroups_file), normalizePath(files[["orthogroups"]]))
        # unassigned genes: one per row, in the column of the species, numbered on
        un <- utils::read.delim(files[["unassigned"]], check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
        expect_identical(nrow(un), nrow(og$unassigned))
        expect_identical(un$Orthogroup[1], sprintf("OG%07d", nrow(og$wide)))
        expect_setequal(unlist(un[ , -1])[nzchar(unlist(un[ , -1]))], og$unassigned$gene)
        counts <- utils::read.delim(files[["gene_count"]], check.names = FALSE)
        expect_identical(sum(counts$Total), nrow(og$groups))
        # pairs.tsv
        p <- data.table::fread(files[["pairs"]], sep = "\t", colClasses = list(character = c("query_id", "subject_id")), quote = "")
        expect_identical(p$query_id, og$pairs$query_id)
        expect_equal(p$bit_score, og$pairs$bit_score)
        expect_identical(p$relation, og$pairs$relation)
        
        # the pairwise tables: import_ortholog_tables_all() and retrieve_core_orthologs()
        tabs <- suppressMessages(import_ortholog_tables_all(file.path(d, "ortholog_tables")))
        expect_setequal(tabs$subject_species, c("B", "C", "D"))
        rbh_a <- og$pairs[og$pairs$query_species == "A" & og$pairs$rbh, ]
        expect_identical(nrow(tabs), nrow(rbh_a))
        core2 <- suppressMessages(retrieve_core_orthologs(tabs, species_order = c("B", "C", "D")))
        # the genes of A with a reciprocal best hit in all three species: the 11 families present everywhere
        expect_identical(length(unique(core2$query_gene_locus_id)), 11L)
        expect_error(write_orthogroups(og, d, ortholog_tables = TRUE, query_species = "Z"), "query_species must be one of")
})

test_that("dNdS_orthogroups(): pairwise alignments give the values of dNdS(); results per pair and per group", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        skip_if_not(perl_available, "Perl is not installed")
        files <- c(thal = thal_cds(), lyra = lyra_cds())
        og <- quiet_og(files)
        ids <- og$wide$orthogroup[1:4]
        res <- suppressMessages(dNdS_orthogroups(og, groups = ids, aa_aln_tool = "pairwise"))
        expect_identical(names(res), c("pairs", "summary", "failed"))
        expect_identical(names(res$pairs), c("orthogroup", "query_species", "query_id", "subject_species", "subject_id",
                                             "method", "dN", "dS", "dNdS", "note"))
        expect_identical(res$pairs$orthogroup, ids)
        expect_true(all(res$pairs$query_species == "thal"))
        ref <- suppressMessages(suppressWarnings(utils::capture.output(
                d <- dNdS(thal_cds(), lyra_cds(), delete_corrupt_cds = TRUE, print_citation = FALSE))))
        m <- merge(as.data.frame(res$pairs), as.data.frame(d)[ , c("query_id", "subject_id", "dN", "dS", "dNdS")],
                   by = c("query_id", "subject_id"), suffixes = c("", ".ref"))
        expect_identical(nrow(m), 4L)
        expect_equal(m$dN, m$dN.ref, tolerance = 1e-12)
        expect_equal(m$dS, m$dS.ref, tolerance = 1e-12)
        expect_equal(m$dNdS, m$dNdS.ref, tolerance = 1e-12)
        s <- res$summary
        expect_identical(s$orthogroup, ids)
        expect_true(all(s$n_genes == 2 & s$n_species == 2 & s$n_pairs == 1 & s$n_ok == 1 & s$n_na == 0))
        expect_equal(s$dNdS_median, res$pairs$dNdS)
        expect_equal(s$dS_q1, res$pairs$dS)
        expect_identical(nrow(res$failed), 0L)
        
        # argument checks
        expect_error(dNdS_orthogroups(og, engine = "paml"), "engine = 'paml' is not supported")
        expect_error(dNdS_orthogroups(og, dnds_est.method = "XX"), "not supported with engine")
        expect_error(dNdS_orthogroups(og, aa_aln_tool = "NW"), "aa_aln_tool = 'NW' is not supported")
        expect_error(dNdS_orthogroups(og, groups = "OG9999999"), "no such orthogroup")
        expect_error(dNdS_orthogroups(og, eval = 1), "Unused argument")
        ogp <- quiet_og(c(thal = thal_aa(), lyra = lyra_aa()), seq_type = "protein")
        expect_error(dNdS_orthogroups(ogp), "computed from proteins")
        expect_error(dNdS_orthogroups(ogp, cds_files = c(thal = thal_cds())), "one file with coding sequences for each")
        expect_error(dNdS_orthogroups(ogp, cds_files = c(a = thal_cds(), b = lyra_cds())), "must be the species names")
})

test_that("dNdS_orthogroups(): one codon alignment per group, all between-species pairs, failures are NA rows", {
        skip_if_not(diamond_available, "DIAMOND is not installed")
        skip_if_not(perl_available, "Perl is not installed")
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        toy <- make_toy_species()
        og <- quiet_og(toy$files, inparalogs = TRUE)
        fam <- function(k) unique(og$groups$orthogroup[is.element(og$groups$gene, toy$truth$gene[toy$truth$family == k])])
        ids <- c(fam(1), fam(11), fam(12), fam(13))
        res <- suppressMessages(dNdS_orthogroups(og, groups = ids, aa_aln_tool = "clustalo"))
        p <- res$pairs
        s <- res$summary[match(ids, res$summary$orthogroup), ]
        # family 1: 4 genes, 6 pairs; family 11: 5 genes, the two copies of C are not compared (4*3/2 + 3 = 9 pairs);
        # family 12: 3 genes, 3 pairs; family 13: 2 genes, 1 pair
        expect_identical(s$n_pairs, c(6L, 9L, 3L, 1L))
        expect_identical(s$n_genes, c(4L, 5L, 3L, 2L))
        expect_identical(s$n_species, c(4L, 4L, 3L, 2L))
        expect_true(all(p$query_species != p$subject_species))
        expect_true(all(match(p$query_species, og$species) < match(p$subject_species, og$species)))
        expect_false(anyNA(p$dN))
        expect_true(all(is.na(p$note)))
        expect_true(all(p$dN > 0 & p$dN < 0.2))
        expect_equal(s$dN_median[1], stats::median(p$dN[p$orthogroup == ids[1]]))
        expect_equal(s$dNdS_mean[2], mean(p$dNdS[p$orthogroup == ids[2]], na.rm = TRUE))
        expect_equal(s$dS_q3[1], unname(stats::quantile(p$dS[p$orthogroup == ids[1]], 0.75, na.rm = TRUE)))
        expect_identical(s$n_ok + s$n_na, s$n_pairs)
        # the IDs are those of the input, whatever their form
        expect_true(all(is.element(c(p$query_id, p$subject_id), toy$truth$gene)))
        
        # close to the values of pairwise alignments (the sequences differ by substitutions only: the same alignment)
        pw <- suppressMessages(dNdS_orthogroups(og, groups = ids[c(1, 4)], aa_aln_tool = "pairwise"))
        m <- merge(as.data.frame(p), as.data.frame(pw$pairs), by = c("query_id", "subject_id"), suffixes = c("", ".pw"))
        expect_identical(nrow(m), 7L)
        expect_equal(m$dN, m$dN.pw, tolerance = 1e-10)
        expect_equal(m$dS, m$dS.pw, tolerance = 1e-10)
        
        # another method, and two processes: the same values
        li1 <- suppressMessages(dNdS_orthogroups(og, groups = ids, aa_aln_tool = "clustalo", dnds_est.method = "Li"))
        li2 <- suppressWarnings(suppressMessages(dNdS_orthogroups(og, groups = ids, aa_aln_tool = "clustalo", dnds_est.method = "Li",
                                                                  comp_cores = 2)))
        expect_identical(li1$pairs, li2$pairs)
        expect_true(all(li1$pairs$method == "Li"))
        expect_gt(stats::cor(li1$pairs$dN, p$dN), 0.95)
        
        # groups = "single_copy" (default) and "all" with max_group_size
        suppressMessages(expect_message(all4 <- dNdS_orthogroups(og, groups = "all", max_group_size = 4, aa_aln_tool = "clustalo"),
                                        "could not be used"))
        expect_identical(all4$failed$orthogroup, fam(11))
        expect_match(all4$failed$reason, "above max_group_size = 4")
        expect_identical(nrow(all4$summary), nrow(og$wide) - 1L)
        
        # a gene without usable coding sequence is left out and listed; a file from which a gene is missing as well
        d <- new_test_dir()
        cds <- lapply(toy$files, function(f) orthologr:::read_fasta_vector(f))
        g1 <- og$groups[og$groups$orthogroup == ids[1], ]
        bad_b <- g1$gene[g1$species == "B"]
        cds$B[[bad_b]] <- paste0(substr(cds$B[[bad_b]], 1, 30), "NNN", substring(cds$B[[bad_b]], 34))
        cds$D <- cds$D[names(cds$D) != g1$gene[g1$species == "D"]]
        cds_files <- vapply(names(cds), function(k) {
                f <- file.path(d, paste0(k, "_cds.fa"))
                writeLines(as.vector(rbind(paste0(">", names(cds[[k]])), cds[[k]])), f)
                f
        }, character(1))
        part <- suppressMessages(dNdS_orthogroups(og, cds_files = cds_files, groups = ids[1], aa_aln_tool = "clustalo"))
        expect_identical(nrow(part$pairs), 1L)
        expect_identical(c(part$pairs$query_species, part$pairs$subject_species), c("A", "C"))
        expect_setequal(part$failed$gene, c(bad_b, g1$gene[g1$species == "D"]))
        expect_match(part$failed$reason[part$failed$gene == bad_b], "other than A, C, G, T")
        expect_match(part$failed$reason[part$failed$gene != bad_b], "no coding sequence with this ID")
        expect_identical(part$summary$n_genes, 2L)
        
        # a pair that cannot be estimated is a row of NA with the reason; the other pairs are computed
        res_na <- suppressWarnings(suppressMessages(testthat::with_mocked_bindings(
                dNdS_orthogroups(og, groups = ids[3], aa_aln_tool = "clustalo"),
                substitutionrate = function(file, ...) {
                        if (grepl("_p2\\.fasta$", file)) stop("no estimate for this pair")
                        data.table::data.table(query_id = "q", subject_id = "s", dN = 0.1, dS = 0.5, dNdS = 0.2)
                }, .package = "orthologr")))
        expect_identical(is.na(res_na$pairs$dNdS), c(FALSE, TRUE, FALSE))
        expect_match(res_na$pairs$note[2], "no estimate for this pair")
        expect_identical(c(res_na$summary$n_ok, res_na$summary$n_na), c(2L, 1L))
        expect_equal(res_na$summary$dNdS_median, 0.2)
        expect_warning(suppressMessages(testthat::with_mocked_bindings(
                dNdS_orthogroups(og, groups = ids[3], aa_aln_tool = "clustalo"),
                multi_aln = function(...) stop("aligner broke"), .package = "orthologr")),
                "3 of 3 gene pairs could not be processed")
})

test_that("dNdS_orthogroups() runs orthogroups() when it is given files; orthofinder2() recommends orthogroups()", {
        expect_warning(expect_error(orthofinder2(file.path(tempdir(), "nowhere"))), "use orthogroups\\(\\) instead")
        expect_error(orthologs(thal_aa(), lyra_aa(), ortho_detection = "Orthofinder2"), "use orthogroups\\(\\)")
        skip_if_not(diamond_available, "DIAMOND is not installed")
        skip_if_not(perl_available, "Perl is not installed")
        skip_if_not(clustalo_available, "Clustal Omega is not installed")
        toy <- make_toy_species()
        res <- suppressMessages(dNdS_orthogroups(toy$files[c("A", "B")], aa_aln_tool = "clustalo", min_identity = 50))
        expect_identical(nrow(res$pairs), 13L)
        expect_true(all(res$summary$n_pairs == 1))
        expect_error(dNdS_orthogroups(toy$files, cds_files = toy$files), "'cds_files' must be NULL")
})
