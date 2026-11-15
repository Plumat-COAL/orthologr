# Which external programs are there? (asked once for all test files)
has_program <- function(x) nzchar(Sys.which(x))

blast_available   <- has_program("blastp") && has_program("makeblastdb")
diamond_available <- has_program("diamond")
perl_available    <- has_program("perl")
kaks_available    <- has_program("KaKs_Calculator")
clustalo_available <- has_program("clustalo")
clustalw_available <- has_program("clustalw2")
mafft_available   <- has_program("mafft")
muscle_available  <- has_program("muscle")
tcoffee_available <- has_program("t_coffee")

# Is the package loaded from its source folder (devtools::load_all()) rather than installed? Then new R processes
# (PSOCK workers) cannot load it. Asked without naming a package that is not a dependency: a package can only
# be loaded that way when the namespace of that tool is loaded.
loaded_with_load_all <- function() {
        dev_tool <- paste0("pkg", "load")
        if (!is.element(dev_tool, loadedNamespaces()))
                return(FALSE)
        isTRUE(tryCatch(get("is_dev_package", envir = asNamespace(dev_tool))("orthologr"), error = function(e) FALSE))
}

thal_cds <- function() system.file("seqs/ortho_thal_cds.fasta", package = "orthologr")
lyra_cds <- function() system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr")
thal_aa  <- function() system.file("seqs/ortho_thal_aa.fasta", package = "orthologr")
lyra_aa  <- function() system.file("seqs/ortho_lyra_aa.fasta", package = "orthologr")

# A folder for a restricted PATH that has the shell utilities R itself starts: parallel::detectCores() runs
# 'grep ... | wc -l' on Linux, and without them the test log showed "sh: 1: wc: not found". The programs a test
# wants to be missing stay missing.
with_coreutils <- function(dir) {
        for (p in c("grep", "wc")) {
                w <- Sys.which(p)
                if (nzchar(w) && !file.exists(file.path(dir, p)))
                        file.symlink(w, file.path(dir, p))
        }
        dir
}

# A new empty folder inside tempdir(); `name` may hold spaces.
new_test_dir <- function(name = "orthologr_test_") {
        d <- tempfile(name)
        dir.create(d, recursive = TRUE)
        d
}

# Synthetic annotation (NCBI-style GFF3 and Ensembl-style GTF) and a proteome with a second, shorter isoform for
# the first five genes, built from one of the example proteomes. Returns the paths and the ids.
make_test_annotation <- function(aa_file, species, gene_id_fun, dir) {
        aa <- Biostrings::readAAStringSet(aa_file)
        ids <- sub(" .*", "", names(aa))
        gid <- gene_id_fun(ids)
        iso <- Biostrings::subseq(aa[1:5], 1, pmax(10, Biostrings::width(aa[1:5]) %/% 2))
        iso_ids <- paste0(ids[1:5], "_iso2")
        all <- c(aa, iso)
        names(all) <- c(ids, iso_ids)
        gids <- c(gid, gid[1:5])
        tids <- names(all)
        proteome <- file.path(dir, paste0(species, ".faa"))
        Biostrings::writeXStringSet(all, proteome)
        st <- seq(1000, by = 5000, length.out = length(all))
        en <- st + Biostrings::width(all) * 3 + 2
        u <- !duplicated(gids)
        gff <- file.path(dir, paste0(species, ".gff"))
        writeLines(c(
                "##gff-version 3",
                sprintf("chr1\tRefSeq\tgene\t%d\t%d\t.\t+\t.\tID=gene-%s;Name=%s;gene=%s;gene_biotype=protein_coding", st[u], en[u], gids[u], gids[u], gids[u]),
                sprintf("chr1\tRefSeq\tmRNA\t%d\t%d\t.\t+\t.\tID=%s;Parent=gene-%s;Name=%s;gene=%s", st, en, tids, gids, gids, gids),
                sprintf("chr1\tRefSeq\tCDS\t%d\t%d\t.\t+\t0\tID=cds-%s;Parent=%s;Name=%s;gene=%s;protein_id=%s", st, en, tids, tids, gids, gids, tids)
        ), gff)
        gtf <- file.path(dir, paste0(species, ".gtf"))
        writeLines(c(
                sprintf("chr1\tensembl\tgene\t%d\t%d\t.\t+\t.\tgene_id \"%s\"; gene_biotype \"protein_coding\";", st[u], en[u], gids[u]),
                sprintf("chr1\tensembl\ttranscript\t%d\t%d\t.\t+\t.\tgene_id \"%s\"; transcript_id \"%s\"; gene_biotype \"protein_coding\";", st, en, gids, tids),
                sprintf("chr1\tensembl\tCDS\t%d\t%d\t.\t+\t0\tgene_id \"%s\"; transcript_id \"%s\"; gene_biotype \"protein_coding\";", st, en, gids, tids)
        ), gtf)
        list(proteome = proteome, gff = gff, gtf = gtf, ids = ids, gene_ids = gid, all_ids = tids, n_genes = sum(u))
}

# Both example species annotated in `dir`, and a dN/dS-like table that pairs the i-th genes (all isoforms).
make_test_ortholog_input <- function(dir = new_test_dir()) {
        q <- suppressWarnings(make_test_annotation(thal_aa(), "Athaliana", function(x) sub("\\..*", "", x), dir))
        s <- suppressWarnings(make_test_annotation(lyra_aa(), "Alyrata", function(x) paste0("G", sub("\\|.*", "", x)), dir))
        n <- min(length(q$all_ids), length(s$all_ids))
        # the columns of a dNdS() result
        dnds <- tibble::tibble(
                query_id      = q$all_ids[seq_len(n)],
                subject_id    = s$all_ids[seq_len(n)],
                dN            = seq(0.01, 0.2, length.out = n),
                dS            = seq(0.1, 0.9, length.out = n),
                dNdS          = seq(0.01, 0.2, length.out = n) / seq(0.1, 0.9, length.out = n),
                perc_identity = 90,
                num_ident_matches = 90L,
                alig_length   = rev(seq_len(n)) + 100L,
                mismatches    = 10L,
                gap_openings  = 0L,
                n_gaps        = 0L,
                pos_match     = 95L,
                ppos          = 95,
                q_start       = 1L,
                q_end         = 100L,
                q_len         = 100L,
                qcov          = 100,
                qcovhsp       = 100,
                s_start       = 1L,
                s_end         = 100L,
                s_len         = 100L,
                evalue        = 10^-(seq_len(n) + 5),
                bit_score     = 100 + seq_len(n),
                score_raw     = 200 + seq_len(n)
        )
        list(query = q, subject = s, dnds = dnds, dir = dir)
}
