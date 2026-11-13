# engine = "coal" of dNdS_orthogroups(): the gene pairs of all groups are handed to COAL
# (https://github.com/Plumat-COAL/COAL) through the R package 'coalr', as dNdS_clusters() does. The two functions
# here are the whole interface of the engine (see dnds_groups_orthologr() for the other one).

# Check the arguments of the COAL engine; returns the method(s) to compute.
dnds_groups_coal_check <- function(dnds_est.method, aa_aln_tool) {
        if (is.null(dnds_est.method))
                dnds_est.method <- "YN"
        check_choice(aa_aln_tool, c("clustalo", "mafft", "pairwise"), "aa_aln_tool",
                     hint = "(with engine = \"coal\"; 'muscle', 'clustalw' and 't_coffee' need engine = \"orthologr\").")
        coal_backend(need = "coalr")
        coal_check_methods(dnds_est.method, "engine = \"coal\"")
        dnds_est.method
}

# dN/dS of the gene pairs `pairs` (data.table: orthogroup, query_id, subject_id) with COAL.
#   aa_aln_tool = "pairwise"            : coalr::coal_dnds_pairs() - every pair aligned on its own
#   aa_aln_tool = "clustalo" / "mafft"  : coalr::coal_dnds_clusters() - one multiple alignment per group, as
#                                         dNdS_clusters(); the pairs of genes of the same species it also computes
#                                         are left out
# Returns the table dnds_groups_orthologr() returns: orthogroup, query_id, subject_id, method, dN, dS, dNdS, note.
#   aa_aln_scoring : the scoring of the pairwise protein alignments COAL makes itself ("pwalign": that of dNdS();
#                    "needle": COAL's own defaults) - the same argument, with the same default, as in
#                    dNdS(engine = "coal"), dNdS_coal() and dNdS_clusters()
dnds_groups_coal <- function(pairs, cds, dnds_est.method, aa_aln_tool, aa_aln_path = NULL, comp_cores = 1, quiet = TRUE,
                             aa_aln_scoring = "pwalign") {
        used <- unique(c(pairs$query_id, pairs$subject_id))
        if (aa_aln_tool == "pairwise") {
                res <- do.call(coalr::coal_dnds_pairs, c(list(
                        query_cds   = cds[unique(pairs$query_id)],
                        subject_cds = cds[unique(pairs$subject_id)],
                        pairs       = data.frame(query_id = pairs$query_id, subject_id = pairs$subject_id, stringsAsFactors = FALSE),
                        methods     = dnds_est.method,
                        threads     = comp_cores,
                        preset      = aa_aln_scoring
                ), coal_flavour_args(coalr::coal_dnds_pairs)))
        } else {
                exe <- if (is.null(aa_aln_path)) NULL else {
                        if (dir.exists(aa_aln_path)) file.path(aa_aln_path, aa_aln_tool) else aa_aln_path
                }
                msa_aligner <- tryCatch(
                        coalr::coal_msa_aligner(aa_aln_tool, path = exe, threads = 1),
                        error = function(e) stop(conditionMessage(e),
                                                 "\nOr use aa_aln_tool = \"pairwise\" (no external program needed).",
                                                 call. = FALSE)
                )
                # the groups as clusters (representative: the first gene of the group)
                members <- unique(data.frame(orthogroup = c(pairs$orthogroup, pairs$orthogroup),
                                             member_id = c(pairs$query_id, pairs$subject_id), stringsAsFactors = FALSE))
                clusters <- data.frame(representative_id = members$member_id[match(members$orthogroup, members$orthogroup)],
                                       member_id = members$member_id, stringsAsFactors = FALSE)
                out <- do.call(coalr::coal_dnds_clusters, c(list(
                        cds              = cds[used],
                        clusters         = clusters,
                        aligner          = msa_aligner,
                        methods          = dnds_est.method,
                        threads          = comp_cores,
                        max_cluster_size = max(table(members$orthogroup), 2L),
                        large_clusters   = "skip",
                        # groups of two genes are aligned pairwise by COAL
                        preset           = aa_aln_scoring
                ), coal_flavour_args(coalr::coal_dnds_clusters)))
                res <- out$pairs
        }
        # Comeron and Li as orthologr's own methods return them (saturation rule, 4 digits, Li's markers)
        res <- data.table::as.data.table(coal_orthologr_values(res))
        # COAL names the two genes seq1 and seq2, in an order of its own: both orders are looked up
        key <- function(a, b) paste(a, b, sep = "\r")
        at <- match(key(res$seq1, res$seq2), key(pairs$query_id, pairs$subject_id))
        swapped <- is.na(at)
        at[swapped] <- match(key(res$seq2, res$seq1), key(pairs$query_id, pairs$subject_id))[swapped]
        res <- res[!is.na(at)]
        at <- at[!is.na(at)]
        note <- ifelse(res$status == "ok", NA_character_,
                       ifelse(is.na(res$note) | !nzchar(res$note), as.character(res$status), paste0(res$status, ": ", res$note)))
        got <- data.table::data.table(
                orthogroup = pairs$orthogroup[at], query_id = pairs$query_id[at], subject_id = pairs$subject_id[at],
                method = as.character(res$method), dN = as.numeric(res$dN), dS = as.numeric(res$dS),
                dNdS = as.numeric(res$dNdS), note = note
        )
        # pairs COAL returned no row for (e.g. a group it could not align)
        lost <- setdiff(seq_len(nrow(pairs)), at)
        if (length(lost) > 0)
                got <- rbind(got, data.table::data.table(
                        orthogroup = pairs$orthogroup[lost], query_id = pairs$query_id[lost], subject_id = pairs$subject_id[lost],
                        method = dnds_est.method[1], dN = NA_real_, dS = NA_real_, dNdS = NA_real_,
                        note = "COAL returned no result for this pair"))
        got
}

# The versions of the programs engine = "coal" used.
dnds_groups_coal_versions <- function(aa_aln_tool, aa_aln_path = NULL) {
        c(if (aa_aln_tool != "pairwise") run_tool_versions(aa_aln_tool = aa_aln_tool, aa_aln_path = aa_aln_path) else
                list(Biostrings = package_version_text("Biostrings"), R = paste(R.version$major, R.version$minor, sep = ".")),
          coal_versions())
}
