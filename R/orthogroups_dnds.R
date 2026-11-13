#' @title dN/dS within orthogroups
#' @description Estimates dN, dS and dN/dS for every pair of genes of different species within each orthogroup and
#' summarises each group. For every group the proteins are aligned (one multiple alignment per group, or pairwise),
#' the codons are threaded onto the protein alignment (codon alignment) and the substitution rates are estimated for
#' the pairs of the group.
#' @param x the orthogroups: an object of class \code{"orthogroups"} (see \code{\link{orthogroups}}), or the files of
#' the coding sequences of the species (as the argument \code{files} of \code{orthogroups()}); in the second case
#' \code{orthogroups(x, seq_type = "cds", ...)} is run first.
#' @param cds_files the files with the coding sequences (FASTA) of the species of \code{x}: a character vector or
#' list named with the species names (or in the order of \code{x$species}). The sequence IDs must be those of the
#' genes of \code{x}. Default is \code{NULL}: the files \code{x} was computed from, when they hold coding sequences
#' (\code{seq_type = "cds"}).
#' @param groups \code{"single_copy"} (Default): the groups with exactly one gene in every species;
#' \code{"all"}: all groups with at most \code{max_group_size} genes (pairs of genes of the same species are not
#' computed); or a character vector with the IDs of the groups to compute (e.g. \code{c("OG0000000", "OG0000007")}).
#' @param dnds_est.method the dN/dS estimation method, see \code{\link{dNdS}} and \code{\link{substitutionrate}}.
#' Default (\code{NULL}) is \code{"Comeron"} with \code{engine = "orthologr"}.
#' @param engine \code{"orthologr"} (Default): the programs \code{\link{dNdS}} uses - \code{\link{multi_aln}} or
#' \code{\link{pairwise_aln}}, PAL2NAL (\code{\link{codon_aln}}) and \code{\link{substitutionrate}}.
#' \code{"coal"}: COAL (\url{https://github.com/Plumat-COAL/COAL}), which computes all groups in one call and is much
#' faster; only in versions of orthologr with the COAL integration (\code{dNdS_clusters()}).
#' @param aa_aln_tool how the proteins of a group are aligned: \code{"clustalo"} (Default), \code{"mafft"},
#' \code{"muscle"}, \code{"clustalw"} or \code{"t_coffee"} - one multiple alignment per group with that program,
#' from which the alignment of every pair is taken; or \code{"pairwise"} - every pair is aligned on its own
#' (global alignment, as \code{dNdS()} does by default; no external alignment program is needed).
#' @param aa_aln_path the folder that holds the alignment program (in case it is not on the PATH).
#' @param aa_aln_params parameters for the alignment program (see \code{\link{multi_aln}}).
#' @param max_group_size groups with more genes than this are not computed (\code{groups = "all"}); they are listed
#' in \code{failed}. Default is \code{50}.
#' @param dS_saturation pairs with a dS above this value are counted as saturated in the summary. Default is \code{2}.
#' @param kaks_calc_path the folder that holds KaKs_Calculator (for its methods, see \code{\link{substitutionrate}}).
#' @param delete_corrupt_cds a logical value: leave out coding sequences whose length is not a multiple of 3.
#' Default is \code{TRUE} (such sequences cannot be used in a codon alignment).
#' @param comp_cores the number of R processes (\code{engine = "orthologr"}) or threads.
#' @param quiet a logical value: do not show the messages of the alignment and estimation steps.
#' @param ... arguments for \code{\link{orthogroups}} when \code{x} holds files.
#' @details A gene, pair or group that cannot be used never ends the run. A gene without a usable coding sequence
#' (missing in \code{cds_files}, length not a multiple of 3, characters other than A, C, G, T) is left out and listed
#' in \code{failed}; a pair whose alignment or estimation fails has \code{NA} values and the reason in \code{note}.
#' The proteins are translated from the coding sequences, so that protein and coding sequence of a gene always fit.
#'
#' With a multiple alignment, the codon alignment of a pair is the two rows of the codon alignment of the group
#' without the columns in which both have a gap.
#' @return A list of three tibbles:
#' \itemize{
#' \item \code{pairs}: one row per pair of genes of different species in a group (and per method):
#' \code{orthogroup}, \code{query_species}, \code{query_id}, \code{subject_species}, \code{subject_id}, \code{method},
#' \code{dN}, \code{dS}, \code{dNdS}, \code{note}. The query is the gene of the species that comes first in \code{x$species};
#' \item \code{summary}: one row per group (and method): \code{orthogroup}, \code{method}, \code{n_genes},
#' \code{n_species}, \code{n_pairs}, \code{n_ok} (pairs with a dN/dS value), \code{n_na}, \code{n_saturated}
#' (dS above \code{dS_saturation}), and median, mean, first and third quartile of dN, dS and dN/dS over the pairs
#' with values;
#' \item \code{failed}: genes and groups that could not be used (\code{orthogroup}, \code{gene}, \code{reason}).
#' }
#' @author Chloe Newton
#' @examples \dontrun{
#' files <- c(
#'   Athaliana = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   Alyrata   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' og <- orthogroups(files)
#' res <- dNdS_orthogroups(og, aa_aln_tool = "pairwise")
#' res$pairs
#' res$summary
#'
#' # in one step, all groups, one Clustal Omega alignment per group
#' res <- dNdS_orthogroups(files, groups = "all", aa_aln_tool = "clustalo", comp_cores = 2)
#' }
#' @seealso \code{\link{orthogroups}}, \code{\link{dNdS}}, \code{\link{substitutionrate}}
#' @export
dNdS_orthogroups <- function(x,
                             cds_files          = NULL,
                             groups             = c("single_copy", "all"),
                             dnds_est.method    = NULL,
                             engine             = c("orthologr", "coal"),
                             aa_aln_tool        = "clustalo",
                             aa_aln_path        = NULL,
                             aa_aln_params      = NULL,
                             max_group_size     = 50,
                             dS_saturation      = 2,
                             kaks_calc_path     = NULL,
                             delete_corrupt_cds = TRUE,
                             comp_cores         = 1,
                             quiet              = TRUE,
                             ...) {
        
        group_ids <- NULL
        if (is.character(groups) && length(groups) >= 1 && !anyNA(groups) && !is.element(groups[1], c("single_copy", "all"))) {
                group_ids <- unique(groups)
                groups <- "all"
        } else {
                groups <- check_choice(groups[1], c("single_copy", "all"), "groups")
        }
        engine <- check_choice(engine[1], c("orthologr", "coal"), "engine")
        check_count(max_group_size, "max_group_size")
        check_flag(quiet, "quiet")
        if (!is.numeric(dS_saturation) || length(dS_saturation) != 1 || is.na(dS_saturation) || dS_saturation <= 0)
                stop("dS_saturation must be one positive number.", call. = FALSE)
        if (engine == "orthologr") {
                if (is.null(dnds_est.method))
                        dnds_est.method <- "Comeron"
                if (!is.character(dnds_est.method) || length(dnds_est.method) != 1 || !is.dnds_est_method(dnds_est.method))
                        stop("dnds_est.method = '", paste(dnds_est.method, collapse = "', '"), "' is not supported with engine = \"orthologr\". ",
                             "Please choose one of: 'Comeron', 'Li', 'NG', 'LWL', 'LPB', 'MLWL', 'YN', 'MYN', 'GY', 'MA', 'MS', 'GMYN', 'ALL'.",
                             call. = FALSE)
                check_choice(aa_aln_tool, c("clustalo", "mafft", "muscle", "clustalw", "t_coffee", "pairwise"), "aa_aln_tool")
        } else {
                dnds_est.method <- dnds_groups_coal_check(dnds_est.method, aa_aln_tool)
        }
        comp_cores <- check_comp_cores(comp_cores)
        
        if (!inherits(x, "orthogroups")) {
                if (!is.null(cds_files))
                        stop("When 'x' holds the files of the coding sequences, 'cds_files' must be NULL.", call. = FALSE)
                x <- orthogroups(x, seq_type = "cds", comp_cores = comp_cores, quiet = quiet, ...)
        } else if (...length() > 0) {
                stop("Unused argument", if (...length() > 1) "s", ": ", paste(names(list(...)), collapse = ", "),
                     " ('...' is handed to orthogroups() when 'x' holds files).", call. = FALSE)
        }
        cds_files <- dnds_groups_cds_files(x, cds_files)
        
        # the groups to compute and their genes
        orthogroup <- gene <- n <- query_species <- query_id <- subject_species <- subject_id <- method <- NULL
        members <- data.table::as.data.table(x$groups)
        failed <- list()
        if (!is.null(group_ids)) {
                unknown <- setdiff(group_ids, members$orthogroup)
                if (length(unknown) > 0)
                        stop("groups = '", paste(utils::head(unknown, 5), collapse = "', '"), "'", if (length(unknown) > 5) ", ..." else "",
                             ": no such orthogroup in 'x'. Please give \"single_copy\", \"all\" or IDs of x$wide$orthogroup.", call. = FALSE)
                members <- members[is.element(orthogroup, group_ids)]
        }
        if (groups == "single_copy") {
                members <- members[is.element(orthogroup, x$single_copy$orthogroup)]
        } else {
                size <- members[ , list(n = .N), by = orthogroup]
                large <- size[n > max_group_size]
                if (nrow(large) > 0) {
                        failed[[1]] <- data.table::data.table(
                                orthogroup = large$orthogroup, gene = NA_character_,
                                reason = paste0("group of ", large$n, " genes is above max_group_size = ", max_group_size))
                        members <- members[!is.element(orthogroup, large$orthogroup)]
                }
        }
        if (nrow(members) == 0)
                stop("There are no groups to compute (groups = \"", groups, "\": ",
                     if (groups == "single_copy") "no group has exactly one gene in every species"
                     else paste0("no group with at most ", max_group_size, " genes"), ").", call. = FALSE)
        
        # the coding sequences of these genes
        cds <- dnds_groups_read_cds(x, cds_files, delete_corrupt_cds = delete_corrupt_cds)
        usable <- is.element(members$gene, names(cds$cds))
        if (!any(usable))
                stop("None of the ", nrow(members), " genes of the groups has a coding sequence in 'cds_files'. ",
                     "The sequence IDs of the coding sequences must be the IDs of the genes (e.g. '", members$gene[1],
                     "'; first ID read: '", cds$first_id, "').", call. = FALSE)
        if (any(!usable)) {
                why <- cds$unusable[members$gene[!usable]]
                why[is.na(why)] <- "no coding sequence with this ID in cds_files"
                failed[[length(failed) + 1L]] <- data.table::data.table(
                        orthogroup = members$orthogroup[!usable], gene = members$gene[!usable], reason = unname(why))
        }
        
        # all pairs of genes of different species within a group; the query is the gene of the earlier species
        sp <- species <- sp.q <- sp.s <- gene.q <- gene.s <- NULL
        m <- members[usable]
        m[ , sp := match(species, x$species)]
        pairs <- merge(m, m, by = "orthogroup", allow.cartesian = TRUE, suffixes = c(".q", ".s"))
        pairs <- pairs[sp.q < sp.s]
        data.table::setorder(pairs, orthogroup, sp.q, gene.q, sp.s, gene.s)
        pair_tbl <- data.table::data.table(orthogroup = pairs$orthogroup, query_id = pairs$gene.q, subject_id = pairs$gene.s)
        
        message("dNdS_orthogroups(): ", nrow(pair_tbl), " gene pairs in ", length(unique(pair_tbl$orthogroup)), " groups (",
                dnds_est.method[1], ", engine = \"", engine, "\", aa_aln_tool = \"", aa_aln_tool, "\") ...")
        
        if (nrow(pair_tbl) == 0) {
                rates <- data.table::data.table(orthogroup = character(0), query_id = character(0), subject_id = character(0),
                                                method = character(0), dN = numeric(0), dS = numeric(0), dNdS = numeric(0),
                                                note = character(0))
        } else if (engine == "orthologr") {
                rates <- dnds_groups_orthologr(
                        pairs = pair_tbl, cds = cds$cds, dnds_est.method = dnds_est.method, aa_aln_tool = aa_aln_tool,
                        aa_aln_path = aa_aln_path, aa_aln_params = aa_aln_params, kaks_calc_path = kaks_calc_path,
                        comp_cores = comp_cores, quiet = quiet
                )
        } else {
                rates <- dnds_groups_coal(
                        pairs = pair_tbl, cds = cds$cds, dnds_est.method = dnds_est.method, aa_aln_tool = aa_aln_tool,
                        aa_aln_path = aa_aln_path, comp_cores = comp_cores, quiet = quiet
                )
        }
        
        gene_species <- stats::setNames(m$species, m$gene)
        res_pairs <- data.table::data.table(
                orthogroup      = rates$orthogroup,
                query_species   = unname(gene_species[rates$query_id]),
                query_id        = rates$query_id,
                subject_species = unname(gene_species[rates$subject_id]),
                subject_id      = rates$subject_id,
                method          = rates$method,
                dN              = as.numeric(rates$dN),
                dS              = as.numeric(rates$dS),
                dNdS            = as.numeric(rates$dNdS),
                note            = rates$note
        )
        data.table::setorder(res_pairs, orthogroup, query_species, query_id, subject_species, subject_id, method)
        
        n_bad <- sum(!is.na(res_pairs$note))
        if (n_bad > 0)
                warning(n_bad, " of ", nrow(res_pairs), " gene pairs could not be processed; their dN, dS and dNdS values are NA ",
                        "and the reason is in the column 'note' of the table 'pairs'. First one: '", res_pairs$query_id[!is.na(res_pairs$note)][1],
                        "' / '", res_pairs$subject_id[!is.na(res_pairs$note)][1], "': ",
                        strsplit(res_pairs$note[!is.na(res_pairs$note)][1], "\n", fixed = TRUE)[[1]][1], call. = FALSE)
        
        failed <- if (length(failed) > 0) data.table::rbindlist(failed) else
                data.table::data.table(orthogroup = character(0), gene = character(0), reason = character(0))
        if (nrow(failed) > 0)
                message(nrow(failed), " gene(s) or group(s) could not be used; see the 'failed' table of the result.")
        
        list(
                pairs   = tibble::as_tibble(res_pairs),
                summary = tibble::as_tibble(dnds_groups_summary(res_pairs, m, dS_saturation = dS_saturation)),
                failed  = tibble::as_tibble(failed)
        )
}

# ---- input ----------------------------------------------------------------------------------------------------

# The CDS file of every species of x: a character vector named with the species.
dnds_groups_cds_files <- function(x, cds_files) {
        sp <- x$species
        if (is.null(cds_files)) {
                if (!identical(x$params$seq_type, "cds"))
                        stop("The orthogroups were computed from proteins: please name the files with the coding sequences ",
                             "of the species in 'cds_files' (", paste(sp, collapse = ", "), ").", call. = FALSE)
                cds_files <- x$params$files
        }
        if (is.list(cds_files))
                cds_files <- unlist(cds_files)
        if (!is.character(cds_files) || length(cds_files) != length(sp))
                stop("'cds_files' must name one file with coding sequences for each of the ", length(sp), " species (",
                     paste(sp, collapse = ", "), ").", call. = FALSE)
        if (is.null(names(cds_files)))
                names(cds_files) <- sp
        if (!setequal(names(cds_files), sp))
                stop("The names of 'cds_files' must be the species names: '", paste(sp, collapse = "', '"), "' (given: '",
                     paste(names(cds_files), collapse = "', '"), "').", call. = FALSE)
        cds_files <- cds_files[sp]
        missing_files <- cds_files[!file.exists(cds_files)]
        if (length(missing_files) > 0)
                stop("The following input file(s) were not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
        cds_files
}

# The coding sequences of all species, named with the gene names of x. Returns
#   cds      : named character vector (upper case) of the usable sequences
#   unusable : named character vector gene -> reason
dnds_groups_read_cds <- function(x, cds_files, delete_corrupt_cds = TRUE) {
        geneids <- seqs <- NULL
        genes <- x$genes
        all <- list()
        unusable <- list()
        first_id <- NA_character_
        for (k in x$species) {
                dt <- read.cds(file = cds_files[[k]], format = "fasta", delete_corrupt_cds = delete_corrupt_cds)
                ids <- as.character(dt[ , geneids])
                check_duplicate_ids(ids, what = "coding sequences", file = cds_files[[k]])
                if (is.na(first_id) && length(ids) > 0)
                        first_id <- ids[1]
                seq <- toupper(as.character(dt[ , seqs]))
                gk <- genes[genes$species == k, , drop = FALSE]
                at <- match(gk$id, ids)
                seq <- seq[at]
                names(seq) <- gk$gene
                seq <- seq[!is.na(at)]
                bad_len <- nchar(seq) %% 3 != 0 | nchar(seq) == 0
                bad_chr <- grepl("[^ACGT]", seq)
                unusable[[k]] <- c(
                        stats::setNames(rep("the length of the coding sequence is not a multiple of 3", sum(bad_len)), names(seq)[bad_len]),
                        stats::setNames(rep("the coding sequence has characters other than A, C, G, T", sum(bad_chr & !bad_len)),
                                        names(seq)[bad_chr & !bad_len])
                )
                all[[k]] <- seq[!bad_len & !bad_chr]
        }
        list(cds = unlist(unname(all)), unusable = unlist(unname(unusable)), first_id = first_id)
}

# ---- engine "orthologr" ---------------------------------------------------------------------------------------

# dN/dS of the gene pairs `pairs` (data.table: orthogroup, query_id, subject_id) with the programs of dNdS().
#   cds : named character vector of the coding sequences (names: gene names)
# Returns a data.table: orthogroup, query_id, subject_id, method, dN, dS, dNdS, note (NA, or why the values are NA).
# This is the interface of an engine; dnds_groups_coal() has the same one.
dnds_groups_orthologr <- function(pairs, cds, dnds_est.method, aa_aln_tool, aa_aln_path = NULL, aa_aln_params = NULL,
                                  kaks_calc_path = NULL, comp_cores = 1, quiet = TRUE) {
        
        end_tmpdir <- orthologr_tmpdir_start("orthologr_dNdS_orthogroups_")
        on.exit(end_tmpdir(), add = TRUE)
        
        used <- unique(c(pairs$query_id, pairs$subject_id))
        aa <- stats::setNames(toupper(translate_cds(cds[used])), used)
        
        if (aa_aln_tool == "pairwise") {
                # every pair on its own: the pipeline of dNdS() (global protein alignment, PAL2NAL, estimation)
                tbl <- data.table::data.table(
                        query_id = pairs$query_id, subject_id = pairs$subject_id,
                        query_cds = unname(cds[pairs$query_id]), subject_cds = unname(cds[pairs$subject_id]),
                        query_aa = unname(aa[pairs$query_id]), subject_aa = unname(aa[pairs$subject_id])
                )
                run <- function(expr) if (quiet) suppressMessages(expr) else expr
                res <- run(withCallingHandlers(
                        compute_dnds(complete_tbl = tbl, aa_aln_type = "pairwise", aa_aln_tool = "NW",
                                     codon_aln_tool = "pal2nal", dnds_est.method = dnds_est.method,
                                     kaks_calc_path = kaks_calc_path, quiet = TRUE, comp_cores = comp_cores),
                        # (the pairs that failed are reported by dNdS_orthogroups())
                        warning = function(w) {
                                if (grepl("gene pairs could not be processed", conditionMessage(w), fixed = TRUE))
                                        invokeRestart("muffleWarning")
                        }
                ))
                fp <- attr(res, "failed_pairs")
                res <- data.table::as.data.table(res)
                if (!is.element("method", names(res)))
                        res[["method"]] <- dnds_est.method
                res[["note"]] <- if (is.null(fp) || nrow(fp) == 0) NA_character_ else
                        fp$reason[match(paste(res$query_id, res$subject_id, sep = "\r"), paste(fp$query_id, fp$subject_id, sep = "\r"))]
                res[["orthogroup"]] <- pairs$orthogroup[match(paste(res$query_id, res$subject_id, sep = "\r"),
                                                              paste(pairs$query_id, pairs$subject_id, sep = "\r"))]
                return(res[ , c("orthogroup", "query_id", "subject_id", "method", "dN", "dS", "dNdS", "note"), with = FALSE])
        }
        
        # one multiple alignment and one codon alignment per group
        for (d in c("AA", "CDS", "pairs")) {
                dir.create(file.path(orthologr_tmpdir(), "_alignment", "orthogroups", d), recursive = TRUE, showWarnings = FALSE)
        }
        by_group <- split(seq_len(nrow(pairs)), factor(pairs$orthogroup, levels = unique(pairs$orthogroup)))
        one_group <- function(g) {
                rows <- by_group[[g]]
                dnds_groups_one_msa(
                        name = paste0("g", g), orthogroup = pairs$orthogroup[rows[1]], query_id = pairs$query_id[rows],
                        subject_id = pairs$subject_id[rows], cds = cds, aa = aa, dnds_est.method = dnds_est.method,
                        aa_aln_tool = aa_aln_tool, aa_aln_path = aa_aln_path, aa_aln_params = aa_aln_params,
                        kaks_calc_path = kaks_calc_path
                )
        }
        run <- function(expr) if (quiet) suppressMessages(expr) else expr
        if (comp_cores > 1 && length(by_group) > 1) {
                cl <- dnds_cluster(min(comp_cores, length(by_group)))
                on.exit(parallel::stopCluster(cl), add = TRUE)
                res <- parallel::parLapply(cl, seq_along(by_group), one_group)
        } else {
                res <- run(lapply(seq_along(by_group), one_group))
        }
        data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
}

# One group: protein alignment, codon alignment, and the rates of its pairs. A failure of the alignment gives NA rows
# for all pairs of the group, a failure of one pair an NA row for that pair.
dnds_groups_one_msa <- function(name, orthogroup, query_id, subject_id, cds, aa, dnds_est.method, aa_aln_tool,
                                aa_aln_path, aa_aln_params, kaks_calc_path) {
        na_rows <- function(q, s, reason) {
                data.table::data.table(orthogroup = orthogroup, query_id = q, subject_id = s, method = dnds_est.method,
                                       dN = NA_real_, dS = NA_real_, dNdS = NA_real_, note = as.character(reason))
        }
        ids <- unique(c(query_id, subject_id))
        # the sequences are called s1, s2, ... in the files: alignment programs and their formats change IDs
        local <- stats::setNames(paste0("s", seq_along(ids)), ids)
        base <- file.path(orthologr_tmpdir(), "_alignment", "orthogroups")
        aa_file  <- file.path(base, "AA", paste0(name, "_aa.fasta"))
        cds_file <- file.path(base, "CDS", paste0(name, "_cds.fasta"))
        writeLines(as.vector(rbind(paste0(">", local), aa[ids])), aa_file)
        writeLines(as.vector(rbind(paste0(">", local), cds[ids])), cds_file)
        
        codon <- tryCatch({
                utils::capture.output(
                        f <- multi_aln(file = aa_file, tool = aa_aln_tool, get_aln = FALSE, multi_aln_name = name,
                                       path = aa_aln_path, params = aa_aln_params, quiet = TRUE)
                )
                codon_file <- codon_aln(file_aln = f, file_nuc = cds_file, format = "fasta", tool = "pal2nal",
                                        codon_aln_name = name, get_aln = FALSE, quiet = TRUE)
                aln <- read_fasta_vector(codon_file)
                if (!all(is.element(local, names(aln))) || length(unique(nchar(aln))) != 1 || nchar(aln[1]) %% 3 != 0)
                        stop("the codon alignment of the group does not hold all its genes in full codons.", call. = FALSE)
                aln
        }, error = function(e) e)
        if (inherits(codon, "error"))
                return(na_rows(query_id, subject_id, paste0("alignment of the group failed: ", conditionMessage(codon))))
        
        n_codons <- nchar(codon[1]) %/% 3
        starts <- seq.int(1, by = 3, length.out = n_codons)
        codons_of <- lapply(stats::setNames(local, local), function(l) substring(codon[[l]], starts, starts + 2L))
        
        res <- lapply(seq_along(query_id), function(i) {
                tryCatch({
                        a <- codons_of[[local[[query_id[i]]]]]
                        b <- codons_of[[local[[subject_id[i]]]]]
                        keep <- !(a == "---" & b == "---")
                        if (!any(keep))
                                stop("the two sequences share no aligned codon.", call. = FALSE)
                        pair_file <- file.path(base, "pairs", paste0(name, "_p", i, ".fasta"))
                        writeLines(c(">q", paste(a[keep], collapse = ""), ">s", paste(b[keep], collapse = "")), pair_file)
                        on.exit(unlink(pair_file), add = TRUE)
                        r <- substitutionrate(file = pair_file, est.method = dnds_est.method, format = "fasta", quiet = TRUE,
                                              kaks_calc_path = kaks_calc_path, subst_name = paste0(name, "_p", i))
                        r <- data.table::as.data.table(r)
                        data.table::data.table(
                                orthogroup = orthogroup, query_id = query_id[i], subject_id = subject_id[i],
                                method = if (is.element("method", names(r))) as.character(r$method) else dnds_est.method,
                                dN = as.numeric(r$dN), dS = as.numeric(r$dS), dNdS = as.numeric(r$dNdS), note = NA_character_
                        )
                }, error = function(e) na_rows(query_id[i], subject_id[i], conditionMessage(e)))
        })
        unlink(c(aa_file, cds_file))
        data.table::rbindlist(res)
}

# The sequences of a FASTA file as a named character vector (name: header up to the first white space).
read_fasta_vector <- function(file) {
        x <- sub("\r$", "", readLines(file, warn = FALSE))
        x <- x[nzchar(x)]
        head_at <- startsWith(x, ">")
        if (!any(head_at))
                return(character(0))
        grp <- cumsum(head_at)
        seqs <- vapply(split(x[!head_at], factor(grp[!head_at], levels = seq_len(sum(head_at)))), paste, character(1), collapse = "")
        stats::setNames(unname(seqs), sub("\\s.*$", "", sub("^>\\s*", "", x[head_at])))
}

# ---- summary --------------------------------------------------------------------------------------------------

# One row per group and method.
dnds_groups_summary <- function(pairs, members, dS_saturation = 2) {
        orthogroup <- method <- dN <- dS <- dNdS <- species <- NULL
        if (nrow(pairs) == 0)
                return(data.table::data.table(orthogroup = character(0), method = character(0)))
        q <- function(v, p) if (all(is.na(v))) NA_real_ else unname(stats::quantile(v, p, na.rm = TRUE, names = FALSE))
        mn <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
        size <- members[ , list(n_genes = .N, n_species = data.table::uniqueN(species)), by = orthogroup]
        s <- pairs[ , list(
                n_pairs     = .N,
                n_ok        = sum(!is.na(dNdS)),
                n_na        = sum(is.na(dNdS)),
                n_saturated = sum(!is.na(dS) & dS > dS_saturation),
                dN_median = q(dN, 0.5), dN_mean = mn(dN), dN_q1 = q(dN, 0.25), dN_q3 = q(dN, 0.75),
                dS_median = q(dS, 0.5), dS_mean = mn(dS), dS_q1 = q(dS, 0.25), dS_q3 = q(dS, 0.75),
                dNdS_median = q(dNdS, 0.5), dNdS_mean = mn(dNdS), dNdS_q1 = q(dNdS, 0.25), dNdS_q3 = q(dNdS, 0.75)
        ), by = list(orthogroup, method)]
        s <- merge(size, s, by = "orthogroup")
        data.table::setcolorder(s, c("orthogroup", "method"))
        data.table::setorder(s, orthogroup, method)
        s
}
