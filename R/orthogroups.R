#' @title Orthogroups of two or more species with DIAMOND
#' @description Infers groups of orthologous genes (orthogroups) for 2 to N species with DIAMOND alone - no
#' OrthoFinder, MCL or tree building is needed. Two methods are available:
#' \itemize{
#' \item \code{method = "rbh"} (Default; recommended for about 2 to 15 species): one all-against-all DIAMOND search of
#' the proteins of all species, reciprocal best hits (RBH) for every pair of species, and the connected components
#' of the graph whose edges are the reciprocal best hits. Optionally, in-paralogs are attached (\code{inparalogs}).
#' \item \code{method = "deepclust"} (many species or very large protein sets): the proteins of all species are
#' clustered with \code{\link{deepclust}} (linear-time clustering); a cluster with genes of at least two species is an
#' orthogroup. With \code{refine = TRUE} clusters that hold more than one gene of a species are split by the RBH graph
#' inside the cluster.
#' }
#' @param files the sequence files (FASTA), one per species: a named character vector or list
#' (\code{c(Athaliana = "a.fa", Alyrata = "b.fa")}; the names are the species names), an unnamed one (the species
#' names are the file names without extension), or the path of one folder (all FASTA files in it are used; species
#' names from the file names). Every sequence is treated as one gene: use one sequence per gene, e.g. the longest
#' isoform (see \code{\link{retrieve_longest_isoforms}}).
#' @param method \code{"rbh"} (Default) or \code{"deepclust"}, see Description.
#' @param seq_type the sequences in \code{files}: \code{"cds"} (Default; coding sequences, translated before the
#' search) or \code{"protein"}.
#' @param best_hit_rule how the best hit of a gene in another species is chosen among its hits: \code{"bitscore"}
#' (Default) or \code{"evalue"}; see \code{\link{diamond_rec}}. The code that chooses is the one \code{diamond_rec()} uses.
#' @param eval the e-value threshold. \code{method = "rbh"}: a hit counts when its e-value \emph{in a search against
#' the proteome of the subject species alone} is not above \code{eval} (see Details); Default \code{"1E-5"}.
#' \code{method = "deepclust"}: handed to \code{\link{deepclust}}.
#' @param sensitivity_mode the DIAMOND sensitivity mode of the search (\code{method = "rbh"}, \code{refine = TRUE}),
#' see \code{\link{diamond}}. Default is \code{"fast"}, the default of \code{\link{diamond_rec}}; for distantly
#' related species use \code{"sensitive"} or higher.
#' @param min_qcov,min_identity optional filters on the hits (Default \code{NULL}: off): hits whose alignment covers
#' less than \code{min_qcov} percent of the query, or with a percent identity below \code{min_identity}, are left out
#' \emph{before} best hits are chosen. As a reciprocal best hit needs a hit in both directions, \code{min_qcov} is
#' in effect a coverage of both proteins.
#' @param inparalogs a logical value. \code{FALSE} (Default): groups are built from reciprocal best hits only.
#' \code{TRUE}: in addition, two genes of the same species are joined when they hit each other with a bit score
#' (the smaller one of both directions) higher than the best between-species hit of either gene - a duplication after
#' the speciation, in the manner of InParanoid. This is a heuristic on bit scores, not a gene-tree reconciliation.
#' @param refine \code{method = "deepclust"} only: split the clusters that hold more than one gene of a species by the
#' reciprocal best hits inside the cluster (one further DIAMOND search, of the members of these clusters only).
#' Default is \code{FALSE}.
#' @param realign \code{method = "deepclust"} only: run \code{\link{deepclust_realign}} to obtain alignment statistics
#' of every cluster member with its representative (reported in \code{pairs}). Default is \code{FALSE}.
#' @param max.target.seqs the largest number of hits DIAMOND reports per query in the all-against-all search
#' (\code{--max-target-seqs}), over all species together. Default is \code{10000}. A warning is given when a query
#' reaches this number, as its best hit in some species may then be missing.
#' @param pairs_max_group_size groups with more genes than this are left out of the table \code{pairs} (the number of
#' pairs grows with the square of the group size). Default is \code{500}.
#' @param delete_corrupt_cds a logical value: leave out coding sequences whose length is not a multiple of 3
#' (\code{seq_type = "cds"}). Default is \code{TRUE}, as in \code{\link{diamond_rec}}.
#' @param diamond_params further parameters for \code{diamond blastp} (one character string), e.g.
#' \code{"--masking 0"}. Default is \code{NULL}.
#' @param comp_cores the number of threads DIAMOND uses. The result does not depend on it.
#' @param path the folder that holds the DIAMOND executable (in case it is not on the PATH).
#' @param quiet a logical value: run DIAMOND quietly and show only the summary messages of this function.
#' @param ... \code{method = "deepclust"}: further arguments for \code{\link{deepclust}}, e.g. \code{approx_id},
#' \code{mutual_cover}, \code{deepclust_params}. Not used with \code{method = "rbh"}.
#' @details
#' \strong{One search instead of N x N searches.} With \code{method = "rbh"} the proteins of all species are written
#' to one file and searched against themselves once. Two things differ between such a search and separate searches of
#' every pair of species, and both are taken care of:
#' \itemize{
#' \item \emph{e-values} grow with the size of the database. The e-value of every hit is therefore converted to the
#' e-value of a search against the subject species alone (multiplied by the share of that species in the residues of
#' all species; DIAMOND's e-value is proportional to the database size), and \code{eval} is applied to the converted
#' value. Bit scores do not depend on the database. The orthologs of two species therefore do not change when further
#' species are added, and \code{orthogroups()} with two species returns the reciprocal best hits of
#' \code{\link{diamond_rec}} (verified in the tests of the package and, on human and mouse, on 41,330 proteins).
#' \item DIAMOND reports at most \code{max.target.seqs} hits per query. A gene with more better hits than that
#' (members of very large families in many species) can lose its best hit in a species; a warning counts such queries.
#' }
#' \strong{What a group is.} An orthogroup holds genes of at least two species. Genes that are in no group are
#' listed in \code{unassigned} (also the members of clusters or components with genes of one species only; OrthoFinder
#' reports those as species-specific orthogroups). With \code{method = "rbh"} two genes are in one group when a chain
#' of reciprocal best hits connects them; a group therefore can hold several genes of a species (co-orthologs) even
#' without \code{inparalogs}. Without \code{inparalogs}, a gene that was duplicated after the speciation is
#' represented by the copy that is the reciprocal best hit; the other copy is unassigned.
#'
#' \strong{Reproducibility.} The proteins are sorted by species name and ID before the search, and groups are numbered
#' by decreasing size and then by their first gene in that order: the result does not depend on the order of
#' \code{files} or on \code{comp_cores}.
#'
#' \strong{Sequence IDs.} The ID of a sequence is its FASTA header up to the first white space. IDs are replaced by
#' numbers during the search, so any ID works. When an ID occurs in more than one species, all genes are named
#' \code{<species>|<ID>} in the result and a message says so (\code{x$genes} holds both names).
#' @return An object of class \code{"orthogroups"}: a list with
#' \itemize{
#' \item \code{groups}: one row per gene in a group (\code{orthogroup}, \code{species}, \code{gene});
#' \item \code{wide}: one row per group and one column per species, genes separated by \code{", "} - the layout of
#' OrthoFinder's \code{Orthogroups.tsv};
#' \item \code{pairs}: one row per ordered pair of genes of different species in the same group
#' (\code{orthogroup}, \code{query_species}, \code{query_id}, \code{subject_species}, \code{subject_id},
#' \code{bit_score}, \code{evalue}, \code{perc_identity}, \code{alig_length}, \code{qcov}, \code{scov}, \code{rbh},
#' \code{relation}). The alignment statistics are those of the hit of the query with the subject (\code{NA} when
#' DIAMOND reported none, e.g. for two genes connected through a third); \code{rbh} says whether the two genes are
#' reciprocal best hits; \code{relation} is \code{"1:1"} when both species have exactly one gene in the group and
#' \code{"co-ortholog"} otherwise. Every pair is listed in both directions;
#' \item \code{single_copy}: the rows of \code{wide} of the groups with exactly one gene in every species;
#' \item \code{species}: the species names, in the order of \code{files};
#' \item \code{unassigned}: the genes in no group (\code{species}, \code{gene});
#' \item \code{genes}: all genes (\code{species}, \code{gene}, \code{id}: the ID in the file, \code{length}: protein
#' length, \code{orthogroup});
#' \item \code{params}: the parameters of the call, the DIAMOND version and the files.
#' }
#' See \code{\link{write_orthogroups}} for files in OrthoFinder's layout and \code{\link{dNdS_orthogroups}} for dN/dS
#' within the groups.
#' @author Chloe Newton
#' @examples \dontrun{
#' files <- c(
#'   Athaliana = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   Alyrata   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' og <- orthogroups(files)
#' og
#' summary(og)
#' og$wide
#' og$pairs
#'
#' # proteins of all species in a folder, DeepClust route, clusters refined by reciprocal best hits
#' og <- orthogroups("proteomes/", method = "deepclust", seq_type = "protein", refine = TRUE, comp_cores = 8)
#'
#' # files in the layout of OrthoFinder
#' write_orthogroups(og, "orthogroups_out")
#' }
#' @seealso \code{\link{write_orthogroups}}, \code{\link{dNdS_orthogroups}}, \code{\link{diamond_rec}},
#' \code{\link{deepclust}}, \code{\link{orthologs}}
#' @export
orthogroups <- function(files,
                        method               = c("rbh", "deepclust"),
                        seq_type             = c("cds", "protein"),
                        best_hit_rule        = getOption("orthologr.best_hit_rule", "bitscore"),
                        eval                 = "1E-5",
                        sensitivity_mode     = "fast",
                        min_qcov             = NULL,
                        min_identity         = NULL,
                        inparalogs           = FALSE,
                        refine               = FALSE,
                        realign              = FALSE,
                        max.target.seqs      = 10000,
                        pairs_max_group_size = 500,
                        delete_corrupt_cds   = TRUE,
                        diamond_params       = NULL,
                        comp_cores           = 1,
                        path                 = NULL,
                        quiet                = TRUE,
                        ...) {

        method   <- check_choice(method[1], c("rbh", "deepclust"), "method")
        seq_type <- check_choice(seq_type[1], c("cds", "protein"), "seq_type")
        best_hit_rule <- orthologr_best_hit_rule(best_hit_rule)
        check_sensitivity_mode(sensitivity_mode)
        eval <- check_evalue(eval)
        check_flag(inparalogs, "inparalogs")
        check_flag(refine, "refine")
        check_flag(realign, "realign")
        check_flag(quiet, "quiet")
        min_qcov     <- check_percent(min_qcov, "min_qcov")
        min_identity <- check_percent(min_identity, "min_identity")
        check_count(max.target.seqs, "max.target.seqs")
        check_count(pairs_max_group_size, "pairs_max_group_size")
        if (method == "rbh") {
                if (...length() > 0)
                        stop("Unused argument", if (...length() > 1) "s", " with method = \"rbh\": ",
                             paste(names(list(...)), collapse = ", "),
                             " ('...' is handed to deepclust() with method = \"deepclust\").", call. = FALSE)
                if (refine || realign)
                        stop("'refine' and 'realign' belong to method = \"deepclust\".", call. = FALSE)
        }
        if (method == "deepclust" && inparalogs && !refine)
                stop("inparalogs = TRUE needs hits between the genes: with method = \"deepclust\" please set refine = TRUE.",
                     call. = FALSE)

        input <- orthogroups_input(files)

        path <- diamond_resolve_path(path)
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        version <- diamond_version(path = path)
        comp_cores <- check_comp_cores(comp_cores)

        # one working folder for this call (removed when it ends), also used by deepclust() and deepclust_realign()
        end_tmpdir <- orthologr_tmpdir_start("orthologr_orthogroups_")
        on.exit(end_tmpdir(), add = TRUE)
        work <- call_dir_start("orthogroups_")

        prot <- orthogroups_read(input, seq_type = seq_type, delete_corrupt_cds = delete_corrupt_cds)
        genes <- prot$genes
        n_sp <- nrow(input)

        message("orthogroups(): ", nrow(genes), " genes of ", n_sp, " species (",
                paste0(input$species, ": ", as.vector(table(factor(genes$species, levels = input$species))), collapse = ", "),
                "), method = \"", method, "\", ", version)

        hits <- NULL
        rbh_edges <- NULL

        if (method == "rbh") {
                hits <- orthogroups_search(
                        genes = genes, aa = prot$aa, work = work, eval = eval, sensitivity_mode = sensitivity_mode,
                        max.target.seqs = max.target.seqs, diamond_params = diamond_params, comp_cores = comp_cores,
                        path = path, quiet = quiet, min_qcov = min_qcov, min_identity = min_identity
                )
                graph <- orthogroups_graph(hits, genes = genes, rule = best_hit_rule, inparalogs = inparalogs)
                label <- graph$label
                rbh_edges <- graph$rbh
        }

        if (method == "deepclust") {
                fasta <- file.path(work, "all_species.fasta")
                orthogroups_write_fasta(seq_len(nrow(genes)), prot$aa, fasta)
                run <- function(expr) if (quiet) suppressMessages(expr) else expr
                clusters <- run(deepclust(input_file = fasta, seq_type = "protein", path = path, comp_cores = comp_cores,
                                          eval = as.numeric(eval), quiet = quiet, direct_input = TRUE, ...))
                rep_id <- as.integer(clusters$representative_id)
                mem_id <- as.integer(clusters$member_id)
                label <- seq_len(nrow(genes))
                label[mem_id] <- rep_id

                if (realign) {
                        re <- run(deepclust_realign(input_file = fasta, clusters = clusters, seq_type = "protein",
                                                    path = path, comp_cores = comp_cores, quiet = quiet))
                        hits <- orthogroups_realign_hits(re)
                }

                if (refine) {
                        refined <- orthogroups_refine(
                                label = label, genes = genes, aa = prot$aa, work = work, eval = eval,
                                sensitivity_mode = sensitivity_mode, max.target.seqs = max.target.seqs,
                                diamond_params = diamond_params, comp_cores = comp_cores, path = path, quiet = quiet,
                                min_qcov = min_qcov, min_identity = min_identity, rule = best_hit_rule,
                                inparalogs = inparalogs
                        )
                        label <- refined$label
                        rbh_edges <- refined$rbh
                        # the hits of the refined clusters come from the search; realign statistics are kept for
                        # the pairs the search does not hold
                        hits <- orthogroups_merge_hits(refined$hits, hits)
                }
        }

        res <- orthogroups_assemble(
                genes = genes, label = label, hits = hits, rbh = rbh_edges, species = input$species,
                pairs_max_group_size = pairs_max_group_size
        )
        res$params <- list(
                method = method, seq_type = seq_type, best_hit_rule = best_hit_rule, eval = eval,
                sensitivity_mode = sensitivity_mode, min_qcov = min_qcov, min_identity = min_identity,
                inparalogs = inparalogs, refine = refine, realign = realign, max.target.seqs = max.target.seqs,
                pairs_max_group_size = pairs_max_group_size, delete_corrupt_cds = delete_corrupt_cds,
                diamond_params = diamond_params, deepclust_args = if (method == "deepclust") list(...) else NULL,
                diamond_version = version, files = stats::setNames(input$file, input$species),
                ids_prefixed = prot$ids_prefixed
        )
        class(res) <- "orthogroups"

        message("orthogroups(): ", nrow(res$wide), " orthogroups with ", nrow(res$groups), " genes (",
                nrow(res$single_copy), " single-copy in all ", n_sp, " species), ", nrow(res$unassigned),
                " genes unassigned.")
        res
}

# ---- argument checks ------------------------------------------------------------------------------------------

check_flag <- function(x, arg) {
        if (!is.logical(x) || length(x) != 1 || is.na(x))
                stop(arg, " must be TRUE or FALSE.", call. = FALSE)
        invisible(x)
}

check_percent <- function(x, arg) {
        if (is.null(x))
                return(NULL)
        if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 100)
                stop(arg, " must be NULL (no filter) or one number between 0 and 100 (percent).", call. = FALSE)
        x
}

check_count <- function(x, arg) {
        if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1)
                stop(arg, " must be one number >= 1.", call. = FALSE)
        invisible(x)
}

# ---- input ----------------------------------------------------------------------------------------------------

# The species and their files: data.frame(species, file), in the order given.
orthogroups_input <- function(files) {
        if (is.list(files)) {
                if (!all(vapply(files, function(f) is.character(f) && length(f) == 1, logical(1))))
                        stop("'files' must hold one file path per species.", call. = FALSE)
                files <- unlist(files)
        }
        if (!is.character(files) || length(files) == 0 || anyNA(files))
                stop("'files' must be a (named) character vector or list of FASTA files, one per species, or the path of a folder.",
                     call. = FALSE)
        if (length(files) == 1 && dir.exists(files)) {
                found <- list.files(files, pattern = "\\.(fa|faa|fas|fasta|fna|ffn|pep|cds|aa)(\\.gz)?$",
                                    full.names = TRUE, ignore.case = TRUE)
                if (length(found) == 0)
                        stop("The folder '", files, "' holds no FASTA files (*.fa, *.faa, *.fasta, *.fna, *.pep, *.cds, optionally .gz).",
                             call. = FALSE)
                files <- sort(found, method = "radix")
        }
        missing_files <- files[!file.exists(files) | dir.exists(files)]
        if (length(missing_files) > 0)
                stop("The following input file(s) were not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
        species <- names(files)
        from_file <- sub("\\.[^.]*$", "", sub("\\.gz$", "", basename(files), ignore.case = TRUE))
        if (is.null(species))
                species <- from_file
        species[is.na(species) | !nzchar(species)] <- from_file[is.na(species) | !nzchar(species)]
        if (length(files) < 2)
                stop("orthogroups() needs the files of at least two species (", length(files), " given).", call. = FALSE)
        if (anyDuplicated(species))
                stop("The species names must be unique, but '", paste(unique(species[duplicated(species)]), collapse = "', '"),
                     "' occur", if (sum(duplicated(species)) == 1) "s", " more than once. ",
                     "Please name the files: files = c(species_1 = \"...\", species_2 = \"...\").", call. = FALSE)
        if (any(grepl("[\t\r\n]", species)))
                stop("Species names must not hold tabs or line breaks.", call. = FALSE)
        data.frame(species = unname(species), file = unname(files), stringsAsFactors = FALSE)
}

# Read the proteins of all species. Returns
#   genes : data.table(species, id, gene, length), sorted by species and id (C locale); the row number of a gene is
#           its number in the search
#   aa    : the protein sequences, in the same order
orthogroups_read <- function(input, seq_type, delete_corrupt_cds = TRUE) {
        geneids <- aa <- species <- id <- NULL
        end_seq_cache <- seq_cache_start()
        on.exit(end_seq_cache(), add = TRUE)
        tabs <- lapply(seq_len(nrow(input)), function(i) {
                dt <- protein_table(file = input$file[i], seq_type = seq_type, delete_corrupt_cds = delete_corrupt_cds)
                dt <- stats::na.omit(dt[ , list(geneids, aa)])
                if (nrow(dt) == 0)
                        stop("No usable sequences were read from '", input$file[i], "' (species '", input$species[i], "').",
                             call. = FALSE)
                data.table::data.table(species = input$species[i], id = as.character(dt$geneids), aa = as.character(dt$aa))
        })
        all <- data.table::rbindlist(tabs)
        data.table::setorder(all, species, id)
        ids_prefixed <- anyDuplicated(all$id) > 0
        if (ids_prefixed) {
                n_dup <- length(unique(all$id[duplicated(all$id)]))
                message(n_dup, " sequence ID", if (n_dup > 1) "s occur" else " occurs", " in more than one species (e.g. '",
                        all$id[duplicated(all$id)][1], "'): all genes are named '<species>|<ID>' in the result.")
        }
        genes <- data.table::data.table(
                species = all$species,
                id      = all$id,
                gene    = if (ids_prefixed) paste0(all$species, "|", all$id) else all$id,
                length  = nchar(all$aa)
        )
        list(genes = genes, aa = all$aa, ids_prefixed = ids_prefixed)
}

orthogroups_write_fasta <- function(idx, aa, file) {
        writeLines(as.vector(rbind(paste0(">", idx), aa[idx])), file)
        invisible(file)
}

# ---- search ---------------------------------------------------------------------------------------------------

# One all-against-all DIAMOND search of the genes `subset` (default: all). Returns the hits as a data.table:
#   q, s (gene numbers), perc_identity, alig_length, qcov, scov, evalue, bit_score
# without self hits. The e-values are those of a search against the subject species alone (see ?orthogroups) and
# hits above `eval` are left out, as are hits below min_qcov / min_identity.
orthogroups_search <- function(genes, aa, work, eval, sensitivity_mode, max.target.seqs, diamond_params, comp_cores,
                               path, quiet, min_qcov = NULL, min_identity = NULL, subset = NULL) {
        q <- s <- evalue <- qcov <- perc_identity <- NULL
        if (is.null(subset))
                subset <- seq_len(nrow(genes))
        dir <- tempfile("search_", tmpdir = work)
        dir.create(dir)
        on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
        fasta  <- file.path(dir, "all_species.fasta")
        output <- file.path(dir, "hits.tsv")
        orthogroups_write_fasta(subset, aa, fasta)
        diamond_makedb(in_file = fasta, db_file = fasta, path = path, comp_cores = comp_cores, quiet = quiet)

        # residues per species in this database: the e-value of a hit is proportional to the size of the database
        sp <- genes$species[subset]
        letters_sp <- tapply(as.numeric(genes$length[subset]), sp, sum)
        letters_all <- sum(letters_sp)
        eval_num <- as.numeric(eval)
        eval_search <- eval_num * letters_all / min(letters_sp)

        message("Searching ", length(subset), " proteins against each other with DIAMOND (", sensitivity_mode, ") ...")
        t0 <- proc.time()[["elapsed"]]
        res <- run_tool(
                "diamond",
                args = list("blastp", "--db", fasta, "--query", fasta,
                            "--evalue", format(eval_search, scientific = TRUE, digits = 15),
                            if (sensitivity_mode != "default") paste0("--", sensitivity_mode),
                            "--max-target-seqs", format(max.target.seqs, scientific = FALSE),
                            "--out", output, "--threads", comp_cores,
                            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qcovhsp", "scovhsp", "evalue",
                            "bitscore",
                            if (!is.null(diamond_params)) I(diamond_params),
                            if (quiet) "--quiet"),
                path = path, tool = "diamond blastp", context = "the all-against-all search of orthogroups()",
                error = FALSE, echo = !quiet
        )
        if (!identical(as.integer(res$status), 0L) || !file.exists(output)) {
                said <- utils::tail(res$stderr[nzchar(trimws(res$stderr))], 5)
                stop("diamond blastp failed in the all-against-all search of orthogroups(): exit status ", res$status, ".",
                     if (length(said) > 0) paste0("\ndiamond said: ", paste(said, collapse = "\n")) else "", call. = FALSE)
        }
        message("DIAMOND search completed in ", format_elapsed(proc.time()[["elapsed"]] - t0))

        cols <- c("q", "s", "perc_identity", "alig_length", "qcov", "scov", "evalue", "bit_score")
        classes <- c("integer", "integer", "double", "integer", "double", "double", "double", "double")
        if (isTRUE(file.size(output) == 0)) {
                hits <- data.table::as.data.table(stats::setNames(lapply(classes, vector, length = 0L), cols))
        } else {
                hits <- data.table::fread(output, sep = "\t", header = FALSE, quote = "", colClasses = classes,
                                          col.names = cols, showProgress = FALSE)
        }

        # queries that reached the limit of reported hits: their best hit in some species may be missing
        if (nrow(hits) > 0) {
                n_per_query <- hits[ , list(n = .N), by = q]
                full <- sum(n_per_query$n >= max.target.seqs)
                if (full > 0)
                        warning(full, " of ", length(subset), " proteins have max.target.seqs = ",
                                format(max.target.seqs, scientific = FALSE), " hits or more: DIAMOND reported only the best ",
                                format(max.target.seqs, scientific = FALSE), ", so that the best hit of such a protein in some species ",
                                "may be missing. Please raise max.target.seqs.", call. = FALSE)
        }

        hits <- hits[q != s]
        # the e-value of a search against the subject species alone, in the three significant digits DIAMOND prints
        scale <- as.vector(letters_sp[genes$species[hits$s]]) / letters_all
        data.table::set(hits, j = "evalue", value = signif(hits$evalue * scale, 3))
        hits <- hits[evalue <= eval_num]
        if (!is.null(min_qcov))
                hits <- hits[qcov >= min_qcov]
        if (!is.null(min_identity))
                hits <- hits[perc_identity >= min_identity]
        hits
}

# ---- graph ----------------------------------------------------------------------------------------------------

# The best hit of every gene in every other species: data.table(q, s), chosen by best_hits() - the code
# diamond_best() and diamond_rec() use - within the hits of a gene in one species.
orthogroups_best_hits <- function(hits, genes, rule) {
        q <- s <- grp <- top <- bit_score <- evalue <- query_id <- subject_id <- perc_identity <- alig_length <- NULL
        sp <- as.integer(factor(genes$species))
        b <- hits[sp[q] != sp[s]]
        if (nrow(b) == 0)
                return(data.table::data.table(q = integer(0), s = integer(0)))
        # one group per query and subject species
        b[ , grp := (as.numeric(q) - 1) * max(sp) + sp[s]]
        # The best hit is among the hits with the best value of the first criterion of the rule: only those are
        # handed to best_hits() (with all criteria and the IDs for the last tie-break).
        if (rule == "bitscore") {
                b[ , top := max(bit_score), by = grp]
                b <- b[bit_score == top]
        } else {
                b[ , top := min(evalue), by = grp]
                b <- b[evalue == top]
        }
        b[ , query_id := genes$gene[q]]
        b[ , subject_id := genes$gene[s]]
        # (the "evalue" rule takes the first row of a table sorted by query_id and subject_id)
        data.table::setorder(b, query_id, subject_id)
        best <- best_hits(as.data.frame(b[ , list(query_id, subject_id, perc_identity, alig_length, evalue, bit_score, grp, q, s)]),
                          group_col = "grp", rule = rule)
        data.table::data.table(q = best$q, s = best$s)
}

# Reciprocal best hits (q < s) from the best hits.
orthogroups_rbh <- function(best) {
        q <- s <- NULL
        if (nrow(best) == 0)
                return(best)
        rev <- data.table::data.table(q = best$s, s = best$q)
        both <- data.table::fintersect(best, rev)
        both <- both[q < s]
        data.table::setorder(both, q, s)
        both
}

# In-paralog edges (q < s): two genes of one species that hit each other with a bit score (the smaller one of both
# directions) above the best between-species bit score of either gene.
orthogroups_inparalogs <- function(hits, genes) {
        q <- s <- bit_score <- rev_score <- NULL
        n <- nrow(genes)
        sp <- as.integer(factor(genes$species))
        between <- hits[sp[q] != sp[s], list(best = max(bit_score)), by = q]
        best_between <- numeric(n)
        best_between[between$q] <- between$best
        within <- hits[sp[q] == sp[s], list(q, s, bit_score)]
        if (nrow(within) == 0)
                return(data.table::data.table(q = integer(0), s = integer(0)))
        rev <- data.table::data.table(q = within$s, s = within$q, rev_score = within$bit_score)
        within <- merge(within, rev, by = c("q", "s"))
        within <- within[q < s & pmin(bit_score, rev_score) > pmax(best_between[q], best_between[s]), list(q, s)]
        data.table::setorder(within, q, s)
        within
}

# Connected components of a graph with the nodes 1..n: the label of a node is the smallest node of its component.
connected_components <- function(n, from, to) {
        node <- m <- NULL
        label <- seq_len(n)
        if (length(from) == 0)
                return(label)
        repeat {
                low <- pmin(label[from], label[to])
                upd <- data.table::data.table(node = c(from, to), m = c(low, low))[ , list(m = min(m)), by = node]
                new <- label
                new[upd$node] <- pmin(new[upd$node], upd$m)
                # (the label of a node is a node of the same component: follow it)
                new <- new[new]
                if (identical(new, label))
                        break
                label <- new
        }
        label
}

# Groups from hits: reciprocal best hits per pair of species (and in-paralogs) as edges, components as groups.
orthogroups_graph <- function(hits, genes, rule, inparalogs = FALSE) {
        best <- orthogroups_best_hits(hits, genes, rule)
        rbh <- orthogroups_rbh(best)
        edges <- rbh
        if (inparalogs)
                edges <- rbind(edges, orthogroups_inparalogs(hits, genes))
        list(label = connected_components(nrow(genes), edges$q, edges$s), rbh = rbh)
}

# ---- DeepClust route ------------------------------------------------------------------------------------------

# The table of deepclust_realign() as hits (representative = query, member = subject, and the reverse).
orthogroups_realign_hits <- function(re) {
        q <- s <- NULL
        a <- data.table::data.table(
                q = as.integer(re$representative_id), s = as.integer(re$member_id),
                perc_identity = as.numeric(re$approx_pident), alig_length = NA_integer_,
                qcov = as.numeric(re$qcovhsp), scov = as.numeric(re$scovhsp),
                evalue = as.numeric(re$evalue), bit_score = as.numeric(re$bitscore)
        )
        a <- a[q != s]
        b <- data.table::copy(a)
        data.table::setnames(b, c("q", "s", "qcov", "scov"), c("s", "q", "scov", "qcov"))
        data.table::setcolorder(b, names(a))
        unique(rbind(a, b), by = c("q", "s"))
}

orthogroups_merge_hits <- function(first, second) {
        if (is.null(second))
                return(first)
        if (is.null(first))
                return(second)
        unique(rbind(first, second), by = c("q", "s"))
}

# Split the clusters that hold more than one gene of a species (and genes of at least two species) by the reciprocal
# best hits inside the cluster. One search of the members of these clusters; hits between clusters are not used.
orthogroups_refine <- function(label, genes, aa, work, rule, inparalogs, ...) {
        q <- s <- lab <- species <- n <- n_species <- max_copies <- NULL
        cl <- data.table::data.table(lab = label, species = genes$species)
        per <- cl[ , list(n = .N), by = list(lab, species)]
        info <- per[ , list(n_species = .N, max_copies = max(n)), by = lab]
        split_labels <- info[n_species >= 2 & max_copies >= 2]$lab
        todo <- which(is.element(label, split_labels))
        if (length(todo) == 0 || length(unique(genes$species[todo])) < 2)
                return(list(label = label, rbh = NULL, hits = NULL))
        message("Refining ", length(split_labels), " clusters with more than one gene of a species (", length(todo),
                " proteins) by reciprocal best hits ...")
        hits <- orthogroups_search(genes = genes, aa = aa, work = work, subset = todo, ...)
        # the search numbers the proteins as the full table does (the FASTA IDs are the row numbers)
        hits <- hits[label[q] == label[s]]
        graph <- orthogroups_graph(hits, genes = genes, rule = rule, inparalogs = inparalogs)
        label[todo] <- graph$label[todo]
        list(label = label, rbh = graph$rbh, hits = hits)
}

# ---- result ---------------------------------------------------------------------------------------------------

# The result tables from the group label of every gene.
orthogroups_assemble <- function(genes, label, hits, rbh, species, pairs_max_group_size = 500) {
        lab <- sp <- idx <- orthogroup <- gene <- n <- n_species <- first <- NULL
        n_sp <- length(species)
        g <- data.table::data.table(idx = seq_len(nrow(genes)), lab = label,
                                    sp = match(genes$species, species), gene = genes$gene)
        info <- g[ , list(n = .N, n_species = data.table::uniqueN(sp), first = min(idx)), by = lab]
        # a group has genes of at least two species
        info <- info[n_species >= 2]
        data.table::setorder(info, -n, first)
        info[ , orthogroup := sprintf("OG%07d", seq_len(.N) - 1L)]
        g[ , orthogroup := info$orthogroup[match(lab, info$lab)]]

        in_group <- g[!is.na(orthogroup)]
        data.table::setorder(in_group, orthogroup, sp, gene)
        groups <- data.table::data.table(orthogroup = in_group$orthogroup, species = species[in_group$sp],
                                         gene = in_group$gene)

        # one row per group, one column per species
        wide <- data.table::data.table(orthogroup = info$orthogroup)
        cells <- in_group[ , list(genes = paste(gene, collapse = ", ")), by = list(orthogroup, sp)]
        counts <- matrix(0L, nrow = nrow(info), ncol = n_sp)
        cnt <- in_group[ , list(n = .N), by = list(orthogroup, sp)]
        counts[cbind(match(cnt$orthogroup, info$orthogroup), cnt$sp)] <- cnt$n
        for (k in seq_len(n_sp)) {
                col <- rep("", nrow(info))
                ck <- cells[sp == k]
                col[match(ck$orthogroup, info$orthogroup)] <- ck$genes
                data.table::set(wide, j = species[k], value = col)
        }
        single <- if (nrow(info) > 0) rowSums(counts == 1L) == n_sp else logical(0)

        unassigned <- g[is.na(orthogroup)]
        data.table::setorder(unassigned, sp, gene)
        unassigned <- data.table::data.table(species = species[unassigned$sp], gene = unassigned$gene)

        pairs <- orthogroups_pairs(in_group, hits = hits, rbh = rbh, species = species, counts = counts,
                                   group_ids = info$orthogroup, pairs_max_group_size = pairs_max_group_size)

        gene_tbl <- data.table::data.table(species = genes$species, gene = genes$gene, id = genes$id,
                                           length = genes$length, orthogroup = g$orthogroup)
        gene_tbl <- gene_tbl[order(match(gene_tbl$species, species), gene_tbl$gene, method = "radix")]

        list(
                groups      = tibble::as_tibble(groups),
                wide        = tibble::as_tibble(wide),
                pairs       = tibble::as_tibble(pairs),
                single_copy = tibble::as_tibble(wide[single]),
                species     = species,
                unassigned  = tibble::as_tibble(unassigned),
                genes       = tibble::as_tibble(gene_tbl)
        )
}

# All ordered pairs of genes of different species within a group, with the statistics of the hit query -> subject.
orthogroups_pairs <- function(in_group, hits, rbh, species, counts, group_ids, pairs_max_group_size) {
        orthogroup <- n <- sp <- sp.q <- sp.s <- idx.q <- idx.s <- gene.q <- gene.s <- idx <- gene <- NULL
        empty <- data.table::data.table(
                orthogroup = character(0), query_species = character(0), query_id = character(0),
                subject_species = character(0), subject_id = character(0), bit_score = numeric(0), evalue = numeric(0),
                perc_identity = numeric(0), alig_length = integer(0), qcov = numeric(0), scov = numeric(0),
                rbh = logical(0), relation = character(0)
        )
        if (nrow(in_group) == 0)
                return(empty)
        size <- in_group[ , list(n = .N), by = orthogroup]
        large <- size[n > pairs_max_group_size]$orthogroup
        if (length(large) > 0)
                message(length(large), " group", if (length(large) > 1) "s" else "", " with more than ", pairs_max_group_size,
                        " genes ", if (length(large) > 1) "are" else "is", " left out of the table 'pairs' (pairs_max_group_size).")
        m <- in_group[!is.element(orthogroup, large), list(orthogroup, idx, sp, gene)]
        if (nrow(m) == 0)
                return(empty)
        p <- merge(m, m, by = "orthogroup", allow.cartesian = TRUE, suffixes = c(".q", ".s"))
        p <- p[sp.q != sp.s]
        if (nrow(p) == 0)
                return(empty)
        stat_cols <- c("bit_score", "evalue", "perc_identity", "alig_length", "qcov", "scov")
        if (!is.null(hits) && nrow(hits) > 0) {
                at <- hits[p, on = c(q = "idx.q", s = "idx.s"), which = TRUE, mult = "first"]
                for (col in stat_cols)
                        data.table::set(p, j = col, value = hits[[col]][at])
        } else {
                for (col in stat_cols)
                        data.table::set(p, j = col, value = if (col == "alig_length") NA_integer_ else NA_real_)
        }
        if (!is.null(rbh) && nrow(rbh) > 0) {
                key <- paste(pmin(p$idx.q, p$idx.s), pmax(p$idx.q, p$idx.s))
                data.table::set(p, j = "rbh", value = is.element(key, paste(rbh$q, rbh$s)))
        } else {
                data.table::set(p, j = "rbh", value = if (is.null(rbh)) NA else FALSE)
        }
        data.table::setorder(p, orthogroup, sp.q, gene.q, sp.s, gene.s)
        one <- counts[cbind(match(p$orthogroup, group_ids), p$sp.q)] == 1L &
                counts[cbind(match(p$orthogroup, group_ids), p$sp.s)] == 1L
        data.table::data.table(
                orthogroup = p$orthogroup, query_species = species[p$sp.q], query_id = p$gene.q,
                subject_species = species[p$sp.s], subject_id = p$gene.s,
                bit_score = p$bit_score, evalue = p$evalue, perc_identity = p$perc_identity,
                alig_length = as.integer(p$alig_length), qcov = p$qcov, scov = p$scov, rbh = p$rbh,
                relation = ifelse(one, "1:1", "co-ortholog")
        )
}
