#' @title Print and summarise orthogroups
#' @description \code{print()} shows the number of groups and genes per species; \code{summary()} returns these
#' numbers as tables: genes per species (in groups, in single-copy groups, unassigned), the distribution of the
#' group sizes and of the number of species per group.
#' @param x,object an object of class \code{"orthogroups"} (see \code{\link{orthogroups}}).
#' @param ... not used.
#' @return \code{print()}: \code{x}, invisibly. \code{summary()}: an object of class \code{"summary.orthogroups"}, a
#' list with \code{n_species}, \code{n_genes}, \code{n_groups}, \code{n_single_copy}, \code{n_unassigned},
#' \code{per_species} (data.frame: \code{species}, \code{genes}, \code{in_groups}, \code{percent_in_groups},
#' \code{groups}, \code{unassigned}), \code{group_sizes} (data.frame: \code{genes}, \code{groups}) and
#' \code{species_per_group} (data.frame: \code{species}, \code{groups}).
#' @author Chloe Newton
#' @seealso \code{\link{orthogroups}}
#' @export
print.orthogroups <- function(x, ...) {
        s <- summary(x)
        cat("Orthogroups of ", s$n_species, " species (method = \"", x$params$method, "\"",
            if (isTRUE(x$params$refine)) ", refined by reciprocal best hits" else "",
            if (isTRUE(x$params$inparalogs)) ", with in-paralogs" else "", ")\n", sep = "")
        cat("  ", s$n_groups, " groups with ", sum(s$per_species$in_groups), " of ", s$n_genes, " genes; ",
            s$n_single_copy, " single-copy in all species; ", s$n_unassigned, " genes unassigned\n", sep = "")
        print(s$per_species, row.names = FALSE)
        cat("  tables: $groups, $wide, $pairs (", nrow(x$pairs), " rows), $single_copy, $unassigned, $genes\n", sep = "")
        invisible(x)
}

#' @rdname print.orthogroups
#' @export
summary.orthogroups <- function(object, ...) {
        x <- object
        sp <- x$species
        genes_sp  <- as.vector(table(factor(x$genes$species, levels = sp)))
        in_groups <- as.vector(table(factor(x$groups$species, levels = sp)))
        groups_sp <- vapply(sp, function(k) sum(nzchar(x$wide[[k]])), numeric(1))
        size <- as.vector(table(x$groups$orthogroup))
        n_species_group <- if (nrow(x$wide) > 0) rowSums(as.matrix(x$wide[ , sp, drop = FALSE]) != "") else integer(0)
        size_tab <- table(size)
        species_tab <- table(n_species_group)
        res <- list(
                n_species     = length(sp),
                n_genes       = nrow(x$genes),
                n_groups      = nrow(x$wide),
                n_single_copy = nrow(x$single_copy),
                n_unassigned  = nrow(x$unassigned),
                per_species   = data.frame(
                        species           = sp,
                        genes             = genes_sp,
                        in_groups         = in_groups,
                        percent_in_groups = round(100 * in_groups / pmax(genes_sp, 1), 1),
                        groups            = as.vector(groups_sp),
                        unassigned        = genes_sp - in_groups,
                        stringsAsFactors  = FALSE
                ),
                group_sizes = data.frame(genes = as.integer(names(size_tab)), groups = as.vector(size_tab)),
                species_per_group = data.frame(species = as.integer(names(species_tab)), groups = as.vector(species_tab)),
                method = x$params$method
        )
        class(res) <- "summary.orthogroups"
        res
}

#' @rdname print.orthogroups
#' @export
print.summary.orthogroups <- function(x, ...) {
        cat("Orthogroups of ", x$n_species, " species: ", x$n_groups, " groups, ", x$n_single_copy,
            " single-copy in all species, ", x$n_unassigned, " of ", x$n_genes, " genes unassigned\n\n", sep = "")
        cat("Genes per species:\n")
        print(x$per_species, row.names = FALSE)
        cat("\nGroups by number of genes:\n")
        gs <- x$group_sizes
        if (nrow(gs) > 12) {
                rest <- gs[gs$genes > gs$genes[11], , drop = FALSE]
                gs <- rbind(data.frame(genes = as.character(gs$genes[1:11]), groups = gs$groups[1:11]),
                            data.frame(genes = paste0(">", x$group_sizes$genes[11], " (max ", max(rest$genes), ")"),
                                       groups = sum(rest$groups)))
        }
        print(gs, row.names = FALSE)
        cat("\nGroups by number of species:\n")
        print(x$species_per_group, row.names = FALSE)
        invisible(x)
}

#' @title Orthogroups as a data.frame
#' @description The long table of an \code{"orthogroups"} object: one row per gene in a group.
#' @param x an object of class \code{"orthogroups"} (see \code{\link{orthogroups}}).
#' @param row.names,optional not used.
#' @param unassigned a logical value: add the genes in no group, with \code{orthogroup = NA}. Default is \code{FALSE}.
#' @param ... not used.
#' @return A data.frame with the columns \code{orthogroup}, \code{species} and \code{gene}.
#' @author Chloe Newton
#' @export
as.data.frame.orthogroups <- function(x, row.names = NULL, optional = FALSE, unassigned = FALSE, ...) {
        res <- as.data.frame(x$groups, stringsAsFactors = FALSE)
        if (unassigned && nrow(x$unassigned) > 0)
                res <- rbind(res, data.frame(orthogroup = NA_character_, species = x$unassigned$species,
                                             gene = x$unassigned$gene, stringsAsFactors = FALSE))
        res
}

#' @title Write orthogroups to files in the layout of OrthoFinder
#' @description Writes the tables of an \code{"orthogroups"} object as tab-separated files with the names and
#' layout of OrthoFinder's result folder, so that programs that read OrthoFinder output - and the importers of this
#' package - can read them.
#' @param x an object of class \code{"orthogroups"} (see \code{\link{orthogroups}}).
#' @param dir the folder to write to (created when it does not exist).
#' @param ortholog_tables a logical value: also write, for the species \code{query_species}, one table per other
#' species with its reciprocal best hits in the format of \code{\link{generate_ortholog_tables_all}}, for
#' \code{\link{import_ortholog_tables_all}} and \code{\link{retrieve_core_orthologs}}. Default is \code{FALSE}.
#' @param query_species the query species of these tables. Default is the first species.
#' @details The files:
#' \itemize{
#' \item \code{<dir>/Orthogroups/Orthogroups.tsv}: one row per group, one column per species, genes separated by
#' \code{", "} (first column \code{Orthogroup});
#' \item \code{<dir>/Orthogroups/Orthogroups_SingleCopyOrthologues.txt}: the IDs of the groups with exactly one gene
#' in every species, one per line;
#' \item \code{<dir>/Orthogroups/Orthogroups_UnassignedGenes.tsv}: the genes in no group, one per row in the column of
#' their species, numbered on from the last group as OrthoFinder does;
#' \item \code{<dir>/Orthogroups/Orthogroups.GeneCount.tsv}: genes per group and species, and their total;
#' \item \code{<dir>/pairs.tsv}: the table \code{x$pairs};
#' \item with \code{ortholog_tables = TRUE}: \code{<dir>/ortholog_tables/orthologs_by_gene_<query>_vs_<subject>.csv}
#' (separated by \code{";"}); dN, dS and the columns DIAMOND was not asked for are \code{NA}.
#' }
#' \code{orthofinder2_retrieve_core_orthologs(file.path(dir, "Orthogroups", "Orthogroups.tsv"),
#' file.path(dir, "Orthogroups", "Orthogroups_SingleCopyOrthologues.txt"))} returns \code{x$single_copy}.
#' @return invisibly, a named character vector with the paths of the files.
#' @author Chloe Newton
#' @examples \dontrun{
#' og <- orthogroups(c(
#'   Athaliana = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   Alyrata   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr')))
#' out <- write_orthogroups(og, file.path(tempdir(), "og"))
#' orthofinder2_retrieve_core_orthologs(out[["orthogroups"]], out[["single_copy"]])
#' }
#' @seealso \code{\link{orthogroups}}, \code{\link{orthofinder2_retrieve_core_orthologs}}
#' @export
write_orthogroups <- function(x, dir, ortholog_tables = FALSE, query_species = NULL) {
        if (!inherits(x, "orthogroups"))
                stop("'x' must be the result of orthogroups().", call. = FALSE)
        if (!is.character(dir) || length(dir) != 1 || is.na(dir) || !nzchar(dir))
                stop("'dir' must be the path of one folder.", call. = FALSE)
        check_flag(ortholog_tables, "ortholog_tables")
        og_dir <- file.path(dir, "Orthogroups")
        if (!dir.exists(og_dir) && !dir.create(og_dir, recursive = TRUE))
                stop("The folder '", og_dir, "' could not be created.", call. = FALSE)
        sp <- x$species
        write_tsv <- function(tbl, file) {
                data.table::fwrite(tbl, file, sep = "\t", quote = FALSE, na = "NA")
                file
        }
        
        wide <- as.data.frame(x$wide, stringsAsFactors = FALSE, check.names = FALSE)
        names(wide)[1] <- "Orthogroup"
        files <- c(orthogroups = write_tsv(wide, file.path(og_dir, "Orthogroups.tsv")))
        
        single_file <- file.path(og_dir, "Orthogroups_SingleCopyOrthologues.txt")
        writeLines(as.character(x$single_copy$orthogroup), single_file)
        files <- c(files, single_copy = single_file)
        
        # unassigned genes: one row each, numbered on from the last group
        un <- x$unassigned
        un_tbl <- data.frame(Orthogroup = sprintf("OG%07d", nrow(wide) + seq_len(nrow(un)) - 1L), stringsAsFactors = FALSE)
        for (k in sp)
                un_tbl[[k]] <- ifelse(un$species == k, un$gene, "")
        files <- c(files, unassigned = write_tsv(un_tbl, file.path(og_dir, "Orthogroups_UnassignedGenes.tsv")))
        
        counts <- data.frame(Orthogroup = wide$Orthogroup, stringsAsFactors = FALSE)
        tab <- table(factor(x$groups$orthogroup, levels = wide$Orthogroup), factor(x$groups$species, levels = sp))
        for (k in sp)
                counts[[k]] <- as.vector(tab[ , k])
        counts$Total <- as.vector(rowSums(tab))
        files <- c(files, gene_count = write_tsv(counts, file.path(og_dir, "Orthogroups.GeneCount.tsv")))
        
        files <- c(files, pairs = write_tsv(as.data.frame(x$pairs), file.path(dir, "pairs.tsv")))
        
        if (ortholog_tables) {
                if (is.null(query_species))
                        query_species <- sp[1]
                if (!is.character(query_species) || length(query_species) != 1 || !is.element(query_species, sp))
                        stop("query_species must be one of the species: '", paste(sp, collapse = "', '"), "'.", call. = FALSE)
                tbl_dir <- file.path(dir, "ortholog_tables")
                if (!dir.exists(tbl_dir))
                        dir.create(tbl_dir, recursive = TRUE)
                p <- as.data.frame(x$pairs, stringsAsFactors = FALSE)
                p <- p[p$query_species == query_species & !is.na(p$rbh) & p$rbh, , drop = FALSE]
                len <- stats::setNames(x$genes$length, x$genes$gene)
                for (k in setdiff(sp, query_species)) {
                        pk <- p[p$subject_species == k, , drop = FALSE]
                        out <- data.frame(
                                query_species = pk$query_species, subject_species = pk$subject_species,
                                query_id = pk$query_id, query_gene_locus_id = pk$query_id,
                                subject_id = pk$subject_id, subject_gene_locus_id = pk$subject_id,
                                dN = NA_real_, dS = NA_real_, dNdS = NA_real_,
                                evalue = pk$evalue, bit_score = pk$bit_score, perc_identity = pk$perc_identity,
                                num_ident_matches = NA_real_, alig_length = pk$alig_length, mismatches = NA_integer_,
                                gap_openings = NA_integer_, n_gaps = NA_real_, pos_match = NA_real_, ppos = NA_real_,
                                q_start = NA_integer_, q_end = NA_integer_, q_len = as.numeric(len[pk$query_id]),
                                qcov = pk$qcov, qcovhsp = pk$qcov, s_start = NA_integer_, s_end = NA_integer_,
                                s_len = as.numeric(len[pk$subject_id]), score_raw = NA_real_,
                                stringsAsFactors = FALSE
                        )
                        f <- file.path(tbl_dir, paste0("orthologs_by_gene_", safe_file_name(query_species), "_vs_",
                                                       safe_file_name(k), ".csv"))
                        data.table::fwrite(out, f, sep = ";", quote = FALSE, na = "NA")
                        files <- c(files, stats::setNames(f, paste0("ortholog_table_", k)))
                }
        }
        invisible(files)
}
