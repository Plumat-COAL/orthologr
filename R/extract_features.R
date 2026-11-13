#' @title Helper function to extract gene loci and splice variant IDs from GFF files
#' @description Extract gene loci and splice variant IDs from GFF files.
#' @param x a data.frame converted from a \code{gff} or \code{gtf} file.
#' @param format either \code{format = "gtf"} or \code{format = "gff"}.
#' @export
extract_features <- function(x, format) {
        type <- NULL

        if (!is.element(format, c("gtf", "gff")))
                stop("Please choose a format that is supported by this function: format = 'gtf' or format = 'gff'.", call. = FALSE)
        
        if (format == "gtf") {
                if (!all(c("gene", "CDS") %in% names(table(x$type)))) {
                        stop("Please make sure that the 'type' column in your gtf file are labeled as 'gene' and 'CDS' to denote genes and CDS.", call. = FALSE)
                }
                        
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "CDS")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$gene_id, length(unique(get_mRNA_id$transcript_id))),
                                              query_id = unique(get_mRNA_id$transcript_id))
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        
        if (format == "gff") {
                if (!all(c("gene", "mRNA") %in% names(table(x$type))))
                        stop("Please make sure that the 'type' column in your gtf file are labeled as  'gene' and 'mRNA' to denote genes and mRNA.", call. = FALSE)
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "mRNA")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$ID, length(unique(get_mRNA_id$ID))),
                                              query_id = unique(get_mRNA_id$ID))
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        
        gene_locus_id <- query_id <- NULL
        res <- dplyr::filter(res, !is.na(gene_locus_id), !is.na(query_id))
        return(res)
}

# extract_features() for all gene loci of an annotation at once: what
# dplyr::ungroup(dplyr::do(dplyr::group_by(x, <group_col>), extract_features(., format))) returns (columns:
# <group_col>, gene_locus_id, query_id; loci in the order of dplyr's groups, transcripts in file order), without
# calling an R function for every locus. Annotations in which a locus has more than one 'gene' row, or whose
# 'type' column is not a factor holding both feature types, take the per-locus path.
extract_features_all <- function(x, format, group_col) {
        if (!is.element(format, c("gtf", "gff")))
                stop("Please choose a format that is supported by this function: format = 'gtf' or format = 'gff'.", call. = FALSE)
        
        t_type <- if (format == "gtf") "CDS" else "mRNA"
        id_col <- if (format == "gtf") "gene_id" else "ID"
        t_col  <- if (format == "gtf") "transcript_id" else "ID"
        
        simple <- nrow(x) > 0 && is.factor(x$type) && all(c("gene", t_type) %in% levels(x$type)) &&
                all(c(group_col, id_col, t_col) %in% names(x))
        if (simple) {
                grp <- x[[group_col]]
                g <- NULL
                keys <- dplyr::group_keys(dplyr::group_by(tibble::tibble(g = grp), g))$g
                gi <- match(grp, keys)
                is_gene <- !is.na(x$type) & x$type == "gene"
                is_t    <- !is.na(x$type) & x$type == t_type
                n_gene <- tabulate(gi[is_gene], nbins = length(keys))
                simple <- all(n_gene <= 1)
        }
        if (!simple)
                return(extract_features_do(x, format, group_col))
        
        # the id of the 'gene' row of every locus (NA: no gene row - the locus has no result)
        gene_id <- rep(NA_character_, length(keys))
        gene_id[gi[is_gene]] <- as.character(x[[id_col]][is_gene])
        
        # the transcripts of every locus: unique, in the order of their first row
        t_g  <- gi[is_t]
        t_id <- as.character(x[[t_col]][is_t])
        keep <- !duplicated(data.frame(t_g, t_id))
        t_g <- t_g[keep]
        t_id <- t_id[keep]
        ord <- order(t_g, method = "radix")
        t_g <- t_g[ord]
        t_id <- t_id[ord]
        
        res <- tibble::tibble(group = keys[t_g], gene_locus_id = gene_id[t_g], query_id = t_id)
        res <- res[!is.na(res$gene_locus_id) & !is.na(res$query_id), ]
        names(res)[1] <- group_col
        res
}

# The implementation used until 0.4.3 (one call of extract_features() per gene locus).
extract_features_do <- function(x, format, group_col) {
        '.' <- NULL
        dplyr::ungroup(dplyr::do(dplyr::group_by(x, .data[[group_col]]), extract_features(. , format = format)))
}
