# Does a sequence file hold several transcripts (splice variants) per gene?
#
# dNdS(), orthologs() and divergence_stratigraphy() treat every sequence of a file as a gene. A file with all
# transcripts of every gene (e.g. the 'cds/all' files of Ensembl, or biomartr::getCDS()) gives several rows per
# gene, with different partners and different dN/dS values, and best hits between splice variants of the same
# pair of genes compete with each other. The check below is a heuristic and only warns; it never changes a
# result. It is deliberately conservative (it prefers to say nothing over a false alarm):
#
#   1. Header fields. When at least 90% of the FASTA headers name their gene - Ensembl 'gene:<ID>', NCBI
#      '[locus_tag=<ID>]' or, without locus tags, '[gene=<symbol>]' - the genes are counted from these fields.
#      Several transcripts per gene are reported when a gene occurs in more than one header.
#   2. ID patterns (only when the headers name no genes). A transcript suffix is removed from the end of every ID:
#      '.<number>' (AT1G01010.2), '.t<number>' / '-t<number>' / '-T<number>' (AUGUSTUS, BRAKER),
#      '-R<letter>' / '-P<letter>' (FlyBase, VectorBase, MAKER), '-mRNA-<number>' (MAKER), '_i<number>'
#      (Trinity), '.p<number>' (TransDecoder). Several transcripts per gene are reported only when at least 5
#      of the remaining stems, and at least 1% of them, occur more than once. (IDs with a version, e.g.
#      'ENST00000373020.9', keep unique stems and are not reported.)
#
# options(orthologr.isoform_warning = FALSE) switches the check off.

# The gene of every sequence of `file` (named character vector: sequence ID -> gene), or NULL when neither the
# headers nor the IDs say anything about genes. attr "source" is "header" or "id".
isoform_map <- function(file) {
        headers <- tryCatch(as.character(Biostrings::fasta.index(file)$desc), error = function(e) NULL)
        if (is.null(headers) || length(headers) < 2)
                return(NULL)
        ids <- fasta_header_id(headers)

        field <- function(pattern) {
                has <- grepl(pattern, headers, perl = TRUE)
                if (mean(has) < 0.9)
                        return(NULL)
                gene <- rep(NA_character_, length(headers))
                gene[has] <- sub(paste0("^.*?", pattern, ".*$"), "\\1", headers[has], perl = TRUE)
                gene
        }
        gene <- field("\\sgene:(\\S+)")
        if (!is.null(gene))
                gene <- sub("\\.[0-9]+$", "", gene)
        if (is.null(gene))
                gene <- field("\\[locus_tag=([^\\]]+)\\]")
        if (is.null(gene))
                gene <- field("\\[gene=([^\\]]+)\\]")
        if (!is.null(gene)) {
                gene[is.na(gene)] <- ids[is.na(gene)]
                names(gene) <- ids
                attr(gene, "source") <- "header"
                return(gene)
        }

        stem <- sub("(\\.[0-9]+|[.-][tT][0-9]+|-R[A-Z]|-P[A-Z]|-mRNA-[0-9]+|_i[0-9]+|\\.p[0-9]+)$", "", ids, perl = TRUE)
        n_multi <- sum(table(stem) > 1)
        if (n_multi < 5 || n_multi < 0.01 * length(unique(stem)))
                return(NULL)
        names(stem) <- ids
        attr(stem, "source") <- "id"
        stem
}

# Warn when `file` seems to hold several transcripts per gene. Returns the ID -> gene vector (invisibly) when it
# does, otherwise NULL.
isoform_warning <- function(file, arg = "file") {
        if (!isTRUE(getOption("orthologr.isoform_warning", TRUE)))
                return(invisible(NULL))
        gene <- tryCatch(isoform_map(file), error = function(e) NULL)
        if (is.null(gene))
                return(invisible(NULL))
        per_gene <- table(gene)
        n_multi <- sum(per_gene > 1)
        if (n_multi == 0)
                return(invisible(NULL))
        ex <- names(per_gene)[which.max(per_gene)]
        warning(
                arg, " '", basename(file), "' seems to hold several transcripts (splice variants) per gene: ",
                length(gene), " sequences for ", length(per_gene), " genes, ", n_multi, " genes with more than one sequence ",
                "(up to ", max(per_gene), ", e.g. ", paste(utils::head(names(gene)[gene == ex], 3), collapse = ", "), "; ",
                if (identical(attr(gene, "source"), "header")) "from the gene names in the FASTA headers" else "judged by the sequence IDs",
                "). Every sequence is treated as a gene: the result will hold several rows per gene, with different ",
                "partners and dN/dS values. Please keep one transcript per gene, e.g. with retrieve_longest_isoforms(). ",
                "(options(orthologr.isoform_warning = FALSE) switches this check off.)",
                call. = FALSE
        )
        invisible(gene)
}

# After the run: how many rows the result holds per query gene (a message; only when isoform_warning() found
# several transcripts per gene in the query file).
isoform_result_note <- function(res, gene) {
        if (is.null(gene) || nrow(res) == 0 || !is.element("query_id", names(res)))
                return(invisible(NULL))
        g <- gene[as.character(res$query_id)]
        g <- g[!is.na(g)]
        if (length(g) == 0)
                return(invisible(NULL))
        per_gene <- table(g)
        if (any(per_gene > 1))
                message("The result holds ", length(g), " rows for ", length(per_gene), " query genes: ",
                        sum(per_gene > 1), " genes have more than one row (up to ", max(per_gene),
                        "), because the query file holds several transcripts per gene.")
        invisible(NULL)
}
