# One transcript -> gene table for GFF3 and GTF annotations.
#
# retrieve_longest_isoforms(), select_orthologs() and generate_ortholog_tables() need to know which gene locus a
# transcript or protein of a FASTA file belongs to. They used to answer this with format-specific code that only
# fitted one kind of file each (GFF3: grouping CDS rows by the gene *symbol*, or rows by 'Name'; GTF: literal
# 'transcript_id' matches on rows that carry 'gene_biotype'). The table built here is the same for all of them:
#
#   GFF3 : the Parent / ID relations are followed from every feature up to its top-level feature (the gene).
#          CDS -> mRNA -> gene (eukaryotes), CDS -> gene (NCBI prokaryotes). The locus is the ID of the top-level
#          feature without the prefixes 'gene-' (NCBI) or 'gene:' (Ensembl); for NCBI files that is the locus_tag.
#          It is never the gene symbol ('gene=' / 'Name='): unnamed genes have none and paralogs can share one.
#   GTF  : 'gene_id' is the locus; 'transcript_id' and 'protein_id' (with 'transcript_version' /
#          'protein_version' when given) name the transcript and the protein.
#
# The annotation is read as text in blocks (connections read .gz files too) and only the feature type and the
# attributes that are needed are kept, so memory use stays far below that of importing the whole file as an
# rtracklayer / GenomicRanges object (see NEWS for numbers).

# Value of attribute `key` in a GFF3 attribute column ('key=value;key=value'); NA where it is missing.
gff_attr <- function(attr, key) {
        m <- regexpr(paste0("(?:^|;)\\s*", key, "=[^;]*"), attr, perl = TRUE)
        res <- rep(NA_character_, length(attr))
        hit <- m > 0
        val <- regmatches(attr, m)
        val <- sub("^;?\\s*[^=]*=", "", val)
        # percent-encoded characters that occur in IDs
        if (any(grepl("%", val, fixed = TRUE)))
                val <- vapply(val, utils::URLdecode, "", USE.NAMES = FALSE)
        res[hit] <- val
        res
}

# Value of attribute `key` in a GTF attribute column ('key "value"; key "value";'); NA where it is missing or "".
gtf_attr <- function(attr, key) {
        m <- regexpr(paste0("(?:^|;)\\s*", key, "\\s+\"[^\"]*\""), attr, perl = TRUE)
        res <- rep(NA_character_, length(attr))
        hit <- m > 0
        val <- regmatches(attr, m)
        val <- sub("\"$", "", sub("^[^\"]*\"", "", val))
        res[hit] <- val
        res[!is.na(res) & res == ""] <- NA_character_
        res
}

# Read the lines of an annotation file in blocks and hand each block (type, attributes) to `fun`; returns the
# list of what `fun` returned. Comment lines are skipped; a '##FASTA' line ends the annotation.
#   skip_types : feature types that are not needed (their lines are dropped before anything is parsed)
annotation_blocks <- function(file, fun, skip_types = character(0), n = 200000L) {
        con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, open = "r") else file(file, open = "r")
        on.exit(close(con), add = TRUE)
        out <- list()
        repeat {
                x <- readLines(con, n = n, warn = FALSE)
                if (length(x) == 0)
                        break
                fasta <- which(startsWith(x, "##FASTA"))
                done <- length(fasta) > 0
                if (done)
                        x <- x[seq_len(fasta[1] - 1L)]
                x <- x[!startsWith(x, "#") & nzchar(x)]
                for (ty in skip_types)
                        x <- x[!grepl(paste0("\t", ty, "\t"), x, fixed = TRUE)]
                if (length(x) > 0) {
                        type <- sub("^(?:[^\t]*\t){2}([^\t]*)\t.*$", "\\1", x, perl = TRUE)
                        attr <- sub("^(?:[^\t]*\t){8}", "", x, perl = TRUE)
                        src  <- sub("^[^\t]*\t([^\t]*)\t.*$", "\\1", x, perl = TRUE)
                        keep <- src != "ERCC"
                        res <- fun(type[keep], attr[keep])
                        if (!is.null(res) && nrow(res) > 0)
                                out[[length(out) + 1L]] <- res
                }
                if (done)
                        break
        }
        out
}

# Feature types that hold no ID relation that is needed (most lines of an annotation).
annotation_skip_types <- function() {
        c("exon", "five_prime_UTR", "three_prime_UTR", "UTR", "start_codon", "stop_codon", "region",
          "biological_region", "chromosome", "scaffold", "supercontig", "intron", "Selenocysteine",
          "sequence_feature", "repeat_region", "mobile_genetic_element", "origin_of_replication", "cDNA_match", "match")
}

# 'gene-b0001' -> 'b0001', 'gene:ENSG01' -> 'ENSG01', 'rna-NM_1.1' -> 'NM_1.1', 'transcript:ENST01' -> 'ENST01',
# 'cds-NP_1.1' -> 'NP_1.1', 'CDS:ENSP01' -> 'ENSP01'
strip_id_prefix <- function(x) {
        sub("^(gene|rna|cds|id|transcript|CDS|mRNA|exon)[-:]", "", x, perl = TRUE)
}

# 'ENST00000373020.9' -> 'ENST00000373020'
strip_id_version <- function(x) {
        sub("\\.[0-9]+$", "", x)
}

# The transcript -> gene table of an annotation file.
# Returns a data.table with one row per transcript / protein:
#   gene_id       : the gene locus (see above)
#   transcript_id : the transcript (NA for a CDS that hangs directly on its gene)
#   protein_id    : the protein (NA for transcripts without CDS)
#   gene_biotype  : biotype of the gene when the annotation gives one, else NA
#   has_gene_row / has_cds : whether the locus has a 'gene' row / the transcript has CDS rows
# and attr "keys": data.table(key, gene_id, kind) - every spelling under which a transcript or protein of the
# table may appear in a FASTA header (kind: "transcript" or "protein"), see annotation_keys().
annotation_table <- function(file, format = c("gff", "gtf")) {
        format <- match.arg(format)
        if (!file.exists(file))
                stop("The annotation file '", file, "' does not exist.", call. = FALSE)
        tbl <- if (format == "gtf") annotation_table_gtf(file) else annotation_table_gff(file)
        if (nrow(tbl) == 0)
                stop("No transcripts or coding sequences were found in the annotation file '", file, "' (read as ",
                     toupper(format), if (format == "gff") "3" else "", "). Please check the file and the format.",
                     call. = FALSE)
        attr(tbl, "keys") <- annotation_keys(tbl)
        tbl
}

annotation_table_gtf <- function(file) {
        gene_id <- transcript_id <- protein_id <- gene_biotype <- type <- has_cds <- NULL
        blocks <- annotation_blocks(file, function(type, attr) {
                keep <- is.element(type, c("gene", "transcript", "mRNA", "CDS"))
                if (!any(keep))
                        return(NULL)
                type <- type[keep]
                attr <- attr[keep]
                tv <- gtf_attr(attr, "transcript_version")
                pv <- gtf_attr(attr, "protein_version")
                unique(data.table::data.table(
                        type               = ifelse(type == "mRNA", "transcript", type),
                        gene_id            = gtf_attr(attr, "gene_id"),
                        transcript_id      = gtf_attr(attr, "transcript_id"),
                        protein_id         = gtf_attr(attr, "protein_id"),
                        transcript_version = tv,
                        protein_version    = pv,
                        gene_biotype       = gtf_attr(attr, "gene_biotype")
                ))
        }, skip_types = annotation_skip_types())
        empty <- data.table::data.table(gene_id = character(0), transcript_id = character(0), protein_id = character(0),
                                        gene_biotype = character(0), has_gene_row = logical(0), has_cds = logical(0),
                                        transcript_version = character(0), protein_version = character(0))
        if (length(blocks) == 0)
                return(empty)
        x <- unique(data.table::rbindlist(blocks))
        x <- x[!is.na(gene_id)]

        genes <- x[type == "gene"]
        biotype <- x[!is.na(gene_biotype), list(gene_biotype = gene_biotype[1]), by = gene_id]

        tx <- x[type != "gene" & !is.na(transcript_id)]
        if (nrow(tx) == 0)
                return(empty)
        first_value <- function(v) {
                v <- v[!is.na(v)]
                if (length(v) == 0) NA_character_ else v[1]
        }
        res <- tx[ , list(
                protein_id         = first_value(protein_id),
                transcript_version = first_value(transcript_version),
                protein_version    = first_value(protein_version),
                has_cds            = any(type == "CDS")
        ), by = list(gene_id, transcript_id)]
        res[ , gene_biotype := biotype$gene_biotype[match(gene_id, biotype$gene_id)]]
        res[ , has_gene_row := is.element(gene_id, genes$gene_id)]
        data.table::setcolorder(res, c("gene_id", "transcript_id", "protein_id", "gene_biotype", "has_gene_row", "has_cds",
                                       "transcript_version", "protein_version"))
        res[]
}

annotation_table_gff <- function(file) {
        ID <- Parent <- type <- protein_id <- transcript_id <- gene_id <- has_cds <- NULL
        blocks <- annotation_blocks(file, function(type, attr) {
                id <- gff_attr(attr, "ID")
                parent <- gff_attr(attr, "Parent")
                keep <- !is.na(id) | !is.na(parent)
                if (!any(keep))
                        return(NULL)
                attr <- attr[keep]
                bt <- gff_attr(attr, "gene_biotype")
                bt2 <- gff_attr(attr, "biotype")
                unique(data.table::data.table(
                        type          = type[keep],
                        ID            = id[keep],
                        # a feature with several parents (rare outside exons): the first one
                        Parent        = sub(",.*$", "", parent[keep]),
                        protein_id    = gff_attr(attr, "protein_id"),
                        transcript_id = gff_attr(attr, "transcript_id"),
                        version       = gff_attr(attr, "version"),
                        Name          = gff_attr(attr, "Name"),
                        biotype       = ifelse(is.na(bt), bt2, bt)
                ))
        }, skip_types = annotation_skip_types())
        empty <- data.table::data.table(gene_id = character(0), transcript_id = character(0), protein_id = character(0),
                                        gene_biotype = character(0), has_gene_row = logical(0), has_cds = logical(0),
                                        transcript_version = character(0), protein_version = character(0),
                                        transcript_raw = character(0), protein_raw = character(0))
        if (length(blocks) == 0)
                return(empty)
        x <- unique(data.table::rbindlist(blocks))

        # every feature with an ID, once (a CDS has one row per exon, all with the same ID and Parent)
        feat <- x[!is.na(ID)]
        feat <- feat[!duplicated(ID)]

        # follow Parent up to the top-level feature
        root_of <- function(start) {
                cur <- start
                for (k in 1:10) {
                        p <- feat$Parent[match(cur, feat$ID)]
                        up <- !is.na(p) & is.element(p, feat$ID)
                        if (!any(up))
                                break
                        cur[up] <- p[up]
                }
                cur
        }

        cds <- x[type == "CDS"]
        cds <- cds[!duplicated(data.frame(ID, Parent, protein_id))]
        cds <- cds[!is.na(Parent) & is.element(Parent, feat$ID)]

        # transcripts: the features between a gene and its CDS (mRNA, transcript, lnc_RNA, V_gene_segment, tRNA, ...),
        # i.e. every feature other than a CDS that has a parent (exons and UTRs were not read)
        tx <- feat[type != "CDS" & !is.na(Parent) & is.element(Parent, feat$ID)]

        rows <- list()
        if (nrow(cds) > 0) {
                parent_is_tx <- is.element(cds$Parent, tx$ID)
                prot <- ifelse(is.na(cds$protein_id), strip_id_prefix(ifelse(is.na(cds$ID), cds$Name, cds$ID)), cds$protein_id)
                ti <- match(cds$Parent, tx$ID)
                rows[[1]] <- data.table::data.table(
                        gene_raw           = root_of(cds$Parent),
                        transcript_raw     = ifelse(parent_is_tx, cds$Parent, NA_character_),
                        transcript_id      = ifelse(parent_is_tx,
                                                    ifelse(is.na(tx$transcript_id[ti]), strip_id_prefix(cds$Parent), tx$transcript_id[ti]),
                                                    NA_character_),
                        transcript_version = ifelse(parent_is_tx, tx$version[ti], NA_character_),
                        protein_raw        = cds$ID,
                        protein_id         = prot,
                        protein_version    = cds$version,
                        has_cds            = TRUE
                )
        }
        no_cds <- tx[!is.element(ID, cds$Parent)]
        if (nrow(no_cds) > 0)
                rows[[length(rows) + 1L]] <- data.table::data.table(
                        gene_raw           = root_of(no_cds$ID),
                        transcript_raw     = no_cds$ID,
                        transcript_id      = ifelse(is.na(no_cds$transcript_id), strip_id_prefix(no_cds$ID), no_cds$transcript_id),
                        transcript_version = no_cds$version,
                        protein_raw        = NA_character_,
                        protein_id         = NA_character_,
                        protein_version    = NA_character_,
                        has_cds            = FALSE
                )
        if (length(rows) == 0)
                return(empty)
        res <- unique(data.table::rbindlist(rows, use.names = TRUE))
        gi <- match(res$gene_raw, feat$ID)
        res[ , gene_id := strip_id_prefix(gene_raw)]
        res[ , gene_biotype := feat$biotype[gi]]
        res[ , has_gene_row := is.element(feat$type[gi], c("gene", "pseudogene", "ncRNA_gene", "transposable_element_gene"))]
        res[ , gene_raw := NULL]
        data.table::setcolorder(res, c("gene_id", "transcript_id", "protein_id", "gene_biotype", "has_gene_row", "has_cds",
                                       "transcript_version", "protein_version", "transcript_raw", "protein_raw"))
        res[]
}

# Every spelling under which a transcript or protein of an annotation table may appear as the ID of a FASTA file:
# the ID itself, the ID with its version ('ENST00000373020' + '.9' when the annotation gives the version in a field
# of its own) and the raw GFF3 ID ('rna-NM_001.1', 'transcript:ENST...').
annotation_keys <- function(tbl) {
        key <- gene_id <- NULL
        add <- function(id, kind) {
                ok <- !is.na(id)
                # (built as a data.frame: 'key' is an argument of data.table())
                data.table::as.data.table(data.frame(key = id[ok], gene_id = tbl$gene_id[ok], kind = rep(kind, sum(ok)),
                                                     stringsAsFactors = FALSE))
        }
        versioned <- function(id, version) ifelse(is.na(id) | is.na(version), NA_character_, paste0(id, ".", version))
        keys <- list(
                add(tbl$transcript_id, "transcript"),
                add(versioned(tbl$transcript_id, tbl$transcript_version), "transcript"),
                add(tbl$protein_id, "protein"),
                add(versioned(tbl$protein_id, tbl$protein_version), "protein")
        )
        if (is.element("transcript_raw", names(tbl)))
                keys <- c(keys, list(add(tbl$transcript_raw, "transcript"), add(tbl$protein_raw, "protein")))
        unique(data.table::rbindlist(keys), by = c("key", "gene_id"))
}

# The gene loci of sequence IDs (FASTA header IDs or the IDs of a dN/dS table).
#   ids  : character vector
#   keys : attr(annotation_table(), "keys"), optionally restricted to one kind
# An ID is looked up as it is; an ID that is not found is looked up again (a) as the protein accession inside an
# NCBI 'cds_from_genomic' ID ('lcl|NC_000913.3_cds_NP_414542.1_1'), and (b) without its version on both sides
# ('ENST00000373020.9' = 'ENST00000373020'; only where the version-less ID belongs to one locus).
# Returns a data.table(id, gene_id, how) with one row per ID and locus (an ID can belong to several loci: identical
# proteins of NCBI prokaryote annotations, 'WP_...'); IDs that were not found are not in it.
match_ids_to_genes <- function(ids, keys) {
        ids <- unique(as.character(ids[!is.na(ids)]))
        key_tbl <- data.frame(k = keys$key, gene_id = keys$gene_id, stringsAsFactors = FALSE)
        # one row per (id, gene) for the IDs whose spelling `k_of_id` is a key
        expand <- function(id, k_of_id, how, tbl = key_tbl) {
                m <- merge(data.frame(id = id, k = k_of_id, stringsAsFactors = FALSE), tbl, by = "k")
                if (nrow(m) == 0)
                        return(NULL)
                unique(data.frame(id = m$id, gene_id = m$gene_id, how = how, stringsAsFactors = FALSE))
        }
        ncbi_accession <- function(x) sub("^lcl\\|.*_(?:cds|prot)_(.+)_[0-9]+$", "\\1", x, perl = TRUE)
        res <- list(expand(ids, ids, "exact"))
        left <- setdiff(ids, res[[1]]$id)
        if (length(left) > 0) {
                acc <- ncbi_accession(left)
                changed <- acc != left
                if (any(changed)) {
                        r <- expand(left[changed], acc[changed], "ncbi_cds_id")
                        res[[length(res) + 1L]] <- r
                        left <- setdiff(left, r$id)
                }
        }
        if (length(left) > 0) {
                # without versions on both sides, where that is unambiguous
                nv <- unique(data.frame(k = strip_id_version(key_tbl$k), gene_id = key_tbl$gene_id, stringsAsFactors = FALSE))
                n_genes <- table(nv$k)
                nv <- nv[is.element(nv$k, names(n_genes)[n_genes == 1]), ]
                res[[length(res) + 1L]] <- expand(left, strip_id_version(strip_id_prefix(ncbi_accession(left))), "without_version", nv)
        }
        res <- data.table::rbindlist(res)
        if (nrow(res) == 0)
                return(data.table::data.table(id = character(0), gene_id = character(0), how = character(0)))
        # in the order of `ids`
        res <- res[order(match(res$id, ids))]
        res
}

# The error for "no ID of the sequence file is in the annotation": examples of both sides, so that the user sees
# at once what does not fit.
annotation_mismatch_error <- function(ids, tbl, seq_what, annotation_file) {
        ex_t <- utils::head(unique(stats::na.omit(tbl$transcript_id)), 3)
        ex_p <- utils::head(unique(stats::na.omit(tbl$protein_id)), 3)
        stop(
                "None of the ", length(unique(ids)), " IDs of ", seq_what, " was found in the annotation file '",
                basename(annotation_file), "'. IDs of ", seq_what, ", e.g.: ", paste(utils::head(unique(ids), 3), collapse = ", "),
                ". Transcript IDs of the annotation, e.g.: ", if (length(ex_t)) paste(ex_t, collapse = ", ") else "(none)",
                "; protein IDs, e.g.: ", if (length(ex_p)) paste(ex_p, collapse = ", ") else "(none)",
                ". (IDs are compared as they are, then without version suffix; NCBI 'lcl|..._cds_<protein>_<n>' IDs by their protein accession.) ",
                "Please check that the sequence file and the annotation file belong to the same genome release and that the annotation format is right.",
                call. = FALSE
        )
}
