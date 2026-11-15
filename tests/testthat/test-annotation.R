# retrieve_longest_isoforms(), select_orthologs() and generate_ortholog_tables() on pieces of real annotations
# (NCBI RefSeq GFF3 / GTF, Ensembl GFF3 / GTF; see fixtures/README). Second audit: N03, N06, N08, N23, N40.

fx <- function(f) test_path("fixtures", "annotation", f)

fasta_headers <- function(f) names(Biostrings::readBStringSet(f))
header_field <- function(h, key) sub(paste0("^.* ", key, ":(\\S+).*$"), "\\1", h)

# what is expected for the Ensembl records: the longest sequence of every gene, from the FASTA headers alone
ensembl_truth <- function(fasta, protein_coding_only) {
        x <- Biostrings::readBStringSet(fasta)
        d <- data.frame(id = sub(" .*", "", names(x)), gene = sub("\\..*", "", header_field(names(x), "gene")),
                        biotype = header_field(names(x), "gene_biotype"), width = Biostrings::width(x), stringsAsFactors = FALSE)
        if (protein_coding_only)
                d <- d[d$biotype == "protein_coding", ]
        unname(vapply(split(d, d$gene), function(g) g$id[which.max(g$width)], ""))
}

mk_dnds <- function(q, s) {
        n <- length(q)
        tibble::tibble(query_id = q, subject_id = s, dN = 0.1, dS = 0.5, dNdS = 0.2, perc_identity = 90, num_ident_matches = 90L,
                       alig_length = 100L + seq_len(n), mismatches = 10L, gap_openings = 0L, n_gaps = 0L, pos_match = 95L, ppos = 95,
                       q_start = 1L, q_end = 100L, q_len = 100L, qcov = 100, qcovhsp = 100, s_start = 1L, s_end = 100L, s_len = 100L,
                       evalue = 10^-(seq_len(n) + 5), bit_score = 100 + seq_len(n), score_raw = 200 + seq_len(n))
}

# ---------------------------------------------------------------- reading annotations

test_that("gff_attr() and gtf_attr() read attribute values ...", {
        a <- c("ID=cds-NP_1.1;Parent=gene-b0001;Name=NP_1.1;gene=thrL;protein_id=NP_1.1", "ID=gene-b0002;Name=x%2Cy", "Parent=rna-1")
        expect_identical(gff_attr(a, "ID"), c("cds-NP_1.1", "gene-b0002", NA))
        expect_identical(gff_attr(a, "Parent"), c("gene-b0001", NA, "rna-1"))
        expect_identical(gff_attr(a, "Name"), c("NP_1.1", "x,y", NA))
        expect_identical(gff_attr(a, "protein_id"), c("NP_1.1", NA, NA))
        # 'gene' is not found inside 'gene_biotype' or 'Parent=gene-...'
        expect_identical(gff_attr("ID=x;gene_biotype=protein_coding", "gene"), NA_character_)
        g <- c('gene_id "b0001"; transcript_id ""; gene "thrL"; gene_biotype "protein_coding";',
               'gene_id "ENSG01"; gene_version "9"; transcript_id "ENST01"; transcript_version "8"; protein_id "ENSP01";')
        expect_identical(gtf_attr(g, "gene_id"), c("b0001", "ENSG01"))
        expect_identical(gtf_attr(g, "transcript_id"), c(NA, "ENST01"))
        expect_identical(gtf_attr(g, "transcript_version"), c(NA, "8"))
        expect_identical(gtf_attr(g, "gene"), c("thrL", NA))
        expect_identical(strip_id_prefix(c("gene-b0001", "gene:ENSG01", "rna-NM_1.1", "transcript:ENST01", "cds-NP_1.1", "CDS:ENSP01", "AT1G01010")),
                         c("b0001", "ENSG01", "NM_1.1", "ENST01", "NP_1.1", "ENSP01", "AT1G01010"))
})

test_that("annotation_table() follows Parent / ID in NCBI and Ensembl GFF3 and reads gene_id / transcript_id / protein_id of GTF ...", {
        sak <- annotation_table(fx("ncbi_ecoli_sakai.gff"), "gff")
        cds <- sak[sak$has_cds, ]
        expect_equal(nrow(cds), 39)
        expect_false(anyDuplicated(cds$gene_id) > 0)
        expect_true(all(grepl("^ECs_[0-9]{4}$", cds$gene_id)))
        expect_true(all(c("ECs_0015", "ECs_0684", "ECs_0016", "ECs_0621") %in% cds$gene_id))     # share a gene symbol
        expect_identical(cds$gene_id[cds$protein_id == "NP_308028.1"], "ECs_0001")
        
        k12 <- annotation_table(fx("ncbi_ecoli_k12.gtf"), "gtf")
        expect_equal(sum(k12$has_cds), 34)
        expect_identical(k12$gene_id[which(k12$protein_id == "NP_414542.1")], "b0001")
        expect_true(all(k12$gene_biotype[k12$has_cds] == "protein_coding"))
        
        gtf <- annotation_table(fx("ensembl_human_chr21.gtf"), "gtf")
        gff <- annotation_table(fx("ensembl_human_chr21.gff3"), "gff")
        expect_setequal(gtf$transcript_id, gff$transcript_id)
        expect_equal(length(unique(gtf$gene_id)), 9)
        expect_setequal(gtf$gene_id, gff$gene_id)
        m1 <- gtf[order(gtf$transcript_id), ]
        m2 <- gff[order(gff$transcript_id), ]
        expect_identical(m1$gene_id, m2$gene_id)
        expect_identical(m1$protein_id, m2$protein_id)
        expect_identical(m1$transcript_version, m2$transcript_version)
        # versioned spellings are keys of both
        v <- paste0(m1$transcript_id, ".", m1$transcript_version)
        expect_true(all(v %in% attr(gtf, "keys")$key))
        expect_true(all(v %in% attr(gff, "keys")$key))
        
        # .gz files are read as well
        gz <- tempfile(fileext = ".gff.gz")
        con <- gzfile(gz, "w"); writeLines(readLines(fx("ncbi_ecoli_sakai.gff")), con); close(con)
        expect_equal(as.data.frame(annotation_table(gz, "gff")), as.data.frame(sak), ignore_attr = TRUE)
        
        # the wrong format, or a file without transcripts, is an error that says so
        expect_error(annotation_table(fx("ncbi_ecoli_k12.gtf"), "gff"), "No transcripts or coding sequences were found")
})

test_that("match_ids_to_genes() matches IDs as they are, by NCBI CDS IDs, then without versions when that is unambiguous ...", {
        keys <- data.frame(key = c("ENST1", "ENST1.4", "NP_1.1", "AT1G01010.1", "AT1G01010.2", "AT1G01020.1", "WP_9.1", "WP_9.1"),
                                       gene_id = c("G1", "G1", "b0001", "AT1G01010", "AT1G01010", "AT1G01020", "b0100", "b0200"),
                                       kind = "transcript")
        m <- match_ids_to_genes(c("ENST1.4", "ENST1.9", "ENST1", "lcl|NC_1.3_cds_NP_1.1_1", "NP_1.2", "AT1G01010.3", "AT1G01020.7", "WP_9.1", "nope"), keys)
        get <- function(id) m$gene_id[m$id == id]
        expect_identical(get("ENST1.4"), "G1")
        expect_identical(get("ENST1.9"), "G1")                 # another version of the same transcript
        expect_identical(m$how[m$id == "ENST1.9"], "without_version")
        expect_identical(get("lcl|NC_1.3_cds_NP_1.1_1"), "b0001")
        expect_identical(get("NP_1.2"), "b0001")
        expect_identical(get("AT1G01010.3"), "AT1G01010")      # both '.1' and '.2' are AT1G01010: unambiguous
        expect_identical(get("AT1G01020.7"), "AT1G01020")
        expect_setequal(get("WP_9.1"), c("b0100", "b0200"))    # one protein, two loci
        expect_length(get("nope"), 0)
})

# ---------------------------------------------------------------- retrieve_longest_isoforms()

test_that("retrieve_longest_isoforms(): NCBI RefSeq GFF3 + protein.faa gives one protein per locus (N03) ...", {
        out <- tempfile(fileext = ".faa")
        msg <- testthat::capture_messages(res <- retrieve_longest_isoforms(fx("ncbi_ecoli_sakai_protein.faa"), fx("ncbi_ecoli_sakai.gff"), out, annotation_format = "gff"))
        ids <- sub(" .*", "", fasta_headers(out))
        expect_length(ids, 39)
        expect_setequal(ids, sub(" .*", "", fasta_headers(fx("ncbi_ecoli_sakai_protein.faa"))))
        expect_equal(attr(res, "n_loci"), 39)
        expect_match(paste(msg, collapse = ""), "Writing 39 unique longest peptide sequences \\(from initially 39 isoforms\\) for 39 gene loci")
        
        # the grouping by gene symbol of earlier versions merged the loci without a symbol and those that share one
        skip_if_not_installed("rtracklayer")
        old <- tempfile(fileext = ".faa")
        suppressMessages(retrieve_longest_isoforms_do(fx("ncbi_ecoli_sakai_protein.faa"), fx("ncbi_ecoli_sakai.gff"), old, annotation_format = "gff"))
        expect_lt(length(fasta_headers(old)), 39)
})

test_that("retrieve_longest_isoforms(): NCBI GTF + protein.faa (N06) ...", {
        out <- tempfile(fileext = ".faa")
        suppressMessages(retrieve_longest_isoforms(fx("ncbi_ecoli_k12_protein.faa"), fx("ncbi_ecoli_k12.gtf"), out, annotation_format = "gtf"))
        expect_setequal(sub(" .*", "", fasta_headers(out)), sub(" .*", "", fasta_headers(fx("ncbi_ecoli_k12_protein.faa"))))
        expect_length(fasta_headers(out), 34)
})

test_that("retrieve_longest_isoforms(): Ensembl GTF and GFF3 with pep and cds FASTA, versioned IDs (N06) ...", {
        pep <- fx("ensembl_human_chr21_pep.fa")
        cds <- fx("ensembl_human_chr21_cds.fa")
        run <- function(fasta, ann, fmt) {
                out <- tempfile(fileext = ".fa")
                suppressMessages(retrieve_longest_isoforms(fasta, ann, out, annotation_format = fmt))
                sub(" .*", "", fasta_headers(out))
        }
        # GTF: genes with gene_biotype 'protein_coding' (as before); 7 of the 8 genes with a protein
        expect_setequal(run(pep, fx("ensembl_human_chr21.gtf"), "gtf"), ensembl_truth(pep, protein_coding_only = TRUE))
        expect_length(run(pep, fx("ensembl_human_chr21.gtf"), "gtf"), 7)
        expect_message(retrieve_longest_isoforms(pep, fx("ensembl_human_chr21.gtf"), tempfile(), annotation_format = "gtf"),
                       "1 sequences of 1 gene loci with another biotype are left out \\(IG_V_gene: 1\\)")
        # GFF3: every gene with a protein (as before for GFF3: no biotype filter)
        expect_setequal(run(pep, fx("ensembl_human_chr21.gff3"), "gff"), ensembl_truth(pep, protein_coding_only = FALSE))
        expect_length(run(pep, fx("ensembl_human_chr21.gff3"), "gff"), 8)
        # coding sequences: the FASTA IDs are transcript IDs
        expect_setequal(run(cds, fx("ensembl_human_chr21.gtf"), "gtf"), ensembl_truth(cds, protein_coding_only = TRUE))
        expect_setequal(run(cds, fx("ensembl_human_chr21.gff3"), "gff"), ensembl_truth(cds, protein_coding_only = FALSE))
        
        # headers without versions work as well
        x <- Biostrings::readAAStringSet(pep)
        names(x) <- sub("\\.[0-9]+( |$)", "\\1", sub(" .*", "", names(x)))
        nov <- tempfile(fileext = ".fa")
        Biostrings::writeXStringSet(x, nov)
        expect_setequal(run(nov, fx("ensembl_human_chr21.gtf"), "gtf"), sub("\\.[0-9]+$", "", ensembl_truth(pep, protein_coding_only = TRUE)))
})

test_that("retrieve_longest_isoforms(): .gz input, sequences that fit no annotation, duplicated IDs ...", {
        gz_ann <- tempfile(fileext = ".gff.gz")
        con <- gzfile(gz_ann, "w"); writeLines(readLines(fx("ncbi_ecoli_sakai.gff")), con); close(con)
        gz_pep <- tempfile(fileext = ".faa.gz")
        con <- gzfile(gz_pep, "w"); writeLines(readLines(fx("ncbi_ecoli_sakai_protein.faa")), con); close(con)
        out <- tempfile(fileext = ".faa")
        suppressMessages(retrieve_longest_isoforms(gz_pep, gz_ann, out, annotation_format = "gff"))
        expect_length(fasta_headers(out), 39)
        
        # a proteome of another genome: an error with examples of both kinds of IDs
        expect_error(suppressMessages(retrieve_longest_isoforms(fx("ncbi_ecoli_k12_protein.faa"), fx("ncbi_ecoli_sakai.gff"), out, annotation_format = "gff")),
                     "None of the 34 IDs of the sequence file 'ncbi_ecoli_k12_protein.faa' was found in the annotation file 'ncbi_ecoli_sakai.gff'.*NP_414542.1.*protein IDs, e.g.: NP_308028.1")
        # some sequences without annotation: left out with a message
        x <- Biostrings::readAAStringSet(fx("ncbi_ecoli_sakai_protein.faa"))
        names(x)[1:2] <- c("unknown_1", "unknown_2")
        part <- tempfile(fileext = ".faa")
        Biostrings::writeXStringSet(x, part)
        expect_message(retrieve_longest_isoforms(part, fx("ncbi_ecoli_sakai.gff"), out, annotation_format = "gff"),
                       "2 of 39 sequences are not in the annotation and are left out, e.g.: unknown_1, unknown_2")
        expect_length(fasta_headers(out), 37)
        names(x)[1:2] <- "same"
        Biostrings::writeXStringSet(x, part)
        expect_error(suppressMessages(retrieve_longest_isoforms(part, fx("ncbi_ecoli_sakai.gff"), out, annotation_format = "gff")), "more than once")
})

test_that("retrieve_longest_isoforms() gives what the 0.4.3 implementation gave where that worked ...", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        for (fmt in c("gtf", "gff")) {
                new <- tempfile(fileext = ".faa"); old <- tempfile(fileext = ".faa")
                suppressMessages(retrieve_longest_isoforms(x$query$proteome, x$query[[fmt]], new, annotation_format = fmt))
                suppressMessages(retrieve_longest_isoforms_do(x$query$proteome, x$query[[fmt]], old, annotation_format = fmt))
                expect_identical(unname(tools::md5sum(new)), unname(tools::md5sum(old)), info = fmt)
        }
})

# ---------------------------------------------------------------- select_orthologs(), generate_ortholog_tables()

test_that("select_orthologs() finds the gene loci of versioned Ensembl IDs in GTF and GFF3 (N08, N23) ...", {
        cds <- fx("ensembl_human_chr21_cds.fa")
        h <- fasta_headers(cds)
        ids <- sub(" .*", "", h)
        gene <- sub("\\..*", "", header_field(h, "gene"))
        dnds <- mk_dnds(ids, rev(ids))
        for (fmt in list(c("gtf", "gtf"), c("gff", "gff"), c("gtf", "gff"))) {
                ann <- c(gtf = fx("ensembl_human_chr21.gtf"), gff = fx("ensembl_human_chr21.gff3"))[fmt]
                sv <- suppressMessages(select_orthologs(dnds, ann[1], ann[2], collapse_by = "splice_variant", format = fmt))
                expect_equal(nrow(sv), length(ids), info = paste(fmt, collapse = "/"))
                expect_identical(sv$query_gene_locus_id, gene[match(sv$query_id, ids)])
                expect_identical(sv$subject_gene_locus_id, gene[match(sv$subject_id, ids)])
                gl <- suppressMessages(select_orthologs(dnds, ann[1], ann[2], collapse_by = "gene_locus", format = fmt))
                expect_false(anyDuplicated(gl$query_gene_locus_id) > 0)
                expect_false(anyDuplicated(gl$subject_gene_locus_id) > 0)
                expect_gt(nrow(gl), 0)
        }
        # IDs without version
        sv <- suppressMessages(select_orthologs(mk_dnds(sub("\\.[0-9]+$", "", ids), sub("\\.[0-9]+$", "", rev(ids))),
                                                fx("ensembl_human_chr21.gtf"), fx("ensembl_human_chr21.gff3"),
                                                collapse_by = "splice_variant", format = c("gtf", "gff")))
        expect_equal(nrow(sv), length(ids))
        # raw GFF3 IDs ('transcript:ENST...') are accepted too
        sv <- suppressMessages(select_orthologs(mk_dnds(paste0("transcript:", sub("\\.[0-9]+$", "", ids)), ids),
                                                fx("ensembl_human_chr21.gff3"), fx("ensembl_human_chr21.gtf"),
                                                collapse_by = "splice_variant", format = c("gff", "gtf")))
        expect_equal(nrow(sv), length(ids))
})

test_that("select_orthologs() finds the loci of NCBI protein IDs and 'cds_from_genomic' IDs (N08) ...", {
        sak <- sub(" .*", "", fasta_headers(fx("ncbi_ecoli_sakai_protein.faa")))
        k12 <- sub(" .*", "", fasta_headers(fx("ncbi_ecoli_k12_protein.faa")))[1:30]
        q_cds_ids <- paste0("lcl|NC_002695.2_cds_", sak[1:30], "_", 1:30)
        res <- suppressMessages(select_orthologs(mk_dnds(q_cds_ids, k12), fx("ncbi_ecoli_sakai.gff"), fx("ncbi_ecoli_k12.gtf"),
                                                 collapse_by = "splice_variant", format = c("gff", "gtf")))
        expect_equal(nrow(res), 30)
        expect_true(all(grepl("^ECs_[0-9]{4}$", res$query_gene_locus_id)))
        expect_true(all(grepl("^b[0-9]{4}$", res$subject_gene_locus_id)))
        expect_identical(res$query_species[1], "ncbi_ecoli_sakai")
        # IDs that are in no annotation: an error with examples, not 'Something went wrong'
        expect_error(suppressMessages(select_orthologs(mk_dnds(paste0("x", 1:3), k12[1:3]), fx("ncbi_ecoli_sakai.gff"), fx("ncbi_ecoli_k12.gtf"),
                                                       format = c("gff", "gtf"))),
                     "None of the 3 IDs of the column 'query_id' of the dN/dS table was found in the annotation file 'ncbi_ecoli_sakai.gff'")
        # some IDs missing: reported
        expect_message(select_orthologs(mk_dnds(c(sak[1:5], "x1"), c(k12[1:5], "y1")), fx("ncbi_ecoli_sakai.gff"), fx("ncbi_ecoli_k12.gtf"),
                                        collapse_by = "splice_variant", format = c("gff", "gtf")),
                       "1 of 6 query_ids are not in the annotation file")
})

test_that("the loci of select_orthologs() are those of the 0.4.3 implementation where that worked ...", {
        skip_if_not_installed("rtracklayer")
        x <- make_test_ortholog_input()
        key <- function(d) { d <- as.data.frame(d); d <- d[order(d$id), ]; rownames(d) <- NULL; d }
        new <- annotation_loci_of_ids(x$dnds$query_id, x$query$gtf, "gtf")
        old <- annotation_loci_of_ids_do(x$dnds$query_id, x$query$gtf, "gtf")
        expect_identical(key(new), key(old))
        # GFF3: the locus was the raw ID of the gene row ('gene-AT1G01010'); it is the ID without the prefix now
        new <- annotation_loci_of_ids(x$dnds$query_id, x$query$gff, "gff")
        old <- annotation_loci_of_ids_do(x$dnds$query_id, x$query$gff, "gff")
        old$gene_locus_id <- strip_id_prefix(old$gene_locus_id)
        expect_identical(key(new), key(old))
})

test_that("generate_ortholog_tables() works on Ensembl files with versioned IDs ...", {
        cds <- fx("ensembl_human_chr21_cds.fa")
        ids <- sub(" .*", "", fasta_headers(cds))
        d <- new_test_dir("ortho_tables_")
        dnds_file <- file.path(d, "dnds.csv")
        utils::write.table(mk_dnds(ids, rev(ids)), dnds_file, sep = ";", row.names = FALSE, quote = FALSE)
        res <- suppressMessages(generate_ortholog_tables(dnds_file, fx("ensembl_human_chr21.gtf"), fx("ensembl_human_chr21.gff3"),
                                                         output_folder = d, format = c("gtf", "gff")))
        expect_s3_class(res, "tbl_df")
        expect_false(anyDuplicated(res$query_gene_locus_id) > 0)
        expect_true(all(grepl("^ENSG[0-9]+$", res$query_gene_locus_id)))
        expect_length(list.files(d, pattern = "orthologs_by_.*csv$", recursive = TRUE), 2)
})

# ---------------------------------------------------------------- check_annotation() (N40)

test_that("check_annotation() finds over-long lines without holding the file in memory ...", {
        ok <- fx("ncbi_ecoli_sakai.gff")
        expect_identical(suppressMessages(check_annotation(ok, remove_annotation_outliers = TRUE)), ok)
        x <- readLines(ok)
        long <- paste0(x[10], ";Note=", strrep("A", 70000))
        f <- file.path(new_test_dir("ann_"), "long.gff")
        writeLines(c(x[1:9], long, x[10:length(x)]), f)
        expect_warning(res <- suppressMessages(check_annotation(f)), "1 line that have more than 65000 characters")
        expect_identical(res, f)
        expect_warning(cleaned <- suppressMessages(check_annotation(f, remove_annotation_outliers = TRUE)),
                       "1 line of more than 65000 characters is not in the copy")
        expect_false(identical(cleaned, f))
        expect_identical(readLines(cleaned), x)
        expect_identical(readLines(f)[10], long)        # the file of the user is not changed
        # small block size: the same answer
        expect_identical(annotation_long_lines(f, n = 7L), 10L)
})
