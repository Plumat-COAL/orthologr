# A toy set of four "species" with known orthology, simulated from the coding sequences of the Arabidopsis thaliana
# example: every species holds a copy of an ancestral gene in which 5% of the codons were replaced at random.
#   families  1-10 : one gene in every species (single-copy)
#   family   11    : duplicated in species C after the speciation (C has a second copy, 4% further codons replaced)
#   family   12    : lost in species D
#   family   13    : in species A and B only
#   family   14    : in species A only; family 15: in species B only (no ortholog)
# Each species has another style of sequence IDs (pipes and hyphens, numbers with leading zeros, long IDs, dots).
# Returns list(files = <named CDS files>, truth = data.frame(family, species, gene, copy), dir).
make_toy_species <- function(dir = new_test_dir("og_toy_"), seed = 2026) {
        old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
        on.exit({
                if (is.null(old_seed)) {
                        if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
                } else {
                        assign(".Random.seed", old_seed, envir = globalenv())
                }
        }, add = TRUE)
        set.seed(seed)
        anc <- Biostrings::readDNAStringSet(thal_cds())
        anc <- toupper(as.character(anc))
        anc <- unname(anc[nchar(anc) %% 3 == 0 & !grepl("[^ACGT]", anc)])[1:15]
        bases <- c("A", "C", "G", "T")
        codons <- apply(expand.grid(bases, bases, bases), 1, paste, collapse = "")
        sense <- setdiff(codons, c("TAA", "TAG", "TGA"))
        mutate <- function(cds, p) {
                n <- nchar(cds) %/% 3
                cod <- substring(cds, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
                # (the first and the last codon are kept)
                change <- which(stats::runif(n) < p & seq_len(n) > 1 & seq_len(n) < n)
                cod[change] <- sample(sense, length(change), replace = TRUE)
                paste(cod, collapse = "")
        }
        id <- list(
                A = function(i) sprintf("sp|A%04d|GENE-%d", i, i),
                B = function(i) sprintf("%03d", i),
                C = function(i) paste0("C_", strrep("longname", 8), "_", i),
                D = function(i) sprintf("d.%d-v1.2", i)
        )
        present <- list(A = 1:14, B = c(1:13, 15), C = 1:12, D = 1:11)
        seqs <- list()
        truth <- list()
        for (sp in names(present)) {
                fam <- present[[sp]]
                s <- vapply(anc[fam], mutate, character(1), p = 0.05)
                names(s) <- id[[sp]](fam)
                truth[[sp]] <- data.frame(family = fam, species = sp, gene = names(s), copy = 1L, stringsAsFactors = FALSE)
                if (sp == "C") {
                        dup <- mutate(s[[11]], 0.04)
                        s <- c(s, stats::setNames(dup, paste0(id$C(11), "_dup")))
                        truth[[sp]] <- rbind(truth[[sp]], data.frame(family = 11L, species = "C", gene = paste0(id$C(11), "_dup"),
                                                                     copy = 2L, stringsAsFactors = FALSE))
                }
                seqs[[sp]] <- s
        }
        files <- vapply(names(seqs), function(sp) {
                f <- file.path(dir, paste0(sp, ".fasta"))
                writeLines(as.vector(rbind(paste0(">", names(seqs[[sp]])), seqs[[sp]])), f)
                f
        }, character(1))
        list(files = files, truth = do.call(rbind, truth), dir = dir)
}

# The groups of an orthogroups object as a sorted character vector, one string per group ("gene gene gene").
og_group_strings <- function(x) {
        g <- split(x$groups$gene, x$groups$orthogroup)
        sort(unname(vapply(g, function(v) paste(sort(v), collapse = " "), character(1))))
}
