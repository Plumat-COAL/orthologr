# Site and difference counts of the LWL family (LWL, LPB, MLWL, MLPB), recomputed from a codon alignment.
#
# kaks_saturated() uses them where the counts KaKs_Calculator prints (6 significant digits) cannot tell the sign
# of 1 - 2P - Q or 1 - 2Q of a class of sites. The method is that of the papers:
#   - Li, Wu & Luo (1985) Mol Biol Evol 2:150: every position of a sense codon is nondegenerate (0-fold: every
#     change is nonsynonymous), fourfold (every change synonymous) or twofold (otherwise; the threefold third
#     position of the isoleucine codons is counted as twofold). L0, L2, L4: the mean number of sites of each class
#     in the two sequences. Differences are transitions (A <-> G, C <-> T) or transversions, counted per class
#     (S_i, V_i); for codons that differ at two or three positions every order of the changes (pathway) that does
#     not pass through a stop codon has the same weight. P_i = S_i / L_i, Q_i = V_i / L_i.
#   - Li (1993) J Mol Evol 36:96 and Pamilo & Bianchi (1993) Mol Biol Evol 10:271 (LPB): the same counts.
#   - Tzeng, Pan & Li (2004) Mol Biol Evol 21:2290 (MLWL, MLPB): the synonymous transversions between the
#     arginine codons CGA <-> AGA and CGG <-> AGG at the first position and between the isoleucine codons
#     ATA <-> ATT, ATC at the third are counted as transitions, the nonsynonymous transition ATA <-> ATG as a
#     transversion (the rule of the method: at twofold sites synonymous changes count as transitions).
# The genetic codes are NCBI's (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), those KaKs_Calculator
# 2.0 offers with '-c'.
#
# The details that the papers leave open were established by running KaKs_Calculator 2.0 on constructed input and
# reading its output (all codes, all codon pairs; see CLEANROOM.md of branch chore/gpl2-cleanroom):
#   - a change between two codons counts half at the class of the position in each of the two codons;
#   - a codon pair with a stop codon, or with a codon that has a character other than A, C, G, T or U (any case),
#     contributes nothing, not even sites; U is read as T for the amino acid but is a letter of its own when
#     changes are counted (C <-> U is a transition, every other change of U a transversion);
#   - the arginine pairs above are special in LWL and LPB too: the half at the first codon counts as a transition,
#     the half at the second as a transversion; the modified methods apply their rule to the same codon pairs in
#     every genetic code;
#   - the counts of a codon pair are the number of halves over its pathways divided by the number of pathways, and
#     are added up codon by codon in alignment order (sequentially, in double precision). This decides the sign
#     of 1 - 2P - Q where it is 0 in exact arithmetic.

# The genetic codes of KaKs_Calculator 2.0 ('-c'), from NCBI's tables: the amino acids of the 64 codons in the
# order TTT, TTC, TTA, TTG, TCT, ..., GGG (bases T, C, A, G), "*" for a stop codon.
kaks_genetic_codes <- function() {
        c("1"  = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "2"  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
          "3"  = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "4"  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "5"  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
          "6"  = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "9"  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
          "10" = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "11" = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "12" = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "13" = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
          "14" = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
          "15" = "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "16" = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "21" = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
          "22" = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "23" = "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
}

# The 64 codons in the order of kaks_genetic_codes().
kaks_codons <- function() {
        b <- c("T", "C", "A", "G")
        paste0(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16))
}

# The translation table of a genetic code of KaKs_Calculator: a character vector of the 64 amino acids (one-letter
# codes, "*" for stop), named by codon (TTT ... GGG); NULL for a code KaKs_Calculator does not have.
kaks_transl_table <- function(genetic_code) {
        if (length(genetic_code) != 1 || is.na(suppressWarnings(as.integer(genetic_code))))
                return(NULL)
        tables <- kaks_genetic_codes()
        key <- as.character(as.integer(genetic_code))
        if (!is.element(key, names(tables)))
                return(NULL)
        stats::setNames(strsplit(tables[[key]], "")[[1]], kaks_codons())
}

# Tables of the LWL family for a genetic code: list(aa, cls, L, modified), or NULL for a code KaKs_Calculator does
# not have.
#   aa       : kaks_transl_table()
#   cls      : 64 x 3 integer matrix, the class of each position of each codon (1: 0-fold, 2: twofold, 3: fourfold;
#              NA for a stop codon), rows in the order of kaks_codons()
#   L        : 64 x 3 matrix, the number of 0-, 2- and 4-fold sites of each codon (0 for a stop codon)
#   modified : TRUE for MLWL and MLPB (Tzeng et al. 2004)
#   key      : the key of the tables in kaks_lwl_store, where they and the counts of the codon pairs met so far
#              (kaks_lwl_pair()) are kept for the session
kaks_lwl_tables <- function(genetic_code, modified) {
        aa <- kaks_transl_table(genetic_code)
        if (is.null(aa))
                return(NULL)
        modified <- isTRUE(modified)
        key <- paste(as.integer(genetic_code), modified)
        tab <- kaks_lwl_store[[key]]
        if (!is.null(tab))
                return(tab)
        codons <- kaks_codons()
        cls <- matrix(NA_integer_, 64, 3, dimnames = list(codons, c("0", "2", "4")))
        for (i in which(aa != "*"))
                for (p in 1:3)
                        cls[i, p] <- kaks_lwl_class(aa, codons[i], p)
        L <- matrix(0, 64, 3, dimnames = list(codons, c("0", "2", "4")))
        for (i in which(aa != "*"))
                L[i, ] <- tabulate(cls[i, ], nbins = 3)
        tab <- list(aa = aa, cls = cls, L = L, modified = modified, key = key)
        assign(key, tab, envir = kaks_lwl_store)
        tab
}
kaks_lwl_store <- new.env(parent = emptyenv())

# The class of position p of a sense codon under translation table aa (1: 0-fold, 2: twofold, 3: fourfold): the
# number of the three other bases at p that keep the amino acid (stop codons count as an amino acid of their own).
kaks_lwl_class <- function(aa, codon, p) {
        b <- strsplit(codon, "")[[1]]
        a <- aa[[codon]]
        syn <- 0L
        for (n in setdiff(c("T", "C", "A", "G"), b[p])) {
                b2 <- b
                b2[p] <- n
                syn <- syn + (aa[[paste(b2, collapse = "")]] == a)
        }
        if (syn == 0L) 1L else if (syn == 3L) 3L else 2L
}

# The counts of one codon pair (x, y: three characters of T, C, A, G, U, upper case): list(L, S, V), each a vector of
# the classes 0-, 2-, 4-fold. Every pathway from x to y (an order of the positions that differ) that passes no stop
# codon has the same weight; each change counts half at the class of its position in the codon before and half
# in the codon after it.
kaks_lwl_pair <- function(tab, x, y) {
        none <- c(0, 0, 0)
        aa <- tab$aa
        sense <- function(codon) aa[[chartr("U", "T", codon)]] != "*"
        if (!sense(x) || !sense(y))
                return(list(L = none, S = none, V = none))
        clsof <- function(codon, p) tab$cls[chartr("U", "T", codon), p]
        L <- (tab$L[chartr("U", "T", x), ] + tab$L[chartr("U", "T", y), ]) / 2
        bx <- strsplit(x, "")[[1]]
        by <- strsplit(y, "")[[1]]
        diff <- which(bx != by)
        if (length(diff) == 0)
                return(list(L = unname(L), S = none, V = none))
        orders <- kaks_permutations(diff)
        nS <- nV <- c(0, 0, 0)  # number of halves (times 0.5) over the pathways, exact
        npath <- 0L
        for (o in orders) {
                u <- bx
                halves <- list()
                ok <- TRUE
                for (p in o) {
                        v <- u
                        v[p] <- by[p]
                        cu <- paste(u, collapse = "")
                        cv <- paste(v, collapse = "")
                        if (!sense(cv)) {
                                ok <- FALSE
                                break
                        }
                        kind <- kaks_lwl_kinds(cu, cv, u[p], v[p], tab$modified)
                        halves[[length(halves) + 1]] <- list(kind[1], clsof(cu, p))
                        halves[[length(halves) + 1]] <- list(kind[2], clsof(cv, p))
                        u <- v
                }
                if (!ok)
                        next
                npath <- npath + 1L
                for (h in halves) {
                        if (h[[1]])
                                nS[h[[2]]] <- nS[h[[2]]] + 0.5
                        else
                                nV[h[[2]]] <- nV[h[[2]]] + 0.5
                }
        }
        if (npath == 0L)
                return(list(L = unname(L), S = none, V = none))
        list(L = unname(L), S = nS / npath, V = nV / npath)
}

# Whether the two halves of the change u -> v (codons as strings; a, b the bases that change) count as transitions
# (TRUE) or transversions (FALSE).
kaks_lwl_kinds <- function(u, v, a, b, modified) {
        pair <- paste(u, v)
        arg <- c("CGA AGA", "AGA CGA", "CGG AGG", "AGG CGG")
        if (modified) {
                if (is.element(pair, c(arg, "ATA ATT", "ATT ATA", "ATA ATC", "ATC ATA")))
                        return(c(TRUE, TRUE))
                if (is.element(pair, c("ATA ATG", "ATG ATA")))
                        return(c(FALSE, FALSE))
        } else if (is.element(pair, arg)) {
                return(c(TRUE, FALSE))
        }
        ts <- is.element(paste0(a, b), c("AG", "GA", "CT", "TC", "CU", "UC"))
        c(ts, ts)
}

# All orders of the elements of x (length 1 to 3), in the order of the positions changed first.
kaks_permutations <- function(x) {
        if (length(x) <= 1)
                return(list(x))
        out <- list()
        for (i in seq_along(x))
                for (rest in kaks_permutations(x[-i]))
                        out[[length(out) + 1]] <- c(x[i], rest)
        out
}

# The counts of the LWL family for two aligned coding sequences: list(L, Si, Vi, P, Q), each a vector of the
# classes 0-, 2- and 4-fold (sites, transitions, transversions, P = Si / L, Q = Vi / L), as KaKs_Calculator
# computes them for 'method' (LWL, LPB, MLWL or MLPB) and genetic code 'genetic_code' ('-c'). NULL when they cannot
# be computed here: another method, a genetic code KaKs_Calculator does not have, or sequences of different
# length or a length that is not a multiple of 3.
kaks_lwl_counts <- function(seq1, seq2, method, genetic_code = 1L) {
        if (length(method) != 1 || !is.element(method, kaks_lwl_methods()))
                return(NULL)
        if (!is.character(seq1) || !is.character(seq2) || length(seq1) != 1 || length(seq2) != 1 ||
            is.na(seq1) || is.na(seq2))
                return(NULL)
        n <- nchar(seq1)
        if (n != nchar(seq2) || n %% 3 != 0)
                return(NULL)
        tab <- kaks_lwl_tables(genetic_code, modified = is.element(method, c("MLWL", "MLPB")))
        if (is.null(tab))
                return(NULL)
        none <- c(0, 0, 0)
        if (n == 0)
                return(list(L = none, Si = none, Vi = none, P = none / none, Q = none / none))
        at <- seq(1, n, by = 3)
        c1 <- substring(toupper(seq1), at, at + 2)
        c2 <- substring(toupper(seq2), at, at + 2)
        # codons with a character other than A, C, G, T, U drop the pair; the rest in alignment order
        valid <- grepl("^[ACGTU]{3}$", c1) & grepl("^[ACGTU]{3}$", c2)
        c1 <- c1[valid]
        c2 <- c2[valid]
        codons <- kaks_codons()
        t1 <- match(chartr("U", "T", c1), codons)
        t2 <- match(chartr("U", "T", c2), codons)
        sense <- tab$aa[t1] != "*" & tab$aa[t2] != "*"
        c1 <- c1[sense]
        c2 <- c2[sense]
        t1 <- t1[sense]
        t2 <- t2[sense]
        # sites: halves of whole numbers, exact in any order
        L <- (colSums(tab$L[t1, , drop = FALSE]) + colSums(tab$L[t2, , drop = FALSE])) / 2
        # differences: one row per codon pair that differs, in alignment order
        d <- which(c1 != c2)
        pairs <- paste(c1[d], c2[d])
        S <- matrix(0, length(d), 3)
        V <- matrix(0, length(d), 3)
        for (u in unique(pairs)) {
                pk <- paste(tab$key, u)
                pc <- kaks_lwl_store[[pk]]
                if (is.null(pc)) {
                        cc <- strsplit(u, " ", fixed = TRUE)[[1]]
                        pc <- kaks_lwl_pair(tab, cc[1], cc[2])
                        assign(pk, pc, envir = kaks_lwl_store)
                }
                k <- pairs == u
                S[k, ] <- rep(pc$S, each = sum(k))
                V[k, ] <- rep(pc$V, each = sum(k))
        }
        # added up one codon pair after the other in double precision (sum() would use extended precision)
        total <- function(m) vapply(1:3, function(i) Reduce(`+`, m[, i], 0), numeric(1))
        Si <- total(S)
        Vi <- total(V)
        L <- unname(L)
        list(L = L, Si = Si, Vi = Vi, P = Si / L, Q = Vi / L)
}
