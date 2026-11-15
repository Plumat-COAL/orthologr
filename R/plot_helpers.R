# Internal checks shared by the plot functions.

# The 'scope' column of an ortholog table: the share of the query that the alignment covers,
# 1 - |q_len - alig_length| / q_len. generate_ortholog_tables_all() and retrieve_core_orthologs() add it; tables
# from import_ortholog_tables_all(), select_orthologs() or dNdS() do not have it, and the plot functions then failed
# with "..2 must be of size ... or 1, not size 0". It is computed here when it is missing.
ensure_scope <- function(ortho_tables, fun) {
        if (is.element("scope", names(ortho_tables)))
                return(ortho_tables)
        if (!all(is.element(c("q_len", "alig_length"), names(ortho_tables))))
                stop(fun, ": 'ortho_tables' has no column 'scope' and no columns 'q_len' and 'alig_length' to compute it from.",
                     call. = FALSE)
        ortho_tables$scope <- 1 - (abs(ortho_tables$q_len - ortho_tables$alig_length) / ortho_tables$q_len)
        ortho_tables
}

# Every species of the table must be named in 'species_order' (exactly): a species that is not there became an
# NA factor level and was plotted at 'NA' without a word.
check_species_order <- function(species, species_order, fun) {
        missing_species <- setdiff(unique(as.character(species)), as.character(species_order))
        if (length(missing_species) > 0)
                stop(fun, ": 'species_order' does not name ", length(missing_species), " species of the table: '",
                     paste(missing_species, collapse = "', '"), "'. 'species_order' must hold the values of the species column ",
                     "exactly (it holds: '", paste(species_order, collapse = "', '"), "').", call. = FALSE)
        invisible(TRUE)
}
