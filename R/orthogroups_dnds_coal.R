# engine = "coal" of dNdS_orthogroups(). COAL is integrated on the branch feat/coal (dNdS_coal(), dNdS_clusters());
# in a version of orthologr without that integration the engine is not available. The two functions here are the
# whole interface: the version with COAL replaces this file.

# Check the arguments of the COAL engine; returns the method(s) to compute.
dnds_groups_coal_check <- function(dnds_est.method, aa_aln_tool) {
        stop("engine = \"coal\" is not available in this version of orthologr (it has no COAL integration). ",
             "Please use engine = \"orthologr\".", call. = FALSE)
}

# Same interface as dnds_groups_orthologr().
dnds_groups_coal <- function(pairs, cds, dnds_est.method, aa_aln_tool, aa_aln_path = NULL, comp_cores = 1, quiet = TRUE) {
        dnds_groups_coal_check(dnds_est.method, aa_aln_tool)
}
