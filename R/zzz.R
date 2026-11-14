# The note shown when the package is attached (U06): a script written for orthologr 0.4.3 runs unchanged on this
# version, but several defaults changed its numbers (and the orthologs found), without any sign at run time. Shown
# once per session (when the package is attached), not while knitr runs (vignettes, reports), and not when
# options(orthologr.quiet_upgrade = TRUE) or any other orthologr.* option is set (a user who sets options knows
# the version). packageStartupMessage(): suppressPackageStartupMessages(library(orthologr)) hides it as well.
upgrade_note <- function(pkgname = "orthologr") {
        paste0(pkgname, " ", utils::packageVersion(pkgname), ": some defaults that change dN/dS values and the ",
               "orthologs found are not those of 0.4.3 (best hits by bit score, DIAMOND's default sensitivity, ",
               "BLOSUM62 protein alignments, NA for saturated estimates, a status column). ",
               "vignette(\"upgrading_from_0.4.3\", package = \"orthologr\") lists them; orthologr_legacy_options() ",
               "gives the results of 0.4.3. options(orthologr.quiet_upgrade = TRUE) hides this note.")
}

show_upgrade_note <- function() {
        if (isTRUE(getOption("knitr.in.progress")))
                return(FALSE)
        quiet <- getOption("orthologr.quiet_upgrade")
        if (!is.null(quiet))
                return(!isTRUE(quiet))
        !any(startsWith(names(options()), "orthologr."))
}

.onAttach <- function(libname, pkgname) {
        if (show_upgrade_note())
                packageStartupMessage(paste(strwrap(upgrade_note(pkgname)), collapse = "\n"))
}
