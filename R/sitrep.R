#' @title Which external programs and packages orthologr can use
#' @description \code{orthologr_sitrep()} ("situation report") looks for every external program and optional
#' package that functions of orthologr start or load, and returns a table that says what was found, where, which
#' version, and what each one is needed for - with an installation hint for what is missing. Run it after
#' installing orthologr, or when a function stops because a program is not found.
#' @details Each program is started once with its version option (at most 60 seconds), as the functions that use it
#' do; a program that is found is remembered for the rest of the R session (as by those functions), one that is not
#' found is looked for again at the next call - after installing a program, call \code{orthologr_sitrep()} again.
#'
#' Covered: DIAMOND (on the \code{PATH}, and the copy installed by \code{\link{diamond_install}}), BLAST+
#' (\code{blastp}, \code{makeblastdb}), Perl (for PAL2NAL), KaKs_Calculator, the multiple alignment programs
#' (Clustal Omega, MAFFT, MUSCLE, T-Coffee, ClustalW), COAL (the R package \pkg{coalr} with the methods it computes,
#' and the program \code{coal-cpp}) and OrthoFinder. A program in another folder than the \code{PATH} is used
#' with the argument for its folder (e.g. \code{aligner_path}, \code{aa_aln_path}, \code{kaks_calc_path}); this
#' table looks on the \code{PATH} (and, for DIAMOND, in the cache folder of \code{diamond_install()};
#' for coal-cpp, \code{options(orthologr.coal_path)} and the environment variable \code{COAL_CPP} first).
#' @return A tibble with one row per component and the columns
#' \describe{
#' \item{\code{component}}{the program, package or feature;}
#' \item{\code{needed_for}}{the functions and arguments that use it;}
#' \item{\code{found}}{\code{TRUE} or \code{FALSE};}
#' \item{\code{path}}{the executable (or package folder) found, \code{NA} when none;}
#' \item{\code{version}}{the version line the program prints, or the version of the package; \code{NA} when not
#' found;}
#' \item{\code{hint}}{how to install what is missing (\code{NA} when it was found).}
#' }
#' @examples
#' \donttest{
#' orthologr_sitrep()
#' }
#' @seealso \code{\link{diamond_install}}, \code{\link{orthologr_options}}
#' @export
orthologr_sitrep <- function() {
        rows <- list()
        add <- function(component, needed_for, found, path = NA_character_, version = NA_character_, hint = NA_character_) {
                rows[[length(rows) + 1L]] <<- list(component = component, needed_for = needed_for, found = isTRUE(found),
                                                   path = if (isTRUE(found)) as.character(path) else NA_character_,
                                                   version = if (isTRUE(found)) as.character(version) else NA_character_,
                                                   hint = if (isTRUE(found)) NA_character_ else hint)
        }
        # (a program that is there but does not run: the hint starts with what went wrong)
        probe <- function(component, needed_for, program, version_args, pattern, hint, path = NULL) {
                p <- sitrep_probe(program, version_args, pattern, path)
                add(component, needed_for, p$found, p$path, p$version, if (p$broken) paste(p$error, hint) else hint)
                p
        }

        diamond_for <- "orthologs(), dNdS() and divergence_stratigraphy() with aligner = \"diamond\" (the default), diamond_*(), deepclust(), orthogroups()"
        probe("DIAMOND", diamond_for, "diamond", "version", "diamond",
              "orthologr::diamond_install(), or 'conda install -c bioconda diamond' (DIAMOND >= 2.1.0).")
        cache <- diamond_cache_dir()
        if (file.exists(file.path(cache, diamond_exe_name()))) {
                probe("DIAMOND (copy of diamond_install())", paste0(diamond_for, " when DIAMOND is not on the PATH"),
                      "diamond", "version", "diamond",
                      "The copy in the cache folder does not run: orthologr::diamond_install(force = TRUE).", path = cache)
        } else {
                add("DIAMOND (copy of diamond_install())", paste0(diamond_for, " when DIAMOND is not on the PATH"), FALSE,
                    hint = paste0("Not installed (", cache, "); orthologr::diamond_install() downloads it."))
        }
        blast_for <- "aligner = \"blast\" (orthologs(), dNdS(), divergence_stratigraphy()), blast*()"
        blast_hint <- "'conda install -c bioconda blast' (BLAST+)."
        bp <- sitrep_probe("blastp", "-version", "blast")
        mk <- sitrep_probe("makeblastdb", "-version", "makeblastdb")
        add("BLAST+ (blastp, makeblastdb)", blast_for, bp$found && mk$found, bp$path, bp$version,
            if (bp$found && !mk$found) paste("blastp was found but makeblastdb was not:", blast_hint) else blast_hint)
        probe("Perl", "codon_aln_tool = \"pal2nal\" (PAL2NAL; the default codon alignments need COAL instead)", "perl", "-v",
              "This is perl", "Install Perl (on Windows Strawberry Perl, https://strawberryperl.com), or use COAL (coalr).")
        kc <- probe("KaKs_Calculator", "dNdS() methods other than Comeron and Li with kaks_tool = \"kaks_calculator\" (the default)",
                    "KaKs_Calculator", "-h", "KaKs_Calculator",
                    "'conda install -c bioconda kakscalculator2', or kaks_tool = \"coal\" (coalr) for NG, LWL, MLWL, LPB, MLPB, YN, MYN.")
        # (KaKs_Calculator prints its version on the line after its name)
        if (kc$found)
                rows[[length(rows)]]$version <- paste("KaKs_Calculator", kaks_version_line(NULL))
        msa_for <- function(tool) paste0("aa_aln_type = \"multiple\", aa_aln_tool = \"", tool, "\" (dNdS()); multi_aln()")
        probe("Clustal Omega", msa_for("clustalo"), "clustalo", "--help", "Clustal", "'conda install -c bioconda clustalo'.")
        probe("MAFFT", msa_for("mafft"), "mafft", "--version", "^v[0-9]|MAFFT", "'conda install -c bioconda mafft'.")
        probe("MUSCLE", msa_for("muscle"), "muscle", "-version", "muscle", "'conda install -c bioconda muscle'.")
        probe("T-Coffee", msa_for("t_coffee"), "t_coffee", "-version", "T-COFFEE", "'conda install -c bioconda t-coffee'.")
        probe("ClustalW", msa_for("clustalw"), "clustalw2", "-help", "CLUSTAL",
              "'conda install -c bioconda clustalw' (installs clustalw2 as well).")

        coalr_hint <- "remotes::install_github(\"Plumat-COAL/COAL\", subdir = \"packaging/r/coalr\")"
        has_coalr <- requireNamespace("coalr", quietly = TRUE)
        info <- if (has_coalr) tryCatch(coalr::coal_info(), error = function(e) NULL) else NULL
        add("coalr (COAL in R)", "engine = \"coal\", codon_aln_tool = \"coal\" (the default when installed), kaks_tool = \"coal\", dNdS_coal(), dNdS_clusters()",
            has_coalr, if (has_coalr) system.file(package = "coalr") else NA, package_version_text("coalr"), coalr_hint)
        features <- as.character(info$features)
        methods <- as.character(info$methods)
        feature_row <- function(component, needed_for, ok, what) {
                add(component, needed_for, has_coalr && ok, if (has_coalr) system.file(package = "coalr") else NA,
                    package_version_text("coalr"),
                    if (has_coalr) paste0("The installed coalr ", package_version_text("coalr"), " does not report ", what,
                                          "; install a newer one: ", coalr_hint) else coalr_hint)
        }
        feature_row("coalr: GY method", "COAL's maximum-likelihood GY method (coalr feature \"gy94_ml\")",
                    is.element("gy94_ml", features), "the feature \"gy94_ml\"")
        feature_row("coalr: Comeron's method", "engine = \"coal\" / dNdS_coal() with dnds_est.method = \"Comeron\"",
                    is.element("Comeron", methods) && is.element("comeron_min_sites_positive", features),
                    "the method \"Comeron\" with the feature \"comeron_min_sites_positive\"")
        feature_row("coalr: Li's method", "engine = \"coal\" / dNdS_coal() with dnds_est.method = \"Li\"",
                    is.element("Li", methods), "the method \"Li\"")
        exe <- coal_binary()
        if (is.null(exe)) {
                add("coal-cpp", "engine = \"coal\" without coalr (options(orthologr.coal_backend = \"binary\"))", FALSE,
                    hint = "Not needed with coalr. A release of COAL (https://github.com/Plumat-COAL/COAL/releases): put coal-cpp on the PATH, or options(orthologr.coal_path = \"<file>\").")
        } else {
                probe("coal-cpp", "engine = \"coal\" without coalr (options(orthologr.coal_backend = \"binary\"))",
                      basename(exe), "--version", "coal", "coal-cpp was found but does not run.", path = dirname(exe))
        }
        probe("OrthoFinder", "orthofinder2(), retrieve_core_orthologs()", "orthofinder", "-h", "OrthoFinder version",
              "'conda install -c bioconda orthofinder'.")

        tibble::tibble(component  = vapply(rows, `[[`, "", "component"),
                       needed_for = vapply(rows, `[[`, "", "needed_for"),
                       found      = vapply(rows, `[[`, TRUE, "found"),
                       path       = vapply(rows, function(r) as.character(r$path)[1], ""),
                       version    = vapply(rows, function(r) as.character(r$version)[1], ""),
                       hint       = vapply(rows, `[[`, "", "hint"))
}

# Is `program` there (tool_available(), at most 60 seconds)? list(found, path, version, broken, error); a program
# that is not found or does not run is found = FALSE (and is looked for again next time); broken = TRUE when it
# is there but does not run or does not print what it should (error: what went wrong).
sitrep_probe <- function(program, version_args, pattern, path = NULL) {
        key <- tool_timeout_key(program, program)
        version <- tryCatch(timeout_cap(key, 60, tool_available(program, version_args = version_args, pattern = pattern,
                                                                  path = path)),
                            orthologr_error = function(e) e)
        if (inherits(version, "error")) {
                exists <- nzchar(if (is.null(path)) Sys.which(program) else
                        Sys.which(tool_command(program, path)))
                broken <- inherits(version, "orthologr_error_tool_failed") ||
                        inherits(version, "orthologr_error_tool_timeout") || exists
                return(list(found = FALSE, path = NA_character_, version = NA_character_, broken = broken,
                            error = if (inherits(version, "orthologr_error_tool_missing"))
                                    paste0("'", program, "' is there but does not print its version (",
                                           paste(version_args, collapse = " "), ").")
                                    else conditionMessage(version)))
        }
        where <- if (is.null(path)) unname(Sys.which(program)) else tool_command(program, path)
        list(found = TRUE, path = if (nzchar(where)) normalizePath(where, mustWork = FALSE) else program,
             version = trimws(version), broken = FALSE, error = NA_character_)
}
