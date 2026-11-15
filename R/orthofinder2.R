#' @title Interface function to OrthoFinder (soft-deprecated)
#' @description Run OrthoFinder from R. \strong{This function is soft-deprecated}: it is a thin, unmaintained
#' wrapper that starts \code{orthofinder -f <proteome_folder>} and gives a warning when it is called. For
#' anything beyond that, run OrthoFinder on the command line and import its tables, e.g. with
#' \code{\link{orthofinder2_retrieve_core_orthologs}}.
#' @param proteome_folder file path to a folder storing the proteome sequences of the species for which orthology inference shall be performed.
#' @param use_existing_output a logical value indicating whether or not an existing \code{OrthoFinder}
#' output folder (\code{<proteome_folder>/OrthoFinder/Results_<name of proteome_folder>}) shall be used for import
#' instead of running OrthoFinder.
#' @param import_type type of \code{OrthoFinder} output that shall be imported after running \code{OrthoFinder}.
#' Options are:
#' \itemize{
#' \item \code{import_type = NULL} (Default): nothing is imported; the path of the results folder is returned invisibly.
#' \item \code{import_type = "orthogroups_core"}: the single-copy orthogroups, see \code{\link{orthofinder2_retrieve_core_orthologs}}.
#' They are only written by a full OrthoFinder run (\code{orthogroup_only = FALSE}).
#' }
#' \code{import_type = "orthogroups_pairwise"} was never implemented and is an error.
#' @param comp_cores number of cores that shall be used for parallel processing. Default is \code{cores = 1}.
#' @param orthogroup_only logical, default \code{TRUE}. When \code{TRUE}, OrthoFinder is run with the
#' \code{-og} flag, which stops the analysis after orthogroup inference and skips the more
#' computationally expensive gene-tree and orthologue inference steps. Set to \code{FALSE} to
#' run the full OrthoFinder pipeline.
#' @param of_path a character string specifying the folder that holds the \code{orthofinder} executable, e.g.
#' \code{of_path = "/opt/miniconda3/bin/"}. The default \code{of_path = NULL} looks for \code{orthofinder} on the
#' \code{PATH} (until version 0.4.3 the path \code{/opt/miniconda3/bin/orthofinder} was hard-coded).
#' @return invisibly, the path of the OrthoFinder results folder; with \code{import_type = "orthogroups_core"} the
#' table of single-copy orthogroups.
#' @author Hajk-Georg Drost
#' @note OrthoFinder and the programs it uses (DIAMOND, MCL, ...) must be installed. Tested with OrthoFinder 3.
#' @examples \dontrun{
#' # specify species names
#' orgs <- c("Arabidopsis lyrata",
#'           "Capsella rubella", "Solanum lycopersicum")
#' # download proteome files for all species
#' biomartr::getProteomeSet(db = "refseq", organisms = orgs, path = "of_proteomes")
#' # download annotation files for all species
#' biomartr::getGFFSet(db = "refseq", organisms = orgs, path = "of_gff")
#' # select longest splice variant per gene locus
#' retrieve_longest_isoforms_all(proteome_folder = "of_proteomes",
#'                               annotation_folder = "of_gff",
#'                               annotation_format = "gff",
#'                               output_folder = "of_proteomes_longest_sv")
#' # run OrthoFinder to infer ortho groups for the specified species
#' orthofinder2(proteome_folder = "of_proteomes_longest_sv", comp_cores = 4)
#' }
#' @export
orthofinder2 <-
        function(proteome_folder,
                 use_existing_output = FALSE,
                 import_type = NULL,
                 comp_cores = 1,
                 orthogroup_only = TRUE,
                 of_path = NULL) {
                
                warning(
                        "orthofinder2() is soft-deprecated: it only starts 'orthofinder -f <proteome_folder>' and is not maintained. ",
                        "Consider running OrthoFinder on the command line and importing its tables ",
                        "(e.g. with orthofinder2_retrieve_core_orthologs()).",
                        call. = FALSE
                )
                
                if (!is.null(of_path))
                        of_path <- normalizePath(of_path, mustWork = FALSE)
                
                if (!is.null(import_type)) {
                        if (identical(import_type, "orthogroups_pairwise"))
                                stop("import_type = 'orthogroups_pairwise' is not implemented. ",
                                     "Please use import_type = NULL and import the OrthoFinder tables you need yourself.",
                                     call. = FALSE)
                        if (!identical(import_type, "orthogroups_core"))
                                stop(
                                        "The specified 'import_type' is not supported by this function. Please choose import_type = NULL or import_type = 'orthogroups_core'.",
                                        call. = FALSE
                                )
                }
                
                if (!file.exists(proteome_folder))
                        stop(
                                "Please provide a valid path to the proteome folder. Your specified path '",
                                proteome_folder,
                                "' seems not to exist.",
                                call. = FALSE
                        )
                
                proteome_folder <- normalizePath(proteome_folder)
                run_name <- basename(proteome_folder)
                results_folder <- file.path(proteome_folder, "OrthoFinder", paste0("Results_", run_name))
                
                if (!use_existing_output) {
                        is_installed_orthofinder(path = of_path)
                        
                        comp_cores <- check_comp_cores(comp_cores)
                        
                        if (file.exists(results_folder))
                                stop("The OrthoFinder results folder '", results_folder, "' exists already. ",
                                     "Please remove it, or use use_existing_output = TRUE to import from it.", call. = FALSE)
                        
                        message("Running OrthoFinder with ", comp_cores, " cores ...")
                        message(
                                "Do your input fasta files store unique sequences per gene locus, e.g. the longest splice variant? If not, then have a look at ?retrieve_longest_isoforms."
                        )
                        
                        # output is stored in <proteome_folder>/OrthoFinder/Results_<run_name>
                        run_tool(
                                "orthofinder",
                                args = list("-f", proteome_folder, "-t", comp_cores, "-a", comp_cores,
                                            "-S", "diamond", "-n", run_name, if (orthogroup_only) "-og"),
                                path = of_path, tool = "OrthoFinder", context = paste0("'", proteome_folder, "'")
                        )
                        
                        if (!file.exists(results_folder))
                                stop("OrthoFinder finished, but its results folder '", results_folder, "' was not found.",
                                     call. = FALSE)
                        
                        message("OrthoFinder finished successfully and stored all results in ", results_folder)
                } else if (!file.exists(results_folder)) {
                        stop("use_existing_output = TRUE, but the OrthoFinder results folder '", results_folder,
                             "' does not exist.", call. = FALSE)
                }
                
                if (identical(import_type, "orthogroups_core")) {
                        files <- orthofinder_core_files(results_folder)
                        return(orthofinder2_retrieve_core_orthologs(
                                orthogroups_file = files$orthogroups_file,
                                single_copy_file = files$single_copy_file
                        ))
                }
                
                invisible(results_folder)
        }

# The orthogroup table and the list of single-copy orthogroups in an OrthoFinder results folder
# (OrthoFinder 2.3-2.5: Orthogroups/Orthogroups.tsv; OrthoFinder 3 with -og: Orthogroups/MCL_Orthogroups.tsv, or
# Phylogenetic_Hierarchical_Orthogroups / Orthogroups/Orthogroups.tsv after a full run).
orthofinder_core_files <- function(results_folder) {
        og_dir <- file.path(results_folder, "Orthogroups")
        tables <- file.path(og_dir, c("Orthogroups.tsv", "MCL_Orthogroups.tsv", "Orthogroups.csv"))
        tables <- tables[file.exists(tables)]
        single <- file.path(og_dir, c("Orthogroups_SingleCopyOrthologues.txt", "SingleCopyOrthogroups.txt"))
        single <- single[file.exists(single)]
        if (length(tables) == 0 || length(single) == 0)
                stop(
                        "The single-copy orthogroups cannot be imported from '", results_folder, "': ",
                        if (length(tables) == 0) "no orthogroup table (Orthogroups/Orthogroups.tsv) " else "",
                        if (length(single) == 0) "no list of single-copy orthogroups (Orthogroups/Orthogroups_SingleCopyOrthologues.txt) " else "",
                        "was found. The list of single-copy orthogroups is only written by a full OrthoFinder run ",
                        "(orthogroup_only = FALSE). Files in '", og_dir, "': ",
                        paste(list.files(og_dir), collapse = ", "), ".",
                        call. = FALSE
                )
        list(orthogroups_file = tables[1], single_copy_file = single[1])
}
