#' @title Compute a Codon Alignment
#' @description This function takes a protein alignment and the corresponding coding sequences 
#' as fasta files and computes the corresponding codon alignment based on the models specified in the PAL2NAL program.
#' @param file_aln a character string specifying the path to the file storing the protein alignment in CLUSTAL or FASTA format.
#' @param file_nuc a character string specifying the path to the file storing the coding sequences in multiple FASTA format.
#' @param format a character string specifying the file format used to store the codon alignment, e.g. "fasta", "clustal".
#' @param tool a character string specifying the program that should be used: \code{"pal2nal"} (Default; needs Perl)
#' or \code{"coal"} (COAL, \url{https://github.com/Plumat-COAL/COAL}: the same codon alignment without Perl, computed inside R
#' when the package \pkg{coalr} is installed and otherwise by the \code{coal-cpp} program). 
#' @param codon_aln_name a character string specifying the name of the stored alignment file. 
#' @param params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string specifying additional parameters that shall be passed to PAL2NAL. For example
#' \code{params} = \code{"-codontable 2"} would specify a codon table: \emph{vertebrate mitochondrial code}
#' instead if the \emph{universal code} (default). The default value is \code{params} = \code{NULL}.
#' Default is \code{codon_aln_name} = \code{NULL} denoting a default name: 'tool_name.aln' .
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between the R language and common codon alignment tools such as "PAL2NAL".
#' Codon alignments can be used to quantify the evolutionary pressure acting on protein sequences based on dNdS estimation. 
#' 
#' The PAL2NAL program is used to convert a multiple sequence alignment or pairwise alignment of proteins and the corresponding DNA (or mRNA) sequences into 
#' a codon-based DNA alignment. The output codon-based DNA alignment can then be used for \code{\link{dNdS}} estimation.
#' 
#' The popular codon alignment tool PAL2NAL is included in this package and is being called when choosing \code{tool} = \code{"pal2nal"}.
#' 
#' \code{tool} = \code{"coal"} builds the codon alignment with COAL, a re-implementation of PAL2NAL's algorithm that does
#' not need Perl and has no limit on the protein length (PAL2NAL fails for proteins of more than about 21,800 residues).
#' COAL is optional: install its R package with
#' \code{remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")}, or put the \code{coal-cpp} program on
#' the \code{PATH} (\code{pip install coal-align}; \code{options(orthologr.coal_path = )} names another location).
#' The PAL2NAL parameters in \code{params} are not available with \code{tool = "coal"}.
#' @note
#' 
#' The PAL2NAL program automatically checks
#' \itemize{
#' \item If you use the same IDs are used in the protein alignment file and DNA (or mRNA) file ->
#' in this case the sequences don't have to be in the same order
#' 
#' \item If you don't use the same IDs are used in the protein alignment file and DNA (or mRNA) file
#' -> you have to rearrange the IDs or sequences so that sequences in the protein alignment file and DNA (or mRNA) file
#' have the same order
#' 
#' }
#' 
#' @examplesIf nzchar(Sys.which("perl"))
#' # performing a codon alignment using PAL2NAL
#' codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#'                        file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'),
#'                        format   = "clustal",
#'                        tool     = "pal2nal",
#'                        get_aln  = TRUE,
#'                        quiet    = TRUE)
#'
#' # running PAL2NAL with additional parameters (e.g. a different codon table)
#' codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#'                        file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'),
#'                        format   = "clustal",
#'                        tool     = "pal2nal",
#'                        get_aln  = TRUE,
#'                        params   = "-codontable 2",
#'                        quiet    = TRUE)
#'                        
#' @return In case \code{get_aln} = \code{TRUE} an seqinr alignment object is returned.
#' @references 
#' 
#' Mikita Suyama, David Torrents, and Peer Bork (2006)
#' PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Res. 34, W609-W612.
#' 
#' https://www.bork.embl.de/services.html
#' 
#' https://github.com/abacus-gene/paml/wiki
#' 
#' @seealso \code{\link{pairwise_aln}}, \code{\link{multi_aln}}, \code{\link{substitutionrate}},
#'  \code{\link{dNdS}}, \code{\link{divergence_stratigraphy}}
#' @export
codon_aln <- function(file_aln, 
                      file_nuc, 
                      format         = "clustal", 
                      tool           = "pal2nal", 
                      params         = NULL, 
                      codon_aln_name = NULL, 
                      get_aln        = FALSE, 
                      quiet          = FALSE){
        
        
        if (!file.exists(file_aln))
                stop("The file you specified in 'file_aln' does not seem to exist: ", file_aln,".", call. = FALSE)
        if (!file.exists(file_nuc))
                stop("The file you specified in 'file_nuc' does not seem to exist: ", file_nuc,".", call. = FALSE)
        
        # the Pal2Nal program is stored as executable within
        # the R package environment: 'exec' folder
        # this is not an elegant version, only motivated by this discussion:
        # http://stackoverflow.com/questions/13463103/inst-and-extdata-folders-in-r-packaging
        # is there a better choice?
        if (!is.element(tool, c("pal2nal", "coal")))
                stop("Please choose a tool that is supported by this function: 'pal2nal' or 'coal'.", call. = FALSE)
        
        # the Pal2Nal program is stored as executable within
        # the R package environment
        if (tool == "pal2nal")
                path <-
                        system.file(
                                file.path("pal2nal", "pal2nal.v14", "pal2nal.pl"),
                                package = "orthologr",
                                mustWork = TRUE
                        )
        
        if (!is.element(format, c("clustal", "fasta")))
                stop("Please choose a format that is supported by this function.", call. = FALSE)
        
        if (!file.exists(file.path(orthologr_tmpdir(), "_alignment", "codon_aln"))) {
                dir.create(file.path(orthologr_tmpdir(), "_alignment", "codon_aln"), recursive = TRUE)
        }
        
        if (is.null(codon_aln_name))
                file.out <-
                file.path(orthologr_tmpdir(), "_alignment", "codon_aln", paste0(tool, ".aln"))
        
        if (!is.null(codon_aln_name))
                file.out <-
                file.path(orthologr_tmpdir(),
                          "_alignment",
                          "codon_aln",
                          paste0(codon_aln_name, "_", tool, ".aln"))
        
        
        if (tool == "pal2nal") {
                
                context <- pair_context(file_nuc, codon_aln_name)
                
                # A file left by an earlier call must never be taken for the result of this one: it is removed,
                # PAL2NAL writes to a file with a name of its own, and that file becomes `file.out` only when
                # PAL2NAL has finished successfully.
                if (file.exists(file.out))
                        unlink(file.out)
                tmp.out <- tempfile("aln_", tmpdir = dirname(file.out), fileext = ".tmp")
                on.exit(unlink(tmp.out), add = TRUE)
                
                # PAL2NAL reports a protein that does not fit its CDS on stderr, exits with status 0 and writes
                # nothing: both the exit status and the output are checked.
                run_tool(
                        "perl",
                        args = list(path, file_aln, file_nuc, "-output", format,
                                    if (!is.null(params)) I(params)),
                        stdout  = tmp.out,
                        tool    = "PAL2NAL",
                        context = context,
                        output  = tmp.out
                )
                publish_file(tmp.out, file.out)
        }
        
        # COAL: the same codon alignment without Perl: in process through the R package 'coalr' when it is
        # installed, otherwise with the coal-cpp program
        if (tool == "coal") {
                if (!is.null(params))
                        stop("'params' are PAL2NAL command-line parameters; they cannot be used with tool = \"coal\".", call. = FALSE)
                coal_codon_aln_file(file_aln = file_aln, file_nuc = file_nuc, file.out = file.out, format = format,
                                    context = pair_context(file_nuc, codon_aln_name))
        }
        
        if (!quiet) {
                message(paste0("Codon Alignment successfully written in ", file.out, "."))
        }
        
        
        if (get_aln) {
                dna_aln <- tryCatch({
                        aln <- suppressWarnings(seqinr::read.alignment(file = file.out, format = format))
                        if (format == "clustal" && (is.null(aln$nb) || aln$nb == 0)) {
                                # seqinr does not read the CLUSTAL format PAL2NAL writes (it returned an empty
                                # alignment)
                                seqs <- read_clustal(file.out)
                                aln <- seqinr::as.alignment(nb = length(seqs), nam = names(seqs),
                                                            seq = as.list(tolower(unname(seqs))), com = NA)
                        }
                        aln
                }, error = function(e) {
                        stop(
                                "Something went wront with ", tool, " .",
                                "\n",
                                file.out,
                                " could not be read properly.", call. = FALSE
                        )
                })
                return(dna_aln)
        }
        
        invisible(file.out)
}
