#' @title Compute Multiple Sequence Alignments
#' @description This function takes a FASTA file containing DNA or amino acid sequences
#' that shall be aligned and computes a multiple alignment using a defined multiple alignment tool.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program that should be used: "clustalw", "t_coffee", "muscle", "clustalo", and "mafft". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @param multi_aln_name a character string specifying the name of the stored alignment file. 
#' Default is \code{multi_aln_name} = \code{NULL} denoting a default name: 'toolname.aln' .
#' @param params (an argument list, not shell text: it is split at white space, quotes group, and no character is interpreted by a shell - see the section "Parameter strings" of \code{\link{orthologr_options}}) a character string listing the input paramters that shall be passed to corresponding alignment tool specified in the \code{tool} argument. 
#' Default is \code{NULL}, implicating that a set of default parameters is used when running the correponding tool, e.g. clustalw. 
#' Example: \code{params} = "-PWMATRIX=BLOSUM -TYPE=PROTEIN" for \code{tool} = "clustalw.
#' For \code{tool} = "t_coffee" the default is T-Coffee's own default mode with CLUSTAL output. Until version 0.4.3
#' the default was \code{-mode accurate}, which needs access to a BLAST server and fails without network; use
#' \code{params = "-mode accurate"} to get it back.
#' @param quiet a logical value specifying whether the output of the corresponding multiple alignment tool shall be printed out to the console.
#' Default is \code{quiet} = \code{FALSE}.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}. 
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function provides an interface between R and common multiple alignment programs
#' such as "clustalw", "t_coffee", "muscle", "clustalo", and "mafft".
#' 
#' \itemize{
#' 
#' \item CLUSTALW : 
#' 
#' 
#' Different operating systems perform different execution calls to the clustalw program:
#' 
#' \itemize{
#' \item MacOS: 'clustalw2' 
#' \item Linux: 'clustalw2'
#' \item  Windows: 'clustalw2.exe'
#' }
#' 
#' In case you use the default path to the clustalw program, depending on your operating system,
#' the following calls to clustalw should work properly on your system:
#' 
#' \itemize{
#' \item MacOS: \code{system("clustalw2 -help")}
#' \item Linux: \code{system("clustalw2 -help")}
#' \item Windows: \code{system("clustalw2.exe -help")}
#' 
#' }
#' 
#' In case these procedures don't work properly, please use the \code{path} argument
#' to specify the 'clustalw' execution path on your system:
#' 
#' \itemize{
#' \item MacOS: \code{system("path/to/clustalw/clustalw2 -help")}
#' \item Linux: \code{system("path/to/clustalw/clustalw2 -help")}
#' \item Windows: \code{system("path/to/clustalw/clustalw2.exe -help")}
#' 
#' }
#' 
#' \item T_COFFEE :
#' 
#' In case you use the default path to the t_coffee program,
#' the following calls to clustalw should work properly on your system:
#' 
#' \code{system("t_coffee -version")}
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 't_coffee' execution path on your system:
#' 
#' \code{system("path/to/t_coffee/t_coffee -version")}
#' 
#'
#'  
#' \item MUSCLE : 
#' 
#' In case you use the default path to the muscle program,
#' the following calls to muscle should work properly on your system:
#' 
#' \code{system("muscle -version")}
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 'muscle' execution path on your system:
#' 
#' \code{system("path/to/muscle/muscle -version")}
#' 
#' Both MUSCLE 3.8 and MUSCLE 5 are supported; the version is detected from \code{muscle -version}.
#' MUSCLE 3.8 is run as \code{muscle -in <file> -out <aln> -clwstrict -quiet} (CLUSTAL output). MUSCLE 5 has
#' another command line and no CLUSTAL output: it is run as \code{muscle -align <file> -output <aln>}, the
#' alignment is written in aligned FASTA format and the sequences are put back into the order of the input
#' file. With MUSCLE 5, \code{params} must hold MUSCLE 5 options (the MUSCLE 3.8 options \code{-diags},
#' \code{-clwstrict}, ... are not accepted by it).
#' 
#' 
#' \item CLUSTALO : 
#' 
#' In case you use the default path to the clustalo program,
#' the following calls to clustalo should work properly on your system:
#' 
#' \code{system("clustalo --help")}
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 'clustalo' execution path on your system:
#' 
#' \code{system("path/to/clustalo/clustalo --help")}
#' 
#' 
#' \item MAFFT : 
#' 
#' In case you use the default path to the mafft program,
#' the following calls to mafft should work properly on your system:
#' 
#' \code{system("mafft --version")}
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 'mafft' execution path on your system:
#' 
#' \code{system("path/to/mafft/mafft --version")}
#' 
#' }
#' @note Note that when using the \code{clustalw.params}, ... parameters, make sure the corresponding alignment tool
#' returns a file in clustal format *.aln. This is only important when \code{get_aln} = \code{TRUE}.
#' @examplesIf nzchar(Sys.which("clustalw2"))
#' # CLUSTALW Example:
#'
#' # in case the default execution path of clustalw runs properly on your system
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "clustalw",
#'           get_aln = TRUE)
#'
#' # running clustalw using additional parameters
#' # details: system("clustalw2 -help")
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "clustalw",
#'           get_aln = TRUE,
#'           params  = "-PWMATRIX=BLOSUM -TYPE=PROTEIN")
#'
#' @examplesIf nzchar(Sys.which("t_coffee"))
#' # T_COFFEE Example:
#'
#' # in case the default execution path of t_coffee runs properly on your system
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "t_coffee",
#'           get_aln = TRUE)
#'
#' # running t_coffee using additional parameters
#' # details: https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "t_coffee",
#'           get_aln = TRUE,
#'           params  = "-output clustalw_aln -quiet")
#'
#' @examplesIf nzchar(Sys.which("muscle"))
#' # (modes such as "-mode expresso" or "-mode accurate" need BLAST and access to
#' # sequence / structure databases)
#'
#'
#'
#' # MUSCLE Example:
#'
#' # in case the default execution path of muscle runs properly on your system
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "muscle",
#'           get_aln = TRUE)
#'
#' # running muscle using additional parameters
#' # details: https://www.drive5.com/muscle/manual/
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "muscle",
#'           get_aln = TRUE,
#'           params  = "-threads 1")          # MUSCLE 5 (what bioconda installs)
#'
#' @examplesIf nzchar(Sys.which("clustalo"))
#' # MUSCLE 3.8 has other options, e.g. params = "-diags -clwstrict";
#' # MUSCLE 5 stops with "Unknown option diags"
#'
#'
#' # CLUSTALO Example:
#'
#' # in case the default execution path of clustalo runs properly on your system
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "clustalo",
#'           get_aln = TRUE)
#'
#' # running clustalo using additional parameters
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "clustalo",
#'           get_aln = TRUE,
#'           params  = "--outfmt clu")
#'
#' @examplesIf nzchar(Sys.which("mafft"))
#' # MAFFT Example:
#'
#' # in case the default execution path of mafft runs properly on your system
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "mafft",
#'           get_aln = TRUE)
#'
#' # running mafft using additional parameters
#' # details: https://mafft.cbrc.jp/alignment/software/manual/manual.html
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "mafft",
#'           get_aln = TRUE,
#'           params  = "--maxiterate 1 --clustalout")
#'
#' \dontrun{
#' # in case the default execution path of clustalw is not set within the default path
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "clustalw",
#'           get_aln = TRUE,
#'           path    = "path/to/clustalw/")
#'
#' # in case the default execution path of t_coffee is not set within the default path
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "t_coffee",
#'           get_aln = TRUE,
#'           path    = "path/to/t_coffee/")
#'
#' # in case the default execution path of muscle is not set within the default path
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "muscle",
#'           get_aln = TRUE,
#'           path    = "path/to/muscle/")
#'
#' # in case the default execution path of clustalo is not set within the default path
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "clustalo",
#'           get_aln = TRUE,
#'           path    = "path/to/clustalo/")
#'
#' # in case the default execution path of mafft is not set within the default path
#' multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool    = "mafft",
#'           get_aln = TRUE,
#'           path    = "path/to/mafft/")
#' }
#' @references 
#' 
#' \itemize{
#' \item CLUSTALW:
#' 
#' Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA, McWilliam H, Valentin F,
#' Wallace IM, Wilm A, Lopez R, Thompson JD, Gibson TJ, Higgins DG. (2007).
#' Clustal W and Clustal X version 2.0. Bioinformatics, 23, 2947-2948.
#' 
#' \url{https://www.clustal.org/clustal2/}
#' 
#' \url{https://www.ebi.ac.uk/jdispatcher/}
#' 
#' 
#' \item T_COFFEE
#' 
#' T-Coffee: A novel method for multiple sequence alignments. 
#' Notredame, Higgins, Heringa, JMB, 302(205-217). 2000.
#' 
#' \url{https://tcoffee.readthedocs.io/en/latest/}
#' 
#' \url{https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html}
#' 
#' 
#' \item MUSCLE:
#' 
#' Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 32(5):1792-1797. 
#' 
#' Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, (5) 113. 
#' 
#' \url{https://www.drive5.com/muscle/}
#' 
#' \url{https://www.drive5.com/muscle/manual/}
#' 
#' 
#' \item CLUSTALO:
#' 
#' Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Soeding J, Thompson JD, Higgins DG (2011). 
#' Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75
#' 
#' \url{https://www.ebi.ac.uk/jdispatcher/msa/clustalo/}
#' 
#' \url{https://www.genome.jp/tools-bin/clustalw}
#' 
#' \item MAFFT : 
#' 
#' Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780)  
#' MAFFT multiple sequence alignment software version 7: improvements in performance and usability. 
#' 
#' \url{https://mafft.cbrc.jp/alignment/software/}
#' 
#' \url{https://mafft.cbrc.jp/alignment/software/manual/manual.html}
#' 
#' \url{https://mafft.cbrc.jp/alignment/software/tips0.html}
#' 
#' }
#' @return In case the argument \code{get_aln} is set \code{TRUE}, an object of class alignment of the seqinr package is returned.
#' @export
multi_aln <- function(file, 
                      tool, 
                      get_aln        = FALSE, 
                      path           = NULL,
                      multi_aln_name = NULL, 
                      params         = NULL, 
                      quiet          = FALSE, 
                      clean_folders  = FALSE){
        
        if (!is.null(path))
                path <- normalizePath(path, mustWork = FALSE)
        
        if(!is.multiple_aln_tool(tool))
                stop("Please choose a tool that is supported by this function.")
        
        if (!file.exists(file))
                stop("The file you specified in 'file' does not seem to exist: ", file, ".", call. = FALSE)
        
        out_dir <- file.path(orthologr_tmpdir(), "_alignment", "multi_aln")
        if (!file.exists(out_dir))
                dir.create(out_dir, recursive = TRUE)
        
        file.out <- file.path(
                out_dir,
                paste0(if (is.null(multi_aln_name)) basename(file) else multi_aln_name, "_", tool, ".aln")
        )
        
        # The program writes to a file with a name of its own; it becomes `file.out` only when the program has
        # finished successfully. An alignment left by an earlier call is removed first, so it can never be
        # returned in place of the result of this call.
        if (file.exists(file.out))
                unlink(file.out)
        tmp.out <- tempfile("aln_", tmpdir = out_dir, fileext = ".tmp")
        on.exit(unlink(tmp.out), add = TRUE)
        
        context <- paste0("'", basename(file), "'")
        
        if (tool == "clustalw") {
                is_installed_clustalw(path = path)
                run_tool(
                        "clustalw2",
                        args = list(paste0("-infile=", file), paste0("-outfile=", tmp.out),
                                    if (is.null(params)) "-quiet" else I(params)),
                        path = path, tool = "ClustalW", context = context, output = tmp.out
                )
        }
        
        # the format the alignment file is written in
        aln_format <- "clustal"
        
        if (tool == "t_coffee") {
                is_installed_tcoffee(path = path)
                # Default: T-Coffee's own default mode, CLUSTAL output only, the guide tree next to the output
                # (not in the working directory). Until 0.4.3 the default was '-mode accurate', which needs a
                # BLAST server (network) and fails without one; use params = "-mode accurate" to get it back.
                # The sequences go in with '-seq': with '-infile', T-Coffee 12 takes a file whose sequences end in
                # '*' (as the translated coding sequences do) for an alignment that is already made and returns
                # it unaligned (rows padded at the end, out of register).
                tree.out <- paste0(tmp.out, ".dnd")
                on.exit(unlink(tree.out), add = TRUE)
                run_tool(
                        "t_coffee",
                        args = if (is.null(params))
                                list("-seq", file, "-output", "clustalw_aln", "-newtree", tree.out, "-quiet",
                                     "-outfile", tmp.out)
                        else
                                list("-seq", file, I(params), "-outfile", tmp.out),
                        path = path, tool = "T-Coffee", context = context, output = tmp.out
                )
        }
        
        if (tool == "muscle") {
                is_installed_muscle(path = path)
                if (muscle_major_version(path = path) >= 5) {
                        # MUSCLE 5: -align / -output, aligned FASTA only, sequences in an order of its own
                        run_tool(
                                "muscle",
                                args = list("-align", file, if (!is.null(params)) I(params), "-output", tmp.out),
                                path = path, tool = "MUSCLE", context = context, output = tmp.out
                        )
                        reorder_fasta(tmp.out, fasta_ids(file))
                        aln_format <- "fasta"
                } else {
                        run_tool(
                                "muscle",
                                args = if (is.null(params))
                                        list("-in", file, "-out", tmp.out, "-clwstrict", "-quiet")
                                else
                                        list("-in", file, I(params), "-out", tmp.out),
                                path = path, tool = "MUSCLE", context = context, output = tmp.out
                        )
                }
        }
        
        if (tool == "clustalo") {
                is_installed_clustalo(path = path)
                run_tool(
                        "clustalo",
                        args = if (is.null(params))
                                list("-i", file, "-o", tmp.out, "--outfmt", "clu", "--force")
                        else
                                list("-i", file, "-o", tmp.out, "--force", I(params)),
                        path = path, tool = "Clustal Omega", context = context, output = tmp.out
                )
        }
        
        if (tool == "mafft") {
                is_installed_mafft(path = path)
                # To avoid error with pal2nal when deleting * during multiple
                # alignment add --anysymbol
                # https://academic.oup.com/mbe/content/early/2013/02/08/molbev.mst010.full
                # MAFFT cuts the sequence names of its CLUSTAL output to 15 characters unless --namelength says
                # otherwise: IDs such as 'ENSMUST00000012345' or 'lcl|NC_000913.3_cds_NP_414542.1_1' came back
                # truncated (and so did the IDs of everything computed from the alignment).
                name_length <- max(c(15L, nchar(fasta_ids(file)) + 1L), na.rm = TRUE)
                run_tool(
                        "mafft",
                        args = if (is.null(params))
                                list("--quiet", "--anysymbol", "--clustalout", "--namelength", name_length, file)
                        else
                                list(I(params), file),
                        stdout = tmp.out,
                        path = path, tool = "MAFFT", context = context, output = tmp.out
                )
        }
        
        publish_file(tmp.out, file.out)
        
        if (!quiet) {
                message("Multiple Alignment successfully written in ", file.out, ".")
        }
        
        if (get_aln) {
                aln <- tryCatch({
                        tryCatch(
                                Biostrings::readAAMultipleAlignment(file.out, format = aln_format),
                                # e.g. the CLUSTAL format written by T-Coffee
                                error = function(e) Biostrings::AAMultipleAlignment(read_clustal(file.out))
                        )
                }, error = function(e) {
                        stop(
                                "Something went wront with ",
                                tool,
                                " .",
                                "\n",
                                file.out,
                                " could not be read properly.",
                                call. = FALSE
                        )
                })
        }
        
        if (clean_folders)
                clean_all_folders(file.path(orthologr_tmpdir(), "_alignment"))
        
        if (get_aln)
                return(aln)
        
        invisible(if (clean_folders) NULL else file.out)
}
