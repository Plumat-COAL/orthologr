## `orthologr` development version

### Changes to pairwise alignments that change dN/dS values (methods change)

`dNdS()`, `divergence_stratigraphy()` and `pairwise_aln()` align each pair of orthologous proteins with
`pwalign::pairwiseAlignment(type = "global")`. Three things about this step have changed. Each of them changes
some alignments and therefore, slightly, the dN, dS and dN/dS values of the affected gene pairs; results computed
with this version are not numerically identical to results from `orthologr` <= 0.4.3. Please report the
`orthologr` version (and, if you changed it, `aa_sub_matrix`) in the methods of analyses based on these values.

- **Amino acid alignments are scored with BLOSUM62.** Up to 0.4.3 no substitution matrix was passed to
  `pairwiseAlignment()`, so every mismatch cost the same, whichever amino acids were exchanged, and the aligner took
  most of the run time of `dNdS()`. The new argument `aa_sub_matrix` (in `dNdS()`, `compute_dnds()` and
  `pairwise_aln()`, default `"BLOSUM62"`) names the substitution matrix; `aa_sub_matrix = NULL` gives the old scoring.
- **The alignment is written full length.** `pwalign::pattern()`/`subject()` return only the aligned range, so the
  overhanging ends of the longer protein were cut off before the codon alignment. The full rows
  (`pwalign::alignedPattern()`/`alignedSubject()`) are written now: overhangs are aligned against gaps.
- **A terminal stop codon is no longer aligned.** The `*` at the end of a translated protein is removed before the
  alignment, so the stop codon is not a column of the codon alignment and is not scored in the dN/dS estimate.

`options(orthologr.legacy_pairwise_aln = TRUE)` restores all of the 0.4.3 behaviour of this step (old scoring,
truncated overhangs, stop codon aligned, score file always written) for reproducing earlier results.

### Speed

- `pairwise_aln()` no longer writes a second copy of every alignment (with its score, via
  `pwalign::writePairwiseAlignments()`) into `tempdir()`; that file is only written when `store_locally = TRUE`.
  Together with the BLOSUM62 scoring this removes most of the run time of `dNdS()` with the default
  `aa_aln_tool = "NW"`.

### Robustness of `dNdS()`

- **One gene pair can no longer end a run.** Each gene pair is processed inside `tryCatch()`. A pair for which a
  step fails is returned with `NA` for dN, dS and dNdS; `dNdS()` gives a warning with the number of such pairs, and
  `attr(result, "failed_pairs")` is a `data.frame` (`query_id`, `subject_id`, `reason`) that names the pair, the
  program and its message. Before, the first failing pair stopped the whole run (and, with `comp_cores > 1`, without
  saying which pair it was).
- **The outcome of every external program on the dN/dS path is checked.** `codon_aln()` (PAL2NAL),
  `substitutionrate()` (KaKs_Calculator) and `blast()` (blastp) now check the exit status, and that an output file
  was written, and stop with a message that names the gene pair and the program and quotes what the program
  reported. PAL2NAL exits with status 0 and an empty output when a protein does not fit its coding sequence; this
  used to surface later as an unrelated error ("Please check the correct path to Comeron ..."). (`diamond()`,
  `set_diamond()` and `set_blast()` already checked their exit status.) Paths are quoted, so folders with spaces work.
- **Each `dNdS()` call works in a folder of its own**, created inside `tempdir()` and removed when the call ends,
  also after an error. Intermediate files of an earlier or failed call can no longer be mistaken for results of the
  current one. `options(orthologr.keep_tmp = TRUE)` keeps the folder and prints its path. `pairwise_aln()`,
  `codon_aln()`, `substitutionrate()`, `diamond()` etc. called on their own still write to `tempdir()` as before.
  `codon_aln()` and `substitutionrate()` delete an output file of the same name before they start.

### Bug fixes

- Gene IDs that contain a hyphen (e.g. `AT1G01010-1`) are no longer split when
  `dnds_est.method` is one of the KaKs_Calculator methods (`"YN"`, `"NG"`, ...). The IDs were recovered by splitting
  KaKs_Calculator's `<query>-<subject>` label at `-`, so `query_id` became `AT1G01010` and `subject_id` `1`, and
  `dNdS()` then lost these pairs (or returned an empty table) when joining with the hit table. The IDs are now
  taken from the codon alignment. The AXT input of KaKs_Calculator is written by R; the Perl script
  `parseFastaIntoAXT.pl` is no longer run (one process less per gene pair).
- `dNdS()` no longer mixes up the two organisms when `query_file` and `subject_file` have the same file name in
  different folders (e.g. `a/cds.fa` and `b/cds.fa`). The proteins of both were read back from files named after
  `basename()` of the inputs, i.e. twice from the same file, so the query proteins were replaced by the subject's
  ("No orthologs could be found!", or codon alignments built from the wrong protein). The proteins are now
  translated from the coding sequences that were read.
- `substitutionrate(est.method = "Li")` (and `dNdS(dnds_est.method = "Li")`) returns numeric `dN`, `dS` and `dNdS`
  columns. They were character columns, and the checks for values that cannot be used (`< 0`, `> 9.9`) compared
  strings: for example a dS of `9.95e-05` was taken for a saturated value and set to `NA`. A dN/dS that is not
  finite (dS = 0) is now `NA`.
- `compute_dnds()` with `comp_cores > 1` and `aa_aln_type = "multiple"` looked for the protein alignment in the
  wrong folder when `store_locally = TRUE`.

### Miscellaneous

- `compute_dnds()` uses one code path for `comp_cores = 1` and `comp_cores > 1`.
- Unit tests for all of the above (`tests/testthat/test-robustness.R`).
- Unit tests for `pairwise_aln()`.

## `orthologr` version 0.4.3

### New Function

- `diamond_install()` downloads the pre-compiled DIAMOND2 binary matching the user's operating system (Linux x86_64, macOS, Windows) into the R user cache. `diamond()`, `diamond_best()`, `diamond_rec()`, `set_diamond()`, `deepclust()` and `deepclust_realign()` use it automatically if `path = NULL` and `diamond` is not on the `PATH`. `diamond_path()` reports which DIAMOND executable is used.
- `deepclust()` enables `diamond deepclust` functionalities for clustering protein sequences. Multiple input `.fasta` files or `.fasta.gz` are supported.
- `deepclust_realign()` enables `diamond deepclust` functionalities for realigning protein sequences from the output of `deepclust()`. Multiple input `.fasta` files or `.fasta.gz` are supported.
- `deepclust_annotate()` enables `diamond deepclust` functionalities for annotating protein sequences from the output of `deepclust()` based on user provided names for each proteome file (or the name of the proteome file if not provided). `deepclust_annotate()` can also be used to annotate the output of `deepclust_realign()`. Multiple input `.fasta` files or `.fasta.gz` are supported.

### Miscellaneous

- Unit tests for `deepclust()`, `deepclust_realign()` and `deepclust_annotate()` have been added to ensure the correct functionality of the new features.
- Unit tests for `dNdS()`, `divergence_stratigraphy()` and `divergence_map()` functions have been created.

## `orthologr` version 0.4.2

### New Function

- `diamond()`, `set_diamond()`, `diamond_best()` and `diamond_rec()` enable a massive speed-up in pairwise sequence alignment functionalities.

### New Features

- a logo is now set for `orthologr`.
- DIAMOND2 is now used by default in `dNdS()` and `divergence_stratigraphy()` unless blast is specified (`aligner = "blast"`).

## `orthologr` version 0.4.1

### New Features

- the `divergence_stratigraphy()` and `divergence_map()` functions now include the parameter `n_quantile`, which enables users to choose the number of quantiles to generate for the `divergence map`. This could allow users to get a higher-resolution `divergence map` if `n_quantile` is greater than 10. Alternatively, this can resolve the issue of empty divergence strata when deciling the dNdS values for closely related organisms with dNdS = 0 for over 10% of the genes.

## `orthologr` version 0.4.0

### New Functions

- new function `check_annotation()` helps to detect corrupt GFF or GTF annotation files and removes such outlier lines

### New Features

- the `generate_ortholog_tables()` and `retrieve_longest_isoforms()` now include the new `check_annotation()` function to capture corrupt GFF or GTF files and fix them 
- adding a new argument `of_path` to `orthofinder2()` to allow users to specify their own path
to their locally installed `orthofinder` executable

- adding new argument `task` to `map_generator_lnc()` and `orthologs_lnc()` to allow users
to use the full `blastn` range provided by [blast_nucleotide_to_nucleotide()](https://drostlab.github.io/metablastr/reference/blast_nucleotide_to_nucleotide.html)
- adding new argument `path` to `map_generator_lnc()` to allow users to specify their local installation path of BLAST

## `orthologr` version  0.3.0

### New Functions

- new function `plot_pairwise_orthologs()` allows users to plot pairwise orthologs tables for multiple pairwise comparisons

- new function `retrieve_core_orthologs()` allows users to retrieve a core set of orthologous gene loci from several pairwise ortholog tables

- new functions `generate_ortholog_tables()` and `generate_ortholog_tables_all()` allow users to generate ortholog tables by gene locus and splice varaint for a set of species

- new function `retrieve_longest_isoforms_all()` allows users to specify folders
and retrieve the longest splice variants for all proteomes stored in a folder

- new functions `translate_cds_to_protein()` and `translate_cds_to_protein_all()` which
translate coding sequences into amino acid sequences for single or multiple files

### New Features

- in `orthologs()` the default value of `delete_corrupt_cds` changed from `delete_corrupt_cds = TRUE` to `delete_corrupt_cds = FALSE` to be consistent with
`dNdS()` and `divergence_stratigraphy()`

- the `divergence_stratigraphy()` function received a new argument `dnds_est.method`
which now allows users to select different dNdS estimation methods when running `divergence_stratigraphy()` (suggested by Momir Futo)

- the `divergence_stratigraphy()` function now allows to change the `eval` argument which wasn't passed down to the `dNdS()` call within the function (Many thanks to Momir Futo) 

- the function `map.generator()` was renamed to `map_generator_dnds()` to be more consistent with the notation of other functions

- the function `map.generator.lnc()` was renamed to `map_generator_lnc()` to be more consistent with the notation of other functions

- the function `DivergenceMap()` was renamed to `divergence_map()` to be more consistent with the notation of other functions

- the function `DivergenceMap()` was renamed to `divergence_map()` to be more consistent with the notation of other functions

- the function `orthologs.lnc` was renamed to `orthologs_lnc` to be more consistent with the notation of other functions

- the function `OF2CoreOrthologs()` was renamed to `orthofinder2_retrieve_core_orthologs()`

### Removed Functions

- the function `advanced_blast()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `advanced_makedb()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `blast.nr()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `delta.blast()` is not supported anymore and thus is not available to users anymore (please consult the [metablastr package](https://github.com/drostlab/metablastr) in case you need this functionality)

- the function `ProteinOrtho()` is not supported anymore and thus is not available to users anymore


## `orthologr` version  0.2.0

### New Features

- new function `retrieve_longest_isoforms()` which enables retrieval of the longest isoforms from a proteome file and save results as fasta file for downstream analyses
- new function `OF2CoreOrthologs()` to retrieve core orthologs across multiple species from Orthofinder2 output
- new function `extract_features()`: Helper function to extract gene loci and splice variant IDs from GFF files
- new function `filter_best_hits()`: Helper function to select best BLAST hit based on minimum evalue
- new function `generate_ortholog_tables()`: Generate ortholog tables by gene locus and splice varaint

## `orthologr` version  0.1.0

### New Features

- `read.cds()` now trimms corrupted CDS (= CDS not divisible by 3) when `delete_corrupt_cds = FALSE` is specified
- the new default value of argument `delete_corrupt_cds` in `dNdS()` is now `FALSE`.
Thus, given the new trimming feature in `read.cds()`, corrupted CDS equences will be trimmed before being translated.


## `orthologr` version  0.0.5

### Bug fixes

- The default setting of the `BLAST` argument `max_target_seqs 1` was removed from `blast_best()` and `blast_rec()` due to
the misunderstood functionality of the `BLAST` argument (See details [here](https://doi.org/10.1093/bioinformatics/bty833) and #9 ; Many thanks to @armish)

## `orthologr` version  0.0.4

### New Features 

- Users can now control via the new `delete_corrupt_cds` argument in `dNdS()` and related downstream functions
whether or not corrupted input coding sequences shall be removed prior to dN/dS inference.
In case corrupted CDS exist, the `dNdS()` function will now generate a separate fasta file which
stores all corrupted CDS so that they can be investigated. See issue #8 for details.

#### Function updates

- `dNdS()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `read.cds()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `cds2aa()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `set_blast()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast_best()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- `blast_rec()` receives new argument `delete_corrupt_cds` to remove corrupted input coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

## `orthologr` version  0.0.3

- Fixing internal path bug that caused that wrong pal2nal paths were generated when using multiple sequence aligners -> see issue https://github.com/drostlab/orthologr/issues/5 (Many thanks to Dr. Mario López-Pérez)

## `orthologr` version  0.0.2

- Fixing a major bug that caused KaKs_Calculator to not be able to correctly 
parse the kaks computation output (Many thanks to Hongyi Li who spotted the bug and found a solution).
