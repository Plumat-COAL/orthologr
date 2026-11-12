## `orthologr` development version

### Options and inputs: documented options work, unsupported ones are rejected with a clear message

`dNdS()`
- Arguments are checked before any program is started, and the error names the argument, the value and the valid
  values. `ortho_detection` accepts `"RBH"` and `"BH"` only: `"PO"` and `"DELTA"` passed the check and then failed
  with `object 'q_cds' not found`; `"Orthofinder2"` was documented but rejected. A mistyped `aligner` gave the same
  `q_cds` error. `seq_type` other than `"cds"` and `format` other than `"fasta"` are rejected (`seq_type = "protein"`
  could hang in `diamond makedb`; `format = "gbk"` was documented but never read).
- `sensitivity_mode` is checked against the DIAMOND modes; the documented default is the real one (`"fast"`).
  `sensitivity_mode = "default"` now runs DIAMOND in its own default mode (it was passed on as `--default`, which
  DIAMOND rejects). Also in `diamond()`, `diamond_best()`, `diamond_rec()`, `divergence_stratigraphy()`.
- **`aa_aln_params` is now passed on** to the multiple alignment program (`multi_aln(params = )`). It was accepted
  and ignored, so results change for calls that set it together with `aa_aln_type = "multiple"`. With
  `aa_aln_type = "pairwise"` it is not used and a warning says so.
- `comp_cores` larger than the number of available cores: the available number is used and a warning is given (was
  an error). The same in `blast()`, `diamond()`, `deepclust()`, `deepclust_realign()`, `compute_dnds()`.
- `cdd.path` and `ortho_path` are deprecated (they were never used); setting them gives a warning.
- `clean_folders = TRUE` leaves nothing in `tempdir()` (all intermediate files are in the folder of the call, which
  is always removed).

Coding sequences that are left out
- Coding sequences with characters other than A, C, G, T (N, IUPAC ambiguity codes) cannot be translated and were
  left out without a word. They are still left out, but a message now names them. The check itself used the pattern
  `[^(A|C|G|T|a|c|g|t)]`, which also let `(`, `)` and `|` pass; it is `[^ACGTacgt]` now.
- `read.cds()` wrote `<file>_corrupted_cds_seqs.fasta` into the **working directory**. The file now goes to
  `tempdir()`, or to the folder set with `options(orthologr.corrupted_cds_dir = "<folder>")`.

Multiple alignment programs
- **MUSCLE 5 is supported** (`multi_aln(tool = "muscle")`, `dNdS(aa_aln_tool = "muscle")`). The version is detected;
  MUSCLE 5 is run as `muscle -align <in> -output <out>`, its aligned FASTA output is put back into the order of the
  input. MUSCLE 3.8 is run as before.
- **T-Coffee: the default mode changed.** `multi_aln(tool = "t_coffee")` ran `t_coffee -mode accurate`, which needs
  a BLAST server and failed (unnoticed) without network. The default is now T-Coffee's own default mode
  (`-output clustalw_aln`, guide tree not written to the working directory). Alignments, and dN/dS values computed
  from them, differ from those of `-mode accurate`; `params = "-mode accurate"` gives the old command.
- `multi_aln(get_aln = TRUE)` reads the CLUSTAL files written by T-Coffee (Biostrings rejected them).
- The MAFFT check runs `mafft --version` (was `mafft -v`, not an option). The MUSCLE, MAFFT and KaKs_Calculator
  checks return `TRUE` (they returned `NULL`).

Other fixes
- `pairwise_aln()` works when called on its own (it needed a folder that only `dNdS()` created; fixed on
  `fix/robustness`, now tested).
- `read.genome(format = "fastq")` works (`object 'genome.dt' not found`).
- `set_blast()` / `set_diamond()` with `seq_type = "dna"`: the sequences are returned as read; `set_blast(makedb =
  TRUE, makedb_type = "nucleotide")` builds a nucleotide database from them (`blastdb_<file>_nucleotide.fasta`).
  Asking for a protein database from DNA is an error that says why (was `column 'aa' not found`).
- `codon_aln(get_aln = TRUE)` returned an empty alignment for `format = "clustal"` and failed for
  `format = "fasta"`; both are read now.
- `promotor_divergence_estimation()` works when all aligned pairs have the same length.
- `check_annotation(remove_annotation_outliers = TRUE)` **no longer overwrites the annotation file**. The file
  without the outlier lines is written to `tempdir()` under the same name (the extension was lost for names with
  more than one dot) and its path is returned. `retrieve_longest_isoforms()` uses that copy.
- `lnc_map_counts()` accepts the species column of the package's ortholog tables (`subject_species`) as well as
  `species`.
- **`orthofinder2()` is soft-deprecated** (a warning on every call). It did not work: the executable path was
  hard-coded to `/opt/miniconda3/bin/orthofinder`, it always returned `NULL` and the imports passed folders as
  files. It now looks for `orthofinder` on the `PATH` (or in `of_path`), checks the exit status, returns the
  results folder invisibly, and `import_type = "orthogroups_core"` reads the tables of OrthoFinder 2.5 and 3 when a
  full run has written them. `import_type = "orthogroups_pairwise"` (never implemented) is an error.
- `orthologs(ortho_detection = "Orthofinder2")` is an error with an explanation (it printed "under development"
  and returned `NULL`).

### External programs are run through one checked helper

- **No more stale results.** `blast()`, `blast_best()`, `blast_rec()`, `diamond()`, `deepclust()`,
  `deepclust_realign()`, `multi_aln()`, `codon_aln()` and `substitutionrate()` start their programs through one
  internal helper (`run_tool()`), which quotes all arguments, checks the exit status, checks that the expected
  output was written and puts the program's own message into the error. Before, a failed BLAST, MSA or
  KaKs_Calculator run could go unnoticed and the result of the *previous* call was returned: e.g.
  `blast(..., blast_params = "-max_target_seqs 1")` (an example of the man page, which BLAST rejects) returned the
  hit table of the call before it, and `substitutionrate(est.method = "kaks_calc")` returned the rows of an earlier
  run. These calls now stop with the program's message.
- **Every call works in files of its own.** `blast()`, `diamond()`, `deepclust()` and `deepclust_realign()` work
  in a folder that is created for the call and removed when it ends (kept with `options(orthologr.keep_tmp = TRUE)`);
  `multi_aln()`, `codon_aln()` and `substitutionrate()` let the program write to a uniquely named file. Results
  that were documented to stay in `tempdir()` still do, under the same names: the hit table
  (`_blast_db/blastresult_<query file>.csv`), alignments (`_alignment/multi_aln`, `_alignment/codon_aln`) and
  databases built with `set_blast()` / `set_diamond()`. A file of an earlier call is removed before the program
  starts and replaced only when the program has succeeded. `substitutionrate()` no longer leaves `*_gestimout`,
  `.axt` and `.axt.kaks` files in `_calculation`.
- **File and folder names with spaces work** in all of the functions above and in `dNdS()` (`makeblastdb`,
  `diamond makedb` and the searches failed before). `set_blast()` / `set_diamond()` replace characters other than
  letters, digits, `.`, `_`, `+`, `-` in the *database* name they return (names without such characters are unchanged).
- **The R session is left alone.** No function changes the working directory any more (`setwd()` is gone from
  `blast()`, `diamond()`, `set_blast()`, `set_diamond()`, `deepclust()`, `deepclust_realign()`), and a `path`
  argument no longer adds that folder to the `PATH` of the session for good: it is used to start the program and
  is on the `PATH` only while the program runs. As a consequence `blast(path = )` now also finds `makeblastdb`
  there, and `multi_aln(path = )` / `substitutionrate(kaks_calc_path = )` work when the program is *only* in that
  folder (the availability check ignored the argument).
- **Programs are looked up once per session.** The "is it installed" checks started the program on every call
  (2-150 ms each, once per gene pair in `dNdS()` with a multiple-alignment tool or a KaKs_Calculator model); the
  answer is now remembered per program and path.
- `select_orthologs()` no longer stores its joined table in `tempdir()` under a name made of the two annotation
  file names. A second call with the same annotation files but another `dnds_tbl` returned the table of the first
  call.
- `save.output` of `blast()`, `diamond()`, `deepclust()` and `deepclust_realign()` overwrites a file of the same
  name (it was silently left as it was).
- `multi_aln()` and `codon_aln()` return the path of the alignment file invisibly when `get_aln = FALSE` (was
  `NULL`). `multi_aln()` stops with a clear message when `file` does not exist.
- `promotor_divergence_estimation()` writes its intermediate FASTA files to a folder of its own and removes it.

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
