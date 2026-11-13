## `orthologr` 0.5.0

Version 0.5.0 instead of 0.4.4 because results and behaviour change in ways a patch release should not bring:
some default results change (see "Changes to pairwise alignments" where that change is merged, the T-Coffee
default mode, `aa_aln_params` being passed on), functions that returned `NULL` or stale results now stop with an
error, and packages moved from Imports to Suggests.

### Sequence IDs, joins, reported drops (second audit)

Fixes for defects found by a second audit of the package on real data (human, mouse, yeast, *E. coli*, plant
genomes). None of the changes in this section alters dN, dS or dN/dS values or which orthologs are found for
input that worked before; what changes is listed under "Changes of results and row counts".

Sequence IDs
- **The IDs of a gene pair are carried through the dN/dS pipeline** instead of being read back from the output
  of the external programs, and `dNdS()` joins its result with the hit table by `query_id` *and* `subject_id`
  (was `query_id`). This fixes, structurally, every case in which a program or a file reader changed an ID on
  the way:
  - `aa_aln_tool = "mafft"`: MAFFT cuts names to 15 characters in CLUSTAL output. Query IDs longer than that
    gave a result with **0 rows without error or warning**; longer subject IDs (e.g. Ensembl mouse
    `ENSMUST00000012345`) gave truncated IDs that are in no input file. `multi_aln(tool = "mafft")` also passes
    `--namelength` now, so the alignment file itself keeps the IDs.
  - IDs that look like numbers: the table written by gestimator (Comeron's method) was read with
    `data.table::fread()` type guessing. Purely numeric IDs (Entrez GeneIDs, JGI protein IDs) stopped `dNdS()`
    in its last join after all pairs had been computed ("incompatible types"), and IDs such as `007`, `1e5`,
    `1.10`, `NA` were rewritten (`7`, `1e+05`, `1.1`, missing) and vanished from the result. All tables written
    by programs are read with fixed column types (IDs as text, no `NA` strings, no quote characters).
  - (IDs with a hyphen and the KaKs_Calculator methods were fixed on `fix/robustness`; that fix is now a
    consequence of the same rule.)
- **Duplicated sequence IDs are an error that names them** (`cds2aa()`, `set_diamond()`, `set_blast()`, `dNdS()`,
  ...; for coding sequences and proteins). They used to end in `Wrong argument type in s2c()` / "could not be
  translated properly".
- The ID of a sequence is its FASTA header up to the first white space; it was cut at the first *space*, so
  `>g1<TAB>description` gave the ID `g1\tdescription` (`read.cds()`, `read.proteome()`, `read.genome()`).

Sequences that are left out are reported
- `dNdS()` returns `attr(result, "dropped_sequences")`: a `data.frame` (`file`, `id`, `reason`) of the input
  sequences that were left out before the search, and a message gives their number per reason. The message
  about coding sequences with characters other than A, C, G, T names the file, the number of sequences read and
  how many of them hold `N` only (5'-incomplete Ensembl CDS are padded with `N`: 9,035 of the 245,535 records of
  the human `cds.all` file of Ensembl 115 are left out for this reason). Empty sequences, which were left out
  silently, are reported too.
- With `delete_corrupt_cds = FALSE` sequences whose length is not a multiple of 3 are trimmed at their 3' end, as
  before; the message now says that a 5'-incomplete sequence (not starting in phase 0) is then translated in the
  wrong reading frame.

Splice variants
- **New warning** when an input file of `dNdS()`, `divergence_stratigraphy()` or `orthologs()` seems to hold
  several transcripts per gene (every sequence is treated as a gene: human chr21 with all isoforms gave 228 RBH
  rows for 168 genes, with dN/dS differing by more than a factor of 2 within 6 genes). The heuristic is
  conservative and documented in `?dNdS` (section "Splice variants"): gene names in Ensembl / NCBI FASTA headers,
  otherwise transcript suffixes of the IDs when at least 5 and at least 1% of the ID stems repeat. `dNdS()` also
  says how many rows per query gene the result holds. `options(orthologr.isoform_warning = FALSE)` switches the
  check off. Results are not changed.

Searches and arguments
- A search without any hit returns a hit table with no rows and the usual columns (`diamond()`, `blast()`,
  `diamond_best()`, `diamond_rec()`, ...) and says "found no hits"; it was reported as "the hit table could not
  be read" with a `data.table` message attached. `dNdS()` then stops with "No orthologs could be found: the
  DIAMOND search of ... returned no reciprocal best hits at e-value 1E-5".
- `substitutionrate(kaks_calc.params = )` **is added to the method instead of replacing it**: `est.method = "YN"`
  with `kaks_calc.params = "-c 1"` ran KaKs_Calculator's default model (MA) and returned its values under the
  heading of the call. `est.method = "kaks_calc"` still takes the methods from `kaks_calc.params`.
- `eval` must be a number (or a number as text); anything else is rejected before a program is started.
- The message shown when DIAMOND is not found recommends `orthologr::diamond_install()` and the arguments `path` /
  `aligner_path` (it named an argument `diamond_exec_path` that no exported function has).
- `map_generator_dnds()` and `map_generator_lnc()` with **one** subject file and the default `progress_bar = TRUE`
  stopped with `object 'pb' not found` after the run.
- `orthologs(seq_type = "cds")` no longer writes `<file>_translated.fasta` copies (they held the untranslated
  DNA) into `tempdir()`; the searches read the input files. Same result.
- `plot_diverse_homology_thresholds()` and `plot_diverse_homology_thresholds_core_orthologs()` compute the
  `scope` column when the table does not have it (tables of `import_ortholog_tables_all()`,
  `select_orthologs()`); `plot_diverse_homology_thresholds()`, `plot_pairwise_orthologs()` and `lnc_map_counts()`
  stop when `species_order` does not name a species of the table (it was plotted at `NA`).

Changes of results and row counts in this section
- `dNdS(aa_aln_tool = "mafft")` with IDs longer than 15 characters: rows that were lost (0 of 8 pairs returned
  in the test of the audit) are returned (8 of 8), and `subject_id` holds the input IDs. No option restores the
  old behaviour (it was a loss of data).
- `dNdS()` with numeric-looking IDs: 6 of 6 instead of an error or 2 of 6 rows in the audit's test. No option.
- `substitutionrate(est.method = <method>, kaks_calc.params = <text without -m>)`: values are now those of
  `<method>`, not of model averaging (MA). The old call is `est.method = "kaks_calc"`,
  `kaks_calc.params = "<text>"`.
- Input files with duplicated IDs, and `eval` values that are not numbers, are now errors.

### Dependencies, documentation, continuous integration

Dependencies
- **Lighter installation.** Imports went from 21 to 11 packages. `ggplot2`, `scales` (the `plot_*()` functions),
  `ape` (`promotor_divergence_estimation()`), `rtracklayer` (`select_orthologs()`, `retrieve_longest_isoforms()`,
  `generate_ortholog_tables()`) and `benchmarkme` (memory check of `check_annotation()`) are now in Suggests: the
  functions say which package to install when it is missing. `fs` (replaced by `file.exists()` / `dir.exists()`),
  `Rdpack` (one reference, now plain text), `R.utils`, `foreach` and `doParallel` are no longer used. Unused
  Suggests (`DBI`, `IRanges`, `RSQLite`, `devtools`, `pkgload`) were dropped; `tidyr` and `biomartr`, which the
  examples use, were added.
- **R >= 4.4 instead of (in effect) R >= 4.5.** `pwalign (>= 1.6.0)` (Bioconductor 3.22) forced R 4.5 although
  the DESCRIPTION said R >= 4.0. The functions used (`pairwiseAlignment()`, `alignedPattern()`,
  `alignedSubject()`, `pattern()`, `subject()`, `writePairwiseAlignments()`) are exported by `pwalign` 1.0.0
  (Bioconductor 3.19, R 4.4), so the version requirement was dropped and `Depends: R (>= 4.4.0)`,
  `Biostrings (>= 2.72.0)` state what is needed. (Checked against the NAMESPACE of pwalign 1.0.0; the package
  itself was run on R 4.5 only - the `oldrel-1` job of the check workflow covers R 4.4.)
- `SystemRequirements` lists Perl and DIAMOND and the optional programs; the `C++11` entry is gone (R CMD check
  NOTE).
- `NAMESPACE` imports three functions of `data.table` instead of the whole of `data.table` and `Biostrings`.

Documentation
- README: installation rewritten (the listed Bioconductor packages `GenomicRanges`, `GenomicFeatures`,
  `Rsamtools` and the GitHub package `rdiamond` are not needed; DIAMOND and Perl, which are, were not mentioned),
  the lists of orthology and dN/dS methods include the defaults (DIAMOND reciprocal best hit, Comeron), the
  example output was regenerated, the Windows section no longer refers to `?build_github_devtools` and
  `install.packages("rtools")`, the function overview covers the DIAMOND, DeepClust and stratigraphy functions.
  The quoted libsequence licence header now reads "version 3" as the files in `src/` do.
- Vignettes: `Install` (R version, BLAST+ 2.2.31 example with a wrong file name and `usr/local/bin`, ClustalW2
  executable name, MUSCLE 3.8/5, Li's method is from `seqinr` not `ape`, KaKs_Calculator 2.0 from bioconda),
  `dNdS_estimation` (the pairwise aligner is Needleman-Wunsch, the default method is Comeron, all methods listed,
  `aligner_params` instead of the non-existing `blast_params`), `sequence_alignments` (`clustalw2`, MUSCLE,
  `mafft --version`, a `pairwise_aln()` example that is run), `orthology_inference` (OrthoFinder), `orthologr`
  (installation). `deepclust.Rmd` used `comp_cores = 4` in a chunk that is run, which stopped `R CMD build` on
  machines with fewer cores.
- Man pages: `blast()` example with `blast_params = "-max_target_seqs 1"` (BLAST stops: the option is set by
  `max.target.seqs`), `format = "gbk"` (never supported), `orthologs()` arguments `outgroup_file` and `cdd.path`
  (deprecated; a warning when set), `diamond_install()` (no checksum, no Linux aarch64 build).
- `environment.yml` replaces `orthologr.yml` (R 3.5.1, python 2.7, macOS build hashes) and
  `orthologr-conda-env.sh` (defunct `biobuilds` channel): R >= 4.4, the R packages and all external programs from
  conda-forge and bioconda; it solves on linux-64 (micromamba, 289 packages, R 4.5.3, DIAMOND 2.2.8, BLAST+ 2.17.0).
- Removed `.travis.yml` and `Read-and-delete-me`.

Continuous integration and tests
- The R-CMD-check workflow installs DIAMOND, BLAST+, Clustal Omega, MAFFT and KaKs_Calculator (bioconda; on
  Windows DIAMOND from its release), so the tests that were skipped without them run.
- A pkgdown workflow and `_pkgdown.yml` (grouped reference index) were added.
- 16 of the 20 empty test files now hold tests (`multi_aln`, `orthologs`, `select_orthologs`,
  `retrieve_longest_isoforms`, `read.genome`, `read.proteome`, `read.dnds.tbl`, `set_blast`, `substitutionrate`,
  `codon_aln`, `cds2aa`, `filter_best_hits`, `extract_features`, `generate_ortholog_table`, `orthofinder2`,
  `orthofinder2_retrieve_core_orthologs`). Still empty: `map_generator`, `map_generator_lnc`, `orthologs_lnc`
  (they need long runs or the GitHub package `metablastr`) and `pairwise_aln` (its tests come with
  `fix/pairwise-aln-speed`).

### Changes to pairwise alignments that change dN/dS values (methods change)

`dNdS()`, `divergence_stratigraphy()` and `pairwise_aln()` align each pair of orthologous proteins with
`pwalign::pairwiseAlignment(type = "global")`. Three things about this step have changed. Each of them changes
some alignments and therefore the dN, dS and dN/dS values of the affected gene pairs; results computed
with this version are not numerically identical to results from `orthologr` <= 0.4.3. Please report the
`orthologr` version (and, if you changed it, `aa_sub_matrix`) in the methods of analyses based on these values.

**How large the change is.** Measured on the 1,000-gene example of the package (`ortho_thal_cds_1000.fasta` vs
`ortho_lyra_cds_1000.fasta`, 984 reciprocal best hit pairs), 0.4.3 against this version:

| | Comeron (default) | YN (KaKs_Calculator) |
|---|---|---|
| pairs with unchanged dN and dS | 66% (651 of 984) | 71% (699 of 984) |
| Pearson r of dN, old vs new | 0.971 | 0.967 |
| Pearson r of dS, old vs new | 0.940 | 0.997 |
| Pearson r of dN/dS, old vs new | 0.988 | 0.984 |
| largest change of dN / dS / dN/dS | 0.19 / 0.41 / 0.31 | 0.19 / 0.48 / 0.44 |
| pairs whose dN/dS changes by more than 0.01 | 22% | 20% |

So about a third of the gene pairs get other values. For most of them the change is small (median change of dN
0.002, of dS 0.004-0.005 among the changed pairs), but for single pairs it is large, and rankings or thresholds
on dN/dS (e.g. `filter_dNdS()`, divergence strata) can move. `options(orthologr.legacy_pairwise_aln = TRUE)`
reproduces the values of 0.4.3.

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
### Faster R code around the external programs (no change of results)

The following changes do not change any result: the old per-gene code is kept as the reference and the test suite
compares both on the example data, including the 1,000-gene example (`tests/testthat/test-vectorise.R`).
R time (CPU of the R process, 1,000 x 1,000 genes; elapsed times depend on DIAMOND):

| | before | after |
|---|---|---|
| `set_diamond()` on CDS | 0.47 s | 0.15 s |
| `diamond_best()` on CDS | 3.3 s | 0.6 s |
| `diamond_rec()` on CDS | 5.0 s | 0.5 s |
| `diamond_rec()` on proteins | 3.7 s | 0.1 s |
| `select_orthologs()` (gtf / gff, 986 pairs) | 14.1 s / 12.3 s | 0.11 s / 0.15 s |
| `retrieve_longest_isoforms()` (gff, without loading rtracklayer) | 1.3 s | 0.2 s |

- **Translation**: coding sequences are translated in one call to `Biostrings::translate()` instead of one
  `seqinr::translate()` call per gene. Only sequences of A, C, G, T with complete codons take the new path, with
  `no.init.codon = TRUE` (Biostrings would otherwise translate an alternative start codon such as TTG or CTG at the
  first position to M, which seqinr does not); every other sequence is translated as before.
- **Each input file is read and translated once** in `dNdS()` (was six times), `diamond_rec()` and `blast_rec()`
  (was four times). As a side effect the messages of `read.cds()` about incomplete codons appear once per file.
- **Best hits** (`diamond_best()`, `blast_best()`, and the splice-variant selection of `select_orthologs()`): one
  sort instead of a `dplyr::do()` call per query. Ties are broken as before: smallest e-value, then longest
  alignment, then the first row of the hit table; the order of the result is that of `dplyr::group_by()`.
- **`select_orthologs()`**: the gene locus / transcript table of an annotation is built in one pass instead of one
  `extract_features()` call per locus (98% of the run time). Annotations in which a locus has more than one `gene`
  row take the old path. One difference: an annotation with a locus that has transcripts but no `gene` row (or
  the other way round) stopped with an obscure `vctrs` type error; such loci are now left out, as intended.
- **`retrieve_longest_isoforms()`**: one sort instead of a `dplyr::do()` call per gene (the first of equally long
  isoforms, as before).
- **`deepclust()` / `deepclust_realign()`** with several protein files read each file once with Biostrings instead
  of reading every sequence with seqinr, writing one merged file and reading it again. The database is the same
  (sequences sorted by ID).
- The protein FASTA files for the searches are written with one `writeLines()` call (same bytes as before).
- **Parallel runs (`comp_cores > 1`) use the `parallel` package only**: forked processes as before where the
  operating system can fork, and otherwise - on Windows, where `comp_cores > 1` failed - `PSOCK` processes that
  load the package and get the options and the working folder of the call
  (`options(orthologr.cluster_type = "PSOCK")` chooses them elsewhere). `foreach` and `doParallel` are no longer
  imported.

New, off by default
- `deepclust(direct_input = TRUE)` hands the protein files to `diamond makedb` without parsing them in R (no
  sequence is held in memory; for very large sets). The IDs are the same, but `diamond deepclust` depends on the
  order of the sequences: with file order instead of ID order the clusters of the 1,988 example proteins differ
  from the default. This is why it is not the default.

Bug fix
- `deepclust()` failed for IDs with `|` and `:` (e.g. the four `example_*_aa.faa` files together): the cluster
  table was read with a guessed separator (`Can't assign 2 names to a 3-column data.table`). It is read as
  tab-separated now.

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
