## `orthologr` 0.6.0 (in development)

### Performance (results unchanged)

The dN/dS step of `dNdS()` no longer spends most of its time on work around the estimates. Every result table was
compared with that of 0.5.0 (all columns, the attributes and the record of the run): identical for the 20-gene
example, the 1,000-gene example and S. cerevisiae / C. glabrata, with the default settings, `engine = "coal"`,
`codon_aln_tool = "pal2nal"`, the methods Comeron, Li and YN (KaKs_Calculator and `kaks_tool = "coal"`), 1 and 2
cores, and under `orthologr_legacy_options()` (which still reproduces the stored 0.4.3 tables).

- Pairwise protein alignments (`aa_aln_type = "pairwise"`, the default): the proteins of up to 100 gene pairs are
  aligned with one `pwalign::pairwiseAlignment()` call (pwalign aligns the i-th pattern with the i-th subject:
  the same alignments, scores and tie-breaks as one call per pair) instead of building an `AAString` object and
  reading a FASTA file per pair. Not used with `store_locally = TRUE`, `aa_sub_matrix = NULL` and the legacy
  pairwise alignment (`options(orthologr.legacy_pairwise_aln = TRUE)`): those pairs are aligned one by one as before,
  and so is a pair the joint call cannot handle (above `options(orthologr.max_aln_cells)`, unusual characters).
- With COAL as codon aligner (the default when `coalr` is installed) the codon alignment is built from these rows
  in memory instead of reading the two files of the pair back.
- Comeron's method (`gestimator()`, C++): the redundancy tables of the genetic code and the Grantham distances
  are built once per R process instead of once per gene pair and once per pair of codons that differ at two or
  three positions: about 10 times faster on divergent genes, with the same arithmetic.
- Per gene pair: no `data.table` sub-assignments (`:=`), the working files written with one `writeLines()` (the
  same bytes as `seqinr::write.fasta()`); the checks of coding sequences for non-ACGT characters use PCRE (the
  same answers, about 50 times faster on long sequences).
- PAL2NAL and KaKs_Calculator: one process for a chunk of up to 100 gene pairs instead of one per pair.
  KaKs_Calculator gets one AXT file with a record per pair (it computes every record on its own) and its rows are
  split by record; PAL2NAL runs through `inst/pal2nal/pal2nal_batch.pl`, which compiles `pal2nal.pl` once and
  computes every pair in a forked child process (a fresh copy of the script's state, the arguments and output file
  of a single call; not on Windows). A pair the batch does not compute - its job failed, its rows cannot be matched
  to its record one to one, the program failed - is computed on its own as before, which also gives its error.
  Not under `orthologr_legacy_options()` (one pair at a time, as in 0.4.3); `options(orthologr.batch_tools =
  FALSE)` switches it off. Not for GMYN: KaKs_Calculator 2.0 (bioconda 2.0.1) decides whether it writes the header
  line of that method with a class member it never initialises (`tempt`), so its output depends on that memory
  (on a CI runner most pairs came back without a result with one process per pair, and other pairs in a batch).
  Using COAL for the codon step of `codon_aln_tool = "pal2nal"` (its codon alignments are the same) was not done:
  that they are the same for a pair can only be known by running PAL2NAL.
  Compared with one process per pair (all columns, attributes, `failed_pairs`, the record of the run), identical
  with 1 and 2 cores: the 20-gene example for every method of KaKs_Calculator (NG, LWL, LPB, MLWL, MLPB, YN, MYN,
  GY, MS, MA, ALL), PAL2NAL with Comeron and YN, and the legacy settings; the 1,000-gene example and
  S. cerevisiae / C. glabrata for NG, LWL, LPB, MLWL, MLPB, YN, MYN (with GY on the 1,000-gene example; yeast:
  the 53 MYN pairs that fail with `nan` included) and PAL2NAL with Comeron and YN; the dN/dS step on
  40 pairs of the 1,000-gene example for MS, MA and ALL, and on yeast pairs for GY (200) and MA (12).
- Measured: the whole `dNdS()` call including DIAMOND, CPU seconds (user + system, with child processes; the
  machine was shared, so elapsed times varied by a factor of 2), 0.5.0 -> 0.6.0:

  | data, settings | 1 core | 2 cores |
  |---|---|---|
  | 1,000-gene example, default (Comeron, COAL codon step) | 29.2 -> 12.1 | 32.5 -> 15.0 |
  | 1,000-gene example, `codon_aln_tool = "pal2nal"` | 43.1 -> 33.9 | |
  | 1,000-gene example, Li | 26.6 -> 14.6 | |
  | 1,000-gene example, YN (KaKs_Calculator) | 33.0 -> 21.1 | |
  | 1,000-gene example, YN with `kaks_tool = "coal"` | 25.5 -> 14.1 | |
  | S. cerevisiae / C. glabrata (4,384 pairs), default | 106.6 -> 53.9 | 122.1 -> 50.6 |
  | S. cerevisiae / C. glabrata, YN / Li / PAL2NAL | | 153 -> 98 / 108 -> 59 / 197 -> 142 |

  The dN/dS step alone (`compute_dnds()`, 1,000-gene example, 1 core): 27.4 s -> 10.4 s before the change of
  `gestimator()`. What is left per pair is mostly the alignments themselves (pwalign), the estimate, and for
  PAL2NAL / KaKs_Calculator the work of the programs themselves. One process per chunk instead of one per pair
  (CPU seconds of the whole `dNdS()` call, mean of 2 runs, before -> after; tables identical):

  | data, settings | CPU s |
  |---|---|
  | 1,000-gene example, YN (KaKs_Calculator), 1 core | 21.0 -> 14.6 |
  | 1,000-gene example, NG, 2 cores | 21.3 -> 14.2 |
  | 1,000-gene example, `codon_aln_tool = "pal2nal"` (Comeron), 1 core | 29.1 -> 19.7 |
  | 1,000-gene example, `codon_aln_tool = "pal2nal"`, YN, 1 core | 38.0 -> 21.7 |
  | S. cerevisiae / C. glabrata, YN, 2 cores | 90.0 -> 64.2 |
  | S. cerevisiae / C. glabrata, `codon_aln_tool = "pal2nal"`, 2 cores | 131.5 -> 90.0 |
  | 20-gene example, MA, 2 cores (the estimates take the time) | 125.8 -> 126.5 |

  The dN/dS step alone (`compute_dnds()`, 1,000-gene example, 1 core): YN 17.7 s -> 10.6 s, PAL2NAL 25.5 s ->
  16.7 s; what is left of PAL2NAL is its own work (about 9 ms per pair).
- `stringr` (and with it `stringi`) is no longer imported: its 27 calls were splitting and matching of file names
  and are base R now.

### Changed results

- Comeron's method (both engines; COAL #76): the rule of 0.5.0 that a site class with fewer than
  `options(orthologr.comeron_min_sites)` (3) sites cannot make dN or dS `NA` now applies only where the value that
  results is greater than 0. Setting such a class's undefined logarithm to 0 could leave dS (in principle dN) at
  exactly 0 although a change was seen (`ATGGCTAAA` / `ATGGCCAAA`, one synonymous change: dS 0), which reads as "no
  synonymous change" (status `"dS_zero"`, left out by `filter_dNdS()`). Such a value is `NA` again (status
  `"saturated"`), as with `options(orthologr.comeron_min_sites = 0)`. No value greater than 0 changes, and nothing
  changes under `orthologr_legacy_options()` (0.4.3's saturation rule). Seeded short pairs of 5 - 20 codons (the test
  `test-substitutionrate.R`, 160 pairs): 2 dS of 0 -> `NA`, the 9 other values the rule gives in place of `NA` unchanged.
  Real data: no value changes - the tables of the 20-gene and 1,000-gene examples and of S. cerevisiae /
  C. glabrata (4,384 pairs, where the rule gives 15 dN, all > 0) are identical to those of 0.5.0 (default
  settings, `engine = "coal"` with a `coalr` that has COAL #76, and the legacy switch).

### Fixed

- A `coalr` from before COAL #76 (whose `coalr::coal_info()$features` lacks `"comeron_min_sites_positive"`) gave
  Comeron values that differ from `engine = "orthologr"` without a warning when it already had the arguments
  `comeron_twosubs` / `comeron_min_sites` (e.g. dS = 0 instead of `NA` for `ATGGCTAAA` / `ATGGCCAAA`). The
  warning of class `orthologr_warning_coal_version` is now given whenever that feature is missing and Comeron is
  computed by COAL (`engine = "coal"`, `dNdS_coal()`, `dNdS_clusters()`, `dNdS_orthogroups(engine = "coal")`).

### Other changes

- `dNdS_orthogroups()`: the summary has the minimum and maximum of dN, dS and dN/dS per group (`dN_min`, `dN_max`,
  ...), as `dNdS_clusters()`, `coalr` and COAL's command line (the statistics median, mean, `_q25`, `_q75`,
  `_min`, `_max`, in that order).
- A `dnds_est.method` that is not a method (`dNdS()`, `substitutionrate()`, `dNdS_orthogroups()`, the COAL engine):
  the error lists the methods of each engine and names the method a mistyped name probably meant (ignoring case,
  or one character off: `"comeron"` -> `"Comeron"`, `"YNN"` -> `"YN"`); method names stay case-sensitive. The COAL
  engine says that a method "is only available with engine = \"orthologr\"" only for a method the native engine
  has (it said so for any name, `"YNN"` included); `"comeron"` got a generic message without the list.
- `library(orthologr)` notes (once per session, as a package startup message) that some defaults that change dN/dS
  values and orthologs are not those of 0.4.3, and points to `vignette("upgrading_from_0.4.3")` and
  `orthologr_legacy_options()`: a 0.4.3 script ran unchanged without any sign that its numbers changed. Not shown
  while knitr runs, when any `orthologr.*` option is set, or with `options(orthologr.quiet_upgrade = TRUE)` (or
  `suppressPackageStartupMessages()`); the option is not part of the record of a run.
- README: the links to pages that the package site (drostlab.github.io/orthologr) does not have yet - the
  articles "Upgrading from 0.4.3" and "Orthogroups" and `?orthologr_options` - point to the files in the repository.

### Internal

- Tests (`test-perf.R`): the batched code paths against the per-pair code (alignments, codon alignments, tables of
  `compute_dnds()` for Comeron, Li and YN with COAL and PAL2NAL), independence of the chunk size and of the number
  of processes, progress messages. `test-perf-tools.R`: PAL2NAL and KaKs_Calculator per chunk against one process
  per pair (every method of KaKs_Calculator that takes seconds, a pair that fails inside a chunk, a batch whose output
  cannot be matched, that fails or whose program is missing, chunk sizes, two processes, legacy mode), the driver
  against one `perl pal2nal.pl` per job.

## `orthologr` 0.5.0

Results change in this version: defaults were changed and several errors that changed numbers were fixed.
`orthologr_legacy_options()` (or `options(orthologr.legacy = TRUE)`) restores the behaviour of 0.4.3 in one line;
the article **"Upgrading from 0.4.3"** (`vignette("upgrading_from_0.4.3")`) gives the measured effect of every
change. Please report the version of orthologr in the methods of analyses based on dN/dS values.

### Upgrading from 0.4.3

Every default that changed, and what restores the behaviour of 0.4.3. All rows marked "switch" at once:
`orthologr_legacy_options()` (or `options(orthologr.legacy = TRUE)`); the stored 0.4.3 tables of the 20-gene
example (Comeron, YN) are reproduced exactly with it. Measurements: the table "Changed defaults" below.

| Default in 0.5.0 | 0.4.3 | Restore | Switch |
|---|---|---|---|
| best hits by bit score | e-value, then alignment length | `best_hit_rule = "evalue"`, `options(orthologr.best_hit_rule = "evalue")` | yes |
| DIAMOND `sensitivity_mode = "default"` | `"fast"` | `sensitivity_mode = "fast"`, `options(orthologr.sensitivity_mode = "fast")` | yes |
| pairwise protein alignments: BLOSUM62, full length, stop codon not aligned | no matrix, aligned range, stop aligned | `options(orthologr.legacy_pairwise_aln = TRUE)` | yes |
| no pwalign alignment above 5e8 cells (with `coalr` installed such a pair is aligned by COAL with pwalign's scoring: the same alignment) | no limit | `options(orthologr.max_aln_cells = Inf)` | yes |
| coding sequences with IUPAC ambiguity codes (N, ...) are kept, codons with them not compared | left out | `options(orthologr.ambiguous_cds = "drop")` | yes |
| Comeron: `NA` for a saturated dN / dS | a finite, too small value | `options(orthologr.legacy_comeron_saturation = TRUE)` | yes |
| Comeron: codons that differ at two positions counted along Comeron's two pathways | libsequence's branches (the estimate depended on which sequence is first) | `options(orthologr.comeron_twosubs = "libsequence")` | yes |
| Comeron: a site class with fewer than 3 sites cannot make dN / dS `NA` (0.6.0: unless the value would be 0) | (0.4.3 had no `NA` for saturation at all) | `options(orthologr.comeron_min_sites = 0)` | (not needed: covered by the saturation switch) |
| KaKs_Calculator methods: dN and dS of 0 are 0 | `NA` | `options(orthologr.legacy_kaks_zero_dn = TRUE)` | yes |
| KaKs_Calculator methods: capped values (99, JC69's cap, empty LWL classes) and `nan` are `NA` (`"saturated"` / `"failed"`) | returned as estimates | `options(orthologr.legacy_kaks_saturation = TRUE)` | yes |
| `filter_dNdS()` / `divergence_stratigraphy()`: pairs with dS = 0 left out | the same | `dS_zero = "keep"` keeps them | (unchanged) |
| dN/dS tables end with the column `status` | no such column | `options(orthologr.status_column = FALSE)` | yes |
| `divergence_stratigraphy(quiet = TRUE)` (progress line in interactive sessions) | `quiet = FALSE` | `quiet = FALSE` | no (messages only) |
| T-Coffee in its default mode | `-mode accurate` | `aa_aln_params = "-mode accurate"` | no (argument) |
| `aa_aln_params` passed to the alignment program | ignored | do not set it | no (argument) |
| with `coalr` installed, codon alignments by COAL (byte-identical to PAL2NAL's) | PAL2NAL | `options(orthologr.codon_aln_tool = "pal2nal")` | no (same result) |

### New

- **`orthogroups()`**: groups of orthologs of 2 to N species from one all-against-all DIAMOND search (reciprocal
  best hits, connected components) or from DIAMOND DeepClust clusters; `print()`, `summary()`,
  `write_orthogroups()`. **`dNdS_orthogroups()`**: dN/dS for the gene pairs within the groups, with a summary per
  group. OrthoFinder is no longer needed for more than two species; `orthofinder2()` is soft-deprecated.
- **Engine COAL** (optional R package `coalr`): `dNdS(engine = "coal")` / **`dNdS_coal()`** hand all gene pairs to
  [COAL](https://github.com/Plumat-COAL/COAL) in one call (protein alignment, codon alignment and estimation in
  memory, on `comp_cores` threads; 1,000-gene example, 1 core: about 6 s against 50 s); **`dNdS_clusters()`** computes
  dN/dS within DeepClust clusters. `codon_aln(tool = "coal")` and `substitutionrate(kaks_tool = "coal")` replace
  PAL2NAL / KaKs_Calculator for single steps. Methods: NG, LWL, MLWL, LPB, MLPB, YN, MYN, and Comeron and Li with
  a `coalr` that has them.
- **Column `status`** (last column of every dN/dS table): `"ok"`, `"dS_zero"`, `"saturated"`, `"no_data"` or
  `"failed"` for each gene pair - the same words with both engines and in the `pairs` tables of
  `dNdS_clusters()` / `dNdS_orthogroups()` (see `?dNdS`). `options(orthologr.status_column = FALSE)` leaves it out.
- **One summary message at the end of a run**: gene pairs by status, median dS, sequences left out before the
  search, and a note when most pairs are saturated or the median dS is above 1.5 (species too distant for
  pairwise dN/dS). With `quiet = FALSE`, a progress message about every tenth of the gene pairs.
- **Record of the run**: `attr(result, "orthologr_run")` on results of `dNdS()`, `dNdS_coal()`, `dNdS_clusters()`,
  `dNdS_orthogroups()` and `divergence_stratigraphy()`: package version, arguments, every `orthologr.*` option in
  force (with its default when it is not set, `orthologr.ambiguous_cds` included), engine (and the COAL backend
  used: `coalr` or the coal-cpp program), program versions (KaKs_Calculator's from the "Version:" line of its
  `-h`), date, size and MD5 of the input files, counts by status; the call as written (`call`) and, when it is of
  another function (`dNdS(engine = "coal")` runs `dNdS_coal()`), its name and all its arguments (`called_as`);
  the DIAMOND / BLAST+ command lines of the run (`commands`); with `engine = "coal"` what orthologr gave COAL
  (`engine_settings`: preset, matrix, tie-break, Comeron / Li flavour, terminal stops) and the equivalent
  `coal pipeline` command line (`engine_settings$coal_command`; `?dNdS_coal`, "The same run with the COAL
  command line").
  **`write_dNdS()`** writes a table as TSV with this record in `<file>.run.json`; **`read_dNdS()`** reads both back
  (needs `jsonlite`). `orthogroups()$params` holds the package version.
- **`orthologr_legacy_options()`** and `options(orthologr.legacy = TRUE)`: all settings of 0.4.3 at once. A test
  compares the result with tables stored from 0.4.3 (20-gene example, Comeron and YN). Only the default engine
  gives the numbers of 0.4.3: the COAL engine cannot make the pairwise alignments of 0.4.3 and with these settings
  warns (`orthologr_warning_legacy_engine`; 20-gene example: 10 of 20 pairs differ from 0.4.3, largest
  difference in dN/dS 0.147, r = 0.984).
- **`?orthologr_options`**: every option and environment variable the package reads.
- **Classed errors and warnings**: the common failures have the class `orthologr_error` plus one for their kind -
  `orthologr_error_tool_missing`, `_package_missing`, `_tool_failed` (with `_alignment_failed` for an alignment
  program), `_no_hits`, `_bad_input`, `_bad_argument` - and the warnings about failed gene pairs and isoforms the
  classes `orthologr_warning_failed_pairs` / `_isoforms` (and `_legacy_engine`, above), so that they can be handled
  with `tryCatch()` without matching message text (`?orthologr_options`, section "Errors and warnings"). Messages
  are unchanged. The argument and input-file checks of `dNdS()`, `dNdS_coal()`, `dNdS_clusters()`, `orthogroups()`,
  `dNdS_orthogroups()`, `divergence_stratigraphy()`, `filter_dNdS()`, the alignment functions, `substitutionrate()`
  and the readers are classed (17 probes of the review: 0 classed before, 17 after); a missing COAL is
  `orthologr_error_package_missing` (with `_tool_missing` when neither `coalr` nor coal-cpp is there); `dNdS()` and
  `dNdS_coal()` check that their input files exist before any program runs.
- `divergence_left_out()`: the gene pairs that are not in a divergence map, with the reason.
- `attr(result, "failed_pairs")` and `attr(result, "dropped_sequences")` on dN/dS tables: gene pairs that could
  not be processed and input sequences that were left out, each with the reason.
- A warning when an input file of `dNdS()`, `divergence_stratigraphy()`, `orthologs()` or `orthogroups()` seems
  to hold several transcripts per gene (`options(orthologr.isoform_warning = FALSE)` switches it off).
- `aa_sub_matrix` (`dNdS()`, `pairwise_aln()`): the substitution matrix of pairwise protein alignments.
- `options(orthologr.ambiguous_cds)`: coding sequences with IUPAC ambiguity codes (N, R, Y, ...) are kept (a codon
  with such a base is translated to X; the default `"keep"`, see "Changed defaults") or left out (`"drop"`, as in
  0.4.3).
- `options(orthologr.sensitivity_mode)`, `options(orthologr.comeron_twosubs)`, `options(orthologr.comeron_min_sites)`,
  `options(orthologr.large_pair_aligner)`: see "Changed defaults" and `?orthologr_options`.
- `best_hit_rule` (`dNdS()`, `orthologs()`, `diamond_best()`, `diamond_rec()`, `blast_best()`, `blast_rec()`,
  `divergence_stratigraphy()`, `select_orthologs()`, `orthogroups()`).
- `dS_zero` (`filter_dNdS()`, `divergence_stratigraphy()`): `"keep"` returns gene pairs with dS = 0 instead of
  leaving them out.
- MUSCLE 5 is supported; `"MLPB"` is an accepted `dnds_est.method`.
- `environment.yml` (conda environment with R, the R packages and the external programs); a pkgdown site
  configuration; articles on orthogroups and on the setup on Windows.

### Changed defaults (and how to get the old behaviour)

| Change | Effect (measured) | Old behaviour |
|---|---|---|
| DIAMOND's sensitivity mode is `"default"` in `dNdS()`, `orthologs()`, `orthogroups()`, `dNdS_coal()`, `divergence_stratigraphy()`, `diamond()`, `diamond_best()`, `diamond_rec()` (was `"fast"`, which DIAMOND designs for hits of more than 90 % identity; DIAMOND's own default mode is for hits of more than 60 %); the hint for distant species now suggests `"sensitive"` | reciprocal best hits as `dNdS()` searches (e-value 1e-5, 2 threads), references OrthoFinder 3.2 (four yeasts) and Ensembl Compara 115 (human / mouse): S. cerevisiae against C. glabrata / K. lactis / A. gossypii 4,165 / 3,846 / 3,695 -> 4,385 / 4,206 / 4,044 pairs (+5.3 to +9.4 %), recall of the 1:1 orthologs 88.9 / 85.7 / 84.7 % -> 93.8 / 94.3 / 93.2 %, precision 99.45 / 99.66 / 99.70 % -> 99.38 / 99.67 / 99.78 %; human / mouse 16,943 -> 17,104 (+0.9 %), recall 97.45 -> 98.14 %; CPU time of the search x 1.3 - 1.6 (human / mouse 81 -> 131 CPU s). Re-measured with this version's `diamond_rec()`: S. cerevisiae / C. glabrata 4,164 -> 4,384 pairs, 9.7 -> 13.6 CPU s; 1,000-gene example 985 -> 987 | `sensitivity_mode = "fast"` or `options(orthologr.sensitivity_mode = "fast")` |
| Best hits are chosen by bit score (was: e-value, then alignment length) | human / mouse: 3.2 % of the best hits change; pairs that contradict Ensembl Compara fall from 2.0 % to 0.4 % | `best_hit_rule = "evalue"` or `options(orthologr.best_hit_rule = "evalue")` |
| Pairwise protein alignments: BLOSUM62, full length, stop codon not aligned | 1,000-gene example: dN or dS changes for a third of the gene pairs, mostly slightly (r of dN/dS old vs new 0.988) | `options(orthologr.legacy_pairwise_aln = TRUE)` |
| Comeron's method: `NA` for a saturated dN or dS (was: a finite, too small value) | 1 dN of 995 pairs (example); 17 dN and 1 dS of 1,500 pairs (human / mouse) | `options(orthologr.legacy_comeron_saturation = TRUE)` |
| KaKs_Calculator methods: `dN = 0` for a pair without nonsynonymous change (was `NA`, with `dNdS = 0`) | 11 of 995 pairs; `filter_dNdS()` keeps them | `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| KaKs_Calculator methods: dN and dS that are 0 are 0 (KaKs_Calculator prints `NA` below 1e-6): identical sequences have dN = dS = 0 and status `"dS_zero"` (was `NA` and `"saturated"`), a pair with nonsynonymous differences only has dS = 0 - the values and status of the COAL engine | E. coli K-12 / O157, YN (3,663 pairs): 230 pairs `"saturated"` -> `"dS_zero"` (the summary no longer counts them as saturated), 12 dN of about 1e-12 -> 0; a file against itself (20 genes): "20 saturated" and the note that the species are probably too distant -> "20 with dS = 0"; 1,000-gene example: 0 pairs | `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| KaKs_Calculator methods: an estimate that KaKs_Calculator replaced with a capped value because the correction for multiple substitutions is not defined is `NA`, status `"saturated"` (was returned as an estimate, status `"ok"`): Ka or Ks = 99, JC69's cap 0.75 log(sites) of YN / MYN, a class of sites counted with 0 substitutions in LWL / MLWL / LPB / MLPB (`?substitutionrate`; where the printed counts are within rounding of 0, the counts are computed again from the codon alignment with KaKs_Calculator's arithmetic, so that the status is what KaKs_Calculator did) | S. cerevisiae / C. glabrata, YN (4,164 pairs): `"ok"` 4,164 -> 605, `"saturated"` 0 -> 3,559 (172 with dS = 99, whose median dN/dS was 0.0028 against 0.081 overall, and 3,387 at JC69's cap, dS 2.56 - 5.99); `filter_dNdS(dnds.threshold = 2)` kept 4,164 rows (the 172 with dS = 99 included) -> 605; E. coli K-12 / O157, YN: 21 of 3,663 pairs (divergent pairs of 35 - 88 % identity, dS 2.8 - 4.8); 1,000-gene example: YN 1 of 995, NG 0; LWL family, S. cerevisiae / C. glabrata: Ks `NA` for 1,925 (LWL, LPB) and 2,100 (MLWL, MLPB) of 4,164 pairs - every pair whose class KaKs_Calculator found undefined, and no other (checked against a KaKs_Calculator that reports its decision, also on S. cerevisiae / K. lactis MLPB and E. coli; with the printed counts alone an intermediate version took 21 pairs the wrong way, among them 2 MLPB Ks of 0.97 and 2.26 that are underestimates) | `options(orthologr.legacy_kaks_saturation = TRUE)` |
| KaKs_Calculator methods: an estimate printed as `nan` (Ka or Ks) is a failed estimate - `NA`, status `"failed"`, the reason with the method in `attr(result, "failed_pairs")` and the failed-pairs warning (was `NaN`, in one pair with a dN, and status `"saturated"`); `substitutionrate()` lists it in `attr(, "kaks_failed")` | S. cerevisiae / C. glabrata, MYN (4,164 pairs): 50 pairs `"saturated"` -> `"failed"` (all of them iterations that did not converge, not saturation: COAL estimates 43 of them, dS 0.12 - 4.1), the dN of one pair 0.0547 -> `NA`; YN, NG and the LWL family on the same pairs, and E. coli K-12 / O157: no `nan` | `options(orthologr.legacy_kaks_saturation = TRUE)` |
| dN/dS tables end with the column `status` | one more column | `options(orthologr.status_column = FALSE)` |
| T-Coffee runs in its default mode (was `-mode accurate`, which needs further programs) | other alignments with `aa_aln_tool = "t_coffee"` | `aa_aln_params = "-mode accurate"` |
| `aa_aln_params` is passed on to the alignment program (was ignored) | only when it is set | do not set it |
| A pairwise protein alignment of more than 5e8 cells (length x length; pwalign needs about 3 bytes per cell, so about 1.5 GB) is not computed by pwalign (was: aligned whatever the memory; `options(orthologr.max_aln_cells)`, `?pairwise_aln`). With `coalr` installed such a pair is aligned by COAL's aligner with its `"pwalign"` preset (pwalign's scoring and tie-break; a message names the pair); without it the pair is returned with status `"failed"` and the reason | titin (35,991 x 35,390 residues, 1.27e9 cells) through `pairwise_aln()`: COAL (coalr) 6.4 - 10.6 s and a peak of 423 MB for the R process against pwalign (`max_aln_cells = Inf`) 16 - 18 s and 3,938 MB, identical rows (and score, 166,805); the alignment itself takes 4.4 s in coalr and 5.3 s with COAL's command line; human 1,000 transcripts / all mouse CDS: 1 of 946 pairs (titin); without `coalr` it fails instead of dS 0.50 (peak memory 4.3 GB -> 0.79 GB); 1,000-gene example: 0 pairs | `options(orthologr.max_aln_cells = Inf)` (pwalign for every pair); `options(orthologr.large_pair_aligner = "none")` keeps failing such pairs with `coalr` |
| `divergence_stratigraphy(quiet = TRUE)`, as `dNdS()` (was `FALSE`); in an interactive session it still shows a progress line about every tenth of the gene pairs (stratigraphy runs are long) | 1,000-gene example: 10 messages instead of 45 (in scripts); same result | `quiet = FALSE` |
| Coding sequences with IUPAC ambiguity codes (N, R, Y, ...) are kept (`options(orthologr.ambiguous_cds = "keep")`): a codon with such a base is translated to X and not compared (was: the sequence was left out); `dNdS()` lists them in `attr(, "ambiguous_sequences")` with the number of such codons and counts their rows in the summary, `dNdS_orthogroups()` (which left them out whatever the option) adds a note to their pairs | Ensembl, one transcript per gene: 87 of 19,746 human and 75 of 21,584 mouse coding sequences contain such characters (median share of ambiguous codons 0.4 %); one N in about 1 % of the codons of every human gene: Comeron's dS within 5 % of the N-free value for 99.1 % of 946 human / mouse pairs (YN 95.3 %); 1,000-gene example: 2 more pairs (with `sensitivity_mode = "default"` 999 pairs instead of 995) | `options(orthologr.ambiguous_cds = "drop")` |
| Comeron's method counts codons that differ at two positions along Comeron's two pathways (libsequence's `TwoSubs`, and so 0.4.3, passes its arguments in another order and counts codon1 - int1, int1 - int2, codon1 - codon2 and codon2 - int2, with which the estimate depends on which sequence is first; libsequence's repository is archived). The COAL engine does the same (`coalr` with `comeron_twosubs`) | swapping the two sequences of 300 random 150-codon pairs changed Ka in 279 and Ks in 190 of them; now in none. dS as `dNdS()` returns it (default settings, only this option switched): 1,000-gene example 451 of the 998 pairs with a dS both ways change, median 0 %, p95 2.9 %, max 12 %, no status changes; S. cerevisiae / C. glabrata (4,384 pairs, most of them saturated) 815 of 856, median 5.9 %, p95 39 %, max 166 %, 317 status changes (311 `"ok"` -> `"saturated"`, 6 back); human / mouse (measured with COAL's gestimator flavour, i.e. this code with the saturation rule of 0.4.3) 813 of 946, median 1.0 %, p95 5.0 %, 4 status changes. The corrected dS is slightly closer to YN's (median difference 7.68 -> 7.53 % on the 1,000-gene example) | `options(orthologr.comeron_twosubs = "libsequence")` |
| Comeron's method: a site class with fewer than 3 sites (the 2-fold transversional class of short alignments) cannot make dN or dS `NA`; its undefined logarithm is 0, as in 0.4.3 | simulated 30-codon pairs: dN `NA` 3.7 % -> 0 %, the recovered values' mean absolute error 8.8 % (all 30-codon pairs 11.2 %); real data: 1,000-gene example 0 pairs, human / mouse 3, yeasts 14. Never changes a value that is not `NA` otherwise. (0.6.0: only where the value is > 0; see 0.6.0, "Changed results") | `options(orthologr.comeron_min_sites = 0)` |
| `orthogroups()` and `dNdS_orthogroups()` trim coding sequences whose length is not a multiple of 3 (`delete_corrupt_cds = FALSE`, as `dNdS()` does) | - (new functions) | `delete_corrupt_cds = TRUE` |

All of these except the T-Coffee mode, `aa_aln_params`, `divergence_stratigraphy(quiet)`, the minimum class size of
Comeron (which does not matter with the 0.4.3 saturation rule) and `delete_corrupt_cds` at once:
`orthologr_legacy_options()`.

Further behaviour that differs from 0.4.3 and cannot be switched back:

- A gene pair whose alignment or estimation fails no longer ends the run: it is a row of `NA` with status
  `"failed"` and its reason in `attr(, "failed_pairs")`.
- A missing program (perl, the alignment program, KaKs_Calculator) and unsupported argument values are errors
  before the search starts.
- Parameter strings (`aa_aln_params`, `aligner_params`, `diamond_params`, `blast_params`, `params`,
  `kaks_calc.params`, ...) are lists of arguments, not shell text: `;`, pipes and redirections inside them no
  longer work.
- `dNdS()` writes nothing to `stdout`. With `quiet = TRUE` (default) it gives a line at the start and the
  summary at the end (31 messages before); the citation note is a message shown once per session. Input sequences
  that are left out or trimmed to a multiple of 3 are counted in the summary (listed in
  `attr(, "dropped_sequences")` / `attr(, "trimmed_sequences")`) instead of a block of messages about each kind:
  1,000-gene example, 12 messages before (8 about the 12 trimmed sequences), 2 after (also with
  `engine = "coal"`); `quiet = FALSE` shows the old messages.
- Coding sequences with characters other than A, C, G, T and the IUPAC ambiguity codes (or, with
  `options(orthologr.ambiguous_cds = "drop")`, with any character other than A, C, G, T) are left out before the
  search, as before, but no longer without a word: they are counted in the summary and listed in
  `attr(, "dropped_sequences")`.
- Intermediate files are kept in a folder of the call inside `tempdir()` and removed afterwards
  (`options(orthologr.keep_tmp = TRUE)` keeps them until the end of the R session, which removes `tempdir()`;
  `options(orthologr.keep_tmp = "<folder>")` keeps them in that folder, also after `Rscript` or a batch job, and
  the message says which applies); no function changes the working directory or writes into it unasked
  (`read.cds()` wrote `<file>_corrupted_cds_seqs.fasta` there; the summary says when that file is in `tempdir()`).
- `ggplot2`, `scales`, `ape`, `rtracklayer`, `benchmarkme`, `tidyr`, `biomartr`, `metablastr` moved from Imports
  to Suggests (Imports: 21 -> 11 packages); R >= 4.4 is enough.
- `cdd.path` and `ortho_path` of `dNdS()` are deprecated (never used); `ortho_detection = "Orthofinder2"`, `"PO"`
  and `"DELTA"` are errors with an explanation.

### Fixed

Results:

- `dNdS()` mixed up the two species when query and subject file had the same file name in different folders.
- Sequence IDs are returned as they are in the input files. IDs with a hyphen, with `|`, numeric-looking IDs
  (`007`, `1e5`), IDs of more than 15 characters (MAFFT) and the same ID in both files lost rows, came back
  changed or gave `NA`. The dN/dS rows are joined with the hit table by gene pair (was: by query only).
- Li's method returned text columns and took small values (`9.95e-05`) for the saturation marker.
- `filter_dNdS()` and `divergence_stratigraphy()` say how many gene pairs they leave out and why
  (`divergence_left_out()`); gene pairs with dS = 0 are left out as before unless `dS_zero = "keep"`.
- `divergence_map()`: strata are tie-safe, `n_quantile = 49` or `98` no longer gives `NA` for the largest value,
  a failed binning no longer loses the finished dN/dS run, and the caller's table is not changed.
- Comeron's method (gestimator): input that is not a codon alignment is an R error (a sequence shorter than one
  codon ended the R session; sequences of unequal length gave numbers).
- Best hits no longer depend on the order of the hit table.
- `retrieve_longest_isoforms()` and `select_orthologs()` work on real NCBI and Ensembl GFF3 / GTF files
  (proteins were lost without notice, annotation lines of more than 65,000 characters were deleted, GFF3 loci
  were grouped by name).
- `promotor_divergence_of_orthologous_genes()` pairs promoters by ortholog pair.
- `substitutionrate(kaks_calc.params = )` is added to the method instead of replacing it.
- `dnds_est.method = "MLPB"` (default engine, KaKs_Calculator) works: `substitutionrate()` returned `NULL` for it
  and `dNdS()` ended after the whole search with an error about a missing column (20-gene example: before an
  error after 4 s; after 20 pairs, the values KaKs_Calculator gives for `-m MLPB`). There is one list of the
  KaKs_Calculator methods for the whole package, and a method that `substitutionrate()` does not handle is an
  error (`orthologr_error_bad_argument`) instead of `NULL`.

Running programs:

- Stale results can no longer be returned: `blast()`, `diamond()`, `deepclust()` and the functions built on them
  check the exit status and the output of every program and work in files of their own.
- File and folder names with spaces work. Programs are looked up once per session.
- Parallel runs use the `parallel` package only; on Windows (PSOCK) arguments given as variables reach the
  worker processes, and worker processes are ended after an interrupt.
- With `R CMD INSTALL --no-byte-compile`, `devtools::load_all()` or `covr` every gene pair of a default
  `dNdS()` run was `NA` in development versions of 0.5.0 (the native routine `gestimator` was registered with five
  instead of six arguments). The registration is now generated by `Rcpp::compileAttributes()`, and a test
  compares it with the R wrappers.
- `multi_aln()`: T-Coffee and MUSCLE examples run with current versions; CLUSTAL files of T-Coffee are read;
  the MAFFT check uses `mafft --version`.
- `multi_aln(tool = "t_coffee")` and `dNdS(aa_aln_tool = "t_coffee")` returned unaligned rows in development
  versions of 0.5.0 (a bug, no option restores it): T-Coffee 12 was given the sequences with `-infile` and took a
  file whose sequences end in `*` (as the translated coding sequences do) for a finished alignment. They are now
  given with `-seq`. 20-gene example, `aa_aln_type = "multiple"`: before 10 of 20 pairs `"saturated"` (dN/dS
  `NA`), 7 of 20 within 0.1 of ClustalW, r = 0.743 (max difference 0.553); after 0 saturated, 20 of 20 within 0.1
  of ClustalW, r = 0.993 (max difference 0.094). 0.4.3 used `-mode accurate`, which fails without network.
- `pairwise_aln()`, `read.genome(format = "fastq")`, `codon_aln(get_aln = TRUE)`,
  `promotor_divergence_estimation()`, `lnc_map_counts()`, `map_generator_dnds()` with one subject file: errors
  on valid input fixed. `check_annotation(remove_annotation_outliers = TRUE)` no longer overwrites the
  annotation file.
- An unknown `aa_sub_matrix` name (`dNdS()`, `pairwise_aln()`) is an error (`orthologr_error_bad_argument`)
  before the search; it ran the whole search and then failed every gene pair (1,000-gene example,
  `aa_sub_matrix = "FOO"`: before 6.4 s and 995 of 995 pairs `NA` with one warning; after an error at once).
  With `engine = "coal"` the name is checked against the matrices of COAL (BLOSUM62, BLOSUM45, BLOSUM80, PAM250).
- A protein file given where coding sequences are expected is an error that says so (was: a DIAMOND error about
  memory).
- `orthogroups()` runs the check for several transcripts per gene; the check no longer warns for a bacterial
  genome with a handful of genes that have two coding sequences.
- FASTA files with a UTF-8 byte-order mark were rejected ("could not be read properly"), and files with old-Mac
  line ends (CR only) reached DIAMOND as one sequence and failed after about 10 s with a memory error. Both are now
  read like the plain file (same result as for the file without the mark / with LF line ends; also gzipped), as
  COAL reads them.
- `dnds_est.method = "kaks_calc"` in `dNdS()` / `dNdS_orthogroups()` is an error before the search (it ran the
  search and failed every gene pair: `kaks_calc.params` exists only in `substitutionrate()`); `"GY"`, `"GMYN"`,
  `"MS"`, `"MA"` and `"ALL"` are documented in `?dNdS` (MA and MS take about 6.5 s per gene pair); with `"ALL"`
  the summary counts gene pairs x methods (20-gene example: "340 rows (20 gene pairs x 17 methods)", was "340 gene
  pairs").
- `dNdS_orthogroups(engine = "coal")` no longer warns that saturated pairs "could not be processed" (4 yeasts, 60
  groups: a warning for 322 of 360 pairs, none of which failed): the warning (`orthologr_warning_failed_pairs`)
  counts pairs with status `"failed"` only; saturated pairs are counted in the summary message, as in `dNdS()`.
- `dNdS_orthogroups()` summary: `n_saturated` counted pairs with dS above `dS_saturation`, while the status
  `"saturated"` means that dN or dS could not be estimated (60 yeast groups: 296 pairs with status "saturated",
  `sum(n_saturated)` 30). `n_saturated` now counts status `"saturated"`, as in `dNdS_clusters()`; the new columns
  `n_dS_zero`, `n_no_data` and `n_failed` count the other statuses, and `n_dS_above` is the old count.
- `dNdS_orthogroups()` summary: the quartile columns are `dN_q25`, `dN_q75`, `dS_q25`, `dS_q75`, `dNdS_q25`,
  `dNdS_q75`, as in `dNdS_clusters()`, coalr and COAL's command line (they were `_q1` / `_q3` in development
  versions of 0.5.0; `dNdS_orthogroups()` is new in 0.5.0, so no release had the old names).
- With COAL (`engine = "coal"`, or `codon_aln(tool = "coal")`), a codon alignment that cannot be built is an error
  of class `orthologr_error_alignment_failed` (with `_tool_failed`), as with PAL2NAL; it was a plain error.
- `codon_aln(tool = "pal2nal")` with selenocysteine (`U`) in a protein: PAL2NAL v14 has no codon for it and failed
  ("inconsistency between pep and nuc"), while COAL accepts any codon there. PAL2NAL now gets the protein alignment
  with `U` as `X` (which it matches with any codon): both tools give the same codon alignment (byte-identical FASTA;
  checked with the input as FASTA and CLUSTAL). Inside `dNdS()` the proteins are orthologr's own translations and
  never hold `U`.
- `dNdS_clusters()` summary: a cluster none of whose pairs has an estimate because they are saturated (or have no
  codons to compare) has the status `"saturated"` (`"no_data"`) and a note with the numbers; it was `"failed"` with
  an empty note and `n_failed = 0` (E. coli K-12 + O157: 37 clusters).
- `read_dNdS()` returns IDs with a double quote as written (`d"q` came back as `d""q`), and non-ASCII text with
  the bytes and the encoding marks it had when it was written (kept in `<file>.run.json`), so that IDs read back
  equal those of `dNdS()` in every locale. (An intermediate version marked all text as UTF-8: in a C / POSIX
  locale a non-ASCII ID then no longer matched the ID of the run - `==`, `match()` and `merge()` failed; and
  `write_dNdS()` wrote such IDs of `attr(, "failed_pairs")` as escapes like `<c3><bc>` in a C locale.)
- An invalid `options(orthologr.max_aln_cells)` (`-1`, `"5e8"` as text, `NA`, several values) is an error
  (`orthologr_error_bad_argument`, naming the value) before the search in `dNdS()` and `dNdS_orthogroups()`; it ran
  the search and then failed every gene pair with the message about the option (20-gene example: 20 of 20 pairs
  "failed"). `options(orthologr.keep_tmp)` of `NA` or of several values is an error as well (was taken for
  `FALSE`), and a folder that cannot be created is `orthologr_error_bad_argument` (was an unclassed error). The
  reason of a pair above `orthologr.max_aln_cells` mentions `engine = "coal"`; a run of one gene pair is summarised
  as "1 gene pair".
- A DIAMOND older than 2.1.0 is used with a warning (`orthologr_warning_tool_version`; it was used without a word);
  a program that is there but crashes when asked for its version is `orthologr_error_tool_failed` ("could not be
  started ... exit status 139"), no longer "not installed"; a failed DIAMOND search is `orthologr_error_tool_failed`
  and `comp_cores` above the number of cores `orthologr_warning_comp_cores` (both were unclassed). The column `qcov`
  of DIAMOND hit tables equals `qcovhsp` (DIAMOND has no per-subject coverage; documented in `?diamond`, kept for
  the column layout of BLAST tables and of 0.4.3).

### Speed (results unchanged)

- Default `dNdS()` on the 1,000-gene example, 1 core: 744 s (0.4.3) -> 50 s; with 2 cores 27 s. Half of the last
  step comes from building the full-length alignment rows directly instead of with `pwalign::alignedPattern()`
  (97 s -> 50 s); the rest from BLOSUM62 scoring, reading and translating every file once, vectorised
  translation and best-hit selection, and not writing a second copy of every alignment.
- `select_orthologs()`, `retrieve_longest_isoforms()` and `deepclust()` read annotations and sequences in one
  pass.
- When the R package `coalr` is installed, `dNdS()` and `divergence_stratigraphy()` build the codon alignments
  with COAL instead of starting PAL2NAL (Perl) once per gene pair (`codon_aln_tool = NULL`, the new default;
  `options(orthologr.codon_aln_tool = "pal2nal")` or `codon_aln_tool = "pal2nal"` keeps PAL2NAL; the record of the
  run names the tool used). The codon alignments are byte-identical: 12,605 of 12,605 gene pairs (1,000-gene
  example 995, E. coli K-12 / O157 3,663, Salmonella LT2 / Shigella flexneri 2,820, S. cerevisiae / C. glabrata
  4,164, human 1,000 transcripts / all mouse CDS 945, a file of edited sequences 18), the 20-gene example with
  each of Clustal Omega, MAFFT, MUSCLE, ClustalW and T-Coffee, and the 1,000-gene example with each of Comeron,
  Li, NG, LWL, MLWL, LPB, MLPB, YN and MYN give identical tables; with and without `coalr` installed the
  default run gives the same table (checked on two of these sets), also under `orthologr_legacy_options()`.
  `codon_aln_tool = "coal"` where COAL is not installed uses PAL2NAL, with a message (so the arguments of a run
  record from a machine with `coalr`, which name `"coal"`, can be used again elsewhere; this was an
  `orthologr_error_package_missing`). `dNdS_orthogroups(engine = "orthologr")` builds its codon alignments with
  PAL2NAL in any case (same alignments).
  Default run, 2 cores (Comeron): 37 -> 17 s (1,000-gene example), 99 -> 69 s (E. coli), 79 -> 56 s
  (Salmonella / Shigella), 156 -> 84 s (yeasts), 61 -> 41 s (human / mouse); the deep test measured 178 -> 94 s
  for the yeasts with YN on 1 core.

### Checked against other programs

`dNdS()` (default engine: gestimator for Comeron, KaKs_Calculator 2.0 for YN) against the COAL engine
(`engine = "coal"`), on 300 sampled reciprocal-best-hit pairs of each of three species pairs (NCBI RefSeq coding
sequences; 2026-11-13). Columns: gene pairs in both results, pairs with dS > 0 from both, share of these whose dS
agrees within 1 %, median relative difference of dS, Pearson r of dS and of dN/dS, pairs with dS from one engine
only.

| Species pair | Method (default engine) | pairs | dS both | within 1 % | median rel. diff. | r dS | r dN/dS | one only |
|---|---|---|---|---|---|---|---|---|
| S. cerevisiae / C. glabrata | Comeron (gestimator) | 283 | 73 | 63.0 % | 0.0018 | 0.952 | 0.989 | 3 |
| S. cerevisiae / C. glabrata | YN (KaKs_Calculator) | 283 | 38 | 65.8 % | 2e-06 | 0.883 | 0.965 | 245 |
| S. cerevisiae / K. lactis | Comeron (gestimator) | 265 | 68 | 69.1 % | 0.0023 | 0.977 | 0.999 | 1 |
| S. cerevisiae / K. lactis | YN (KaKs_Calculator) | 265 | 31 | 71.0 % | 2.3e-06 | 0.965 | 0.999 | 234 |
| Salmonella LT2 / Shigella flexneri | Comeron (gestimator) | 294 | 274 | 95.3 % | 0 | 0.999 | 0.998 | 0 |
| Salmonella LT2 / Shigella flexneri | YN (KaKs_Calculator) | 294 | 278 | 92.4 % | 1.2e-06 | 0.972 | 0.987 | 16 |

For the close pair the engines agree (median difference 0, r >= 0.97). The differences in this table come from
the protein alignments: where several alignments score the same, the COAL of that day placed a gap otherwise
than pwalign. **With a `coalr` that has COAL's pwalign tie-break** (`tie_break`; `coalr` from COAL main since 009bdd5,
pull request #64 of 2026-11-14), `engine = "coal"` computes pwalign's protein alignments row for row and the engines agree
to the digits orthologr prints. All reciprocal-best-hit pairs, native engine against `engine = "coal"` (coalr from
COAL main 009bdd5), before / after (dS where both engines give one, equal to 6 significant digits, 4 for Comeron; pairs with the same status):

| Species pair | Method | pairs | dS equal | dS > 1 % apart | same status |
|---|---|---|---|---|---|
| A. thaliana / A. lyrata | YN | 995 | 677 / 994 -> 994 / 994 | 267 -> 0 | 995 -> 995 |
| A. thaliana / A. lyrata | Comeron | 995 | 692 / 994 -> 994 / 994 | 258 -> 0 | 995 -> 995 |
| A. thaliana / A. lyrata | LWL | 995 | 687 / 994 -> 994 / 994 | 245 -> 0 | 995 -> 995 |
| S. cerevisiae / C. glabrata | YN | 4,164 | 322 / 565 -> 587 / 587 | 206 -> 0 | 4,107 -> 4,146 |
| S. cerevisiae / C. glabrata | Comeron | 4,164 | 563 / 1,120 -> 1,155 / 1,155 | 363 -> 0 | 4,097 -> 4,164 |
| S. cerevisiae / C. glabrata | LWL | 4,164 | 821 / 2,166 -> 2,239 / 2,239 | 873 -> 0 | 4,021 -> 4,164 |

Re-measured with the defaults of the decisions of 2026-11-14 (`sensitivity_mode = "default"`, `ambiguous_cds =
"keep"`, Comeron's two-position pathways and minimum class size in both engines; `coalr` from COAL main 956a458,
pull request #65): A. thaliana / A. lyrata 999 pairs, YN, Comeron and LWL: dS equal 998 / 998, dN equal 999 / 999
(Comeron 998 / 998), same status 999 / 999; S. cerevisiae / C. glabrata 4,384 pairs: Comeron dS 862 / 862, dN
1,246 / 1,246, same status 4,384; LWL dS 2,301 / 2,301, same status 4,384; YN dS 619 / 619, dN 4,351 / 4,352, same
status 4,361 (the 23 others are YN iterations that do not converge, as below); no dS more than 1 % apart.

dN agrees in the same way (one yeast YN pair: dN 2.5e-12 from COAL, 0 from KaKs_Calculator, which prints values
below 1e-6 as 0). The 18 yeast YN pairs with another status are iterations that do not converge: COAL gives no
estimate (`"saturated"`), KaKs_Calculator its last round (`"ok"`). The yeast pairs are saturated at synonymous sites (median dS above 2): there COAL returns `NA` (status "saturated")
for most YN estimates for which KaKs_Calculator still prints a number (dS 4 - 5); where both give a value they
agree. (This table was measured before the default engine returned `NA` for KaKs_Calculator's capped values as
well - see "Changed defaults". Since then, on the same codon alignments, the status of `kaks_tool =
"kaks_calculator"` and of `kaks_tool = "coal"` agrees for 150 of 150 E. coli K-12 / O157 pairs with each of NG,
LWL, MLWL, LPB, MLPB, YN and MYN (before 140 - 141), and for S. cerevisiae / C. glabrata in 150 (YN, before 20),
148 (MYN, before 10; 2 pairs for which KaKs_Calculator prints `nan`), 149 (LWL / LPB, before 105; MLWL / MLPB,
before 101) and 150 (NG) of 150 - on these 150-pair samples. On all 4,164 S. cerevisiae / C. glabrata pairs, YN:
4,138 (18 YN iterations that did not converge, which KaKs_Calculator returns as estimates with status `"ok"`
because it does not report convergence, and 8 saturated pairs that COAL calls not converged); MYN: the `nan`
pairs are `"failed"`, 7 non-converged pairs are `"ok"`; the LWL family: all pairs. `kaks_tool = "coal"` /
`engine = "coal"` report the iterations that do not converge (`?substitutionrate`, section "Saturation
(KaKs_Calculator methods)").) DIAMOND `sensitivity_mode`: "fast" (default) finds 5.0 % (S. cerevisiae /
C. glabrata), 8.6 % (S. cerevisiae / K. lactis) and 1.5 % (Salmonella / Shigella) fewer reciprocal best hits than
"default", at about the same run time on these genomes; the default is unchanged, and `dNdS()` suggests
"default" when the median identity of the best hits is below 70 %.

### Internal

- External programs are started through one helper (`run_tool()`); `compute_dnds()` has one code path for serial
  and parallel runs.
- Tests: 16 of the 20 empty test files now hold tests; pinned dN/dS values of fixed codon alignments, the 20
  known ortholog pairs end to end, a tie-break fixture for both best-hit rules, IDs end to end, stored 0.4.3
  tables for the legacy settings, malformed input files (empty, truncated record, lower case, protein instead of
  CDS, duplicated IDs, Windows line endings, gzip, an N in a sequence). The suite passes without any optional
  program.
- GitHub Actions: R CMD check on Linux, macOS and Windows with the external programs installed, one leg without
  them; triggers for pull requests, the default branch and a manual start. Checked with `--as-cran`.
- Examples that need no external program are run by `R CMD check`; the others are guarded by the programs they
  need.
- `.travis.yml`, `orthologr.yml`, `Read-and-delete-me` and the unused `inst/KaKs_Calc_parser` were removed.

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
