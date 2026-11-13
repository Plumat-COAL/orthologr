## `orthologr` 0.5.0

Version 0.5.0 instead of 0.4.4 because results and behaviour change in ways a patch release should not bring:
some default results change (see "Changes to pairwise alignments" where that change is merged, the T-Coffee
default mode, `aa_aln_params` being passed on), functions that returned `NULL` or stale results now stop with an
error, and packages moved from Imports to Suggests.

**Default change that alters which orthologs are returned: best hits are ranked by bit score.** See the next
section; `options(orthologr.best_hit_rule = "evalue")` restores the rule of 0.4.3.

### New: orthogroups of several species with DIAMOND - `orthogroups()`, `write_orthogroups()`, `dNdS_orthogroups()`

orthologr no longer needs OrthoFinder for groups of orthologs of more than two species (`orthofinder2()` never
worked beyond starting the program and is soft-deprecated; its warning and `?orthologs` now point to
`orthogroups()`). Nothing in this section changes the result of an existing function.

- **`orthogroups(files, method = c("rbh", "deepclust"), ...)`** for 2 to N species (a named vector or list of
  FASTA files, or a folder; coding sequences or proteins).
  - `method = "rbh"` (default; about 2 to 15 species): **one** all-against-all DIAMOND search of the proteins of all
    species, reciprocal best hits for every pair of species - chosen by the code `diamond_rec()` uses, so
    `best_hit_rule` applies -, and the connected components of the graph of reciprocal best hits as groups.
    E-values are converted to those of a search against the subject species alone (DIAMOND's e-value is
    proportional to the database size), so that the result for two species does not depend on which other species
    are in the run. `inparalogs = TRUE` attaches a gene without reciprocal best hit to a gene of its species when
    the two are more similar to each other than either is to any gene of another species (heuristic, off by
    default). Optional hit filters `min_qcov`, `min_identity` (off by default).
  - `method = "deepclust"` (many species, large data): `deepclust()` on the proteins of all species; clusters with
    genes of at least two species are the groups. `refine = TRUE` splits the clusters with more than one gene of a
    species by the reciprocal best hits inside the cluster (one further search, of the members of these clusters
    only); `realign = TRUE` adds `deepclust_realign()` statistics.
  - The result (class `"orthogroups"`) holds `groups` (long), `wide` (one row per group, one column per species,
    as OrthoFinder's `Orthogroups.tsv`), `pairs` (every pair of genes of different species in a group with bit
    score, e-value, identity, alignment length, coverage, `rbh`, and `relation`: `"1:1"` or `"co-ortholog"`),
    `single_copy`, `unassigned`, `genes`, `species`, `params`; `print()`, `summary()`, `as.data.frame()`.
  - The result does not depend on `comp_cores` or on the order of the files. Sequence IDs of any form work (they
    are replaced by numbers during the search); IDs that occur in several species are prefixed with the species
    name, and a message says so.
- **`write_orthogroups(x, dir)`** writes `Orthogroups/Orthogroups.tsv`, `Orthogroups_SingleCopyOrthologues.txt`,
  `Orthogroups_UnassignedGenes.tsv`, `Orthogroups.GeneCount.tsv` and `pairs.tsv` in OrthoFinder's layout;
  `orthofinder2_retrieve_core_orthologs()` reads them. `ortholog_tables = TRUE` adds pairwise tables for
  `import_ortholog_tables_all()` / `retrieve_core_orthologs()`.
- **`dNdS_orthogroups(x, cds_files, groups = c("single_copy", "all"), ...)`**: for every group one protein
  alignment (`aa_aln_tool = "clustalo"`, `"mafft"`, `"muscle"`, `"clustalw"`, `"t_coffee"`; or `"pairwise"`), the
  codon alignment (PAL2NAL), and dN, dS, dN/dS for every pair of genes of different species
  (`substitutionrate()`; default method Comeron), plus a summary per group (numbers of pairs, of pairs without
  value and of saturated pairs; median, mean and quartiles of dN, dS, dN/dS). A gene, pair or group that cannot be
  used gives `NA` rows and a reason, never an error. With `aa_aln_tool = "pairwise"` the values are those of
  `dNdS()` for the same pairs (tested).

Measured (DIAMOND 2.2.8, 2 cores; the data are not part of the package):

VALIDATION_NEWS

### Best hits are ranked by bit score (default change)

Decided by the maintainer of this fork for 0.5.0. `diamond_best()`, `blast_best()`, `diamond_rec()`, `blast_rec()`,
`orthologs()`, `dNdS()`, `divergence_stratigraphy()` (and the choice of the splice variant per locus in
`select_orthologs()`) share one piece of code that chooses the best hit of a query, and that code has a new
argument, **`best_hit_rule`**, at the end of each signature:

- `"bitscore"` (**new default**): the hit with the highest bit score; of several, the one with the smallest
  e-value; then the longest alignment; then the highest percent identity; then the alphabetically first
  `subject_id`. The choice is deterministic and does not depend on the order of the hit table.
- `"evalue"`: the rule of orthologr <= 0.4.3 - smallest e-value, then longest alignment, then the first row of
  the hit table (sorted by `query_id`, `subject_id`). The bit score was never used. This rule reproduces earlier
  results exactly: the test suite compares it with the per-query reference implementation (`best_hits_do()`) on
  the 1,000-gene example, for best hits and reciprocal best hits.
- `options(orthologr.best_hit_rule = "evalue")` sets the old rule for all functions of a session.

Why: DIAMOND and BLAST report an e-value of 0 for every strong hit of a long, conserved protein. Under the old
rule such ties were decided by alignment length, so a paralog with a longer but worse alignment could be taken
for the best hit, and a few by alphabetical order.

Measured effect (the default of every function named above changes accordingly):

| | old rule (`"evalue"`) | new default (`"bitscore"`) |
|---|---|---|
| human vs mouse (19,746 x 21,584 genes, DIAMOND): reciprocal best hits | 16,458 | 16,827 |
| ... pairing a human gene with another mouse gene than any Ensembl Compara ortholog | 330 (2.0%) | 68 (0.4%) |
| ... Compara one-to-one orthologs found (of 16,312 present in both inputs) | 93.0% | 96.8% |
| ... pairs in common / only under this rule | 16,135 / 323 | 16,135 / 692 |
| human -> mouse best hits that change | | 596 of 18,676 (3.2%) |
| 1,000-gene example (`ortho_thal_cds_1000.fasta` vs `ortho_lyra_cds_1000.fasta`): best hits that change | | 12 of 986 (1.2%) |
| ... reciprocal best hits (`diamond_rec()`) | 972 | 983 (970 in common, 2 only old, 13 only new) |
| 20-gene example | 20 | 20 (the same pairs) |

The human / mouse numbers were obtained with the best-hit code of this version on the DIAMOND hit tables of the
second audit; the old rule gives the audit's (and 0.4.3's) best hits exactly, and its numbers are those the audit
published (16,458; 330; 93.0%). The audit's own bit-score ranking (bit score, e-value, row order) gave 16,828
pairs; the documented tie-breaks of this version give 16,827, with the same 68 and 96.8%.

Also
- The column `qcov` of DIAMOND hit tables always equals `qcovhsp`. This is not a parsing error: BLAST+ has a
  per-subject query coverage (`qcovs`), DIAMOND has not, and `diamond()` asks for `qcovhsp` twice to keep the 21
  columns of both kinds of table the same. Documented in `?diamond`; the table is unchanged.
- `blast_best()` / `blast_rec()` add `-best_hit_score_edge 0.05 -best_hit_overhang 0.25` to every BLAST call, as
  before (BLAST's own culling of hits; N38 of the audit): the BLAST and DIAMOND paths therefore start from
  slightly different hit lists before the same rule is applied. Not changed here.

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

### Annotations: longest isoforms and gene loci on real GFF3 / GTF files (second audit)

`retrieve_longest_isoforms()`, `select_orthologs()` and `generate_ortholog_tables()` now read which gene locus a
transcript or protein belongs to from **one transcript -> gene table** that is built the same way for every
annotation source: GFF3 through the `Parent` / `ID` relations (CDS -> mRNA -> gene, or CDS -> gene in NCBI
prokaryote files; the locus is the gene `ID` without `gene-` / `gene:`, i.e. the `locus_tag` for NCBI), GTF through
`gene_id` / `transcript_id` / `protein_id`. IDs of the FASTA file or dN/dS table are matched as they are, then
without version suffix, and NCBI `lcl|..._cds_<protein>_<n>` IDs by their protein accession. Tested on pieces of
real annotations that ship with the tests (`tests/testthat/fixtures`): NCBI RefSeq GFF3 + `protein.faa`, NCBI GTF,
Ensembl GFF3 and Ensembl GTF with `pep` and `cds` FASTA.

- **`retrieve_longest_isoforms(annotation_format = "gff")` lost proteins without saying so.** CDS rows were
  grouped by the `gene` attribute, the gene *symbol*: all loci without a symbol were merged into one, and loci
  sharing a symbol into one. *E. coli* O157:H7 Sakai (NCBI RefSeq): 5,155 proteins of 5,155 loci in, **3,816 out**
  (1,283 unnamed loci -> 1 protein, 57 more lost to shared symbols). Now 5,155 out, one per locus. Loci are never
  grouped by symbol.
- `retrieve_longest_isoforms()` only worked for NCBI GFF3 + `protein.faa` and for Ensembl GTF when the FASTA
  headers were unversioned transcript IDs. Measured on the files of the audit (before -> after):

  | annotation + sequences | before | after |
  |---|---|---|
  | NCBI *E. coli* K-12 GFF3 + protein.faa (4,300 proteins, 4,290 loci) | 4,290 | 4,290 |
  | NCBI *E. coli* O157 Sakai GFF3 + protein.faa (5,155 loci) | 3,816 | 5,155 |
  | NCBI *E. coli* K-12 GTF + protein.faa | error | 4,290 |
  | NCBI GFF3 `.gz` + protein `.gz` | 4,290 | 4,290 |
  | Ensembl yeast GTF + pep (6,600) | error | 6,600 |
  | Ensembl yeast GTF + cds | 6,600 | 6,600 |
  | Ensembl yeast GFF3 + pep | error | 6,600 |
  | Ensembl human chr21 GTF + pep (2,212 proteins, 221 protein-coding genes) | error | 221 (all the longest of their gene) |
  | Ensembl human chr21 GFF3 + pep (222 genes with a protein) | error | 222 |

  With `annotation_format = "gtf"` only genes with `gene_biotype "protein_coding"` are kept, as before (a message
  now says how many sequences of which other biotypes are left out); with `"gff"` every locus with a sequence, as
  before.
- **`select_orthologs(format = "gff")` could not work on real GFF3 files** (rows were grouped by `Name`, which a
  gene and its mRNA never share; NCBI: "Something went wrong during the feature extraction", Ensembl: "missing
  value where TRUE/FALSE needed"). **`select_orthologs(format = "gtf")`** matched IDs to `transcript_id` literally
  and failed for the versioned IDs of Ensembl FASTA files (`ENST00000373020.9`). Both work now, also for protein
  IDs and for raw GFF3 IDs (`rna-NM_...`, `transcript:ENST...`).
- When no ID fits, the error shows example IDs of the sequence file / dN/dS table and of the annotation
  (was: "Something went wrong during the feature extraction ..."). Sequences / IDs that are not in the annotation
  are counted in a message.
- **Memory and time.** The annotation is read as text in blocks and only feature type and the needed attributes
  are kept; `check_annotation()` no longer holds the whole file in memory to look for over-long lines.
  `select_orthologs()` on a human / mouse GTF pair cut down to 436k + 123k lines, 1,500 gene pairs:
  peak RSS 1,004 MB -> 369 MB, 16 s -> 10 s (R with the packages loaded takes about 250 MB of that).
  `rtracklayer` is no longer needed by these functions (it stays in Suggests for the reference implementations
  that the tests compare with).
- `retrieve_longest_isoforms()` returns the path of the new file invisibly (was `NULL`), with the number of loci
  and sequences as attributes, and stops on duplicated sequence IDs.

Changes of results in this section
- `retrieve_longest_isoforms(annotation_format = "gff")`: more sequences for annotations with unnamed loci or
  shared gene symbols (Sakai: 3,816 -> 5,155). No option restores the old grouping (it was a loss of data); the old
  code is kept as the internal `retrieve_longest_isoforms_do()` for the tests.
- `select_orthologs(format = "gff")`: `query_gene_locus_id` / `subject_gene_locus_id` hold the gene ID without
  the prefix `gene-` / `gene:` (e.g. `AT1G01010`, was `gene-AT1G01010` in the synthetic files for which the old
  code worked), so that GFF3 and GTF of one annotation give the same loci.
- `select_orthologs(format = "gtf")`: transcripts of loci without a `gene` row, and IDs that only match without
  version, are now found (rows that were left out before).

### Divergence strata, dS = 0, saturation and input checks of Comeron's method (second audit)

`divergence_map()` / `divergence_stratigraphy()`
- **A finished run is no longer lost.** `divergence_map()` stopped when two quantile thresholds coincided, and
  `divergence_stratigraphy()` called it only after the complete dN/dS run and then returned nothing: *E. coli* K-12
  (1,500 genes) against O157 Sakai ran for 40 minutes and failed, because 28% of the gene pairs have dN/dS = 0.
  Now (1) strata are assigned by a tie-safe rule (below), (2) `attr(result, "dNdS_tbl")` of
  `divergence_stratigraphy()` holds the complete dN/dS table of the run, before filtering and binning, (3) should
  the binning fail all the same, the dN/dS table is returned with a warning, and (4) `n_quantile` and `dS_zero`
  are checked before the run starts.
- **Tie-safe strata.** Stratum `k` holds the values in `(q[k-1], q[k]]` of the `n_quantile` quantiles, as before;
  with coinciding thresholds the same rule is applied through ranks
  (`max(1, ceiling((r - 1) * n / (N - 1)))` with `r` the smallest rank of a value). Without ties the strata are
  exactly those of earlier versions (tested for `n_quantile` 2-30, 49, 50, 98, 100). Equal values always share a
  stratum, the lowest one their ranks reach; the strata such a block spans beyond that stay empty, so fewer than
  `n_quantile` strata are used and they differ in size - a warning names the empty strata. See `?divergence_map`.
- `n_quantile = 49` and `98` gave `NA` for the largest value (`seq(0, 1, 1/n)` does not end at 1).
  `n_quantile` must be a whole number >= 1.
- `divergence_map()` no longer changes the table of the caller (a tibble came back as a re-sorted `data.table`).

dS = 0
- **`filter_dNdS()` (and so `divergence_stratigraphy()`) dropped every gene pair with dS = 0 without a word**,
  because their `dNdS` is `NA` and `NA <= dnds.threshold` is `NA`: genes with identical coding sequences
  (dN = dS = 0) and genes with nonsynonymous changes only (dN > 0, dS = 0). In the *E. coli* run 85 of 1,280 pairs
  disappeared, 47 of them identical sequences. New argument `dS_zero`:
  `"keep"` (**the new default**) keeps these pairs with `dNdS = NA` and marks them in a new logical column
  `dS_zero`; in a divergence map they have the stratum `NA`. `"drop"` leaves them out as before. A message gives
  the numbers (identical / dN > 0) in both cases. `options(orthologr.dS_zero = "drop")` restores the old default.

Comeron's method (`dnds_est.method = "Comeron"`, the default)
- **Saturation is `NA`, not a finite number.** When the argument of a logarithm of the correction is not positive
  (more differences than the model can correct for) the term was set to 0 and a finite, too small dS or dN was
  returned - `ATGGCTAAA` / `ATGGCCAAA` (one synonymous substitution) gave dS = 0. Such a dS / dN, and dN/dS, is now
  `NA`, as the documentation always said ("NA when it cannot be computed"). Values of all other gene pairs are
  unchanged. `options(orthologr.legacy_comeron_saturation = TRUE)` reproduces the old numbers. How many pairs
  change (default pipeline of this branch, same codon alignments in both modes):

  | gene pairs (reciprocal best hits) | dS: value -> `NA` | dN: value -> `NA` | dN/dS: value -> `NA` | other values |
  |---|---|---|---|---|
  | *A. thaliana* / *A. lyrata*, 1,000-gene example of the package: 984 | 0 | 5 | 2 (981 -> 979 with a value) | identical |
  | human / mouse, the 1,500 pairs of the audit (gene families, long genes, random pairs): 1,500 | 2 (old dS 0.71 and 1.23) | 23 | 15 (1,490 -> 1,475 with a value) | identical |

  (Most of the pairs concerned are flagged on the nonsynonymous side: the term for twofold-transversional sites,
  of which a gene has few. The general underestimation of large dS by Comeron's method that the audit describes -
  median 0.87 against 1.96 from PAML yn00 for pairs with yn00 dS > 1.5 - is a property of the method, not of this
  flag, and is unchanged: use a maximum-likelihood method such as `"YN"` for distant pairs.)
- **Input that is not a codon alignment is an R error**, in the C++ code and before it is called: a sequence
  shorter than one codon aborted the R session ("irrecoverable exception"; an unsigned `length - 3`), sequences of
  unequal length gave numbers read from beyond the end of the shorter one (dN 0.50, dS 0.26 in the audit's test),
  a length that is not a multiple of 3 gave a console note. `gestimator()` also reads files with Windows line
  ends (the carriage returns were counted as sequence) and frees its memory (about 1.5 kB were lost per gene
  pair).

`promotor_divergence_of_orthologous_genes()`
- **Promoters are paired by ortholog pair.** Query and subject promoters were selected independently and paired
  by position. When each species lacked the promoter of a different gene, both sets had the same length, the pairs
  in between were shifted, distances were computed between promoters of non-orthologous genes and the final join
  dropped them: 4 valid pairs in, 2 rows out, no message. Now 4 rows; a message says how many pairs lack a promoter
  in which species. Promoter names are matched to gene IDs as `<gene>` or `<gene>_<anything>`; they were cut at
  the first `_`, which broke IDs such as `NM_001234` or `LOC_Os01g01010`.

Changes of results and row counts in this section
- `filter_dNdS()`, `divergence_stratigraphy()`: gene pairs with dS = 0 are kept (more rows; `dNdS` / stratum `NA`;
  `filter_dNdS()` returns the extra column `dS_zero`). *E. coli* K-12 vs Sakai in the audit: the 85 of 1,280 pairs
  that were dropped are kept.
  Restored by `dS_zero = "drop"` or `options(orthologr.dS_zero = "drop")`.
- `divergence_map()`, `divergence_stratigraphy()`: tables that stopped with "duplicate threshold" now give
  strata; tables that worked give the same strata (except `n_quantile` 49 / 98, where the `NA` is gone).
- Comeron dN / dS / dN/dS: `NA` instead of a finite value for saturated pairs (numbers above). Restored by
  `options(orthologr.legacy_comeron_saturation = TRUE)`.
- `promotor_divergence_of_orthologous_genes()`: rows that were lost are returned; no option.
- (N11, N12 of the audit - serial and parallel code of `compute_dnds()` differing for `store_locally = TRUE`, and
  worker processes left behind after an error - were fixed by `fix/robustness` / `perf/vectorise`: one function for
  both paths, `on.exit(parallel::stopCluster())`.)

### Golden-value tests (second audit)

No test asserted a single dN/dS value, a known ortholog pair, or an ID with `|`, digits only or more than 15
characters. `tests/testthat/test-golden_values.R` and `tests/testthat/fixtures/` add:
- **Pinned dN, dS and dN/dS** of four fixed codon alignments for Comeron's and Li's method, cross-checked against
  an independent program (KaKs_Calculator 2.0, method MLPB): Comeron and Li agree with it within 12% on these
  pairs, the test states a tolerance of 15% (or 0.005); hand-made alignments with no, one synonymous and one
  nonsynonymous change.
- The **20 known ortholog pairs** of the example data with their values, end to end through `dNdS()`, for both
  best-hit rules; **reciprocal best hits are symmetric** (query and subject swapped; DIAMOND and BLAST).
- A **tie-break fixture** (a hit table file in which every level of the two best-hit rules decides once).
- **IDs end to end**: UniProt / NCBI `|` IDs, numeric IDs, `007`, IDs with a hyphen, 60-character IDs through
  `dNdS()` with Comeron, Li, YN, MAFFT, Clustal Omega and through `divergence_stratigraphy()`.

Tests that passed for the wrong reason, or could fail without a defect, were repaired:
- `test-deepclust.R` "stricter approx_id yields at least as many rows": `deepclust()` returns one row per sequence,
  so the test compared 20 with 20. It now clusters proteins with near-identical copies and compares the numbers
  of clusters.
- The `save.output` tests of `deepclust()` / `deepclust_realign()` wrote into the shared `tempdir()`, where the
  file of an earlier test existed already; they passed although the file was never overwritten. They use a folder
  of their own, start from a stale file and compare the content.
- `test-filter_dNdS.R` asserted that a gene pair with dN = dS = 0 is dropped (see "dS = 0" above).
- `test-diamond.R` "more-sensitive returns at least as many hits as fast" is a tendency of a heuristic, not a
  guarantee; the test now requires the same top hit for each of the 20 example genes in both modes.
- The `comp_cores = 2` test of `divergence_stratigraphy()` is skipped on machines with one core.

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
