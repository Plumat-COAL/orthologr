# The upstream-v2 branch series

A series of branches on top of drostlab/orthologr master (06c3098, version 0.4.3). Each branch builds on the
previous one and passes R CMD check on its own (GitHub Actions: Linux release / devel / oldrel, Linux without
any optional program, Windows, macOS). The least controversial changes come first; the branches that change
results come after the ones that do not; COAL appears only in the last branch. The tree of the last branch is
the tree of `integration/all` plus this file.

`tests/testthat/test-legacy.R` (from branch 01 on) compares `dNdS()` on the 20-gene example with tables stored
from 0.4.3 (Comeron and YN). Every branch keeps it passing: a branch that changes a default adds the option that
restores it to that test, so the test shows that the old numbers can still be had.

| Branch | What it changes | Results change? | Old behaviour |
|---|---|---|---|
| `upstream-v2/01-ci-and-tests` | GitHub Actions (with and without the external programs, Linux / macOS / Windows); stored 0.4.3 tables and `test-legacy.R`; no package code | no | - |
| `upstream-v2/02-robustness` | a failing gene pair is an `NA` row with a reason instead of ending the run; per-call working folder; same-named input files; hyphenated IDs; Li's method numeric | only fixes (same-named files, hyphenated IDs, Li's columns) | - |
| `upstream-v2/03-tool-calls` | every program through `run_tool()` (exit status and output checked, no `setwd()`); parameter strings are argument lists, not shell text (R07) | no (`;`, pipes in parameter strings no longer work) | - |
| `upstream-v2/04-options-and-inputs` | arguments checked up front, missing programs found before the search (R09), documented options work, MUSCLE 5 / T-Coffee, input fixes | T-Coffee default mode; `aa_aln_params` now reaches the program | `aa_aln_params = "-mode accurate"`; do not set `aa_aln_params` |
| `upstream-v2/05-vectorise` | each file read and translated once, vectorised best hits / isoforms, `parallel` instead of foreach / doParallel; PSOCK workers on Windows (R08, R13; working folder handed over) | no | - |
| `upstream-v2/06-deps-and-docs` | version 0.5.0, optional packages to Suggests, README / vignettes / man pages corrected, `environment.yml`, pkgdown, native routines registered by Rcpp | no | - |
| `upstream-v2/07-ids-and-joins` | IDs returned as in the files, results joined by gene pair, duplicate IDs an error, dropped sequences reported, isoform warning, same ID in both files (R06) | only fixes (lost or changed IDs) | `options(orthologr.isoform_warning = FALSE)` silences the new warning |
| `upstream-v2/08-annotation` | longest isoforms and gene loci on real NCBI / Ensembl GFF3 / GTF; long annotation lines kept (R03) | only fixes | - |
| `upstream-v2/09-pairwise-alignment` | pairwise protein alignments: BLOSUM62, full length, no terminal stop; full-length rows built directly (about twice as fast) | **yes**: dN or dS of about a third of the pairs of the 1,000-gene example, mostly slightly (r = 0.988) | `options(orthologr.legacy_pairwise_aln = TRUE)` |
| `upstream-v2/10-strata-and-saturation` | tie-safe strata, dS = 0 policy, gestimator input checks, Comeron saturation (R17), dN = 0 from KaKs_Calculator (R14), left-out pairs reported (R23) | **yes**: Comeron NA for a saturated dN / dS (1 of 995 pairs); dN = 0 instead of NA (11 of 995) | `options(orthologr.legacy_comeron_saturation = TRUE)`, `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| `upstream-v2/11-best-hit-rule` | best hits by bit score; order-independent tie-breaks (R20) | **yes**: 3.2 % of human / mouse best hits | `best_hit_rule = "evalue"` or `options(orthologr.best_hit_rule = "evalue")` |
| `upstream-v2/12-golden-values` | pinned dN / dS of fixed codon alignments, known ortholog pairs, tie-break fixture; fragile tests repaired | no (tests only) | - |
| `upstream-v2/13-orthogroups` | new `orthogroups()`, `write_orthogroups()`, `dNdS_orthogroups()` (engine "orthologr"); `orthofinder2()` soft-deprecated | no (new functions) | - |
| `upstream-v2/14-run-record-and-legacy` | column `status`, summary message, record of the run, `write_dNdS()` / `read_dNdS()`, `orthologr_legacy_options()`, `?orthologr_options` (R22), upgrade article; classed errors and warnings (`orthologr_error_*`, `orthologr_warning_*`), `options(orthologr.ambiguous_cds = "keep")`, a hint for `sensitivity_mode` with distant species, malformed-input tests | one more column | `options(orthologr.status_column = FALSE)`; all at once: `orthologr_legacy_options()` |
| `upstream-v2/15-coal` | optional COAL engine (`coalr`): `dNdS(engine = "coal")`, `dNdS_coal()`, `dNdS_clusters()`, `codon_aln(tool = "coal")`, `substitutionrate(kaks_tool = "coal")`; maintainer notes; this file | no for calls that do not ask for COAL | do not use `engine = "coal"` |

The fixes of the third review (R01 - R23) are in the branches they correct (named in the table and in the
commit messages); there is no separate fix-up branch.

Differences from `integration/all`: this file, and its line in `.Rbuildignore`.
