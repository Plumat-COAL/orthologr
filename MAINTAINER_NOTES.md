# For the maintainer

Notes that need a decision by the maintainer of `drostlab/orthologr`. This file is not part of the built package
(`.Rbuildignore`).

## 1. Licence: changed to GPL-3 on this fork in 0.6.0 (finding R05 of the third review; KaKs_Calculator finding of the final review)

**What was changed** (branch `chore/gpl3`, decided by the owner of the fork, Chloe Newton, choosing option (a)
below):

| Where | Before | Now |
|---|---|---|
| `DESCRIPTION` | `License: GPL-2` (version 2 **only**) | `License: GPL-3`, `Copyright: See file inst/COPYRIGHTS` (no other field changed) |
| `LICENSE`, `COPYING` | text of the GNU GPL version 2 | the GNU GPL version 3, verbatim: copied from `share/licenses/GPL-3` of the R 4.5.3 installation used for the checks; SHA-256 3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986, the same bytes as `/usr/share/common-licenses/GPL-3` (Debian) and as the `LICENSE` of KaKs_Calculator 2.0 (kullrich/kakscalculator2) |
| `inst/COPYRIGHTS` | the third-party parts, no licence for the package | begins with "orthologr: Copyright (C) Hajk-Georg Drost and the orthologr contributors. GPL-3." and why |
| `README.md`, section "Licenses" | quoted the libsequence header; said (until 0.6.0) that no KaKs_Calculator code is part of the package | says GPL-3 and why; names the KaKs_Calculator-derived functions and PAL2NAL; points to `inst/COPYRIGHTS` |

Not changed: the `@author` lines of the R files, `Authors@R`, and the file headers in `src/` (the libsequence
headers still say that the orthologr changes to them are "free under the terms of GNU General Public License
version 2 of the License, or any later version", which permits GPL-3).

**Why** (the facts that remain true):

| Part | Licence as stated by its authors |
|---|---|
| `src/gestimator.cpp` and the headers `src/Comeron95.h`, `Comparisons.h`, `Grantham.h`, `GranthamWeights.h`, `Kimura80.h`, `RedundancyCom95.h`, `Sequence.h`, `SingleSub.h`, `Sites.h`, `ThreeSubs.h`, `TwoSubs.h` | libsequence (Kevin Thornton): GNU GPL "either **version 3** of the License, or (at your option) any later version" (line 11 of the headers); the orthologr changes to them (line 27): "version 2 of the License, or any later version" |
| `R/substitutionrate.R`: `kaks_transl_table()`, `kaks_lwl_tables()`, `kaks_lwl_counts()` (added on this fork in 0.5.0 / 0.6.0) | derived from **KaKs_Calculator 2.0**: the genetic-code strings are copied from `transl_table` in `base.cpp`; the two other functions translate into R the site and transition / transversion counts of `LWL85.cpp` (`CountSiteAndDiff()`, `getCodonClass()`, `TransitionTransversion()`, `MLWL85::TransitionTransversion()`) and `LPB93.cpp` (`MLPB93::TransitionTransversion()`), step for step, so that the counts are KaKs_Calculator's to the last bit. Also derived: the class test of `kaks_saturated()` (1 - 2P - Q > 0, 1 - 2Q > 0) follows the Kimura two-parameter step of `LWL85::preProcess()`; `kaks_calc_header()` holds the column names of KaKs_Calculator's output (an interface fact, listed for completeness). KaKs_Calculator 2.0: Copyright (C) 2009 Beijing Institute of Genomics of CAS. Licence facts: its source files carry only "Copyright (C) 2009, Beijing Institute of Genomics of CAS / All rights reserved" and **no licence and no "or any later version"** statement; the original `KaKs_Calculator2.0.tar.gz` files, imported from SourceForge into github.com/kullrich/kakscalculator2 (the repository bioconda builds) in commit f646786 (2020-11-12), have **no licence file**; that repository's `LICENSE`, the text of the **GPL version 3**, was added by Kristian K. Ullrich in its initial commit bd24013 (2020-11-12), before the import; the SourceForge project page is said to list GPLv3 (not re-checked). Section 14 of the GPL version 3 lets a program that specifies no version number be used under any published version; whether that applies here (does the added `LICENSE` specify version 3?) is for the maintainers and the copyright holders to decide - these notes treat the code as GPL version 3 |
| `inst/pal2nal/pal2nal.v14/` | PAL2NAL v14: **no licence statement** in any of its files (`README`, `pal2nal.pl`, `for_paml/`, test files; only "Mikita Suyama" in the header of `pal2nal.pl`). Debian's package records `License: GPL-2` and (Debian's record, not the files') "Copyright: 2004-2011 Mikita Suyama, David Torrents, and Peer Bork and their respective employing institute" (`debian/copyright` of pal2nal 14.1-3, https://metadata.ftp-master.debian.org/changelogs/main/p/pal2nal/pal2nal_14.1-3_copyright), after a Debian Med thread "Lack of license for pal2nal" (2017-2018). **Unverified**: the upstream statement (PAL2NAL web site, Suyama, Torrents & Bork 2006, NAR 34:W609) was not checked |
| optional package **coalr** (COAL's R package; `dNdS(engine = "coal")`, `dNdS_coal()`, `dNdS_clusters()`, `dNdS_orthogroups(engine = "coal")`) | **GPL-3** (only). In `Suggests`, loaded only when installed; orthologr does not contain its code |

According to the Free Software Foundation's GPL FAQ ("Is GPLv3 compatible with GPLv2?": no, for code licensed
under GPLv2 only), code under "GPL version 3 or later" and code under "GPL version 3" cannot be combined into a
work licensed "GPL-2 only"; on that reading `License: GPL-2` could not be right. **GPL-3 rather than GPL (>= 3)**: the
KaKs_Calculator-derived functions are treated as GPL version 3 without "or later" (the GPL-3 `LICENSE` that the
repository bioconda builds added; see the row above, including the section 14 question), so the package as a
whole can only be distributed under GPL-3. (Alternatives that were not taken: (b) ask the
KaKs_Calculator authors whether "or any later version" applies, which would allow GPL (>= 3); (c) replace the
three functions by code written from the publications - Li, Wu & Luo 1985; Tzeng, Pan & Li 2004 - and the NCBI
genetic codes and check it against KaKs_Calculator's output; they only tell, for a few values near 0, whether
KaKs_Calculator's printed estimate is saturated. Either would still be possible later.)

**Still open - for the maintainer and the copyright holders:**

- **Consent of the orthologr copyright holders.** The code of orthologr up to 0.4.3 (Hajk-Georg Drost, Sarah
  Scharfenberg, Jaruwatana Sodai Lotharukpong and other contributors) was published under `License: GPL-2`
  (version 2 only; only the libsequence-file changes say "or any later version"). Relicensing that code under
  GPL-3 is the decision of its copyright holders, not of this fork. The owner of the fork is obtaining their
  consent, and this branch is not to be released until it is given (the README, `inst/COPYRIGHTS` and NEWS say
  so where they state GPL-3). Whether, until then, the GPL-3 statement can apply to anything beyond the fork's
  own additions and the third-party parts is a legal question for the copyright holders, not settled here.
  (On the FSF's reading cited above, the GPL-2 statement before this change did not fit the GPL-3 parts either.)
- **PAL2NAL.** If PAL2NAL is GPL-2 *only* (Debian's `License: GPL-2` does not say whether "or later" applies),
  the FSF's GPL FAQ (on "mere aggregation") suggests that shipping it next to a GPL-3 package and running it as a
  separate program (`perl pal2nal.pl`, communicating through files and the command line) is aggregation rather
  than one combined program, whereas `inst/pal2nal/pal2nal_batch.pl` (orthologr's code, GPL-3) loads and compiles
  `pal2nal.pl` into its own Perl process, which by the same FAQ may count as one combined program. These are
  readings for the maintainers and copyright holders to confirm, not conclusions of these notes. Asking the
  PAL2NAL authors for their licence would clarify it (a GPL-2-or-later or more permissive answer would, on the
  FSF's reading, remove the question); the per-pair `perl pal2nal.pl` route (`options(orthologr.batch_tools =
  FALSE)`) does not load PAL2NAL into orthologr's process.
- **Upstream (drostlab/orthologr).** The upstream series (`upstream-v2/*`, `upstream-v3/*`) still says GPL-2;
  whether upstream takes the GPL-3 change is the upstream maintainer's decision.

## 2. Defaults and results that differ from orthologr 0.4.3, and the option that restores each

| Change in 0.5.0 | Restored by |
|---|---|
| Best hits are ranked by bit score (was: e-value, then alignment length) | `options(orthologr.best_hit_rule = "evalue")` or the argument `best_hit_rule = "evalue"` |
| Pairwise protein alignments (BLOSUM62; only on the branch with the faster pairwise alignment) | `options(orthologr.legacy_pairwise_aln = TRUE)` |
| Comeron's method: `NA` instead of a finite value when a correction is undefined in a class that has sites | `options(orthologr.legacy_comeron_saturation = TRUE)` |
| KaKs_Calculator methods: `dN = 0` instead of `dN = NA` for a pair without nonsynonymous change; dN = dS = 0 for identical sequences, dS = 0 for nonsynonymous-only pairs | `options(orthologr.legacy_kaks_zero_dn = TRUE)` |
| KaKs_Calculator methods: `NA` (status `"saturated"`) for KaKs_Calculator's capped values (Ka / Ks = 99, JC69's cap of YN / MYN, a class with 0 substitutions in the LWL family) | `options(orthologr.legacy_kaks_saturation = TRUE)` |
| Li's method returns numeric columns; a value below 9.9 that R compared as text is no longer `NA` | none (bug fix) |
| dN/dS tables end with the column `status` (`"ok"`, `"dS_zero"`, `"saturated"`, `"no_data"`, `"failed"`) | `options(orthologr.status_column = FALSE)` |
| T-Coffee default mode (was `-mode accurate`) | `aa_aln_params = "-mode accurate"` / `params = "-mode accurate"` |
| `aa_aln_params` is passed on to the alignment program (was ignored) | do not set it |
| Parameter strings are argument lists, not shell text (no `;`, pipes, redirections) | none (security fix) |
| A missing program (perl, aligner, KaKs_Calculator) is an error before the search (was: all pairs `NA` + warning on this fork; an error at the first pair in 0.4.3) | none |
| Annotation lines of more than 65,000 characters are read (were deleted) | none (bug fix); `check_annotation(remove_annotation_outliers = TRUE)` still writes a copy without them |

**Not changed any more** (a default introduced earlier on this fork was taken back after the third review):
`filter_dNdS()` and `divergence_stratigraphy()` leave gene pairs with dS = 0 out, as 0.4.3 did
(`dS_zero = "drop"`); new are a message with the numbers, `attr(, "left_out")` / `divergence_left_out()`, and
`dS_zero = "keep"` / `options(orthologr.dS_zero = "keep")` to keep them.

**One switch:** `orthologr_legacy_options()` sets the five options of this table (and `orthologr.dS_zero =
"drop"`); `options(orthologr.legacy = TRUE)` is the shorthand. With them the 20-gene and 1,000-gene examples give
the tables of 0.4.3 (master 06c3098) for Comeron and YN. `tests/testthat/test-legacy.R` compares the 20-gene
result with tables stored from 0.4.3 (`tests/testthat/fixtures/dnds_0.4.3/`, made by `make_fixture.R` there with
the 0.4.3 package installed) on every run, for both ways of setting the switch.

## 3. Things that are set for the fork and must change when this is merged

| Where | Now | Change to |
|---|---|---|
| `DESCRIPTION`, `Remotes:` | `coalr=Plumat-COAL/COAL/packaging/r/coalr@03612da73e41a044e2385acb4ac9d599a8a49e8b` - pinned to a commit of COAL's default branch (main 03612da, with COAL #76), which holds the R package coalr since COAL's pull request #54 was merged (syntax `package=owner/repo/subdir[@ref]`; the `coalr=` prefix is needed because the package is not named like its repository - without it pak, which GitHub Actions uses, reports "Can't find package called coalr") | a release tag (`...coalr@v<version>`) once COAL has one with the Comeron / Li methods |
| install command of orthologr in `README.md` (two places), `environment.yml`, `vignettes/orthologr.Rmd` | `remotes::install_github("Plumat-COAL/orthologr@integration/all")`, so that the command installs the version the page describes | `remotes::install_github("drostlab/orthologr")` |
| `.github/workflows/R-CMD-check.yaml`, `on: push: branches` | also `integration/**`, `upstream/**`, `upstream-v2/**`, `fix/**`, `chore/**` (the fork's working branches) and `workflow_dispatch` | the fork branches can be removed; `workflow_dispatch` is worth keeping |

GitHub Actions runs on the fork. All six legs (Linux release / devel / oldrel-1, Linux without the optional
programs, Windows, macOS; the programs installed with micromamba) passed on `upstream-v2/15-coal` (run
43507540285) and on `integration/all` (run 43505256122, macOS leg re-run after a cancellation), and on
`upstream-v2/01-ci-and-tests` and `upstream-v2/08-annotation`; the runs of the other upstream-v2 branches were
cancelled because the fork's runners start about one job at a time (see `STACK.md` on `upstream-v2/15-coal`). The
CI installs coalr through the `Remotes` field; if the COAL repository cannot be read by the workflow's token, or the ref
disappears, the dependency step fails on every leg - remove the `Remotes` line (coalr is in `Suggests` and every
use is guarded) rather than leave it dangling.

## 4. Proposals that were deliberately not made in 0.5.0, and what 0.5.0 does about them

- **`sensitivity_mode = "default"` as the default of `dNdS()`, `orthologs()`, `orthogroups()`** (now `"fast"`,
  DIAMOND's mode "for hits of more than 90 % identity"). Measured with `diamond_rec()`, 2 cores:

  | reciprocal best hits | `"fast"` | `"default"` | `"sensitive"` |
  |---|---|---|---|
  | *A. thaliana* / *A. lyrata*, 1,000-gene example | 983 | 985 | 985 |
  | *S. cerevisiae* / *C. glabrata* | 4,164 (7.3 s) | 4,384, +5.3 % (8.2 s) | 4,502, +8.1 % (24.4 s) |

  The pairs found only with `"default"` have a median identity of 37 % and a median length of 298 amino acids
  (not only short proteins); the search is about 10 s of a 200 s `dNdS()` run. It would be the third change of
  the ortholog set in one release, so 0.5.0 only documents it (`?dNdS`, the dN/dS vignette). Suggested for 0.6.0.
- **Keeping coding sequences with `N`** (translate the codon to `X` and let the estimator skip it) instead of
  leaving the sequence out. The loss is 0 to 1 sequence in RefSeq sets but 2.5 % (human chr21) to 4.7 % (mouse
  chr10/16/17) of Ensembl transcripts, typically one leading `N`. The default is unchanged ("drop"); 0.5.0
  documents the behaviour, reports the sequences, and offers it as an opt-in: `options(orthologr.ambiguous_cds =
  "keep")` (recorded in `attr(, "orthologr_run")$options`). Making "keep" the default is a decision for 0.6.0.

**Done in 0.5.0** (were listed here as not made):

- **Classed conditions**: `orthologr_error_*` / `orthologr_warning_*` (`?orthologr_options`, "Errors and
  warnings"). The argument and input-file checks of the main functions (`dNdS()`, `dNdS_coal()`,
  `dNdS_clusters()`, `orthogroups()`, `dNdS_orthogroups()`, `divergence_stratigraphy()`, `filter_dNdS()`, the
  alignment functions, `substitutionrate()`, the readers) and every check for a program or package are classed;
  105 calls are classed; about 190 `stop()` calls, mostly in other, older functions and inside single gene
  pairs, are still plain errors (`?orthologr_options` says which functions are classed).
- **Shared code of `dNdS()` and `dNdS_coal()`**, first part: the search step is one function, `dnds_search()`.
  Still duplicated: the argument checks, the start message and the set-up of the call (working folder, sequence
  cache, collection of left-out sequences) - a second step for 0.6.0.

## 5. Changes of this pass that touch existing behaviour

- `orthogroups()` and `dNdS_orthogroups()` have `delete_corrupt_cds = FALSE` (trim), the default of `dNdS()`;
  they had `TRUE`. `dNdS()` itself is unchanged. (`diamond_rec()`, `diamond_best()`, `deepclust()` called on
  their own keep `TRUE`.) For two species `orthogroups()` therefore returns the reciprocal best hits of
  `dNdS()` / `diamond_rec(delete_corrupt_cds = FALSE)`: 995 instead of 983 pairs on the 1,000-gene example.
- `dNdS()` output: a new last column `status`; one start line and one summary instead of 31 messages (also on
  real data: sequences that are left out or trimmed to a multiple of 3 are counted in the summary and listed in
  `attr(, "dropped_sequences")` / `attr(, "trimmed_sequences")`; 1,000-gene example 2 messages, `quiet = FALSE`
  36); nothing on `stdout`. `divergence_stratigraphy()` has `quiet = TRUE` too (was `FALSE`; 45 -> 10 messages).
  Scripts that parse the messages or take columns by position from the end are affected.
- `dNdS(quiet = FALSE)` no longer prints one line per alignment and estimation call; it reports progress.
- The detailed, review-by-review log that used to be the NEWS entry is `MAINTAINER_CHANGELOG.md`.

The first Actions run on the fork failed before any R code ran: micromamba received `"diamond>=2.1"` with its
double quotes (inside a folded YAML block no shell removes them) and could not parse the version, and pak could
not match the `Remotes` entry to coalr. Both are corrected (`diamond>=2.1` unquoted in the two workflows; the
`coalr=` prefix). The later runs pass on all legs (section 3).
