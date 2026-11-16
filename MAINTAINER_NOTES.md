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
| `inst/pal2nal/pal2nal_batch.pl` (new on this fork, 0.6.0; orthologr's own code, written for orthologr, not derived from PAL2NAL - it loads `pal2nal.pl` and runs it once per job) | GPL version 3: the licence header ("Copyright (C) 2026 Chloe Newton and Keanu Murati (orthologr contributors). This file is part of orthologr ... GNU General Public License version 3") was added with the release commit of 0.6.0 (REVIEW-OL-3 had left it out until the package licence was decided). Its authors are the fork's maintainers, so this header does not wait for the consent below |

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

**When the copyright holders have agreed: the one commit to make before the tag `v0.6.0`.** Release 0.6.0
(branch `release/0.6.0`) is prepared with the consent-pending wording as audited (AUDIT-OLHYG-2), and is not to be
published before the consent is given. Once it is given (in writing; keep the record outside the repository), make
one commit on `release/0.6.0` that changes exactly these sentences, and nothing else (all four places carry the
same sentence, word for word; it is the only text that says the release waits):

| File | Sentence now | Change |
|---|---|---|
| `README.md`, section "Licenses", first paragraph | "The change from GPL-2 is pending the agreement of the copyright holders of the code that orthologr published under GPL-2 only (versions up to 0.4.3); this version is not to be released until they have agreed." | replace by "The copyright holders of the code that orthologr published under GPL-2 only (versions up to 0.4.3) agreed to the change." (or remove it) |
| `inst/COPYRIGHTS`, the paragraph under "orthologr: Copyright (C) Hajk-Georg Drost and the orthologr contributors. GPL-3." (its last sentence, lines 12-14) | the same sentence, wrapped over three lines | the same replacement (or remove it), keeping the four-space indent |
| `NEWS.md`, 0.6.0, the item "**Licence: GPL-3** (was GPL-2)" | its last sentence, the same sentence | the same replacement (or remove it) |
| `NEWS.md`, the same item | "`MAINTAINER_NOTES.md`, section 1, records what is still open (the consent of the orthologr copyright holders for their code that was published under GPL-2, and PAL2NAL's licence)." | "`MAINTAINER_NOTES.md`, section 1, records what is still open (PAL2NAL's licence)." |
| `MAINTAINER_NOTES.md` (this file, not in the package), section 1, "Still open", the bullet "**Consent of the orthologr copyright holders.**" | the whole bullet | replace by one line: who agreed and when (no private correspondence in the repository), and remove this paragraph and its table |

Then check that nothing else says it: `grep -rn -i "pending\|not to be released\|until they have" README.md NEWS.md
inst/ MAINTAINER_NOTES.md` should find nothing about the licence. `Date:` in `DESCRIPTION` may be moved to the day of
that commit. No code, test or other file changes, so the identity, `R CMD check` and CI results of the release
commit stay valid; re-run `python3 -I .github/scripts/release_tools.py version v0.6.0` and let CI run on the new head
before tagging. The draft release notes (`RELEASE_NOTES_orthologr_0.6.0.md`, kept outside the repository) say
GPL-3 without a consent sentence: they are written for publication after the consent.

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
| `DESCRIPTION`, `Remotes:` | `coalr=Plumat-COAL/COAL/packaging/r/coalr@v0.6.1` - pinned to COAL's release tag v0.6.1 (since the release commit of 0.6.0, so that the coalr the release is checked with is the coalr of the COAL release that carries `coalr_0.6.1.tar.gz`; before that v0.6.0: Comeron / Li, COAL #76's refinements and GY, `coalr::coal_info()$features` `"gy94_ml"`; the pin was main 03612da before, which has no GY, so the GY tests of the COAL engine were skipped in CI) (syntax `package=owner/repo/subdir[@ref]`; the `coalr=` prefix is needed because the package is not named like its repository - without it pak, which GitHub Actions uses, reports "Can't find package called coalr"). The install hints of `R/coal.R` (`coalr_install_hint()`) name the same tag | the tag of a newer COAL release when one is needed (move `coalr_install_hint()` with it) |
| `DESCRIPTION`, `Remotes:` | `drostlab/metablastr@02dad6daf9f4a0675317177ff434ae9bbcec535f` - the head of metablastr's default branch on 2026-11-16, pinned because the release workflow refuses an unpinned `Remotes` entry (section 6) | a newer commit when one is needed |
| install command of orthologr in `README.md` (two places), `environment.yml`, `vignettes/orthologr.Rmd` | `remotes::install_github("Plumat-COAL/orthologr@v0.6.0")` (the release tag; it was `@integration/all` before the release commit), so that the command installs the version the page describes | `remotes::install_github("drostlab/orthologr")` |
| `.github/workflows/R-CMD-check.yaml`, `on: push: branches` | also `integration/**`, `upstream/**`, `upstream-v2/**`, `fix/**`, `chore/**`, `release/**` (the fork's working and release branches; `release/**` since the release commit of 0.6.0, so that the release branches are checked on push: `workflow_dispatch` cannot be used while the workflow with it is not on `master`) and `workflow_dispatch` | the fork branches can be removed; `workflow_dispatch` is worth keeping |

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

## 6. Releases of the fork, and the COAL + orthologr bundle

Two workflows, modelled on COAL's `release.yml`; the checks they share are `.github/scripts/release_tools.py`
(standard library Python, with its tests in `.github/scripts/test_release_tools.py`, which both workflows run first).
Every action is pinned to a full commit SHA, the token is read-only except in the one job that writes, and the
inputs reach the shell through environment variables only. Neither workflow ever writes the notes of an orthologr
release, and neither runs anything it downloads from a release.

**Where the workflows must be for each way of starting them.** The default branch of this fork is `master`, the
upstream 0.4.3 line, which the team does not push to.

- The **release event** (publishing a release `v...`) runs the `release.yml` of the tagged commit: it needs no
  change to `master`, only `release.yml` and `release_tools.py` in the commit that gets tagged.
- **"Run workflow"** (and `gh workflow run` / the dispatch API) is offered only for a workflow file that is on the
  default branch; the run then uses the workflow file and the script of the branch chosen in "Use workflow from".
  So the manual runs of `release.yml` (dry run, re-upload) and **every run of `bundle.yml`, which is started by
  hand only, cannot happen until `bundle.yml` (and `release.yml`) are on the default branch**; dispatching them
  now answers HTTP 404. Either merge `.github/workflows/{release,bundle}.yml` and `.github/scripts/` into `master`,
  or (without touching `master`) switch the default branch (Settings > General > Default branch) to the release
  branch, e.g. `integration/all`, and pick that branch in "Use workflow from". Switching the default branch also
  changes what the repository page and a plain `git clone` show.

**`release.yml`** (on a published release with a tag `v<version>`, or by hand). Job `version`: the tag must be `v`
+ the `Version` of `DESCRIPTION` (not a development version, `.9000`) and the first heading of `NEWS.md` must name
that version; it checks out the tag as `refs/tags/<tag>` (the release event: the event's commit), never a branch of
the same name; it refuses a `Remotes` entry, or a git / GitHub reference in a `Config/Needs/*` field, that is not
pinned to `@<40-digit commit SHA>` (a `@v<version>` tag is accepted for `Plumat-COAL/*` repositories only), and any
`Additional_repositories` field (a dry run only warns). Job `build`: `R CMD build` of that commit, and the SHA-256 of the tarball taken right after the
build, passed on as a job output. Job `check`: `R CMD check --as-cran --no-manual` of a copy of that tarball,
compared with the hash first (Linux, R release, with DIAMOND, BLAST+, Clustal Omega, MAFFT and KaKs_Calculator from
bioconda as in `R-CMD-check.yaml`, and coalr and metablastr from the `Remotes` field); an ERROR or a WARNING fails
the run, NOTEs do not. Job `files` takes the build job's tarball and the check log, compares the tarball with the
build job's hash again and writes `SHA256SUMS.txt` (its own hash is a job output); job `upload` checks the hash of
`SHA256SUMS.txt`, the files against it, the tarball against the build hash and that there are exactly these three
files, then attaches `orthologr_<version>.tar.gz`, `orthologr_<version>-check.log` (the check's `00check.log`) and
`SHA256SUMS.txt`, by name, with `gh release upload --clobber`. A release whose tag does not start with `v` (a
bundle) is ignored. By hand: `tag` empty with `dry_run` on checks the chosen branch; `tag` set with `dry_run` off
re-attaches the files of an existing release.

Where third-party code runs, and what it can reach:

- Job `build` installs only what `R CMD build` needs: R and pandoc (the pinned r-lib actions), the hard
  dependencies of `DESCRIPTION` (Depends, Imports, LinkingTo, from CRAN / Bioconductor through pak) and knitr and
  rmarkdown; no Suggests, no `Remotes` (coalr, metablastr) and no bioconda programs. This code runs before the
  tarball exists and is hashed. The job uses no package cache (`cache: false`): setup-r-dependencies restores a
  cache by a key prefix, so on a miss it could restore a library saved by the check job or by `R-CMD-check.yaml`,
  which ran untrusted code; it installs from CRAN / Bioconductor every time instead. The vignettes run their examples only when DIAMOND or BLAST+ is on the PATH, so they
  are built as CRAN would build them, as text and code without results (R CMD build of `integration/all` passes
  without the programs, coalr and metablastr: tested); none of their evaluated chunks uses coalr, metablastr or
  another Suggests.
- Job `check` runs much more: the bioconda packages, all Suggests and the `Remotes` (which a release must pin, as
  above), the package's tests; it may use the dependency cache. Nothing it writes is published except its log; it cannot change the tarball the
  release gets (the `files` and `upload` jobs compare it with the build job's hash).
- Jobs `version`, `files` and `upload` run only the repository's own `release_tools.py`, coreutils, `gh` and
  GitHub's checkout and artifact actions.

To publish orthologr X.Y.Z:

1. `release.yml` and `.github/scripts/` must be in the commit that gets tagged: merge `ci/release` into the release
   branch first (and see above for the manual runs).
2. In the release commit: set `Version: X.Y.Z` (and `Date:`) in `DESCRIPTION`; make the first heading of `NEWS.md`
   name X.Y.Z; move the `Remotes` pin of coalr (`coalr=Plumat-COAL/COAL/packaging/r/coalr@v...`, and
   `coalr_install_hint()` in `R/coal.R` with it, section 3) to the COAL release the bundle will carry (e.g. `v0.6.1`),
   so that the coalr the release is checked with is the coalr the bundle ships; pin `drostlab/metablastr` to a
   commit SHA (`drostlab/metablastr@<sha>`; a tag is accepted for Plumat-COAL repositories only): the version job
   refuses a release with an unpinned `Remotes` entry. Commit, push.
   Optionally: Actions > Release > Run workflow on that branch, `tag` empty, `dry_run` on.
3. On GitHub: Releases > Draft a new release > tag `vX.Y.Z` (create it on that commit) > title and notes written
   by hand > Publish. The workflow starts by itself and attaches the three files (about 1 to 2 hours; most of it is
   the check).
4. If it fails: fix, and either delete the release and the tag and start again from 2, or (same commit) Actions >
   Release > Run workflow with `tag` `vX.Y.Z` and `dry_run` off. Such a re-run builds the tarball again, and
   `R CMD build` does not make a byte-identical tarball twice (the `Packaged:` time it writes into
   `DESCRIPTION`, the file times in the archive), so the files and `SHA256SUMS.txt` of the release change. A bundle made from the earlier
   files is not checked again by anything: it simply goes stale (its tarball and checksums no longer match the
   orthologr release). Make the bundle again under a new bundle tag; the old tag cannot be reused while its tag and
   release exist (bundle.yml refuses a taken tag), so delete the stale bundle release and tag by hand if wanted.

**`bundle.yml`** (by hand only; see above: not before it is on the default branch). Inputs `coal_tag` (a COAL
release, e.g. `v0.6.1`), `orthologr_tag` (a release of this repository, e.g. `v0.6.0`), `bundle_tag`
(`bundle-<name>`, e.g. `bundle-2026.11`: never `v...`, so it cannot collide with a tag of either project) and
`dry_run` (on by default). It downloads every file of both releases (`gh release download`), checks each against
its release's `SHA256SUMS.txt` and against GitHub's own record of the file (SHA-256 digest and size), and stops if a
release is missing, a draft or a prerelease, if a file is missing, extra or unlisted, or if any checksum differs.
The bundle tag must not exist yet: only a 404 from GitHub for both the tag and the release counts as free, any other
API error stops the run (checked again just before the release is created). The COAL release must contain
`coalr_<version>.tar.gz` (COAL releases do from 0.6.1 on; 0.6.0 has none) and the orthologr tarball must export
`orthologr_sitrep()`, which `INSTALL.md` runs. It then creates, in this repository, the release `bundle_tag` (tag
on the commit of `orthologr_tag`; not a prerelease; not marked "Latest", so that `releases/latest` stays
orthologr's own release) with `coalr_<v>.tar.gz` and `orthologr_<v>.tar.gz` (unchanged), `INSTALL.md` (install
coalr, then orthologr, then `library(orthologr); orthologr_sitrep()`) and `SHA256SUMS.txt` of those three; its
notes (`--notes-file`) list the files and their SHA-256 and end with the authors line. A dry run does all of this
except creating the release; the files and the notes are the workflow artifact `bundle`.

To publish a bundle: Actions > Bundle > Run workflow with the three tags and `dry_run` on; read the artifact
(`INSTALL.md`, `notes.md`); run it again with `dry_run` off. From a shell:
`gh workflow run bundle.yml -R Plumat-COAL/orthologr -f coal_tag=v0.6.1 -f orthologr_tag=v0.6.0 -f bundle_tag=bundle-2026.11 -f dry_run=true`
(add `--ref <branch>` for "Use workflow from").
