# Clean-room log: `kaks_transl_table()`, `kaks_lwl_tables()`, `kaks_lwl_counts()` (branch chore/gpl2-cleanroom)

Why: the three internal functions in `R/substitutionrate.R` (fork commits ad753f2, 0339835) were written by
translating KaKs_Calculator 2.0's GPL-3 C++ source, but orthologr is GPL-2. This branch deletes them and has a new
implementation in `R/lwl_counts.R`. It was written only from published methods and public data, and checked against
the *output* of the KaKs_Calculator 2.0 binary (`/home/claude/tools-g/bin/KaKs_Calculator`, "Version: 2.0, June.
2009"), which was used as a black box.

Agent: OL-CLEANROOM (Claude), 15–16 Nov 2026. Base: integration/all 7d25e35.

## What was and was not read

Not read: KaKs_Calculator's source (any copy), the bodies of the three old functions, their git history, and
MAINTAINER_NOTES.md section 1.

Read:
- Callers and tests: `kaks_saturated()`, `kaks_counts_agree()`, `kaks_columns()`, `kaks_genetic_code()` (all in
  R/substitutionrate.R), and the test "the saturation test of the LWL family follows KaKs_Calculator's own
  numbers ..." in tests/testthat/test-kaks_values.R. From these I took the contract. `kaks_lwl_counts(seq1, seq2,
  method, genetic_code)` returns `list(L, Si, Vi, P, Q)`, each an unnamed length-3 vector for the 0-, 2- and 4-fold
  classes. It returns NULL when it cannot compute (for example `genetic_code = NA`). `kaks_lwl_tables(code,
  modified)$L` is a 64 x 3 matrix with codons in TCAG order (TTT, TTC, TTA, ..., GGG).
- **Disclosure, accidental exposure.** To find the three functions' signatures I ran an `awk` filter for
  `<- function(` lines over the region of the old definitions. Besides the three signatures it printed four
  one-line headers of local helpers inside the old `kaks_lwl_tables()` and `kaks_lwl_counts()`:
  - a codon-index helper (`b[1] * 16 + b[2] * 4 + b[3] + 1`);
  - a by-name variant of it;
  - a helper named `is_arg` that tests membership in a set `arg_pairs`;
  - the start of a helper `tt <- function(c1, c2, pos, sv) {`;
  - in `kaks_lwl_counts()`, `code <- function(s) {` and `add <- function(x) Reduce(`+`, x, 0)`.

  I read no bodies. What this exposure revealed is limited:
  - TCAG codon order. The test needs this anyway.
  - That some special handling of arginine codon pairs exists. The task brief had already said there was
    "special-case handling of particular codon pairs", and Tzeng et al. (2004) describe it.
  - Sequential `Reduce` addition. R's `sum()` uses extended precision, so any bit-exact port needs a sequential sum.

  Every rule in the new code was established independently by the experiments below. In particular, which pairs
  are special, how they are counted in each method, and the summation scheme were each chosen out of competing
  candidates by black-box evidence (section 3).
  When I deleted the old definitions I located them by line types only (blank/comment/top-level). I printed only
  the *name* of the one other top-level object in that span (`kaks_lwl_cache`, a cache used only there) to check
  it had no other users.

## 1. Sources for each rule

| Rule | Source |
|---|---|
| Positions are 0-fold (every change nonsynonymous), 4-fold (every change synonymous) or 2-fold (otherwise; the threefold Ile third position counts as twofold). A change to a stop codon is nonsynonymous. L_i is the mean over the two sequences. | Li, Wu & Luo 1985 (abstract: classes by degeneracy, transitions vs transversions); Tzeng, Pan & Li 2004 (definitions of L0, L2, L4; "threefold … treated as twofold"). Confirmed by KC for all 64 codons x 17 codes (exp 1). |
| Differences are transitions (A<->G, C<->T) or transversions, counted per class (S_i, V_i); P_i = S_i/L_i, Q_i = V_i/L_i. | Li et al. 1985; Li 1993; Pamilo & Bianchi 1993. |
| Codons that differ at 2 or 3 positions: every order of the changes (pathway) has equal weight; pathways through a stop codon are left out. | Li 1993 / Tzeng et al. 2004 ("equal weights (0.5) for the two pathways"). Pathway exclusion and equal weighting were confirmed by KC (exp 2). |
| MLWL/MLPB: CGA<->AGA and CGG<->AGG (first position) count as S2, not V2. ATA's third position is the other exception (synonymous changes are transitions, nonsynonymous ones transversions). | Tzeng et al. 2004, section "Arginine corrections" ("every synonymous difference … at a twofold degenerate site is counted as a transition"). The exact set of codon pairs and how it applies in other codes came from exp 2. |
| Genetic codes 1–6, 9–16, 21–23 (those KC offers with `-c`, from `KaKs_Calculator -h`) | NCBI gc.prt (`https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt`, `ncbieaa` strings, downloaded 15 Nov 2026). These agree with KC's behaviour for every codon and every codon pair in all 17 codes (exp 1, 2). |

## 2. Behavioural details established by black-box experiments

All experiments use the KaKs_Calculator binary with `-m LWL -m MLWL -m LPB -m MLPB` on AXT files I built. They read
the printed "Fold-Sites(0:2:4)", "Fold-S-Substitutions(0:2:4)" (transitions) and "Fold-N-Substitutions(0:2:4)"
(transversions). Scripts are in /home/claude/w-olclean-scratch (kc.py, model.py, exp_*.py, an*.py).

**Exp 1: per-codon classes.** Each codon c (all 64, all 17 codes) is placed in front of a fixed background. The
pair is (c + BG1, c + BG2), and the background's own counts are subtracted.
- 0 of 64 x 17 x 4 rows differ from the NCBI-derived classes.
- A stop codon contributes 0 sites.

**Exp 2: every ordered codon pair (x, y), x != y.** 4032 pairs x 17 codes x 4 methods = 274,176 rows, built the same
way. The model below gives the printed values for all 274,176 rows (`an3.py`, tolerance at printing precision).
- **Pairs that contribute nothing.** A pair where either codon is a stop codon contributes nothing, not even
  sites. If every pathway passes through a stop codon, the sites still count but there are no differences.
- **Half-and-half classes.** A single change u -> v at position p counts half at the class of p in u and half at
  the class of p in v. This was seen where the classes differ (for example in code 13). The two halves are counted
  separately, so a change can be part transition and part transversion.
- **LWL and LPB special case.** For the four ordered pairs CGA->AGA, AGA->CGA, CGG->AGG and AGG->CGG, the half at
  the first codon counts as a transition and the half at the second as a transversion. This happens in every
  genetic code, even where these codons code for different amino acids (codes 5, 9, 13). Evidence:
  - standard code: S2 0.5 / V2 0.5;
  - code 13 (CGA = R, AGA = G): CGA->AGA gives S0 0.5 + V2 0.5, and AGA->CGA gives S2 0.5 + V0 0.5.
- **MLWL and MLPB special cases.** These apply in every genetic code, by nucleotide, whatever the translation:
  - the same four arginine pairs count as transitions on both halves;
  - ATA<->ATT and ATA<->ATC count as transitions;
  - ATA<->ATG counts as a transversion on both halves. This holds even in codes where ATA and ATG are both Met
    (codes 2, 3, 5, 13, 21).

  No other pair is special. In particular AGA<->GGA and CGG<->TGG count by chemistry: every other pair matched the
  plain rule.

**Exp 3: input characters** (`exp_odd.py`, `exp_u.py`).
- A codon pair is dropped (no sites, and Length shrinks by 3) when either codon contains a character other than
  A, C, G, T or U. This covers `-`, `N`, `R`, `X` and `?`, also at a single position.
- Lowercase counts as uppercase.
- U is translated like T, but when changes are counted it is a letter of its own. T<->U is a counted difference
  and a transversion; C<->U is a transition; U<->A and U<->G are transversions. The model with this rule matched
  9012 U-containing pairs (LWL, MLWL) with 0 mismatches. This rule was implemented first (see section 4a) and
  then replaced by "return NULL for sequences with U".

**Exp 4: floating-point summation** (`zero_gen.py`, `schemes.py`, `an_zero.py`, `an_zero2.py`). These decide the
sign of 1 - 2P - Q when it is exactly 0 in rational arithmetic.
- **The cases.** I built 2,370 sequence pairs (codes 1, 2, 5, 13, 23; targets: the 0-, 2- and 4-fold classes; counts
  from the LWL and MLWL models), where one class has 2*S_i + V_i = L_i exactly. Fractional contributions come from
  multi-position codon pairs (1/4, 1/6, 1/12, ...), and codons are shuffled.
- **The observable.** KC's Ka/Ks shows its decision. If the class is saturated, KC uses 0 for it; if not, a term
  0.5*ln(1/a) appears with a ~ 1e-17, so the estimate is > 10. For LWL, the printed estimate gives ln(1/a) to
  about 5 digits. This reveals the magnitude of KC's residual a, not only its sign.
- **The candidates.** Eight summation schemes, from per-half global addition to per-pathway and per-pair sums.
- **The winner.** "For each codon pair: number of halves x 0.5, summed over the valid pathways (exact); divided by
  the number of pathways; added to the running totals one codon pair after another in alignment order." This
  scheme reproduced KC on every case:
  - LWL: 1920/1920, both the decision and ln(1/a) to 1e-3 relative;
  - LPB: 1920/1920 decisions;
  - MLPB: 1170/1170 decisions.

  The best alternative scheme failed on 87 cases. In R the totals are summed with `Reduce(`+`, x, 0)`; `sum()` uses
  extended precision and would differ.
- **Side observations about KC's estimator** (not part of the counts). These only refine my black-box analysis:
  - a saturated class contributes A_i = B_i = 0;
  - a negative A_i is set to 0 in LWL.

## 3. Verification of the R implementation against KC (black box)

Comparison scripts: /home/claude/w-olclean-scratch/cmp_kc.R and cmp_zero.R. They source only `R/lwl_counts.R`.

- **Random and simulated pairs** (`gen_rand.py`): 2 x 10,200 rows (seeds 11 and 12), covering 17 codes x 150 pairs x
  4 methods.
  - Pair lengths are 1–400 codons, with divergence from 2 % to 100 %.
  - 15 % of pairs have random `-`, `N`, `n`, `R` or `U`; 10 % are lowercase; 5 % have T->U.
  - `kaks_lwl_counts()` formatted with `%.6g` (KC's printing) is **character-identical to KC's three printed
    fields in 20,400/20,400 rows**.
- **Boundary cases** (exp 4 set, 8,099 method rows). In R, the sign and magnitude of 1 - 2P - Q (and the 1 - 2Q
  test) give KC's decision in all rows.
- **Existing tests.** The testthat tests pass unchanged (section 4).

## 4. Identity of the package with integration/all

Merged tree 1cd0470 (this branch + integration/all bf8d61b) against bf8d61b. Both were installed from `git archive`
and run with `/home/claude/w-olclean-scratch/run_id.R`, `dNdS()`, `comp_cores = 2`, using the KaKs_Calculator 2.0
binary of tools-g and coalr 0.6.x.

Configurations:
- d20 and d1000 (the package's example data): `dnds_est.method` LWL, MLWL, LPB and MLPB through KaKs_Calculator;
  the default route; `engine = "coal"`; and `engine = "coal"` with LWL (d1000);
- S. cerevisiae / C. glabrata with LWL, and S. cerevisiae / K. lactis with MLPB. These are the two genome pairs
  that the exact-zero fixtures come from. A tracer counted how often `kaks_saturated()` recomputed the counts:
  82 times for ScCg and 90 times for ScKl. d20 and d1000 never reach that path.

Results (`cmp_id.R` / `cmp_cols.R`): in all 15 configurations every column of the result table is identical with
`identical()` and no normalisation. The table attributes differ only in installed-library paths, temp paths and
dates; after normalising those they are `identical()` too.

### 4a. Function-level check after the freeze, and the one change it led to

After the implementation was frozen (commit a1d0712), I ran the old `kaks_lwl_counts()` (installed from
integration/all) and the new one on the same 28,499 inputs: the random sets plus the boundary set. Only their
return values were compared; no code was read. Any mismatch was to be settled by KaKs_Calculator only.
- On every input without `U`, the results were `identical()` (bitwise, all five fields): 25,059 of 25,059. These
  include gaps, N, ambiguity codes, lowercase and all exact-zero boundary cases.
- On all 3,440 inputs with `U`, the old function returned NULL; the new one returned KaKs_Calculator's counts
  (verified against KC's printed output in section 3).

NULL means "the printed counts decide" in `kaks_saturated()`. To keep package results identical to integration/all,
the new function now also returns NULL for sequences containing U/u. This is commit 2 on this branch. It is a
contract-level choice (decline RNA letters), not a port of any code. All KaKs_Calculator comparisons were re-run
after it, and the U rows were dropped from the comparison set.

## 5. Proposed licence wording (README "Licenses" / inst/COPYRIGHTS)

Neither file has the KaKs_Calculator-derived statement on integration/all; that statement is added on
chore/ol-hygiene. If chore/ol-hygiene and this branch both merge, replace its KaKs_Calculator paragraph with this:

> orthologr is free software under the GNU General Public License version 2 (GPL-2). It calls the external programs
> KaKs_Calculator, PAL2NAL, DIAMOND, BLAST, MAFFT and others, which are distributed separately under their own
> licences; no code of those programs is included. The counts of the LWL family that orthologr recomputes in
> `R/lwl_counts.R` (sites and differences per degeneracy class, used to read KaKs_Calculator's saturation
> decisions) were implemented from the published methods (Li, Wu & Luo 1985; Li 1993; Pamilo & Bianchi 1993;
> Tzeng, Pan & Li 2004) and NCBI's genetic code tables, and checked against KaKs_Calculator's output only. They
> contain no code derived from KaKs_Calculator.

For inst/COPYRIGHTS: drop the entry that attributes the three functions to KaKs_Calculator 2.0 (Zhang et al., GPL-3).
If wanted, replace it with: "R/lwl_counts.R: independent implementation (clean room, see CLEANROOM.md on branch
chore/gpl2-cleanroom); the genetic code tables are NCBI data (public domain, US Government work)."

## 6. CI and R CMD check

- **Fork CI** (R-CMD-check.yaml, workflow_dispatch):
  - 7/7 green on a1d0712 (run 43788724616) and 7/7 green on the merged 1cd0470 (run 43796527172);
  - the jobs: ubuntu release / devel / oldrel-1 / 4.4 / no optional tools, macOS, Windows;
  - on ubuntu release (bioconda kakscalculator2 2.0.1) the live black-box test in test-lwl_counts.R ran and passed:
    FAIL 0, 5273 PASS, and no skip for it.
- **Local `R CMD check --no-manual --as-cran`** of 1cd0470: Status 6 NOTEs, 0 errors, 0 warnings. These are the same
  six NOTEs as integration/all (future file timestamps, Rd line widths, Rd cross-references, compilation flags,
  examples timing, and CRAN incoming: no internet + `Remotes` field). None comes from this change.
