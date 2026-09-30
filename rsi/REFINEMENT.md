# rsi — refinement walk (SINEderella code + bat SINE understanding)

A running log of his review of the *Rhinolophus sinicus* (rsi) page, subfamily by subfamily from the
bottom of the list, and what each observation led to. Each entry: **his observation → what was
measured → cause → proposal (for SINEderella) → status.** Nothing marked *proposed* is implemented yet.

Run: therioserver `~/rhin/rsi/run_add_20260927_180847` (bank `rsi_peel10_bank.fa + MEG.resolved.fa`),
republished 2026-09-28 on the fast step8a (SINEderella 332f9fe).

---

## 0. Page defects (all 25 bat pages) — fixed 2026-09-28

- **"Genome: ? · Consensus: ?"** in every page header. Cause: `--add` runs write only `SOURCE_RUN` /
  `ADD_FILE` to `manifest.txt`; step6 read `GENOME_IN` / `CONS_IN`. Fixed in SINEderella `6aba0c1`
  (`run_inputs()` walks back through `SOURCE_RUN`); existing pages filled from the runs.
- **Navigation and the Tribes section missing, no species name.** The republish regenerates
  `report.html`, dropping the blocks Tal injects. Re-injected; the nav block now also carries
  species · family · accession (from `chiroptera/genomes.tsv`). Run `chiroptera/after_pull.sh` after
  every pull of republished reports.
- The header fix first did not reach the pages: `publish_run.sh` ran a copy of step6 frozen in each
  run dir. It now always publishes with the installed step6 (SINEderella `d730836`). All 25 pages
  rebuilt 2026-09-28 21:33 MSK (Tal `aa2d292`): seed 42 shown, plain-language comments, headers filled.

## 1. r9_15seqs — the reference case

**His call:** a textbook SINE and presentation: all three alignments clean, subfam clean and
homogeneous, descriptions correct. One quirk: the alignment edges are rugged although every copy has
enough flank to cut evenly.

**Measured** (`rsi_r9_15seqs_top100.aln.fa`, 100 copies): 5′ edge already even (rows start at
columns 0–3). 3′ edge ragged: every copy has ≥ 118 bp of 3′ flank, but the flank length runs
118–183 bp (≥ 10 % of rows at 160–183), so rows end anywhere between column 323 and 388.

**Cause:** the flank is counted from where each copy's search hit ends, which varies with the length
of its A tail (e.g. `AAAAAAAaaaaaaaaaaaaa…` vs `AAAAAAAaacaacaac…`). `trim_display_flanks.py` keeps a
column while ≥ 10 rows still have a base there, so the longest-flank tenth sets the right edge.

**Proposal (P1, display only):** cut both flanks at a common column — the edge where ≥ 90 % of the
rows still have flank (here 3′ ≈ 120 bp past the element) — so every row ends together. Rows that
genuinely run out earlier (contig end, tandem, see §2b) stay short and stand out. *Proposed.*

## 2. r8_83seqs — good, with two findings

### 2a. A distinct set inside r8: "group B" (521 copies, 11 %)

**His call:** the subfam alignment contains a set of somewhat different copies; indicate it and
analyse it separately, compared with other known groups.

**Measured.** In `rsi_r8_83seqs_subfam.aln.fa` (92 chunk consensuses of 50 copies), 12 chunk rows —
`input_081` … `input_092`, contiguous in SubFam's similarity order — have 71–76 % identity to the
r8 row, against a median of 99 % for the rest; they are also shorter (113–151 bp vs 181).
Their majority consensus (123 bp, `analysis/rsi_r8_groupB_consensus.fa`) scored against all 4 591
r8-assigned copies (ssearch36 `-z 11`; default statistics collapse on an all-homolog library):
**521 copies match group B better by > 10 bits** — a clean split (median copy favours r8 by 117
bits, top 5 % favour group B by 66). Members: `analysis/rsi_r8_groupB_members.txt`; locus length
p10 122 / median 137 / p90 174 bp.

**What group B is** (60 best copies with 50 bp / 100 bp flanks + both consensuses,
`analysis/rsi_r8_groupB_top60.aln.fa`):
- **Shared fixed differences from r8** along the 5′ part, e.g. `CCGACTGCCTCAGCC` vs r8
  `CCGGTTGGCTCAGTT`, `GGAGCGCAAGGCTCATAATACCA` vs `AGAGCGCAGTGCTCATAACACCA` — the same in all 60
  copies: a diagnostic pattern, not decay (MANUAL §6.1.6 subfamily definition).
- **A 20 bp deletion that removes the B box** (`GGTCGCCGGTTCGATTCCCA` in r8), in every copy.
- **A different 3′ end:** after `…CGGCTTGAACCGGAGTG` the copies end in a G run and a G/A (purine)
  tail of ~40–60 bp; r8's 3′ region `GATGGGTCCTGG…CCCCAATAAAAAAAA` is absent.

**Compared with known groups:** its best match in the rsi bank (r1–r10, MEG, Rhin-1) and among the
rebuilt consensuses of all 25 bats (`chiroptera/recreated/recreated.fa`) is the 5′ 136 bp of
r8 / Rhin-1 at only 77 % (21 gaps); r5, r7, r10 match 71–77 bp at 73–79 %. **No known consensus
represents it.** Open for him: a separate subfamily / family (B-box-less, purine-tailed), or a
known structure in Rhinolophidae.

**Proposal (P2):** in the subfam plate, flag contiguous runs of chunk rows whose identity to the
original is well below the rest (here < 0.8 vs median 0.99) as a *variant block*; on the page show
the block and its size; build its consensus, count members by the score split above, write its own
plate, and compare it with the bank and the rebuilt consensuses — i.e. automate this section.
*Proposed.*

### 2b. Four rand100 copies much shorter at the 3′ end — tandem loci

**His call:** in the 100 random copies at least 4 are much shorter at the right end — why?

**Measured** (`rsi_r8_83seqs_rand100.aln.fa`): exactly 4 rows hold 201–239 letters (median 402) and
cover only 86–124 of the 181 element positions, stopping 151–395 columns before r8's 3′ end.
top100 and r9 rand100 have none. Every consensus hit in ±300 bp of each locus, in genome order:

| locus | hits along the locus |
|---|---|
| NC_142501.1:113316563-113316928(−) | r8 1–132, then r8 1–132 again 14 bp later |
| NC_142504.1:30250275-30250638(−) | r5 1–127, r10 6–100, r8 1–153, r6 106–224 |
| NC_142508.1:31735777-31736083(+) | r10 1–79, then r6 1–179 36 bp later, r6 134–225 |
| NC_142513.1:31101796-31102151(+) | r5 1–125, then r6 1–225 19 bp later |
| (normal row) NC_142517.1:96342939-96343335(+) | r8 2–175 only |

**Cause:** each is a **head-to-tail tandem**: a 3′-truncated first copy (cut at consensus 79–132)
immediately followed by the head of another copy. SINEderella treats the whole tandem as one r8
copy. In the plate only the truncated first copy aligns with r8; the second copy's head cannot align
to r8's 3′ half, so MAFFT pushes it into its own gap block — the row looks short, and the whole
rand100 spreads out (the 181 bp original occupies 530 columns vs 186 in top100). top100 is free of
them because single full copies score higher.

**Proposal (P3):** mark such rows `[tandem]` when a second consensus hit starts inside the locus
(the same machinery as `[array]`); leave them out of the continuation decision and put them last,
or split the locus at the junction. *Proposed.*

## 3. r8 group B, second pass — a candidate new group; the consensus lacks its left part

**His call:** a candidate for a new group, but it needs more work: the copies share their left
flank, and the consensus is not properly extended to the left.

**Measured** (60 best group-B copies, 600 bp upstream, column-majority agreement read outward):
0.78–0.90 for **~300 bp to the left** of the current consensus start, then 0.40–0.50 (unrelated
flank) from ~330 bp. So the shared "left flank" is part of the element, and there is proper unique
flank beyond it, well inside step8a's +600 bp limit.

**Left-extended consensus** (`analysis/rsi_r8_groupB_extended_consensus.fa`, 391 bp, majority of the
60 copies with 310 bp upstream):

| positions | what it is |
|---|---|
| 8–146 | **r10** (r10_19seqs 1–145, 99.3 % identity) |
| 147–251 | ~105 bp GC-rich middle — no hit in the bank |
| 252–374 | group B (the 123 bp consensus, 100 %) |
| 375–391 | G/A (purine) tail |

**r10 and group B are one element:** 85 % of r10 copies (3 198 of 3 761 checked) have group B
0–250 bp downstream (mostly 50–150 bp, matching the ~105 bp middle); 76 % of group-B copies
(395 / 521) have r10 upstream. The bank holds this ~390 bp element as two pieces: r10 (its left
part) and group B (its right part, hidden inside r8). *Next:* rebuild it as one consensus, re-run
assignment with it, and see what r10 and r8 become.

## 4. r7 — almost perfect (his call). Nothing to fix.

## 5. r6 — left-side extensions in subfam and top100: tandems with an r5 head

**His call:** subfam contains individual left-side extensions; in top100 the right flank is fine,
but several copies extend further left — tandems? some are nearly identical in flanks. Do these
extended copies exceed the discriminator's maximum length, and if not, do they have proper flank
further away?

**Measured** (`rsi_r6_210seqs_top100.aln.fa`): 22 of 100 copies carry ~236 bp before the element
(the rest 101 bp); their loci start ~135 bp further upstream, and that extra stretch agrees between
copies (0.64–0.75, vs ~0.45 for unrelated flank).
- **What it is:** in all 22 the extra segment is an **r5 head, r5 positions 1–~127** (85–92 %). So
  these are head-to-tail tandems: a 5′ r5 fragment + r6 — the "nearly identical flanks" are SINE
  sequence. (The new Flank context text on the page says it: top 100 left side, 34 of 100 copies
  share flanking DNA, largest group 34.)
- **Length:** the shared part is ~135 bp — far inside step8a's +600 bp continuation limit (the only
  length cap in the pipeline; the discriminator has no separate maximum element length).
- **Proper flank beyond:** yes. Over the next 800 bp upstream agreement is 0.42–0.50 (unrelated) and
  no two far flanks are ≥ 0.8 identical (highest pair 0.59) — independent loci, not segmental
  duplicates.

## 6. r5 — two subgroups; the second is the left half of an r5–r5 dimer

**His call:** even from subfam, r5 has at least two subgroups; what is the other one?

**Measured** (`rsi_r5_27seqs_subfam.aln.fa`): ~40 chunk rows are full length (~178 bp, 0.98–0.99 to
r5); ~36 rows (`input_040`–`074`) are **132–134 bp at 0.83–0.91**. Their consensus
(`r5_subgroupB`, 133 bp) is r5's first ~110 bp followed by a different end,
`…GATTGAAGACAATGAGCTGCCGCTGAGCTTCCG` (r5: `…GATTGAAGGACAACGACTTGGAGCTGATGG…`).
Split of the 4 062 r5 members by score (> 10 bits): **subgroup B 1 546**, full-length 2 045,
unclear 471. Directly downstream (≤ 300 bp) of 200 random copies of each:
- subgroup B: **182 / 200 have another SINE starting at 0–25 bp**, 160 of them r5 from position 1;
- full-length r5: 1 / 200.

So subgroup B is a **5′ unit (~133 bp, own 3′ end) immediately followed by another copy**, which
SINEderella splits into two loci. **Correction (§11):** the partner is mostly **r6**, not r5 — this first
test took the first hit by position, and r5 and r6 share their 5′ part; scored by best match, r5
left halves are followed by r6 (827), r5 (146) or r3 (72). The r6 tandems (§5) start
with an r5 head that is slightly closer to this subgroup than to r5 (91.2 vs 89.3 %) — the same left
unit may pair with several right units; not yet shown.

**Proposal (P3, now broader):** a locus followed within ~30 bp by the head of another copy is part
of a compound element (tandem or dimer). Mark such rows (`[tandem]`), and report on the page how
many copies of each subfamily are such left or right halves; a subfamily that is mostly left halves
(r5 subgroup B, r10) should be rebuilt together with its partner.

## 7. Verdict comments rewritten (Flank context, Overall) — done 2026-09-28

One fixed layout, no internal codes (DISC `8124186`):
- **Flank context:** "Are the copies independent insertions? Checked by comparing the DNA just outside
  each copy." · per plate and side: "left side: 34 of 100 copies share flanking DNA with another copy
  (largest group: 34 copies)" or "all 100 copies have their own flanking DNA" · one-sentence reading.
- **Overall:** "Score 96 of 100: SINE." · "The score is held down because: …" (was
  "capped: SMALL_CORE") · "Evidence against:" · "Also noted (does not change the call): …".
- Fixed on the way: a "group" of one copy counted as shared when few copies were measured
  (cse MEG-TR rand100: 6 unique copies was flagged "Shared context").

## 8. The three questions from the first pass

1. **`continuation.tsv` and `[array]` rows.** step8a asks "do the copies still agree where their
   flank runs out? if so, extract more" — but only over copies that are *not* tandem-array units,
   because array units share flank for their whole length and would make every plate extend to the
   cap. The published `continuation.tsv` / page note re-measures the same thing **over all rows,
   array units included**. On array-heavy plates (cse MEG-RS rand100: 55 of 100 rows are array
   units) the array units alone make it say "unresolved", although step8a, looking at the
   independent copies, correctly decided there was nothing to extend. Proposal: measure the
   published status on independent copies too, and say how many array rows were set aside.
   *Awaiting his decision.*
2. **Subfam plates in `continuation.tsv` — done** (his call: leave them out; DISC `02370dc`, 84 rows
   removed from the 25 published tsv files). The subfam plate is not copies with flanks: each row is
   the consensus of a chunk of 50 copies, element only. "Do copies stay similar past the element's
   end" has no meaning there — its rows simply stop where the element stops, so it reads as
   "unresolved, 0 bp". Proposal: leave subfam plates out of `continuation.tsv`. *Awaiting his
   decision.*
3. **rand100 seed — done** (his call: fix and indicate). step8a draws rand100 with a fixed seed
   (`RAND_SEED`, default 42; SINEderella `785895d`); the page column reads "100 random copies
   (seed 42)" and the button's rule says the same copies come back on every rebuild.

## 9. Array rows and the published continuation status — concrete examples (for his decision)

What is compared, for each published top100 / rand100 plate that has `[array]` rows (all 25 bats):
the published status (all rows) against the same measurement with the `[array]` rows set aside
(what step8a uses to decide on extension). The status is read at the element's edge: `none` = the
copies stop agreeing right at the edge; `ends N` = they keep agreeing for N bp, then stop while
≥ 50 % of copies still have sequence; `unresolved` = the copies run out of sequence (cover < 50 %)
while still agreeing. **18 plate sides change** (of the plates with array rows); in bp the
differences are small (0–25 bp) — what flips is mostly the *cover* at the edge, i.e. how many rows
still have sequence there.

| plate | side | array rows | all rows (published) | independent copies only |
|---|---|---|---|---|
| tbr MEG-RS rand100 | 3′ | 73 / 100 | unresolved 10 bp, cover 0.49 | none, cover 0.74 |
| nth MEG-RS rand100 | 5′ | 34 / 100 | unresolved 4 bp, cover 0.39 | none, cover 0.79 |
| tni MEG-RS rand100 | 5′ | 63 / 100 | unresolved 1 bp, cover 0.43 | none, cover 0.70 |
| mgi MEG-TR top100 | 3′ | 14 / 39 | unresolved 6 bp, cover 0.41 | none, cover 0.76 |
| fho MEG-RS top100 | 5′ | 41 / 100 | unresolved 1 bp, cover 0.49 | none, cover 0.15 |
| nth MEG-TR top100 | 5′ + 3′ | 35 / 75 | unresolved 0 bp | none |
| rsi MEG-TR rand100 | 3′ | 76 / 100 | unresolved 0 bp, cover 0.07 | none, cover 0.00 |
| vmu MEG-TR top100 | 3′ | 18 / 90 | unresolved 10 bp, cover 0.48 | none, cover 0.60 |
| mev MEG-RS rand100 | 5′ | 42 / 100 | ends 8 bp | none |
| rmi MEG-RS top100 | 3′ | 33 / 100 | ends 4 bp | none |
| vmu MEG-RS rand100 | 5′ | 74 / 100 | unresolved 0 bp, cover 0.24 | **unresolved 25 bp**, cover 0.46 |
| vmu MEG-RS rand100 | 3′ | 74 / 100 | unresolved 0 bp, cover 0.38 | none, cover 0.19 |
| rsi MEG-TR top100 | 5′ | 73 / 100 | none, cover 0.03 | unresolved 0 bp, cover 0.11 |
| hla MEG-TR rand100 | 5′ | 72 / 100 | none, cover 0.04 | unresolved 0 bp, cover 0.14 |
| rmi MEG-TR rand100 | 3′ | 49 / 87 | none, cover 0.24 | unresolved 0 bp, cover 0.37 |
| lly MEG-RS rand100 | 3′ | 3 / 100 | none, cover 0.44 | unresolved 0 bp, cover 0.42 |
| nle MEG-RS rand100 | 3′ | 3 / 100 | none, cover 0.24 | unresolved 0 bp, cover 0.24 |

Three kinds, one example each to look at (rows marked `[array]` in the viewer):
1. **Array rows create the "unresolved"** — tbr MEG-RS rand100 3′: with the 73 array units in, only
   49 % of rows still have sequence at the edge, so "the copies ran out while still agreeing" is
   triggered; the 27 independent copies still have sequence there (74 %) and stop agreeing right at
   the edge (`none`). Why the array units run out there is not yet measured.
   <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Ftbr%2Falignments%2Ftbr_MEG-RS_rand100.aln.fa&title=tbr_MEG-RS_rand100>
   (also nth MEG-RS rand100 5′: <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fnth%2Falignments%2Fnth_MEG-RS_rand100.aln.fa&title=nth_MEG-RS_rand100>)
2. **Array rows hide a real extension** — vmu MEG-RS rand100 5′: the independent copies agree for
   25 bp past the edge; with the 74 array rows mixed in, it reads 0 bp.
   <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Fvmu%2Falignments%2Fvmu_MEG-RS_rand100.aln.fa&title=vmu_MEG-RS_rand100>
3. **Too few independent copies left** — rsi MEG-TR top100 5′: 73 of 100 rows are array units;
   the 27 left barely reach the edge (cover 0.11), so either reading is weak.
   <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_MEG-TR_top100.aln.fa&title=rsi_MEG-TR_top100>

All 18 are MEG-RS / MEG-TR plates (the families with tandem arrays). *Awaiting his decision.*

## 10. r10 + group B rebuilt as one consensus and reassigned (2026-09-28)

`r10_groupB` (391 bp, `analysis/rsi_r8_groupB_extended_consensus.fa`) added to the rsi run with
`SINEderella --add` (incremental: an old 10/10 call is re-voted only where the new consensus hits).
New run: therioserver `~/rhin/rsi_gB/run_add_20260928_214348` (not published; the rsi page is unchanged).

| subfamily | firm before | firm after |
|---|---|---|
| r10_groupB | — | **4 181** (+669 soft) |
| r10_19seqs | 5 188 | **1 072** |
| r8_83seqs | 3 951 | 3 992 |
| r5_27seqs | 3 940 | 3 472 |
| r1_9seqs | 3 951 | 3 573 |
| others | | within ±1 % |

- ~4 100 of the 5 188 "r10" copies are the compound element: r10 was its left part.
- Of r8's 521 group-B copies, **505 now lie in an r10_groupB locus**, 11 stay r8, 5 elsewhere.
- r8's total nevertheless stays ~4 000: it lost the group-B copies and gained about as many from
  elsewhere — *not yet traced*.
- r5 (−468) and r1 (−378) also gave copies to r10_groupB — *not yet examined* (r5's dimer halves?).
- r10_groupB mean sim_ratio 0.59 (vs 0.74 for r10 before): many copies are shorter than 391 bp or
  diverged; its plates are the next thing to look at.

## 11. Composite copies, tested genome-wide (2026-09-29)

**Tool:** SINEderella `tools/composite_scan.py` (analysis only, not in the pipeline). For every
assigned copy of the run (66 501 in `~/rhin/rsi_gB/run_add_20260928_214348`, i.e. with
r10_groupB), the locus ± 400 bp is searched with every consensus (ssearch36, both strands);
non-overlapping hits are kept best-first and written out in genome order as a **chain**, each unit
tagged `h` (starts at its consensus 5′ end), `e` (reaches its 3′ end) or `mid`; links `+` ≤ 30 bp,
`~` 30–200 bp, `..` further. A `mid` unit right after another is one element matched piecewise by
two consensuses, not a new copy. Random expectation of a neighbouring copy within 30 bp on a given
side: 0.4–1.3 % (copy density × element length). Validated on the cases found by hand (r5 left
halves, r6 with an r5 head, r9 clean). Outputs: `analysis/composite_chains.tsv`,
`composite_summary.tsv`, `composite_loci.tsv.gz`, `composite_scan_report.txt`.

**Most common layouts** (share of that subfamily's copies):

| subfamily | layout | share |
|---|---|---|
| r9 | `r9[he]` single | 91 % |
| r10_groupB | `r10_groupB[he]` single | 83 % |
| r7 | single | 99 % |
| **r6** (13 375) | `r6[he]` single | 51 % |
| | **`r5[h] + r6[he]`** — r5 head (cut at ~130) + full r6 | **33 %** |
| **r5** (3 472) | `r5[he]` single | 48 % |
| | `r5[h] + r6[he]` (the same structure, seen from the r5 side) | 24 % |
| | `r5[h] + r5[he]` | 4 % |
| **r3** (22 352) | **`r1[he] ~ r3[he]`** — full r1, spacer, full r3 | **31.5 %** |
| | `r1[he] ~ r3[h] + r3[e]` — same, r3 with an internal repeat | 15.5 % |
| | `r10[h] + r3[e]` — r3 whose head matches r10 better (piecewise) | 9.5 % |
| | `r5[h] + r3[he]` | 4.6 % |
| | `r1[he] ~ r2[he] + r3[he]` (+ variants) | 6 % |
| | `r3[he]` alone | only 2.3 % |
| **r1** (3 573) | `r1[he] ~ r3…` | ~50 % |
| | `r1[he]` alone | 22 % |
| **r2** (386) | `r1[he] ~ r2[he] + r3…` | ~50 % |

**The r1–r3 spacer is fixed:** 13 848 r3 copies have r1 30–200 bp upstream; the gap is 39 bp at the
median (p25 39, p75 42; 8 551 in 30–39 bp) — a linker, not chance proximity. A second mode at
170–180 bp is the variant with r2 between them.

**What this says about rsi's families:**
1. **r1 + 39 bp + r3 (~395 bp)** is the main element behind r3 (≥ 47 % of r3 copies, 13.8 k with r1
   nearby) — r3 as searched is its right part, r1 its left part; a variant carries r2 in between.
2. **r5 head (~130 bp) + r6 (~360 bp)** is ~⅓ of r6; the ~130 bp cut is fixed (a left unit, like the
   r8 / r5 cuts at 128–132 seen before).
3. **r10 + 105 bp + group B** — already rebuilt as r10_groupB (§10); 83 % of its copies now read as
   one full unit, so that rebuild worked.
4. **r3 carries an internal repeat** in ~20 % of copies (`r3[h] + r3[e]`: r3 1–163 then 115–201,
   ~50 bp duplicated) that the consensus lacks.
5. Clean single families: r9, r7, r10_groupB (and r8 at 89 %).

**Next (proposed):** rebuild r1+linker+r3 and r5head+r6 as full consensuses the way r10_groupB was
built (60 best copies of the layout, majority), add them, re-run assignment, re-scan; wire a short
form of the scan into the report (per subfamily: % single, top layouts) and mark composite rows
on the plates.

## 12. r1_r3 and r5h_r6 rebuilt and reassigned (2026-09-29)

Built with `tools/build_composite.py` (60 best copies of one exact layout, cut first unit to last
unit, L-INS-i, majority): **r1_r3** 396 bp from `r1[he] ~ r3[he]` (copies 388/396/397 bp,
p10/median/p90); **r5h_r6** 360 bp from `r5[h] + r6[he]` (359/360/361 bp). Added to the
r10_groupB run: therioserver `~/rhin/rsi_comp/run_add_20260928_222535` (not published).

| subfamily | firm before | firm after | read as one full unit after |
|---|---|---|---|
| r1_r3 | — | **18 329** | 49.5 % (`r1_r3[he]`); 22.5 % `r1_r3[h] + r3[e]` |
| r5h_r6 | — | **7 294** | **84.8 %** |
| r3_58seqs | 22 352 | 4 985 | mixed (r1_r3 pieces, r10 + r3) |
| r6_210seqs | 13 375 | 7 734 | 86 % (98.5 % single) |
| r5_27seqs | 3 472 | 2 220 | 98 % single |
| r1_9seqs | 3 573 | 1 674 | mostly r1_r3 pieces |
| r2_3seqs | 386 | 44 | |

- **r5h_r6 works:** one full unit for 85 % of its copies, and what is left of r6 and r5 is clean.
- **r1_r3 half works:** a quarter of its copies have an r3 part that r1_r3 does not represent —
  the r3 variant with the ~50 bp internal repeat (§11 point 4). *Next:* a second consensus from the
  `r1[he] ~ r3[h] + r3[e]` layout.

## 13. What the parts derive from — tRNA, 5S, 7SL? (2026-09-29)

Queries: every consensus of the current run plus the parts of the three rebuilt elements (part
boundaries from the piece-by-piece matches, ± a few bp). Three independent searches:
**Infernal cmscan against all of Rfam** (8 454 models: tRNA, 5S, 7SL, U6, 7SK, vault, Y RNA …;
E ≤ 1e-3, relaxed to E ≤ 1 for the parts without a hit), **tRNAscan-SE -E**, and **Dfam**
(transposon families incl. tRNA-related SINE heads; curated, `--cut_ga`, taxa *Homo sapiens* and
*Myotis lucifugus* — identical results; "Chiroptera" is rejected by the service). Files: therioserver
`~/refs/smallrna/q/` (`cm.tbl`, `trnascan.out`, `dfam_*.json`).

**Controls behave as their names say** (Gogolevsky et al. 2009: R = rRNA-related, T = tRNA-related):
MEG-RL, MEG-RS → **5S rRNA** (1–118; Dfam E ≈ 1e-36, Rfam 1e-24); MEG-T2 → **tRNA-Val** (Dfam
tRNA-Val-GTA, tRNAscan Val, not pseudo); MEG-TR → **tRNA-Val (1–56) + 5S (65–179)**.

| query | Dfam (human / *Myotis*, same) | Rfam | tRNAscan-SE |
|---|---|---|---|
| r9 | tRNA-Ile 5–79, E 1e-23 | tRNA, E 3e-11 | Ile-AAT 65.4, **not pseudo** |
| r7, r8 | tRNA-Ile 5–79, E 1e-18 – 1e-19 | tRNA | Lys / Met, pseudo |
| r5, r6 | tRNA-Ile 5–77, E 1e-10 – 1e-11 | — | Ser / Leu, pseudo |
| r10 | tRNA-Ile 2–74 (8e-11) **and MamSINE1** 8–104 (1e-09) | tRNA | Arg, pseudo |
| r1, r2, r3, r4 | tRNA-Ile 2–75, E 1e-6 – 1e-7 | — | r1 Ser (pseudo); r2–r4 none |
| **r1_r3** | **tRNA-Ile at 2–75 and at 195–262** | — | 1–75 only |
| **r5h_r6** | **tRNA-Ile at 5–77 and at 141–210** | — | 5–77 and 139–223 |
| r10_groupB | MamSINE1 15–106 (the r10 part) | — | 8–77 (Undet) |
| linker (r1_r3 157–195), middle (r10_groupB 147–251), group-B part | **nothing** in any of the three | | |

- **Every rsi head is tRNA-derived** (Dfam's model and tRNAscan's isotype model both say Ile; the
  anticodon calls of degenerate heads — Ser, Lys, Met, Arg, Leu — are not reliable; r9 is the one
  clean tRNA(Ile)). **No 7SL, 5S or other small-RNA part in any rsi r-family** (the 5S in the bank
  is only in MEG-R*/MEG-TR).
- **Correction to the first version of this section:** r2, r3 and r4 *do* have a tRNA head — Rfam and
  tRNAscan missed it, Dfam's more sensitive tRNA model finds it. So **r1_r3 and r5h_r6 are both
  tRNA–tRNA dimers**: two tRNA-derived units head to tail (r1_r3 with a 39 bp linker, r5h_r6 with
  the first unit cut at ~130). Whether r3 and r6 also exist as independent monomers is what the
  residual r3 (4 985) and r6 (7 734, 98.5 % single) copies show — r6 does, clearly.
- **r10's head resembles MamSINE1** (Dfam), a family not in the rsi bank; group B (no B box) and the
  105 bp middle match nothing known.

## 14. Round 2 of manual inspection — page with the rebuilt consensuses (2026-09-29)

<https://toki-bio.github.io/Tal/rsi_v2/report.html> — the run with r10_groupB, r1_r3 and r5h_r6 added
(therioserver `~/rhin/rsi_comp/run_add_20260928_222535`), published on the current code (seed 42,
plain comments). The original page stays at `rsi/`. Not yet on it: composite marks on plate rows,
the r1_r3 variant with r3's internal repeat.

## 15. Composite candidates found automatically by flankscan (2026-09-30)

flankscan stages 1-6 (SINEderella/flankscan) on the original run (`~/rhin/rsi/run_add_20260927_180847`,
no hand rebuilds in the bank), output therioserver `~/tmp/fs_rsi4`: 50 junction peaks -> 14 distinct
candidate consensuses (`rsi_v2/composites/`: `*.aln.fa` = 60 copies of the layout, `*.fa` = candidate,
`candidates.tsv`, `reassign.tsv`). Stage 6 adds all candidates once and re-assigns; "accept" = >= 70 % of the
peak's elements now read as one full unit of the candidate.

| candidate | type | peak copies | full-unit share | verdict | = hand result |
|---|---|---|---|---|---|
| r1__r3_P1 | composite (39 bp linker) | 9 158 | 71.3 % | accept | r1_r3 (99.75 %) |
| r5__r6_P26 | composite (r5 cut ~127) | 3 518 | 89.8 % | accept | r5h_r6 (100 %) |
| r10__r8_P18 | piecewise (r10 + group B) | 90 | 98.3 % | accept | ~r10_groupB (97.2 %, ~30 bp indel) |
| r1__r3_P2 | piecewise (r1 1-79 + r3 from 29) | 3 427 | 91.9 % | accept | prototype "r10[h]+r3[e]" |
| r1__r3_P21 | piecewise | 1 454 | 76.4 % | accept | new |
| r2__r3_P6 | composite (gap 0) | 5 725 | 19.4 % | check | r1~r2+r3 variant |
| r1__r2_P40 | composite (gap 31) | 236 | 61.0 % | check | r1~r2 part |
| r3__r3_P13 | homodimer | 270 | 53.7 % | check | |
| r8__r8_P43 | homodimer | 139 | 54.5 % | check | |
| r3__r3_P9 | piecewise (r3 1-161 + 114-201) | 4 300 | 2.1 % | check | r3 internal repeat; its copies go to r1_r3 |
| r1__r2_P39 | piecewise | 387 | 2.0 % | check | |
| r10__r6_P41 | piecewise | 101 | 7.9 % | check | |
| r7__r3_P20 | piecewise | 329 | 0 % | check | artifact: r7's real 3' tail TAAATAA(A)TAAAAGTT + A run beyond the r7 consensus end (boundary note for r7) |
| r9__r8_P33 | composite | 55 | 0 % | check | |

Open: the 3-part element r1 + 39 bp + r3-with-internal-repeat needs a second round (stages 3-4 again with
the kept candidates in the bank). All calls on these alignments are his.

**Update 2026-09-30 (stage 5 fixed):** candidate consensuses are now built from copies cut with 100 bp flanks and
extended while >= 60 % of copies agree (the element cut had stopped short in r3's simple-repeat tail and before
GGGCC at the r1 5' end). Rerun (therioserver `~/tmp/fs_rsi5`, Tal `rsi_v2/composites/v2/`): 17 candidates,
7 accept - P1 r1+r3 82.3 % (was 71.3), P26 87.1, P2 89.7, P34 r5+r3 83.6 (new), P12 r3 homodimer 89.0 (new),
P18 96.6; P21 fell to 63.3 and P40 to 9.0 (check). P1 = 431 bp.

**r8 dimers (2026-09-30):** the "10 % r8 homodimers" are ONE element seen through two junction splits:
P42 r8(1-132) + 9 bp + r8 and P43 r8(1-83) + 35 bp + r8 are 99.0 % identical over 325 bp; P38 (r8 head + unit
assigned r6) = P43 100 %. The element (P43, 300 bp) is ~88 % to r5h_r6 - related, not the same. Stage 5 folded
P42 into P26 (r5h_r6, 90.6 %) because P26 came from the larger peak; the fold should go to the MOST similar
candidate (P43, 99 %). Stage 6: P43 54.5 % one full unit (check). Alignments: rsi_v2/composites/v2/cand/ P43,
P42, P37, P38.

## 16. Singletons: do r1, r2, r3, r10 exist alone? (2026-09-30, flankscan stage 8)

Per family three groups aligned with 250 bp flanks (S = singles, B = copies in composites, A = r9 singles as
the clean control); ends from column conservation (ViewAlign auto mode), one TSD per copy at those ends
(ViewAlign detector, 3' slack 25 bp, minimum calibrated on shuffled pairs). Alignments: Tal rsi_v3/singletons/.
- Controls: r9 singles end at the unit, TSD 19-42 % vs 3.5-5 % chance -> rsi SINEs make TSDs. Composite members'
  ends move to the partner by themselves (r3 5' +197 = r1 + 39 bp; r1 3' +225; r2 +192 / +223; r10 3' +232).
- r3 singles (529): end at the unit, 3' end ~23 bp beyond the consensus + tail; TSD 7.5 % (chance 3 %, r9 19.5 %).
- r1 singles (142): end at the unit; TSD 15.5 % (chance 3.5 %, r9 30.5 %); 56 % full-length; relaxed partner 3' 30 %.
- r2: 3 singles only, and they are composites -> r2 never occurs alone.
- r10 "singles": 3' end +217 like the composites; partner 3' 82 % -> they are r10 + group B.
Open (his call on the alignments): are the r3 and r1 singles partly standalone (TSD above chance, below r9)?

## Viewer links

- group B, 60 best copies: <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Fanalysis%2Frsi_r8_groupB_top60.aln.fa&title=rsi%20r8%20group%20B%20top60>
- r8 subfam (chunk rows 081–092 are group B): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r8_83seqs_subfam.aln.fa&title=rsi%20r8%20subfam>
- r8 rand100 (the 4 tandem rows): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r8_83seqs_rand100.aln.fa&title=rsi%20r8%20rand100>
- r9 top100 (ragged 3′ edge): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r9_15seqs_top100.aln.fa&title=rsi%20r9%20top100>

Analysis scripts: therioserver `~/tmp/r8short/` (tandem loci), `~/tmp/r8grp/` (group B).
