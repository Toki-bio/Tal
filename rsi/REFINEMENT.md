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

So subgroup B is a **5′ unit (~133 bp, own 3′ end) immediately followed by a full r5**: the left half
of a head-to-tail r5–r5 dimer, which SINEderella splits into two loci. The r6 tandems (§5) start
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

## Viewer links

- group B, 60 best copies: <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Fanalysis%2Frsi_r8_groupB_top60.aln.fa&title=rsi%20r8%20group%20B%20top60>
- r8 subfam (chunk rows 081–092 are group B): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r8_83seqs_subfam.aln.fa&title=rsi%20r8%20subfam>
- r8 rand100 (the 4 tandem rows): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r8_83seqs_rand100.aln.fa&title=rsi%20r8%20rand100>
- r9 top100 (ragged 3′ edge): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r9_15seqs_top100.aln.fa&title=rsi%20r9%20top100>

Analysis scripts: therioserver `~/tmp/r8short/` (tandem loci), `~/tmp/r8grp/` (group B).
