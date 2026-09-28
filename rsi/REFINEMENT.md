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

## Open questions from the republish (2026-09-28)

1. `continuation.tsv` counts `[array]` rows, which step8a's extension decision skips — array-heavy
   plates (cse MEG-RS rand100, 55 % array) are reported `unresolved`. Count independent copies only?
2. Subfam plates (chunk consensuses, no flanks) get meaningless `unresolved 0 bp` entries — drop them?
3. rand100 is drawn unseeded (`shuf`), so it changes every republish (cse MEG-RS rand100 3′ went
   `ends 647` → `unresolved 652`). Fix the seed?

## Viewer links

- group B, 60 best copies: <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Fanalysis%2Frsi_r8_groupB_top60.aln.fa&title=rsi%20r8%20group%20B%20top60>
- r8 subfam (chunk rows 081–092 are group B): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r8_83seqs_subfam.aln.fa&title=rsi%20r8%20subfam>
- r8 rand100 (the 4 tandem rows): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r8_83seqs_rand100.aln.fa&title=rsi%20r8%20rand100>
- r9 top100 (ragged 3′ edge): <https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2Fmain%2Frsi%2Falignments%2Frsi_r9_15seqs_top100.aln.fa&title=rsi%20r9%20top100>

Analysis scripts: therioserver `~/tmp/r8short/` (tandem loci), `~/tmp/r8grp/` (group B).
