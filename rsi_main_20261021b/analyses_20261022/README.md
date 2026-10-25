# rsi_main_20261021b analyses (2026-10-22)

Alignments asked for on 2026-10-22, made from the rsi run `rsi_main_20261021b` (therioserver
`/home/toki/rhin/rsi_main_20261021b/run_20261021_192808`; genome *Rhinolophus sinicus* GCF_057876585.1, RefSeq annotation
GCF_057876585.1_RSI-T2T). Work directory: therioserver `/home/toki/tmp/rsi_aln_20261022/`. Scripts in `scripts/` and
`trna_genes_downstream/`. MAFFT 7.526; ssearch36 36.3.8g; tRNAscan-SE 2.0.13.

## `trna_heads/`: tRNA heads with the best and the close second tRNA genes (task C)

Source: the run's tRNA-origin table (`step2_output/trna/`, tools/trna_heads.sh: ssearch36 -3 of every consensus and
composite stretch against GtRNAdb hg19, 419 genes; 20 dinucleotide shuffles as the null). A query has a tRNA head when its
best tRNA beats every shuffle (the report's rule). The report itself keeps only one runner-up, the best gene of a *different
amino acid* (`2nd_isotype`), and has no "close" rule. Here **close = bits >= 90 % of the best**, counted over distinct gene
sequences (identical paralogues collapsed), up to 3 per query, other anticodons and amino acids first. `trna_heads_table.tsv`
lists, per query, the best gene, the close genes with their share of the best score, and the best other amino acid.

- `rsi_trna_heads_Ile.aln.fa`: 31 Ile-headed queries (r-units, composites and composite parts; head span + 30 bp) and 15
  human genes (best and close). Ile-AAT and Ile-GAT genes are within 90-100 % of each other for every r-head, so the
  isotype within Ile is not resolved by this test.
- `rsi_trna_heads_Val.aln.fa`: MEG-T2, MEG-TR and three Val-TAC genes.
- `rsi_trna_heads_Lys.aln.fa`: the r8 part of r10_r8_P18 (weak, 26.6 bits against a 24.7-bit shuffle).
- MEG-RL, MEG-RS and five composite stretches have no tRNA head (best tRNA <= the shuffles).
MAFFT `--localpair --maxiterate 1000 --threadit 0`, input order. Script `scripts/c_heads.py`.

## `trna_genes_downstream/`: why the tRNA similarity runs past the tRNA end (task D)

**Genes.** The 47 Ile tRNA genes of the RefSeq annotation (tRNA features, product tRNA-Ile: anticodons AAU, GAU, UAU).
10 lie inside a SINE copy of the run (a hit longer than 110 bp; `ile_in_sine.bed`) and were left out; 37 remain. Each was
extracted with 100 bp downstream (strand-aware; gene upper case, downstream lower case) and compared with the 17 Rhin
consensuses (ssearch36 -3 -m 8, E <= 0.001), against a null with the downstream 100 bp shuffled (5 shuffles per gene).
`extension.tsv`: per SINE x gene, how many bp past the tRNA 3' end the local alignment reaches on the genomic sequence.

**Classes** (`gene_classes.tsv`, `d_class.py`): a genuine Pol III tRNA gene ends in a run of >= 4 T (terminator) soon
after the gene; a SINE copy ends in its own body and an A-rich tail.
- genuine (T-run 2-27 bp after the end): **9 of 37**
- SINE-like (A-rich tail within 40 bp, no T-run): **21**
- neither: **7**

**Answer.** The similarity past the tRNA end comes from annotated "tRNA genes" that are r9 SINE copies. All 28 non-genuine
loci align best with r9 (104-128 bits) and run 24-27 bp past the tRNA end, which is the whole of r9 after its tRNA part
(r9 = 105 bp: tRNA-Ile-like 74 bp, then CACTGTGAGCTGCGCCCTCC... and the A-tail). Example: RefSeq TRNAI-GAU,
NC_142501.1:54295524-54295597, is 97 % identical to r9 over 101 nt including the stretch after the tRNA and the A-tail.
Median distance past the tRNA end, all SINE x gene pairs: SINE-like loci 22 bp (r9 27, r7 36, r8 36, r6 37), genuine genes
**3 bp** (max 35), shuffled downstream 4-8 bp (median of the per-pair maxima). So on genuine Ile genes the SINEs do not
continue past the tRNA end beyond chance, and no single genomic tRNA locus stands out as a source gene that matches the
SINE body further. Of the 9 genuine genes, 6 are Ile-AAU/GAU and 3 are Ile-UAU, which the SINEs barely match (two without a hit at E <= 0.001, one at 32 bits).
The 110-bp rule kept r9-length copies in the gene set (r9 is 105 bp), which is how they showed up; RefSeq's annotation counts
them as tRNA genes.

`rsi_ile_tRNA_genes_dn100_vs_SINEs.aln.fa`: the 17 Rhin consensuses, then the 37 genes + 100 bp, named
`GENUINE|...`, `noTrun_noAtail|...` or `SINElike|...` (MAFFT `--localpair --maxiterate 1000 --threadit 0 --preservecase`).
Scripts: `d_trna.sh`, `d_class.py`.

## `r10_all/`: every r10_19seqs copy (task E)

`rsi_r10_19seqs_all_SINE_grey.aln.fa`: the consensus as row 1, then all **593** copies with verdict SINE and all **168**
grey copies (marked `[grey]`, owner rule 2026-10-18). The **18** copies with verdict not SINE (shared flanks) are left out,
as copies with shared flanks are never shown. Elements aligned with MAFFT `--auto --threadit 0`; each copy's own 50 bp left
and 70 bp right flank (lower case) attached packed against the element (5' right-justified, 3' left-justified), so flank
columns hold no gaps between flank bases. `counts.txt`: verdict counts. Script `scripts/e_r10.sh`.

## `unit_pieces/`: all pieces, including units never found alone (task F)

From `composite_units.tsv`: every stand-alone r-unit consensus (`unit@unit`) and every stretch of the composites assigned to
a unit (`unit@composite:start-end`); unlabelled stretches of >= 30 bp as `unassigned@...`; MEG excluded (no composites).
- `rsi_all_unit_pieces.aln.fa`: all 31 pieces, ordered by unit (r1 ... r10, then unassigned).
- `by_unit/<unit>.aln.fa`: one alignment per unit with at least two pieces (r1, r2, r3, r4, r5, r6, r8, r10, unassigned;
  r7 and r9 occur only alone).
MAFFT `--localpair --maxiterate 1000 --threadit 0`, input order. Script `scripts/f_pieces.sh`.

### `unit_pieces/rsi_unit_pieces_r1module.aln.fa` (+ `.bed`): r1's first B-box module set apart (2026-10-23, owner's layout)
The same 31 pieces. r1 carries the B-box module twice (27-75 and ~76-133; 29 of 42 columns identical), every other unit
once. Here all other pieces were aligned (MAFFT --localpair --maxiterate 1000) against r1 without its first module
(r1 1-26 + 76-155), so their module sits under r1's second one; r1's first module (49 bp) was then put back as its own
column block, gaps in all other rows. The other three r1 pieces were added with `mafft --add --keeplength` (one base of
r1_r3_P1:6-161 falls outside the columns and is dropped; named in its header). Track: r1's 5' part, module 1, module 2.
Script `scripts/shift_r1_module.py` (run on therioserver `/home/toki/tmp/shift/`).


### `unit_pieces/rsi_unit_pieces_r1_two_modules.aln.fa` (+ `.bed`): each piece under the r1 module it is like (2026-10-25)
The same 31 pieces. r1 1-80 is an r10-type head and r1 81-155 an r2/r4-type element (`../../rsi_main_20261022/README.md`,
"Which r1 module is r2's?"), so the two modules are two column blocks and every piece goes under the one it scores
higher against (local alignment to r1 27-80 and to r1 81-133):
- module 1 (r10 type, unshifted): r5, r7, r8, r9, r10 and their composite pieces;
- module 2 (r2/r4 type, shifted right): r2, r3, r4 and their pieces; r6 scores the same against both (its module is half
  one type, half the other) and is put here, marked "tie";
- no module (both scores below 18): five short tail pieces and P18's middle piece; they sit only in the tail columns.
The 5' parts of all pieces and the tails of all pieces share columns (realigned together); the r3, r4 and r6 rows carry
an extra ~45 bp after their module (`CTGCGCCCCCTGCAACTA...GAGCTGAGCTG`) before the tail shared with the r10-type rows.
Track: r1's 5' part, module 1, module 2. Script `scripts/r1_two_modules.py` (run on KIT `/data/V/toki/temp/r1two/`,
therioserver was unreachable).
