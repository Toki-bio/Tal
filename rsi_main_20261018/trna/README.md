# rsi: tRNA / 5S / 7SL origin of the consensuses and of the composites' parts (2026-10-18)

Queries (`queries.fa`, 42): the 21 consensuses of this run and every stretch of the 7 composites as read by
`composite_units.tsv` (units; "-" = no unit). Script `trna_rsi.sh` (therioserver `/home/toki/rhin/rsi_main_20261018/trna/`),
summary `trna_summary.tsv` (`summary.py`). Method chosen after the deep research "tRNA similarity scans for SINEs"
(knowlege `methods/trna-similarity-scans-for-sines.md`):

1. **ssearch36 against the 419 human tRNA genes** (GtRNAdb hg19): best tRNA and its identity / length; SINEBase's rule
   (>= 60 % identity over >= 60 nt); the best hit of a *different* amino acid (2nd isotype) for the margin.
2. **Null:** 20 dinucleotide-preserving shuffles of every query against the same tRNAs; the best shuffled bit score.
3. **Infernal cmscan, all of Rfam** (E <= 1) and **tRNAscan-SE -E -H --detail -X 10** (isotype-model scores, pseudo /
   IPD flags).

## Result
- **Every r-unit has a tRNA-Ile-derived head** at its 5' end: 73-93 % identity over 66-75 nt to human tRNA-Ile (AAT/GAT),
  bit score 39-82 against at most 29 for 20 shuffles of the same sequence; the next amino acid (Lys, Thr, Asn) is 8-37
  bits lower, so the isotype call is clear. tRNAscan's isotype model says Ile too wherever it scores (anticodon calls -
  Ser, Arg, Lys, Met, Leu - disagree: "IPD", the mark of a decayed tRNA-derived repeat). r9 is the best preserved (93 %,
  tRNAscan Ile-AAT 65.4, not pseudo); r7 / r8 86.5 %; the rest 73-79 %.
- **Detectors differ, as the research predicted:** Rfam misses r2, r3, r4 and r6; tRNAscan misses r3 and r4; the
  pairwise search against tRNA genes finds all of them (Dfam was needed for this in the 2026-09-29 run).
- **Composites carry one tRNA head per unit:**
  - r1_r3_P1: Ile head at 6-80 (r1) and at 201-267 (r3); the 39 bp between (162-200) and the 3' end (402-448) match no RNA.
  - r1_r3_P2: Ile head at 6-80 (r1); its r3 part matches only 34 bp here because the unit boundary from
    composite_units cuts the r3 head at 86 (the head lies at 81-150 in P2).
  - r1_r2_r3_r3_C11: three Ile heads, at 6 (r1), 193 (r2) and 339 (r4); the last part (472-589, "r3") has no head.
  - r5_r3_P34: Ile heads at 6 (r5) and at 141 (r3). r5_r6_P26 and r5_r5_P48: Ile heads at 6 (r5) and 142 (r6).
  - r10_r8_P18: Ile head at 6-79 (r10) only; the 141 bp middle and the r8-like part (78 %) have no tRNA head.
- **MEG controls as expected** (Gogolevsky et al. 2009): MEG-RL / MEG-RS = 5S rRNA (Rfam 1-119, E 1e-24);
  MEG-T2 = tRNA-Val (94.5 %, tRNAscan Val-TAC 70.3, not pseudo); MEG-TR = tRNA-Val (1-55) + 5S (65-180).
- No 7SL-related part anywhere (Rfam).
